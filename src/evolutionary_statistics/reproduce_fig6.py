#!/usr/bin/env python
"""Reproduce Figure 6 (+ Supplementary Figure S6) and their Source Data.

Faithful transcription of the evolutionary-statistics analysis contributed by
Jose Alberto de la Paz (Morcos lab, UT Dallas), `Figure6_Analysis.ipynb`, kept
verbatim alongside this script for provenance. Direct-coupling analysis (mean
field) of the Pfam PF13354 family constrains the combinatorial TEM-1 variants:
panel 6A is the family sequence logo at the nine mutated positions present in the
family; 6B the contextual effective alphabet; 6C/6G/6H the family Hamiltonian
distributions of the top performers; S6 the Hamiltonian-versus-fitness relation.

Runs on our shipped AUC-fitness parquets (`data/processed/{amp,azt}_auc_wide_df.parquet`)
and the committed PF13354 MSA; writes the plotted values to `source_data/`.
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.stats import spearmanr

WT_SEQ = "LQMERMAGERTRN"
# Nine positions that vary within PF13354 (subset of the 13 mutated in the library),
# expressed on the TEM-1 263-aa reference before back-mapping to MSA columns.
FAMILY_POSITIONS_PDB = np.array([44, 79, 139, 157, 212, 213, 214, 218, 238])
LABEL_POSITIONS = ["69", "104", "164", "182", "237", "238", "240", "244", "265"]
WT_REF9 = "MERMAGERT"  # WT residues at the nine family positions (Ambler order above)
INTENDED_MUT = {
    0: ["L", "V"], 1: ["K"], 2: ["S", "H", "N"], 3: ["T"], 4: ["T"],
    5: ["S"], 6: ["K"], 7: ["S", "C"], 8: ["M"],
}


def resolve_msa(msa_dir: Path) -> Path:
    """Return the PF13354 family MSA, decompressing a shipped .gz once if needed."""
    fasta = msa_dir / "PF13354_ga.fasta"
    if fasta.exists():
        return fasta
    gz = msa_dir / "PF13354_ga.fasta.gz"
    if gz.exists():
        with gzip.open(gz, "rb") as fi, open(fasta, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        return fasta
    raise FileNotFoundError(f"PF13354_ga.fasta[.gz] not found in {msa_dir}")


def build_average(parquet: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """Per-concentration median AUC-fitness for one drug (Alberto's compute_median)."""
    groups: dict[str, list[str]] = {}
    for col in parquet.columns:
        if not col.startswith(prefix):
            continue
        groups.setdefault(col.split()[1], []).append(col)
    medians = {c: parquet[cols].median(axis=1, skipna=True) for c, cols in groups.items()}
    n_mut = parquet["mut_profile_masked"].apply(lambda p: sum(ch.isalpha() for ch in p))
    df = pd.DataFrame({"Code": parquet[""], "Profile": parquet["mut_profile_masked"],
                       "N. Mut.": n_mut, **medians})
    df["Reduced_Seq"] = df["Profile"].apply(
        lambda p: "".join(WT_SEQ[i] if p[i] == "." else p[i] for i in range(len(WT_SEQ))))
    return df


def hmm_backmap(msa_dir: Path) -> dict[int, int]:
    """PDB->MSA column map from an hmmscan of the TEM-1 reference against PF13354.

    Regenerates the scan with HMMER when available; otherwise uses the committed
    TEM1_263AA.scan (deterministic — the reference-to-profile alignment is fixed).
    """
    from dca import helper_functions

    hmm_file, scan_file = msa_dir / "PF13354.hmm", msa_dir / "TEM1_263AA.scan"
    pdb_fasta, align_file = msa_dir / "TEM1_263AA.fasta", msa_dir / "TEM1_263AA.align"
    if shutil.which("hmmscan") and shutil.which("hmmpress"):
        subprocess.run(["hmmpress", "-f", str(hmm_file)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run(["hmmscan", "-o", str(scan_file), "--notextw",
                        str(hmm_file), str(pdb_fasta)], check=True)
    elif not scan_file.exists():
        raise RuntimeError("HMMER not on PATH and no cached TEM1_263AA.scan")

    important, counter, record = [], 0, False
    for line in scan_file.read_text().splitlines():
        if "domain" in line and "score" in line and "E-value" in line:
            record = True
        elif record and counter < 4:
            try:
                int(line.split()[1]); important.append(line.rstrip())
            except (ValueError, IndexError):
                pass
            counter += 1
    with open(align_file, "w") as out:
        for alignment in important[:2]:
            out.write("\n".join(alignment.split()) + "\n\n")
    fwd = helper_functions.backmap_alignment(str(align_file))
    return {v: k for k, v in fwd.items()}


def main() -> None:
    here = Path(__file__).resolve().parent
    repo = here.parents[1]  # src/evolutionary_statistics -> src -> repo root
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", default=repo / "data" / "processed", type=Path,
                    help="dir with amp/azt_auc_wide_df.parquet + intended CSV")
    ap.add_argument("--intended", default=None, type=Path,
                    help="TEM1-combinatorial-mutagenesis-intended.csv (default: --data/)")
    ap.add_argument("--msa", default=here / "MSAs", type=Path)
    ap.add_argument("--figures", default=here / "figures", type=Path)
    ap.add_argument("--source-data", default=here / "source_data", type=Path)
    args = ap.parse_args()
    # Large per-genotype tables (6C/6G/6H, S6 scatter) go under source_data/full/ — git-ignored,
    # regenerable, packaged into the release Source Data zip. Small summaries stay at top level.
    full = args.source_data / "full"
    for d in (args.figures, args.source_data, full):
        d.mkdir(parents=True, exist_ok=True)

    import logomaker as lm
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    from dca import helper_functions
    from dca.dca_class import dca

    intended_csv = args.intended or (args.data / "TEM1-combinatorial-mutagenesis-intended.csv")

    # ---- data: per-concentration median fitness, restricted to intended designs ----
    amp = build_average(pd.read_parquet(args.data / "amp_auc_wide_df.parquet"), "Ampicillin")
    azt = build_average(pd.read_parquet(args.data / "azt_auc_wide_df.parquet"), "Aztreonam")
    conc_amp = amp.columns[3:9].tolist()
    conc_azt = azt.columns[3:11].tolist()
    intended = pd.read_csv(intended_csv)["mut_profile_masked"]
    amp = amp[amp["Profile"].isin(intended)].copy()
    azt = azt[azt["Profile"].isin(intended)].copy()
    print(f"Cleaned selected AMP={len(amp)}  AZT={len(azt)}")

    # ---- MSA filter (<=5% gaps) + PDB->MSA back-map + nine family positions ----
    pfam = resolve_msa(args.msa)
    msa_len = len(next(SeqIO.parse(pfam, "fasta")))
    max_gaps = round(5 * msa_len / 100)
    pfam_filt = pfam.with_name(pfam.stem + "_filtered_5pct.fasta")
    helper_functions.filter_pfam(str(pfam), max_gaps, str(pfam_filt))
    back = hmm_backmap(args.msa)
    positions = np.array([back[p] for p in FAMILY_POSITIONS_PDB])  # MSA columns (1-based)

    # ---- 6A: family sequence logo at the nine positions ----
    family_red = ["".join(str(rec.seq)[i - 1] for i in positions)
                  for rec in SeqIO.parse(pfam_filt, "fasta")]
    counts = lm.alignment_to_matrix(family_red, to_type="information", pseudocount=0.01)
    fig, ax = plt.subplots(figsize=(10, 10))
    logo = lm.Logo(counts, ax=ax)
    pos_rank = {i: counts.loc[i].sort_values().index.get_loc(WT_REF9[i]) for i in range(9)}
    mut_rank = {i: [counts.loc[i].sort_values().index.get_loc(m) for m in INTENDED_MUT[i]]
                for i in range(9)}
    ncol = counts.shape[1]
    for k, patch in enumerate(logo.ax.patches):
        pos, rank = k // ncol, k % ncol
        patch.set_facecolor("black" if rank == pos_rank[pos] else "#808080")
        if rank in mut_rank[pos]:
            patch.set_facecolor("#FF0000")
    logo.style_spines(visible=False)
    logo.style_spines(spines=["left", "bottom"], visible=True)
    ax.set_title("PF13354 Logo"); ax.set_ylim(0, np.log2(20))
    logo.ax.set_xticklabels(LABEL_POSITIONS)
    ax.hlines(np.log2(20), -0.5, 9, color="black", ls="dashed")
    fig.savefig(args.figures / "FamilyLogo_6a.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    counts.insert(0, "position", LABEL_POSITIONS)
    counts.to_csv(args.source_data / "fig6a_family_logo_information.csv", index=False)

    # ---- DCA mean field + Hamiltonian / effective alphabet ----
    model = dca(str(pfam_filt)); model.mean_field()
    print(f"M={model.M} N={model.N} Meff={model.Meff:.6f}")
    ham, headers = model.compute_Hamiltonian(str(pfam_filt))
    ham = np.asarray(ham)
    tem1_idx = next(i for i, n in enumerate(headers) if "P62593" in n)
    tem1_ham, tem1_name = float(ham[tem1_idx]), headers[tem1_idx]
    tem1_record = next(r for r in SeqIO.parse(pfam_filt, "fasta") if tem1_name in r.description)
    single = args.msa / "TEM1_single_MSA.fasta"
    SeqIO.write([tem1_record], single, "fasta")
    tem1_eff = np.asarray(model.compute_EffAlphabet(str(single)))[0]
    fam_eff = np.asarray(model.compute_EffAlphabet(str(pfam_filt)))  # (M, N)
    n_fam = fam_eff.shape[0]

    # ---- 6B: contextual effective alphabet at the nine positions ----
    fam_mean = fam_eff[:, positions - 1].mean(axis=0)
    fam_std = fam_eff[:, positions - 1].std(axis=0)
    tem1_target = tem1_eff[positions - 1]
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.errorbar(range(9), fam_mean, yerr=fam_std, fmt="o", ms=10,
                label=f"Family Average (n={n_fam})")
    ax.bar(LABEL_POSITIONS, tem1_target, label="Target Positions", color="#f4a582")
    ax.set(title="Contextual Mutability", xlabel="Position",
           ylabel="Effective Alphabet", ylim=(0, 4.5))
    ax.legend(loc="upper left")
    fig.savefig(args.figures / "EffAlphabet_6b.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame({"position": LABEL_POSITIONS, "family_mean_eff_alphabet": fam_mean,
                  "family_sd": fam_std, "tem1_target_eff_alphabet": tem1_target,
                  "n_family_sequences": n_fam}).to_csv(
        args.source_data / "fig6b_effective_alphabet.csv", index=False)

    # ---- synthetic-variant Hamiltonians against the family model ----
    def seq_wrt(profile: str) -> str:
        base = list(str(tem1_record.seq))
        for pos, ch in zip(positions, profile[2:]):
            if ch != ".":
                base[pos - 1] = ch
        return "".join(base)

    for df, tag, key in ((amp, "AMP", "781.0"), (azt, "AZT", "36.0")):
        df["Seq_wrtPF13354"] = df["Profile"].apply(seq_wrt)
        recs = [SeqRecord(Seq(s), id=str(c), description="")
                for c, s in zip(df["Code"], df["Seq_wrtPF13354"])]
        syn = args.msa / f"{tag}_MSA_wrtPF13354.fasta"
        SeqIO.write(recs, syn, "fasta")
        df["H_fam"] = np.asarray(model.compute_Hamiltonian(str(syn))[0])

    # ---- 6G / 6H: Hamiltonian of all vs top-2.5k performers ----
    bins = np.linspace(-1870, -1740, 51)
    for df, tag, key, drug in ((amp, "AMP", "781.0", "Ampicillin"),
                               (azt, "AZT", "36.0", "Aztreonam")):
        top = df.nlargest(2500, key).index
        df["top2500"] = df.index.isin(top)
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.histplot(df["H_fam"], bins=bins, color="#000000", label="All", stat="density", ax=ax)
        sns.histplot(df.loc[df["top2500"], "H_fam"], bins=bins, color="#808080",
                     label=f"Top 2.5k {tag}", stat="density", ax=ax)
        ax.set(xlabel="Hamiltonian", ylabel="Relative Frequency",
               title=f"Best 2500 performers at {drug}"); ax.legend()
        fig.savefig(args.figures / f"Hfam_top2.5k_{tag}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        df[["Code", "Profile", "N. Mut.", key, "H_fam", "top2500"]].to_csv(
            full / f"fig6{'g' if tag == 'AMP' else 'h'}_hamiltonian_{tag.lower()}.csv",
            index=False)

    # ---- 6C: Hamiltonian by dataset (all / top-10k common / top-2.5k common) ----
    keys = ["Profile", "Reduced_Seq", "N. Mut.", "H_fam"]
    top10 = pd.merge(amp.nlargest(10000, "781.0"), azt.nlargest(10000, "36.0"),
                     on=keys, suffixes=("_Amp", "_Azt"))
    top25 = pd.merge(amp.nlargest(2500, "781.0"), azt.nlargest(2500, "36.0"),
                     on=keys, suffixes=("_Amp", "_Azt"))
    print(f"Common top10k={len(top10)}  top2.5k={len(top25)}")
    all_v = amp[["Profile", "H_fam"]].assign(Dataset="All Variants")
    combined = pd.concat([all_v, top10[["Profile", "H_fam"]].assign(Dataset="Top 10k Common"),
                          top25[["Profile", "H_fam"]].assign(Dataset="Top 2.5k Common")],
                         ignore_index=True)
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    for ax, lab in zip(axes, ["Top 2.5k Common", "Top 10k Common", "All Variants"]):
        sns.histplot(combined[combined["Dataset"] == lab]["H_fam"], bins=bins,
                     stat="count", ax=ax, label=lab); ax.legend(); ax.set(xlim=(-1870, -1740))
    axes[-1].set_xlabel("Hamiltonian")
    fig.savefig(args.figures / "Hfam_by_dataset_6c.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    combined.to_csv(full / "fig6c_hamiltonian_by_dataset.csv", index=False)

    # ---- S6: Hamiltonian vs fitness, per concentration, with Spearman rho ----
    s6_stats = []
    for df, tag, concs in ((amp, "AMP", conc_amp), (azt, "AZT", conc_azt)):
        cols = ["Profile", "N. Mut.", "H_fam", *concs]
        df[cols].to_csv(full / f"figS6_hamiltonian_vs_fitness_{tag.lower()}.csv",
                        index=False)
        for conc in concs:
            valid = df[["H_fam", conc]].dropna()
            rho, p = spearmanr(valid["H_fam"], valid[conc])
            s6_stats.append({"drug": tag, "concentration": conc, "spearman_rho": rho, "p_value": p})
    pd.DataFrame(s6_stats).to_csv(args.source_data / "figS6_spearman_stats.csv", index=False)

    json.dump({"M_filtered_family_sequences": int(model.M), "N_positions": int(model.N),
               "Meff": float(model.Meff), "TEM1_Hamiltonian": tem1_ham,
               "n_family_for_6B": int(n_fam), "TEM1_name": tem1_name,
               "family_positions_msa_cols": positions.tolist()},
              open(args.source_data / "fig6_dca_stats.json", "w"), indent=2)
    print(f"TEM1_Hamiltonian={tem1_ham:.6f}  n_family_for_6B={n_fam}")
    print(f"Source Data -> {args.source_data}")


if __name__ == "__main__":
    main()

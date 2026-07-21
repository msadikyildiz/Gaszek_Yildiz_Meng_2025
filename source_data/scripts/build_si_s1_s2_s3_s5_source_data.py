"""Regenerate the Source Data CSVs for Supplementary Figures S1(B), S2, S3 and S5.

Writes tables into ``source_data/derived/`` (git-ignored; assembled into
``source_data.xlsx`` by ``build_source_data.py``). Each reproduces exactly the
values plotted by the corresponding figure generator:

  S1B minimum library sequence-logo counts -- per-position amino-acid counts across
      all 55,296 library genotypes (the panel-B logo of
      ``TEM1-combinatorial-mutagenesis-intended.csv``). Panel A (the clinical
      mutational spectrum) is an external NCBI beta-lactamase analysis described
      in Data Availability and is not regenerated in this repository.
  S2  minimum inhibitory concentrations for the ten strains plotted by
      ``si_figures/s02_mic/CML_MIC_estimations.ipynb`` (its ``strain_order``,
      PASS status; the dose-range-adjusted re-run rows are excluded exactly as
      the figure excludes them).
  S3  the highlighted dose-response trajectories of
      ``si_figures/s03_dose_response/analysis.py`` -- the top-18 highest-AUC
      genotypes plus TEM-1 WT and the dead control, per-replicate AUC (n=3) at
      every concentration, per drug.
  S5  the context-dependent triple-mutant panel of
      ``src/05_epistasis_figures.ipynb`` -- mean AUC-fitness +/- s.d. for the
      eight E104K/R164S/E240K combinations on four backgrounds, at the two
      main-text concentrations, read from ``Epistasis_Combined.parquet``.

(Supplementary Figure S4 graph node/edge tables are produced separately by
``src/graph_analysis/graph_builder/extract_figS4_source_data.py``.)
"""
from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def _repo_root(p: Path) -> Path:
    for a in (p, *p.parents):
        if (a / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return a
    raise FileNotFoundError("repo root not found from " + str(p))


REPO = _repo_root(HERE)
OUT = HERE.parent / "derived"
OUT.mkdir(exist_ok=True)

WT = "." * 13
DEAD = "X" * 13
AMBLER = [21, 39, 69, 104, 164, 182, 237, 238, 240, 244, 265, 275, 276]


def build_s1() -> None:
    """S1 panel B -- per-position amino-acid counts (the library sequence logo)."""
    src = REPO / "data" / "processed" / "TEM1-combinatorial-mutagenesis-intended.csv"
    prof = pd.read_csv(src)["mut_profile"].astype(str)
    rows = []
    for i, amb in enumerate(AMBLER):
        for aa, n in prof.str[i].value_counts().items():
            rows.append({"position_ambler": amb, "amino_acid": aa, "count": int(n)})
    out = pd.DataFrame(rows).sort_values(["position_ambler", "count"], ascending=[True, False])
    out.to_csv(OUT / "figS1B_library_logo_counts.csv", index=False)
    print(f"  S1B: {len(out)} rows -> figS1B_library_logo_counts.csv")


def build_s2() -> None:
    """S2 -- MIC panel: the ten plotted strains, PASS rows only."""
    strain_order = [
        "NEB10B", "TEM-1", "TEM-1-CML",
        "E104K-R164N-M182T-E240K", "L21P-E104K-R164N-E240K-T265M",
        "Q39K-E104K-R164S-E240K-T265M", "E104K-R164N-E240K-T265M",
        "Q39K-E104K-R164N-E240K", "Q39K-E104K-R164N-M182T-E240K",
        "L21P-E104K-R164N-M182T-E240K",
    ]
    keep = ["Strain", "Antibiotic", "MIC", "MIC_ci_lower", "MIC_ci_upper",
            "IC50", "IC50_ci_lower", "IC50_ci_upper", "insufficient_drug"]
    src = REPO / "si_figures" / "s02_mic" / "data" / "growth_features.csv"
    df = pd.read_csv(src)
    sel = df[df["Strain"].isin(strain_order) & (df["Status"] == "PASS") & df["MIC"].notna()].copy()
    sel["_o"] = sel["Strain"].map({s: i for i, s in enumerate(strain_order)})
    sel = sel.sort_values(["_o", "Antibiotic"])[keep]
    sel.to_csv(OUT / "figS2_mic_values.csv", index=False)
    print(f"  S2: {len(sel)} plotted MIC points -> figS2_mic_values.csv")


def build_s3() -> None:
    """S3 -- dose-response: top-18 + WT + dead, per-replicate AUC per concentration."""
    n_top = 18
    drug_files = {
        "Ampicillin": REPO / "data" / "raw" / "Ampicillin_auc_per_genotype.csv",
        "Aztreonam": REPO / "data" / "raw" / "Aztreonam_auc_per_genotype.csv",
    }
    rows = []
    for drug, path in drug_files.items():
        df = pd.read_csv(path, index_col=0)  # already sorted descending by total AUC
        top = df.index[:n_top].tolist()
        wt = df.index[df["mut_profile_masked"] == WT].tolist()
        dead = df.index[df["mut_profile_masked"] == DEAD].tolist()
        rank = {g: i + 1 for i, g in enumerate(df.loc[top, "mut_profile_masked"])}
        value_cols = [c for c in df.columns if c.startswith(drug + " ")]
        for ridx in top + wt + dead:
            geno = df.loc[ridx, "mut_profile_masked"]
            role = "WT" if geno == WT else ("dead" if geno == DEAD else f"top{rank.get(geno, '')}")
            for col in value_cols:
                _, conc, rep = col.split(" ")
                rows.append({"drug": drug, "genotype": geno, "role": role,
                             "concentration": float(conc), "replicate": int(rep),
                             "auc": df.loc[ridx, col]})
    pd.DataFrame(rows).to_csv(OUT / "figS3_dose_response_highlighted.csv", index=False)
    print(f"  S3: {len(rows)} rows -> figS3_dose_response_highlighted.csv")


def build_s5() -> None:
    """S5 -- context-dependent E104K/R164S/E240K on four backgrounds (AMP 781, AZT 36)."""
    mut = ["E104K", "R164S", "E240K"]
    backgrounds = {"WT": "LQMERMAGERTRN", "G238S": "LQMERMASERTRN",
                   "Q39K": "LKMERMAGERTRN", "A237T": "LQMERMTGERTRN"}
    drugs = {"AMP": 781.0, "AZT": 36.0}

    def gen_combos(baseline: str):
        out = [{"Mutant": "Baseline", "Epistatic Order": 0, "Genotype": baseline}]
        for r in range(1, len(mut) + 1):
            for combo in itertools.combinations(mut, r):
                pos = [np.searchsorted(AMBLER, int("".join(filter(str.isdigit, m)))) for m in combo]
                if len(set(pos)) < len(pos):
                    continue
                g = list(baseline)
                for m, idx in zip(combo, pos):
                    g[idx] = m[-1]
                out.append({"Mutant": " + ".join(combo), "Epistatic Order": r, "Genotype": "".join(g)})
        return out

    df = pd.read_parquet(REPO / "data" / "processed" / "Epistasis_Combined.parquet")
    rows, missing = [], []
    for drug, conc in drugs.items():
        look = df[(df.Drug == drug) & (df.Concentration == conc)].drop_duplicates("Genotype").set_index("Genotype")
        for bgname, bgseq in backgrounds.items():
            for c in gen_combos(bgseq):
                g = c["Genotype"]
                if g in look.index:
                    rows.append({"drug": drug, "concentration": conc, "background": bgname,
                                 "mutations": c["Mutant"], "epistatic_order": c["Epistatic Order"],
                                 "genotype": g, "fitness_mean": round(float(look.loc[g, "Fitness"]), 6),
                                 "fitness_sd": round(float(look.loc[g, "Error"]), 6)})
                else:
                    missing.append((drug, bgname, c["Mutant"]))
    pd.DataFrame(rows).to_csv(OUT / "figS5_context_dependent_triple.csv", index=False)
    print(f"  S5: {len(rows)} rows -> figS5_context_dependent_triple.csv" +
          (f"  (WARNING {len(missing)} missing)" if missing else ""))


if __name__ == "__main__":
    build_s1()
    build_s2()
    build_s3()
    build_s5()

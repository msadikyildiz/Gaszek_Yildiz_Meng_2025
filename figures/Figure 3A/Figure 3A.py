import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D

# ── 1) GLOBAL STYLING ─────────────────────────────────────────────
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         24,
    "axes.titlesize":    23,
    "axes.labelsize":    24,
    "xtick.labelsize":   24,
    "ytick.labelsize":   24,
    "legend.fontsize":   24,
    "axes.linewidth":    1.2,
    "lines.linewidth":   1.5,
    "grid.color":        "0.5",
    "grid.linestyle":    "--",
    "grid.linewidth":    0.5,
    "grid.alpha":        0.3,
})

PROJECT_PATH = Path(__file__).resolve().parents[2]
amp_concs = [195.0, 781.0]
azt_concs = [36.0, 108.0, 324.0]

df        = pd.read_parquet(PROJECT_PATH/"data/processed"/"Epistasis_Combined.parquet")
amp_stats = pd.read_parquet(PROJECT_PATH/"data/processed"/"amp_auc_long_df.parquet")
azt_stats = pd.read_parquet(PROJECT_PATH/"data/processed"/"azt_auc_long_df.parquet")

dead_seq = 'XXXXXXXXXXXXX'
wt_seq   = '.............'

colors = {
    'dead':     '#0072B2',
    'wt':       '#000000',
    'tested':   '#FF0000',
    'clinical': '#FFFF00',
}

legend_handles = [
    Line2D([0],[0], marker='o', ls='none',
           color=colors['dead'],   label='Dead variant',
           markerfacecolor='none', ms=14, mew=3.5),
    Line2D([0],[0], marker='o', ls='none',
           color=colors['wt'],     label='Wild-type',
           markerfacecolor='none', ms=14, mew=3.5),
    Line2D([0],[0], marker='D', ls='none',
           color=colors['tested'], label='Tested variants', ms=12, mec='k', mew=1),
    Line2D([0],[0], marker='^', ls='none',
           color=colors['clinical'],label='Clinical isolates',   ms=13, mec='k', mew=1),
]

extra_variants = {
  'c1.1':'LQMKNTAGKRTRN','c1.2':'PQMKNMAGKRMRN','c1.3':'LKMKSMAGKRMRN',
  'c2.1':'LQMKNMAGKRMRN','c2.2':'LKMKNMAGKRTRN','c2.3':'LKMKNTAGKRTRN',
  'c3.1':'PQMKNTAGKRTRN' 
}

try:
    clin = pd.read_csv(
      PROJECT_PATH/"figures"/"known_variants"/"encoded_variants.csv"
    )['Encoded_Sequence'].dropna().unique()
except FileNotFoundError:
    clin = []

# ── 2) SETUP FIGURE ────────────────────────────────────────────────
fig, axes = plt.subplots(
    nrows=2, ncols=3,
    figsize=(20,14),
    constrained_layout=True,
)

# ── 3) LOOP & PLOT ─────────────────────────────────────────────────
for i, a_conc in enumerate(amp_concs):
    for j, z_conc in enumerate(azt_concs):
        ax = axes[i, j]

        # 3a) 2D histogram, with a bit of transparency
        amp_fit = df.query("Drug=='AMP'  & Concentration==@a_conc").Fitness
        azt_fit = df.query("Drug=='AZT'  & Concentration==@z_conc").Fitness
        hb = ax.hist2d(
            amp_fit, azt_fit,
            bins=[np.linspace(0,4.5,50), np.linspace(.95,6.1,50)],
            cmap='BuPu', norm=LogNorm(),
            alpha=0.75
        )

        # 3b) plot dead–WT reference line
        wt_a   = amp_stats.query("mutant_profile==@wt_seq and concentration==@a_conc")['median'].iat[0]
        dead_a = amp_stats.query("mutant_profile==@dead_seq and concentration==@a_conc")['median'].iat[0]
        wt_z   = azt_stats.query("mutant_profile==@wt_seq and concentration==@z_conc")['median'].iat[0]
        dead_z = azt_stats.query("mutant_profile==@dead_seq and concentration==@z_conc")['median'].iat[0]

        # ax.plot([dead_a, wt_a], [dead_z, wt_z],
        #         'k--', lw=2)
        ax.plot(dead_a, dead_z,
                'o', color=colors['dead'],
                ms=12, mew=3, mfc='none', zorder=200)
        ax.plot(wt_a, wt_z,
                'o', color=colors['wt'],
                ms=12, mew=3, mfc='none', zorder=200)

        # 3c) engineered variants
        ex_a, ex_z = [], []
        for seq in extra_variants.values():
            a = df.query("Drug=='AMP' & Concentration==@a_conc & Genotype==@seq").Fitness
            z = df.query("Drug=='AZT' & Concentration==@z_conc & Genotype==@seq").Fitness
            if not a.empty and not z.empty:
                ex_a.append(a.iat[0])
                ex_z.append(z.iat[0])
        ax.plot(ex_a, ex_z,
                'D', color=colors['tested'],
                mec='k', ms=10, alpha=0.8)

        # 3d) clinical variants
        cl_a, cl_z = [], []
        for seq in clin:
            a = df.query("Drug=='AMP' & Concentration==@a_conc & Genotype==@seq").Fitness
            z = df.query("Drug=='AZT' & Concentration==@z_conc & Genotype==@seq").Fitness
            if not a.empty and not z.empty:
                cl_a.append(a.iat[0])
                cl_z.append(z.iat[0])
        ax.plot(cl_a, cl_z,
                '^', color=colors['clinical'],
                mec='k', ms=11, alpha=0.8)

        # 3e) axes, labels, grid
        ax.set_xlim(0,4.6)
        ax.set_ylim(.95,6.1)
        ax.set_xlabel(rf"Fitness in AMP {a_conc:.0f}$\mathrm{{\mu g/ml}}$")
        ax.set_ylabel(rf"Fitness in AZT {z_conc:.0f}$\mathrm{{\mu g/ml}}$")
        ax.grid(True)

        if j!=0:
            ax.set_yticklabels([])
        if i!=0:
            ax.set_xticklabels([])
        if i==0 and j==0:
            ax.legend(handles=legend_handles, loc='upper left')

# Add shared colorbar
cbar_ax = fig.add_axes([1.01, 0.4, 0.03, 0.3], zorder=100)
cbar = fig.colorbar(hb[3], cax=cbar_ax)
cbar.ax.set_ylabel('Mutant Count', labelpad=0)

# ── 4) SAVE ───────────────────────────────────────────────────────
fig.savefig(PROJECT_PATH/"figures"/"Figure 3A"/"Figure 3A. AMP vs. AZT fitness with clinical and tested variants.png", dpi=600, bbox_inches='tight')

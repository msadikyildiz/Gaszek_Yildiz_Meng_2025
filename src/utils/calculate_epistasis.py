#!/usr/bin/env python3

"""
Epistasis calculation module for TEM1 combinatorial mutagenesis.
Implements methods to calculate and analyze various forms of epistasis 
in TEM1 enzyme variants with multiple mutations.
"""

import polars as pl
import numpy as np
import pandas as pd
from itertools import product
import matplotlib.pyplot as plt
import torch
from pathlib import Path
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from math import ceil
from matplotlib.colors import LogNorm
from tqdm import tqdm
import time

# Import shared configuration and utility functions
from device_utils import get_device, clear_memory
from config import PROJECT_PATH, CACHE_DIR, RUN_ID, INTENDED, ALL_MUTATIONS, get_mutation_options

# Get the best available device (CPU or GPU)
device = get_device()

# Define mutation positions in TEM1 enzyme
sign_pos = [19, 37, 67, 102, 162, 180, 235, 236, 237, 241, 261, 271, 272]

def generate_binary_string(mutant, intended):
    """Convert a mutant string to a binary representation based on intended mutations."""
    return ''.join(str(intended[pos].index(mutant[i])) if mutant[i] != "." else "0" 
                  for i, pos in enumerate(sorted(intended.keys())))

def calculate_epistasis_biochemical(W, data, data_error):
    """Calculate epistasis using biochemical definition with weight matrix W."""
    data = torch.tensor(data, device=device, dtype=torch.float16)
    data_error = torch.tensor(data_error, device=device, dtype=torch.float16)
    
    epistasis = torch.matmul(W, data)
    epistasis_error = torch.sqrt(torch.matmul(W * W, data_error * data_error))
    
    return epistasis.cpu().numpy(), epistasis_error.cpu().numpy()

def calculate_Ensemble_Averaged_epistasis(W_ea, test_data, test_data_error):
    """Calculate epistasis using Ensemble Averaging method."""
    start_time = time.time()
    
    test_data = torch.tensor(test_data, device=device, dtype=torch.float16)
    test_data_error = torch.tensor(test_data_error, device=device, dtype=torch.float16)
    
    epistasis_ea = torch.matmul(W_ea, test_data)
    epistasis_ea_error = torch.sqrt(torch.matmul(W_ea * W_ea, test_data_error * test_data_error))
    
    return epistasis_ea.cpu().numpy(), epistasis_ea_error.cpu().numpy()

def calculate_epistasis_linear_regression(W_linregress, Test_data, Test_data_error):
    """Calculate epistasis using linear regression method."""
    Test_data = torch.tensor(Test_data, device=device, dtype=torch.float16)
    Test_data_error = torch.tensor(Test_data_error, device=device, dtype=torch.float16)
    
    epistasis_LG = torch.matmul(W_linregress, Test_data)
    epistasis_LG_error = torch.sqrt(torch.matmul(W_linregress * W_linregress, Test_data_error * Test_data_error))
    
    return epistasis_LG.cpu().numpy(), epistasis_LG_error.cpu().numpy()

def calculate_W_linregress(V, X, S, H):
    """Calculate linear regression weight matrix from components V, X, S, H."""
    return torch.matmul(V, torch.matmul(X.T, torch.matmul(S, H)))

def Epistasis_LG_partial_calculate(E_slice, V, X, H, Epistasis_Combined):
    """
    Calculate epistasis for different orders using partial linear regression.
    Progressively adds interaction terms up to each order.
    """
    epistatic_order = pd.to_numeric(Epistasis_Combined.loc[E_slice,'Epistatic Order']).to_numpy()
    epistatic_order = torch.tensor(epistatic_order, device=device, dtype=torch.float16)
    model_fitness = Epistasis_Combined.loc[E_slice, 'Fitness'].to_numpy()
    model_error = Epistasis_Combined.loc[E_slice,'Error'].to_numpy()

    max_order = int(epistatic_order.max())
    progress_bar = tqdm(range(1, max_order + 1), desc="Processing epistatic orders")
    
    for order in progress_bar:
        progress_bar.set_description(f"Processing order {order}/{max_order}")
        
        # Create mask for terms up to current order
        Q_partial = torch.eye(X.shape[0], dtype=torch.float16, device=device)
        
        mask = torch.nonzero(epistatic_order > order)
        Q_partial[mask, mask] = 0
        
        Q_partial = torch.matmul(Q_partial, Q_partial.T)
        
        # Calculate weights and epistasis for this order
        W_linregress_partial = calculate_W_linregress(V, X, Q_partial, H)
        
        epistasis_LG, epistasis_LG_error = calculate_epistasis_linear_regression(W_linregress_partial, model_fitness, model_error)
        
        # Predict fitness using partial model
        X_partial = torch.matmul(X, Q_partial)
        Fitness_predicted = torch.matmul(X_partial, torch.tensor(epistasis_LG.flatten(), device=device))
        
        # Store results in the dataframe
        Epistasis_Combined.loc[E_slice, f'Epistasis_LG for order {order}'] = epistasis_LG
        Epistasis_Combined.loc[E_slice, f'Error_LG for order {order}'] = epistasis_LG_error
        Epistasis_Combined.loc[E_slice, f'Fitness_predicted for order {order}'] = Fitness_predicted.cpu().numpy()
        
        del Q_partial, X_partial, W_linregress_partial
        
        # Clear memory to prevent GPU/CPU overflow
        clear_memory(device)

def plot_fitness_vs_predicted(df, Drug, gridsize=34, bin_range=(-0.5, 8)):
    """
    Plot heatmaps of measured vs predicted fitness for each concentration and epistatic order.
    Uses 2D histograms to show density of points.
    """
    slice = df['Drug'] == Drug
    concentrations = sorted(df.loc[slice, 'Concentration'].unique())
    max_order = df.loc[slice, 'Epistatic Order'].max()
    
    if pd.isna(max_order) or max_order <= 0:
        print(f"No valid data for {Drug}")
        return

    fig, axes = plt.subplots(len(concentrations), max_order, figsize=(max_order * 6, len(concentrations) * 6))

    x_min, x_max = bin_range
    y_min, y_max = bin_range
    xlim = ylim = (x_min, x_max)

    for i, conc in enumerate(concentrations):
        for order in range(1, max_order + 1):
            ax = axes[i, order - 1] if len(concentrations) > 1 else axes[order - 1]
            conc_slice = slice & (df['Concentration'] == conc)
            fitness_measured = df.loc[conc_slice, 'Fitness']
            fitness_predicted = df.loc[conc_slice, f'Fitness_predicted for order {order}']
            
            # Create 2D histogram
            bins = np.linspace(bin_range[0], bin_range[1], gridsize)
            hist, xedges, yedges, img = ax.hist2d(fitness_measured, fitness_predicted, bins=[bins, bins], cmap='Reds', norm=LogNorm())

            # Format plot
            if i == 0:
                ax.set_title(f'Order {order}', fontsize=18, fontweight='bold')
            if i == len(concentrations) - 1:
                ax.set_xlabel('Actual Fitness', fontsize=18, labelpad=18, fontweight='bold')
                ax.set_xticks(range(int(x_min), ceil(x_max)))
                ax.set_xticklabels(ax.get_xticks(), rotation=45, fontsize=18)
            else:
                ax.set_xticklabels([])
            if order == 1:
                ax.set_ylabel(f'Fitness$_{{predicted}}$ in [{Drug}] = {conc} µg/mL', fontsize=18, labelpad=18, fontweight='bold')
                ax.set_yticklabels([f'{int(tick)}' for tick in ax.get_yticks()], rotation=0, fontsize=18)
            else:
                ax.set_yticklabels([])
            ax.plot([x_min, x_max], [y_min, y_max], 'k--', lw=1, alpha=0.5)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_aspect('equal')
            ax.grid(True)

    plt.tight_layout(rect=[0, 0, 1, 1], h_pad=0.3, w_pad=0.3)
    plt.show()

def plot_epistasis_vs_num_mutations(Epistasis_Combined, drug_name):
    """
    Plot epistasis values against mutation count for each drug concentration.
    Uses jitter to better visualize overlapping points.
    """
    doses = sorted(Epistasis_Combined['Concentration'].unique())
    fig, axes = plt.subplots(len(doses), 1, figsize=(10, 20), sharex=True)

    for i, dose in enumerate(doses):
        ax = axes[i]
        plot_data = Epistasis_Combined[Epistasis_Combined['Concentration'] == dose]
        
        # Add jitter to better visualize overlapping points
        order = plot_data['Epistatic Order']
        jitter = np.where(order > 0, np.random.normal(0, 0.08, size=len(order)), 0)
        ax.scatter(plot_data['Epistatic Order'] + jitter, plot_data['Ensemble Averaging'], s=6, alpha=0.2, label='Ensemble Averaging')

        ax.text(0.01, 0.95, f'[{drug_name}] = {dose} µg/ml', transform=ax.transAxes, ha='left', va='top')
        ax.set_ylabel('Epistasis')
        ax.set_ylim(-25, 25)
        ax.set_xlim(-0.1, 13.3)
        ax.set_xticks(range(14))

        if i == 0:
            ax.set_title(f'Epistasis vs Number of Mutations for Different {drug_name} Doses')

    axes[-1].set_xlabel('Number of Mutations')
    plt.tight_layout()
    plt.show()

def main():
    """
    Main function to calculate epistasis for AMP and AZT drugs across different concentrations.
    Loads precomputed data, calculates epistasis using multiple methods, and saves results.
    """
    # Load drug response data
    amp_stats_df = pl.read_parquet(PROJECT_PATH / "data/processed/amp_auc_long_df.parquet")
    azt_stats_df = pl.read_parquet(PROJECT_PATH / "data/processed/azt_auc_long_df.parquet")

    # Generate all possible mutant combinations
    intended_mutants_size = np.prod(get_mutation_options(INTENDED))
    intended_mutants_all = [''.join(mut) for mut in product(*[INTENDED[pos] for pos in sorted(INTENDED.keys())])]
    all_mutants = [''.join(mut) for mut in product(*[ALL_MUTATIONS[pos] for pos in sorted(ALL_MUTATIONS.keys())])]

    # Create DataFrame with all genotypes and their properties
    F = pd.DataFrame({'Genotype': intended_mutants_all, 'Genotype Masked': all_mutants})
    F['Epistatic Term'] = F['Genotype'].apply(lambda x: generate_binary_string(x, INTENDED))
    F['N Mutations'] = F['Epistatic Term'].apply(lambda x: sum(1 for char in x if char != '0'))
    F = F.sort_values(['Epistatic Term']).reset_index(drop=True)
    F['Epistatic Order'] = F['Epistatic Term'].apply(lambda x: sum(1 for char in x if char != '0'))

    # Load precomputed matrices for epistasis calculations
    V = torch.load(CACHE_DIR / 'V.pt', weights_only=True, map_location=device)
    X = torch.load(CACHE_DIR / 'X.pt', weights_only=True, map_location=device)
    H = torch.load(CACHE_DIR / 'H.pt', weights_only=True, map_location=device)
    W_ea = torch.load(CACHE_DIR / 'W_ea.pt', weights_only=True, map_location=device)
    W_bioch = torch.load(CACHE_DIR / 'W_bioch.pt', weights_only=True, map_location=device)

    # Function to add fitness data for each drug and concentration
    def add_fitness_data(F, drug_stats_df, drug_name):
        """Join fitness data to genotype dataframe for each concentration."""
        Fs = []
        for conct in np.sort(drug_stats_df.select('concentration').unique().to_numpy().flatten()):
            F_pl = pl.DataFrame(F)
            const_conc = drug_stats_df.filter(pl.col('concentration') == conct)
            F_pl = F_pl.join(const_conc[['mutant_profile', 'median', 'std']], left_on='Genotype Masked', right_on='mutant_profile', how='left')
            F_pl = F_pl.rename({"median": "Fitness", "std": "Error"})
            F_pl = F_pl.with_columns([
                pl.lit(drug_name).alias('Drug').first(),
                pl.lit(conct).alias('concentration').first()
            ])
            F_pl = F_pl.sort('Epistatic Term')
            Fs.append(F_pl)
        return Fs

    # Calculate fitness data for each drug
    F_AZT = add_fitness_data(F, azt_stats_df, 'AZT')
    F_AMP = add_fitness_data(F, amp_stats_df, 'AMP')

    # Calculate epistasis for all drugs and concentrations
    Epistasis_Combined = []
    for F in F_AZT + F_AMP:
        model_fitness = F['Fitness'].to_numpy().copy()
        model_error = F['Error'].to_numpy().copy()

        # Handle NaN/Inf values by replacing with minimum values
        min_fitness = np.nanmin(model_fitness[np.isfinite(model_fitness)])
        min_error = np.nanmin(model_error[np.isfinite(model_error)])

        nan_inf_indices = ~(np.isfinite(model_fitness) & np.isfinite(model_error))
        model_fitness[nan_inf_indices] = min_fitness
        model_error[nan_inf_indices] = min_error

        # Calculate epistasis using different methods
        epistasis, epistasis_error = calculate_epistasis_biochemical(W_bioch, model_fitness, model_error)
        epistasis_ea, epistasis_ea_error = calculate_Ensemble_Averaged_epistasis(W_ea, model_fitness, model_error)

        Epistasis_Combined.append(
            pd.DataFrame({
            'Genotype': F['Genotype'],
            'Epistatic Term': F['Epistatic Term'],
            'Epistatic Order': F['Epistatic Order'],
            'Fitness': model_fitness,
            'Error': model_error,
            'Biochemical Definition': epistasis,
            'Error_Bioch': epistasis_error,
            'Ensemble Averaging': epistasis_ea,
            'Error_EA': epistasis_ea_error,
            'Drug': F['Drug'][0],
            'Concentration': F['concentration'][0]
            })
        )
    Epistasis_Combined = pd.concat(Epistasis_Combined).reset_index(drop=True)

    # Free up memory
    del W_bioch, W_ea
    clear_memory(device)
    
    # Run linear regression for every drug and concentration
    for drug in Epistasis_Combined.Drug.unique():
        concs = np.array(sorted(Epistasis_Combined.query(f'Drug == "{drug}"').Concentration.unique()))
        for conc in concs:
            print(f'Calculating Epistasis for {drug} at {conc} µg/ml:')
            E_slice = (Epistasis_Combined['Drug'] == drug) & (Epistasis_Combined['Concentration'] == conc)
            Epistasis_LG_partial_calculate(E_slice, V, X, H, Epistasis_Combined)

    # Convert float16 columns to float32 to avoid ArrowNotImplementedError
    for col in Epistasis_Combined.select_dtypes(include=['float16']).columns:
        Epistasis_Combined[col] = Epistasis_Combined[col].astype('float32')

    # Save results to the project directory
    Epistasis_Combined.to_parquet(PROJECT_PATH / "data/processed/Epistasis_Combined.parquet")

    # Optional: Generate plots
    # plot_fitness_vs_predicted(Epistasis_Combined, Drug='AMP')
    # plot_fitness_vs_predicted(Epistasis_Combined, Drug='AZT')
    # plot_epistasis_vs_num_mutations(Epistasis_Combined, drug_name='AMP')
    # plot_epistasis_vs_num_mutations(Epistasis_Combined, drug_name='AZT')

if __name__ == "__main__":
    main()

# TEM1 Combinatorial Mutagenesis Data

This repository contains data from experiments involving TEM1 combinatorial mutagenesis under ampicillin and aztreonam selection.

## Data Files

### Metadata

#### Experimental Metadata
`data/raw/metadata.csv`: Contains information about each experimental sample with the following columns:
- `Sample Name`: Unique sequencing sample identifier for each sample (format: [Drug]_[Concentration]_[Replicate]_[Timepoint])
- `Drug`: Antibiotic used (Ampicillin, Aztreonam, or BLT1)
- `Concentration`: Drug concentration in μg/ml
- `Timepoint`: Time at which the sample was taken (0h, 3h, 6h, 9h, 12h, or 8h for BLT1)
- `Replicate`: Replicate number (1, 2, or 3)
- `OD600`: Optical density measurement at 600nm, indicating bacterial growth

### Raw Data Tables

#### Ampicillin Read Counts
`data/raw/Ampicillin_read_counts_per_genotype.csv`: Contains read counts for each genotype under ampicillin selection. Columns include genotype identifiers and read counts for various ampicillin concentrations. The sample names in this file match the `Sample Name` column in the metadata sheet.

#### Aztreonam Read Counts
`data/raw/Aztreonam_read_counts_per_genotype.csv`: Contains read counts for each genotype under aztreonam selection, with columns structured similarly to the ampicillin read counts file. The sample names correspond to those in the metadata sheet's `Sample Name` column.

#### Ampicillin AUC Measurements
`data/raw/Ampicillin_auc_per_genotype.csv`: Contains Area Under the Curve (AUC) measurements, representing antibiotic resistance levels, for each genotype under ampicillin selection. Sample names are consistent with the `Sample Name` column in the metadata sheet.

#### Aztreonam AUC Measurements
`data/raw/Aztreonam_auc_per_genotype.csv`: Contains AUC measurements for each genotype under aztreonam selection, structured similarly to the ampicillin AUC file. Sample names are consistent with the `Sample Name` column in the metadata sheet.

The processed data is available in both wide and long formats:

### Wide Format Files
- `data/processed/amp_auc_wide_df.parquet`: Ampicillin AUC data in wide format
- `data/processed/azt_auc_wide_df.parquet`: Aztreonam AUC data in wide format

These files contain columns:
- `mut_profile_masked`: Mutant profile string (13 characters representing mutations at specific positions)
- Concentration-replicate columns (e.g., "Ampicillin 0.0 1"): AUC values for each concentration and replicate

### Long Format Files
- `data/processed/amp_auc_long_df.parquet`: Ampicillin AUC data in long format with statistics
- `data/processed/azt_auc_long_df.parquet`: Aztreonam AUC data in long format with statistics

These files contain columns:
- `mutant_profile`: Mutant profile string
- `concentration`: Drug concentration value
- `replicate1`, `replicate2`, `replicate3`: AUC values for each replicate
- `mean`: Mean AUC across replicates
- `median`: Median AUC across replicates
- `std`: Standard deviation of AUC across replicates
- `cv_std`: Coefficient of variation (std/mean)

### Epistasis Analysis File
- `data/processed/Epistasis_Combined.parquet`: Epistasis analysis data for both drugs at all concentrations

This file contains columns:
- `Genotype`: Mutant profile string
- `Epistatic Term`: Binary representation of mutations
- `Epistatic Order`: Number of mutations in the genotype
- `Fitness`: Median AUC value (from long format files)
- `Error`: Standard deviation of AUC
- `Biochemical Definition`: Epistasis values calculated using biochemical definition
- `Error_Bioch`: Error of biochemical epistasis values
- `Ensemble Averaging`: Epistasis values calculated using ensemble averaging method
- `Error_EA`: Error of ensemble averaging epistasis values
- `Drug`: Drug name (AMP or AZT)
- `Concentration`: Drug concentration value
- `Epistasis_LG for order X`: Epistasis values calculated using linear regression for order X
- `Error_LG for order X`: Error of linear regression epistasis values for order X
- `Fitness_predicted for order X`: Predicted fitness using linear regression model for order X

## Mutation Key

The mutant profile string represents mutations at positions:
19, 37, 67, 102, 162, 180, 235, 236, 237, 241, 261, 271, 272 (Ambler numbering: 21, 39, 69, 104, 164, 182, 237, 238, 240, 244, 265, 275, 276)

Where:
- `.` indicates wild-type
- Letter indicates mutation to that amino acid
- `X` indicates a non-functional (dead) mutant

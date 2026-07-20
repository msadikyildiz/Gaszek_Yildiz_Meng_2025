import pandas as pd
import json
import math
import warnings
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
from functools import partial
from scipy.optimize import least_squares
from scipy.optimize import fsolve

# Suppress specific runtime warnings that are expected in normal operation
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=FutureWarning)


def log_logistic_4pl_(params, x):
    bottom, top, ec50, hill_slope = params
    return bottom + (top - bottom) / (1 + (x / ec50)**hill_slope)

def log_logistic_4pl_solve(y_target, params, x_guess):
    # Equation to solve for x, where four_pl(x) = y_target
    return fsolve(lambda x: log_logistic_4pl_(params, x) - y_target, x0=x_guess)

def log_logistic_3pl_(params, x):
    top, ec50, hill_slope = params
    return top / (1 + (x / ec50)**hill_slope)

def log_logistic_3pl_solve(y_target, params, x_guess):
    # Equation to solve for x, where three_pl(x) = y_target
    return fsolve(lambda x: log_logistic_3pl_(params, x) - y_target, x0=x_guess)

def log_logistic_2pl_(params, x):
    ec50, hill_slope = params
    return 1 / (1 + (x / ec50)**hill_slope)

def log_logistic_2pl_solve(y_target, params, x_guess):
    # Equation to solve for x, where four_pl(x) = y_target
    return fsolve(lambda x: log_logistic_2pl_(params, x) - y_target, x0=x_guess)


class ufpf:
    """UFPF is class definition for a file format that includes:
          - A path to source excel file
          - Parsed dataframes for plate information and OD values
          - Methods to read and process the data
    """        
    def __init__(self, path, timepoint_in_sheetnames=False):
        self.path = path
        self.number_of_plates = None
        self.plate_info = None
        self.OD = None
        self.read_data(timepoint_in_sheetnames)
    
    def read_data(self, timepoint_in_sheetnames=False):
        # Get the number of sheets in the excel file
        self.number_of_plates = len(pd.ExcelFile(self.path).sheet_names)
        # For loop over each sheet read them into a dataframe
        for sheet_id, sheet_name in enumerate(pd.ExcelFile(self.path).sheet_names):
            # Read as-is
            raw = pd.read_excel(self.path, sheet_name=sheet_name, header=None)

            # Locate 96-well plate formatted entries in the excel sheet
            plate_entries = []
            for ix, row in raw.iterrows():
                if np.all(np.array(row.values[1:13]) == np.arange(1,13)):
                        plate_entries.append(ix)

            # 1. Collect OD values
            od_wide = raw.iloc[plate_entries[0]+1:plate_entries[0]+9, 1:13]
            od_wide['Row'] = np.arange(1,9)
            od_long = od_wide.melt(id_vars=['Row'], var_name='Column', value_name='OD')
            od_long['Plate_ID'] = sheet_id
            well_rows = {1:'A', 2:'B', 3:'C', 4:'D', 5:'E', 6:'F', 7:'G', 8:'H'}
            od_long['Well'] = od_long['Row'].map(well_rows) + od_long['Column'].astype(str)

            # 2. Collect plate info
            infos = dict()
            for block_name, block_id in zip(['Antibiotic', 'Dose', 'Strain'], range(1,4)):
                info_wide = raw.iloc[plate_entries[block_id]+1:plate_entries[block_id]+9, 1:13]
                info_wide['Row'] = np.arange(1,9)
                info_long = info_wide.melt(id_vars=['Row'], var_name='Column', value_name=block_name)
                info_long['Plate_ID'] = sheet_id
                well_rows = {1:'A', 2:'B', 3:'C', 4:'D', 5:'E', 6:'F', 7:'G', 8:'H'}
                info_long['Well'] = info_long['Row'].map(well_rows) + info_long['Column'].astype(str)
                infos[block_name] = info_long
  
            # 3. Merge all plate info into single table
            plate_info = pd.merge(infos['Antibiotic'],infos['Dose'], on=['Plate_ID', 'Well', 'Row', 'Column'])
            plate_info = pd.merge(plate_info, infos['Strain'], on=['Plate_ID', 'Well', 'Row', 'Column'])
            
            # If there is timepoint data provided in sheetnames, append this to strain names
            if timepoint_in_sheetnames:
                plate_info['Strain'] = plate_info['Strain'].astype(str) + '_' + sheet_name.split('_')[-1]
                plate_info.loc[plate_info['Strain'].str.contains('nan'), 'Strain'] = 'nan'

            # 4. Stack all plates into one table
            if self.plate_info is None:
                self.plate_info = plate_info
            else:
                self.plate_info = pd.concat([self.plate_info, plate_info])
            
            if self.OD is None:
                self.OD = od_long
            else:
                self.OD = pd.concat([self.OD, od_long])

class static:
    
    def load_plate_info(file_path):
        '''
        This function assumes generalized plate information format specified as follows:

        - Input is an multi-sheet xlsx file
        - Each sheet corresponds to a single plate map (strain, replicate, antibiotic dose etc.)
        - Each sheet name is like "Plate {number} - Layout
        - Each sheet contain following columns: Well, Strain, Replicate, Antibiotic, Dose, Unit (optional)
        - Blank wells should be annotated "Media Only"
        '''
        # Get all sheet names
        sheet_names = pd.ExcelFile(file_path).sheet_names
        # If there are multiple sheets in the excel combine
        if len(sheet_names) > 1:
            dfs = []
            # Iterate through each sheet
            for sheet_name in sheet_names:
                # Extract Plate ID from sheet name
                plate_id = sheet_name.split('-')[0].strip().split(' ')[1].strip()
                # Load the data from the sheet
                df = pd.read_excel(file_path, sheet_name=sheet_name)
                # Add Plate_ID column
                df['Plate_ID'] = int(plate_id)
                dfs.append(df)
            # Concatenate all DataFrames
            df_plate_info = pd.concat(dfs, ignore_index=True)
        # If there is only one sheet
        else:
            df_plate_info = pd.read_excel(file_path)
            if not (('Plate_ID' in df_plate_info.columns) or ('Plate' in df_plate_info.columns)):
                # Add Plate_ID column
                df_plate_info['Plate_ID'] = 1
            # Column name if only plate be graceful
            if 'Plate' in df_plate_info.columns:
                df_plate_info.rename(columns={'Plate':'Plate_ID'}, inplace=True)
                # Multiple experiments combine them
                # Group df_plate_info by experiment and plate_id
                grouped_df = df_plate_info.groupby(['Experiment', 'Plate_ID'])
                # Assign new incrementing plate_ids for each plate in each experiment
                df_plate_info['combined_id'] = grouped_df.ngroup() + 1
                # Merge the combined_ids back into the original dataframe
                df_plate_info = df_plate_info.merge(df_plate_info[['Experiment', 'Plate_ID', 'combined_id']],
                                                    on=['Experiment', 'Plate_ID'], suffixes=['', '_combined'])
                
                df_plate_info['Plate_ID'] = df_plate_info['combined_id']
                # Drop per experiment plate IDs
                df_plate_info.drop(columns=['Experiment','combined_id','combined_id_combined'], inplace=True)

        if 'Drug' in df_plate_info.columns:
            df_plate_info.rename(columns={'Drug':'Antibiotic'}, inplace=True)

        return df_plate_info
    
    def get_unique_plate_ids(df_plate_info):
        """ Extract and return unique Plate_IDs and their numerical parts from df_plate_info """
        # Extract unique Plate_IDs
        unique_plates = df_plate_info['Plate_ID'].unique()
        if not pd.api.types.is_integer_dtype(df_plate_info['Plate_ID']):
            # Extract numbers from Plate_ID strings
            plate_numbers = [int(str(plate).split(' ')[1]) for plate in unique_plates]
        else:
            plate_numbers = list(unique_plates)
        return unique_plates, sorted(plate_numbers)

    def read_biotek_OD_data(path):
        ## Read data from input path
        df_OD_raw = pd.read_excel(path)
        # Identify how many plates exist in the BioTek output file
        found_plate_index = df_OD_raw[df_OD_raw.iloc[:,0] == 'Results'].index
        ## Iterate and collect them into long-format table
        results = []
        for plate_id, i in enumerate(found_plate_index):
            start_row = i+4
            end_row = i+12
            # Extract the block
            data_block = df_OD_raw.iloc[start_row:end_row + 1, 2:14]

            # Map each entry in the data block to the correct well in df
            for idx, row in enumerate(data_block.values):  # iterate through rows in the block
                for jdx, value in enumerate(row):  # iterate through columns in the row
                    row_letter = chr(65 + idx)  # ASCII 'A'
                    well_position = f"{row_letter}{jdx + 1}"

                    # Append a new row to results list as a DataFrame
                    results.append({
                        'Plate_ID': plate_id+1,
                        'Well': well_position,
                        'OD': value,
                        'Row': idx + 1,
                        'Column': jdx + 1
                    })

        # Concatenate all dataframes in the results list
        df_results = pd.DataFrame(results)
        return df_results
    
    def convert_OD_plate_to_long(df_OD_raw,  df_plate_info):
        """ Extract data blocks for each plate and map to the wells along with corresponding Plate_IDs 
        Blocks of data are organized in stacked 8 by 12 arrays as input, they are converted to long format. 
        This is a necessary step.
        """
        
        unique_plates, plate_numbers = static.get_unique_plate_ids(df_plate_info)
        # Initialize an empty list to collect DataFrame slices
        results = []

        # Iterate over each plate number to extract and map the data block
        for n in plate_numbers:
            # Calculate 0-based row indices for the 1-based formula
            start_row = 10 * n - 9 - 1
            end_row = 10 * n - 2 - 1

            # Extract the block
            data_block = df_OD_raw.iloc[start_row:end_row + 1, 1:13]

            # Map each entry in the data block to the correct well in df
            for idx, row in enumerate(data_block.values):  # iterate through rows in the block
                for jdx, value in enumerate(row):  # iterate through columns in the row
                    row_letter = chr(65 + idx)  # ASCII 'A'
                    well_position = f"{row_letter}{jdx + 1}"

                    # Append a new row to results list as a DataFrame
                    results.append(pd.DataFrame({
                        'Plate_ID': n,
                        'Well': [well_position],
                        'OD': [value],
                        'Row': [idx + 1],
                        'Column': [jdx + 1]
                    }))

        # Concatenate all dataframes in the results list
        df_results = pd.concat(results, ignore_index=True)
        return df_results

    def calc_median_background_all_plates(df, df_plate_info, specific_wells=None, excluded_plates=None, plot=False):
        """
        Calculate the common median background OD value for wells labeled as 'Media Only' in specific wells,
        and plot histogram of these OD values.
        - Optionally excluding plates and specific wells is implemented.
        """
       
        # Filter out specific conditions
        if excluded_plates is not None:
            exclude_condition = (~df_plate_info['Plate_ID'].isin(excluded_plates))
        else:
            exclude_condition = True
        if specific_wells is not None:
            specific_wells_condition = (df_plate_info['Well'].isin(specific_wells))
        else:
            specific_wells_condition = True

        condition = (
            ((df_plate_info['Strain'].str.lower() == 'media only') |
             (df_plate_info['Antibiotic'].str.lower() == 'media only')) &
             exclude_condition &
            specific_wells_condition
        )
        media_only = df_plate_info[condition]

        # Select OD values corresponding to 'Media Only' wells meeting the specific conditions
        media_only_values = media_only.merge(df, on=['Plate_ID', 'Well'], how='left')['OD']
        
        # Calculate and return the median of these OD values
        median_value = media_only_values.median()
        # print("Common Background Median OD:", median_value)
        if plot:
            # Plotting the histogram of OD values
            plt.figure(figsize=(12, 1.5),dpi=90)
            plt.hist(media_only_values.dropna(), bins=np.arange(0,1,0.01), color='black')
            plt.title('OD Values for Media Only Wells')
            plt.xlabel('OD Values')
            plt.ylabel('Frequency')
            plt.grid(True)
            plt.show()

        return median_value

    def apply_background_correction(df, df_plate_info, plate_id, background_common):
        """Apply background correction based on 'Media Only' wells for a specific plate.
        If no 'Media Only' wells are present, use a predefined common background value."""
        
        # Filter the current plate information
        df_plate_current = df_plate_info[df_plate_info['Plate_ID'] == plate_id]
        
        # Filter 'Media Only' wells for the current plate
        media_only_wells = df_plate_current[(df_plate_current['Strain'].str.lower() == 'media only')  |
                                            (df_plate_current['Antibiotic'].str.lower() == 'media only')]['Well'].tolist()

        if media_only_wells:
            # Get the data for 'Media Only' wells from df and calculate the median of the first 9 time points
            media_values = df[(df['Plate_ID'] == plate_id) & (df['Well'].isin(media_only_wells))]['OD'].iloc[1:10]
            median_value = media_values.median()
        else:
            # Use predefined common background value if no 'Media Only' wells are present
            median_value = background_common

        # Create a new DataFrame for background corrected values
        df_bc = df.copy()

        # Filter all wells related to the current plate
        all_wells = df_plate_current['Well'].tolist()

        # Apply correction to each well in the current plate
        for well in all_wells:
            well_mask = (df_bc['Plate_ID'] == plate_id) & (df_bc['Well'] == well)
            corrected_values = df_bc.loc[well_mask, 'OD'] - median_value
            # Set any negative values to 10^-6
            corrected_values[corrected_values <= 0] = 10**-6
            df_bc.loc[well_mask, 'OD_final'] = corrected_values

        return df_bc[df_bc['Plate_ID'] == plate_id]

    def plot_single_plate_media_only_wells(df, df_plate_info, plate_id, column_name):
        """ Plot histogram of 'Media Only' values for a specific plate """
        if plate_id is None:
            raise ValueError("No Plate_ID provided")

        # Filter 'Media Only' wells for a specific Plate_ID
        media_only_wells = df_plate_info[((df_plate_info['Strain'].str.lower() == 'media only') |
                                          (df_plate_info['Antibiotic'].str.lower() == 'media only'))
                                         & (df_plate_info['Plate_ID'] == plate_id)]['Well'].tolist()

        if not media_only_wells:
            raise ValueError("No media only wells found for the specified Plate_ID.")

        # Filter df for the selected Plate_ID and wells in the 'Media Only' list
        media_values = df[(df['Plate_ID'] == plate_id) & (df['Well'].isin(media_only_wells))]

        if media_values.empty:
            raise ValueError("No data found for media only wells in the DataFrame.")

        # Extract OD values for plotting
        all_media_values = media_values[column_name]

        # Calculate the median
        median_value = all_media_values.median()

        # Plotting the histogram of 'Media Only' values
        plt.figure(figsize=(3, 3))
        plt.hist(all_media_values, color='skyblue', alpha=0.7, label='Media Values')
        plt.axvline(median_value, color='red', linestyle='dashed', linewidth=2, label=f'Median: {median_value:.3f}')
        plt.title(f'Media Only Values for {plate_id}', fontsize =10)
        plt.xlabel('OD Value')
        plt.ylabel('Frequency')
        plt.legend()
        plt.grid(True)
        plt.show()

        return median_value

    def plot_all_plates_media_only_wells(df, df_plate_info, column_name):
        """ Generate subplots for histograms of 'Media Only' values for all plates """
        # Determine the layout of the subplots grid
        unique_plate_ids = df_plate_info.Plate_ID.unique()
        num_plates = len(unique_plate_ids)
        cols = 4  # Define the number of columns in the subplot grid
        rows = (num_plates + cols - 1) // cols  # Calculate required rows to accommodate all plates

        # Create a figure with subplots
        fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(cols * 5, rows * 5))
        fig.subplots_adjust(hspace=0.4, wspace=0.4)  # Adjust spacing between plots

        # Flatten axes array if multiple rows
        if num_plates > 1:
            axes = axes.flatten()
        
        # Loop through each plate and plot its histogram
        for idx, plate_id in enumerate(unique_plate_ids):
            ax = axes[idx] if num_plates > 1 else axes  # Select subplot axis
            try:
                media_only_wells = df_plate_info[((df_plate_info['Strain'].str.lower() == 'media only') |
                                                  (df_plate_info['Antibiotic'].str.lower() == 'media only'))
                                                 & (df_plate_info['Plate_ID'] == plate_id)]['Well'].tolist()
                media_values = df[(df['Plate_ID'] == plate_id) & (df['Well'].isin(media_only_wells))][column_name]
                
                if media_values.empty:
                    ax.text(0.5, 0.5, 'No data found', transform=ax.transAxes, ha='center', va='center')
                    continue

                median_value = media_values.median()
                ax.hist(media_values, bins=10, color='skyblue', alpha=0.7, label='Media Values')
                ax.axvline(median_value, color='red', linestyle='dashed', linewidth=2, label=f'Median: {median_value:.4f}')
                ax.set_title(f'{plate_id}', fontsize=10)
                ax.set_xlabel('OD Value')
                ax.set_ylabel('Frequency')
                ax.legend()
                ax.grid(True)
            
            except Exception as e:
                ax.text(0.5, 0.5, str(e), transform=ax.transAxes, ha='center', va='center')

        # Hide any unused axes if there are empty subplots
        for ax in axes[num_plates:]:
            ax.axis('off')

        plt.show()

    def plot_dose_response_curve_errorbar(df_analysis, strain, antibiotic, strain_colors, ax):
        
        # Set y-axis limits and ticks for OD_final
        ax.set_ylim(-0.03, 1.1)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
        
        # Ensure every unique strain has a color, assign random colors if not predefined
        unique_strains = df_analysis['Strain'].unique()
        for unique_strain in unique_strains:
            if unique_strain not in strain_colors:
                strain_colors[unique_strain] = 'gray'  # default color
                    # Use predefined colors or assign a random color if the strain is not listed
        
        # Filter data for the specific strain and antibiotic
        filtered_data = df_analysis[
            (df_analysis['Strain'].str.lower() == strain.lower()) &
            (df_analysis['Antibiotic'].str.lower() == antibiotic.lower())
        ]
        
        # Group by Dose and calculate mean and standard deviation
        stats = filtered_data.groupby('Dose')['OD_final'].agg(['mean', 'std']).reset_index()

        # Replace 0 with half of the smallest non-zero dose for log-scale plotting
        non_zero_doses = stats['Dose'] > 0
        min_nonzero_dose = stats.loc[non_zero_doses, 'Dose'].min() if non_zero_doses.any() else 0.1
        stats['Dose'].replace(0, min_nonzero_dose / 2, inplace=True)

        # Plotting mean values with error bars for standard deviation
        ax.errorbar(stats['Dose'], stats['mean'], yerr=stats['std'], fmt='-o', color=strain_colors[strain],
                    ecolor=strain_colors[strain], elinewidth=2, capsize=3, markerfacecolor='none')

        ax.set_xlabel(f'{antibiotic} (μg/mL)', fontsize=10)
        ax.set_ylabel('OD_final', fontsize=10)
        ax.set_title(f'Response for {strain}', fontsize=10)
        ax.text(0.02, 0.96, strain, transform=ax.transAxes, fontsize=10, verticalalignment='top')
        ax.grid(True, linestyle=':')
        ax.set_xscale('log')

        # Set x-ticks to include the actual zero value if necessary
        all_doses = np.sort(np.unique(stats['Dose']))
        ax.set_xticks(all_doses)
        ax.set_xticklabels([f'{dose:.2f}' if dose != min_nonzero_dose / 2 else '0' for dose in all_doses], fontsize=8, fontweight='normal', rotation=60)

    def plot_od_final_for_selected_antibiotic(df_analysis, plot_function, selected_antibiotic, strain_colors=None):
        # Filter the DataFrame for the selected antibiotic
        df_filtered = df_analysis[df_analysis['Antibiotic'] == selected_antibiotic]
        
        # Get unique strains for the selected antibiotic
        unique_strains = df_filtered['Strain'].unique()
        unique_strains = list(set(unique_strains) - set(['Media Only','Cells Only']))
        
        # Calculate the number of rows and columns for the subplots
        num_columns = 4
        num_rows = math.ceil(len(unique_strains) / float(num_columns))

        # Create a figure with calculated rows and columns
        fig, axs = plt.subplots(num_rows, num_columns, figsize=(4 * num_columns, 3 * num_rows), dpi=150)
        
        # Iterate over each strain and plot in the appropriate subplot
        for i, strain in enumerate(unique_strains):
            plot_function(df_filtered, strain, selected_antibiotic, strain_colors, axs.flat[i])

        # Handle any unused subplots by hiding them
        total_plots = num_rows * num_columns
        for i in range(len(unique_strains), total_plots):
            axs.flat[i].axis('off')  # Turn off unused subplots
        
        # Set a super title for the figure, adjust using `y` for vertical placement
        fig.suptitle(f'Antibiotic: {selected_antibiotic}', fontsize=14, y=0.97)
        
        # Adjust layout to prevent overlap and give space for super title
        plt.tight_layout()

        # Show the plot
        plt.show()

    def find_dose_for_effect(effect_level, y_fit, x_fit):
        """ Helper function to find the dose corresponding to a given effect level. """
        try:
            # Check if the effect level ever reaches the target in the fitted range
            if effect_level < np.min(y_fit) or effect_level > np.max(y_fit):
                return f"> {x_fit.max():.2f}"
            else:
                #return np.interp(effect_level, y_fit[::-1], x_fit[::-1])  # Ensure it's sorted correctly for interp
                return 
        except Exception as e:
            return f"Error: {e}"  # Return a string if the interpolation fails

    def prep_valid_combinations(df_analysis, multiplex=['Antibiotic', 'Replicate', 'Strain'], ic50_threshold=0.5, mic_threshold=0.05):

        # Filter out 'media only' strains and select relevant columns including 'Plate_ID'
        valid_combinations = df_analysis[
            (df_analysis['Strain'].str.lower() != 'media only') &
            (df_analysis['Strain'].str.lower() != 'cells only') &
            (df_analysis['Antibiotic'].str.lower() != 'media only') &
            (df_analysis['Antibiotic'].str.lower() != 'cells only')][['Antibiotic', 'Replicate', 'Strain', 'Plate_ID']]

        # Drop duplicates based on 'Antibiotic', 'Replicate', 'Strain' while keeping the first 'Plate_ID'
        valid_combinations = valid_combinations.drop_duplicates(subset=multiplex)

        # Eliminate rows where any of 'Antibiotic', 'Replicate', 'Strain', 'Plate_ID' does not exist
        valid_combinations.dropna(subset=['Antibiotic', 'Replicate', 'Strain', 'Plate_ID'], inplace=True)

        # Reset the index after dropping rows
        valid_combinations.reset_index(drop=True, inplace=True)

        return valid_combinations[multiplex]

    def dose_response_data_selector(df_analysis, strain, antibiotic, replicate=""):
        # Extracting dose-response data for the specified conditions
        dose_response_data = df_analysis[((df_analysis['Antibiotic'] == antibiotic) |
                                          (df_analysis['Antibiotic'].str.lower() == 'cells only')) &
                                         (df_analysis['Strain'] == strain)]
        if replicate != "":
            dose_response_data = dose_response_data[dose_response_data['Replicate'] == replicate]
        # Use psuedo-log transformation to handle zeros
        nonzero_min = dose_response_data['Dose'][dose_response_data['Dose'] > 0].min()
        log_min_dose = (nonzero_min)/10
        dose_response_data.loc[:, 'Dose'] = dose_response_data['Dose'].fillna(log_min_dose)
        dose_response_data.loc[dose_response_data['Dose']<nonzero_min, 'Dose'] = nonzero_min
        return dose_response_data
    
    # Define the residuals function
    def residuals(params, x, y):
        return y - log_logistic_2pl_(params, x)
    
    # Define the residuals function for 3PL
    def residuals_3pl(params, x, y):
        return y - log_logistic_3pl_(params, x)

    # Function to fit model and return EC50 estimate
    def fit_model(x_data, y_data, p0):
        # Robust fitting with better parameter bounds and scaling
        # Set reasonable bounds for EC50 and hill slope
        bounds = ([x_data.min() * 0.001, 0.1], [x_data.max() * 1000, 10])
        
        try:
            # Try with scaled parameters first
            return least_squares(static.residuals, p0, args=(x_data, y_data), 
                               bounds=bounds, method='trf', loss='huber', 
                               ftol=1e-8, xtol=1e-8, max_nfev=1000)
        except:
            # Fallback to simpler bounds
            try:
                return least_squares(static.residuals, p0, args=(x_data, y_data), 
                                   method='lm', ftol=1e-6, xtol=1e-6)
            except:
                # Final fallback - return parameters without optimization
                from scipy.optimize import OptimizeResult
                return OptimizeResult(x=p0, success=False)
    
    # Function to fit 3PL model and return parameters
    def fit_model_3pl(x_data, y_data, p0):
        # Robust fitting for 3PL with parameter bounds
        # Bounds: [top, ec50, hill_slope]
        bounds = ([0.01, x_data.min() * 0.001, 0.1], 
                 [y_data.max() * 10, x_data.max() * 1000, 10])
        
        try:
            return least_squares(static.residuals_3pl, p0, args=(x_data, y_data), 
                               bounds=bounds, method='trf', loss='huber',
                               ftol=1e-8, xtol=1e-8, max_nfev=1000)
        except:
            # Fallback method
            try:
                return least_squares(static.residuals_3pl, p0, args=(x_data, y_data), 
                                   method='lm', ftol=1e-6, xtol=1e-6)
            except:
                # Final fallback
                from scipy.optimize import OptimizeResult
                return OptimizeResult(x=p0, success=False)
    
    def bootstrap_resample(x_data, y_data, p0, seed=42, mic_percentage=0.05):
        np.random.seed(seed)
        x_resampled = []
        y_resampled = []
        unique_doses = np.unique(x_data)
        # Vectorized resampling for each dose
        resampled_indices = np.hstack(
            [np.random.choice(np.where(x_data == dose)[0], size=np.sum(x_data == dose), replace=True) 
                for dose in unique_doses])
        
        # Resample both x and y based on the indices
        x_resampled = x_data[resampled_indices]
        y_resampled = y_data[resampled_indices]

        try:
            # Fit the model on the resampled data and store EC50 estimate
            optimizer = static.fit_model(x_resampled, y_resampled, p0)
            ic50, hill_coeff = optimizer.x

            # Calculate IC50 and MIC using the interpolated values
            mic = log_logistic_2pl_solve(mic_percentage, optimizer.x, ic50)[0]

            return ic50, hill_coeff, mic
        except Exception as e:
            return np.nan, np.nan
    
    def bootstrap_resample_3pl(x_data, y_data, p0_3pl, seed=42, mic_percentage=0.05):
        """Bootstrap resampling for 3PL curve fitting to estimate max_growth_scale robustly"""
        np.random.seed(seed)
        unique_doses = np.unique(x_data)
        # Vectorized resampling for each dose
        resampled_indices = np.hstack(
            [np.random.choice(np.where(x_data == dose)[0], size=np.sum(x_data == dose), replace=True) 
                for dose in unique_doses])
        
        # Resample both x and y based on the indices
        x_resampled = x_data[resampled_indices]
        y_resampled = y_data[resampled_indices]

        try:
            # Fit the 3PL model on the resampled data
            optimizer_3pl = static.fit_model_3pl(x_resampled, y_resampled, p0_3pl)
            top, ec50, hill_coeff = optimizer_3pl.x
            
            return top, ec50, hill_coeff
        except Exception as e:
            return np.nan, np.nan, np.nan
        
    
    def robust_phenotyper(df_analysis, antibiotic, strain, mic_percentage=0.05, n_bootstrap=100, threads=8):

        dose_response_data = static.dose_response_data_selector(df_analysis, strain, antibiotic)
        x_data = np.asarray(dose_response_data['Dose'].to_numpy(), dtype=float)
        y_data = np.asarray(dose_response_data['OD_final'].to_numpy(), dtype=float)
        
        # Data preprocessing to avoid numerical issues
        # Remove any NaN or infinite values
        valid_mask = np.isfinite(x_data) & np.isfinite(y_data) & (x_data > 0) & (y_data >= 0)
        x_data = x_data[valid_mask]
        y_data = y_data[valid_mask]
        
        # Check if we have enough data points
        if len(x_data) < 3:
            return {'Status': "FAIL - Insufficient data points",
                    'IC50': np.nan, 'MIC': np.nan, 'IC50_ci_lower': np.nan,
                    'IC50_ci_upper': np.nan, 'MIC_ci_lower': np.nan, 'MIC_ci_upper': np.nan,
                    'max_growth': np.nan, 'max_growth_ci_lower': np.nan, 'max_growth_ci_upper': np.nan,
                    'hill_coeff': np.nan, 'ic50_threshold': np.nan, 'mic_threshold': np.nan,
                    'x_fit': '[]', 'y_fit': '[]', 'ic50_bootstrap': '[]'}
        
        # Ensure x_data spans reasonable range
        if np.max(x_data) / np.min(x_data) < 2:
            return {'Status': "FAIL - Insufficient dose range",
                    'IC50': np.nan, 'MIC': np.nan, 'IC50_ci_lower': np.nan,
                    'IC50_ci_upper': np.nan, 'MIC_ci_lower': np.nan, 'MIC_ci_upper': np.nan,
                    'max_growth': np.nan, 'max_growth_ci_lower': np.nan, 'max_growth_ci_upper': np.nan,
                    'hill_coeff': np.nan, 'ic50_threshold': np.nan, 'mic_threshold': np.nan,
                    'x_fit': '[]', 'y_fit': '[]', 'ic50_bootstrap': '[]'}
        
        # First, fit a 3PL curve to estimate max_growth_scale empirically
        try:
            # Initial parameter estimates for 3PL: [top, ec50, hill_slope]
            top_initial = np.max(y_data)
            ec50_initial = np.exp((np.log(np.min(x_data))+np.log(np.max(x_data)))/2)
            hill_initial = 1
            p0_3pl = [top_initial, ec50_initial, hill_initial]
            
            # Fit 3PL model to estimate max_growth_scale
            optimizer_3pl = static.fit_model_3pl(x_data, y_data, p0_3pl)
            max_growth_scale = optimizer_3pl.x[0]  # The 'top' parameter is our max_growth_scale
            
            # Bootstrap estimation of max_growth_scale for confidence intervals
            max_growth_bootstrap = []
            with mp.Pool(processes=threads) as pool:
                bootstrap_3pl_results = pool.starmap(static.bootstrap_resample_3pl,
                                                   [(x_data, y_data, p0_3pl, seed) for seed in range(n_bootstrap)])
            
            # Collect max_growth (top parameter) from bootstrap results
            for top_, ec50_, hill_ in bootstrap_3pl_results:
                if not np.isnan(top_):
                    max_growth_bootstrap.append(top_)
            
            max_growth_bootstrap = np.array(max_growth_bootstrap)
            max_growth_ci_lower = np.percentile(max_growth_bootstrap, 2.5) if len(max_growth_bootstrap) > 0 else np.nan
            max_growth_ci_upper = np.percentile(max_growth_bootstrap, 97.5) if len(max_growth_bootstrap) > 0 else np.nan
            
        except Exception as e:
            # Fallback to original method if 3PL fitting fails
            max_growth_scale = dose_response_data.loc[dose_response_data['Dose'] == dose_response_data['Dose'].min(), 'OD_final'].median()
            max_growth_ci_lower = np.nan
            max_growth_ci_upper = np.nan
        
        # Normalize y_data using the empirically estimated max_growth_scale
        y_data = y_data / max_growth_scale
        
        try:
            # Better parameter initialization for 2PL fitting
            # Find approximate EC50 as dose where response is ~50% of max
            sorted_idx = np.argsort(x_data)
            x_sorted = x_data[sorted_idx]
            y_sorted = y_data[sorted_idx]
            
            # Find EC50 estimate using interpolation
            target_response = 0.5  # 50% inhibition
            if np.max(y_sorted) > target_response > np.min(y_sorted):
                ec50_initial = np.interp(target_response, y_sorted[::-1], x_sorted[::-1])
            else:
                ec50_initial = np.exp((np.log(np.min(x_data))+np.log(np.max(x_data)))/2)
            
            # Estimate hill slope from data steepness
            if len(x_data) > 3:
                # Calculate slope between 25% and 75% response points
                y25 = np.percentile(y_data, 25)
                y75 = np.percentile(y_data, 75)
                if y75 != y25:
                    hill_initial = abs(np.log(3) / np.log(np.max(x_data)/np.min(x_data)))
                else:
                    hill_initial = 1.0
            else:
                hill_initial = 1.0
            
            # Ensure hill slope is reasonable
            hill_initial = np.clip(hill_initial, 0.1, 5.0)
            
            p0 = [ec50_initial, hill_initial]

            # Calculate IC50 using all data points
            optimizer = static.fit_model(x_data, y_data, p0)
            ic50, hill_coeff = optimizer.x

            # Store bootstrap estimates for EC50
            ic50_bootstrap = []
            hill_coeff_bootstrap = []
            mic_bootstrap = []

            # Use multiprocessing Pool to parallelize the bootstrap resampling with starmap
            with mp.Pool(processes=threads) as pool:
                results = pool.starmap(partial(static.bootstrap_resample, mic_percentage=mic_percentage),
                                       [(x_data, y_data, p0, seed) for seed in range(n_bootstrap)])

            # Collect results
            for ic50_, hill_coeff_, mic_ in results:
                ic50_bootstrap.append(ic50_)
                hill_coeff_bootstrap.append(hill_coeff_)
                mic_bootstrap.append(mic_)

            ic50_bootstrap = np.array(ic50_bootstrap)
            hill_coeff_bootstrap = np.array(hill_coeff_bootstrap)
            mic_bootstrap = np.array(mic_bootstrap)

            # get rid of floating point overflow errors
            filtered_ic50_bootstrap = ic50_bootstrap[abs(ic50_bootstrap)<1e20]
            filtered_hill_coeff_bootstrap = hill_coeff_bootstrap[abs(ic50_bootstrap)<1e20]
            filtered_mic_bootstrap = mic_bootstrap[abs(ic50_bootstrap)<1e20]

            # Compute the confidence interval (e.g., 95% confidence interval, 2.5th and 97.5th percentiles)
            ic50_ci_lower = np.percentile(filtered_ic50_bootstrap, 2.5)
            ic50_ci_upper = np.percentile(filtered_ic50_bootstrap, 97.5)

            # Now same for mic
            mic_ci_lower = np.percentile(filtered_mic_bootstrap, 2.5)
            mic_ci_upper = np.percentile(filtered_mic_bootstrap, 97.5)

            # Generate a fine range of dose values for interpolation
            x_fit = np.logspace(np.log10(x_data.min()/2), np.log10(x_data.max()*2), 100)
            y_fit = log_logistic_2pl_(optimizer.x, x_fit) * max_growth_scale

            # Calculate IC50 and MIC using the interpolated values
            mic = log_logistic_2pl_solve(mic_percentage, optimizer.x, ic50)[0]

            # Prepare output
            if (not np.isnan(ic50) and not np.isnan(mic)):
                return {"Status": "PASS", 
                        'IC50': ic50,
                        'MIC': mic,
                        'IC50_ci_lower': ic50_ci_lower,
                        'IC50_ci_upper': ic50_ci_upper,
                        'MIC_ci_lower': mic_ci_lower,
                        'MIC_ci_upper': mic_ci_upper,
                        'max_growth': max_growth_scale,
                        'max_growth_ci_lower': max_growth_ci_lower,
                        'max_growth_ci_upper': max_growth_ci_upper,
                        'hill_coeff': hill_coeff,
                        'ic50_threshold': log_logistic_2pl_(optimizer.x, ic50) * max_growth_scale,
                        'mic_threshold': log_logistic_2pl_(optimizer.x, mic) * max_growth_scale,
                        'x_fit': json.dumps(x_fit.tolist()) if len(x_fit)>1 else '[]',
                        'y_fit': json.dumps(y_fit.tolist()) if len(y_fit)>1 else '[]',
                        'ic50_bootstrap': json.dumps(filtered_ic50_bootstrap.tolist())
                        }
            else:
                return {"Status": "FAIL", 
                        'IC50': ic50,
                        'MIC': mic,
                        'IC50_ci_lower': ic50_ci_lower,
                        'IC50_ci_upper': ic50_ci_upper,
                        'MIC_ci_lower': mic_ci_lower,
                        'MIC_ci_upper': mic_ci_upper,
                        'max_growth': max_growth_scale,
                        'max_growth_ci_lower': max_growth_ci_lower,
                        'max_growth_ci_upper': max_growth_ci_upper,
                        'hill_coeff': hill_coeff,
                        'ic50_threshold': np.nan,
                        'mic_threshold': np.nan,
                        'x_fit': json.dumps(x_fit.tolist()) if len(x_fit)>1 else '[]',
                        'y_fit': json.dumps(y_fit.tolist()) if len(y_fit)>1 else '[]',
                        'ic50_bootstrap': json.dumps(filtered_ic50_bootstrap.tolist())
                        }
        except Exception as e:
            return {'Status': "FAIL",
                    'IC50': np.nan,
                    'MIC': np.nan,
                    'IC50_ci_lower': np.nan,
                    'IC50_ci_upper': np.nan,
                    'MIC_ci_lower': np.nan,
                    'MIC_ci_upper': np.nan,
                    'max_growth': np.nan,
                    'max_growth_ci_lower': np.nan,
                    'max_growth_ci_upper': np.nan,
                    'hill_coeff': np.nan,
                    'ic50_threshold': np.nan,
                    'mic_threshold': np.nan,
                    'x_fit': '[]',
                    'y_fit': '[]',
                    'ic50_bootstrap': '[]',
                    'error:': str(e)}

    def apply_phenotyper(df_analysis, valid_combinations, metric='OD_final'):

        # Loop over each row in the valid combinations DataFrame
        df_phenotyped = valid_combinations.copy()
        for ix, row in valid_combinations.iterrows():
            #phenotyped = static.phenotyper(df_analysis, row['Antibiotic'], row['Replicate'], row['Strain'])
            phenotyped = static.robust_phenotyper(df_analysis, row['Antibiotic'], row['Strain'])
            for key, val in phenotyped.items():
                df_phenotyped.loc[ix, key] = val

        return df_phenotyped

    def cap_growth_features_within_experiment_range(growth_features):
        # Introduce insufficient drug flag for handling IC50 estimations not covered by experimental range
        capped = growth_features.copy()
        capped['insufficient_drug'] = False
        for ix, row in capped.iterrows():
            if len(row['x_fit']) > 2:
                x_fit = json.loads(row['x_fit'])
                median_fold_change_in_drug_range = np.median(np.exp(np.diff(np.log(sorted(x_fit)))))
                allowed_ic50_cap = (max(x_fit) * median_fold_change_in_drug_range**2)
                if allowed_ic50_cap < row['IC50']:
                    capped.loc[ix, 'IC50'] = allowed_ic50_cap
                    capped.loc[ix, 'MIC'] = allowed_ic50_cap
                    capped.loc[ix, 'IC50_ci_upper'] = allowed_ic50_cap
                    capped.loc[ix, 'IC50_ci_lower'] = allowed_ic50_cap
                    capped.loc[ix, 'MIC_ci_upper'] = allowed_ic50_cap
                    capped.loc[ix, 'MIC_ci_lower'] = allowed_ic50_cap
                    capped.loc[ix, 'insufficient_drug'] = True
        return capped
    
    def plot_dose_response_curve_fit(df_analysis, growth_features, strain_colors, strain, antibiotic, replicate="", ax=None):

        def display_value(value):
            """Check if the value can be converted to float for plotting purposes."""
            try:
                float_value = float(value)
                return f"{float_value:.1e} μg/mL", float_value, True
            except ValueError:
                return f"{value} μg/mL", None, False
        # Define marker style and strain colors
        marker = 'x'
        color = strain_colors.get(strain.lower(), 'gray')

        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4), dpi=100)
        # Set plot parameters
        ax.set_ylim(-0.01, 1.01)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_yticklabels([f'{ytick:.2f}' for ytick in [0, 0.25, 0.5, 0.75, 1]], fontsize=8)

        # Filter and prepare data
        filtered_data = static.dose_response_data_selector(df_analysis, strain, antibiotic)
        display_text = f"{strain}:{replicate}"
        all_doses = np.array([])

        if not filtered_data.empty:
            filtered_data = filtered_data.sort_values(by='Dose')
            non_zero_doses = filtered_data['Dose'][filtered_data['Dose'] > 0]
            if non_zero_doses.any():
                min_nonzero_dose = non_zero_doses.min()
                adjusted_zero_dose = min_nonzero_dose / 2 if min_nonzero_dose > 0 else 0.1
                temp_doses = filtered_data['Dose'].replace(0, adjusted_zero_dose)
                ax.semilogx(temp_doses, filtered_data['OD_final'], marker=marker, markersize=4, linestyle='', color=color, markerfacecolor='none', label=f'Replicate {replicate}')
                all_doses = np.append(all_doses, non_zero_doses)

        # Results data handling
        filtered_results = growth_features[(growth_features['Strain'].str.lower() == strain.lower()) &
                           (growth_features['Antibiotic'].str.lower() == antibiotic.lower())]
        if replicate != "":
            filtered_results = filtered_results[filtered_results['Replicate'] == replicate]
        result_data = filtered_results.iloc[0] if not filtered_results.empty else {}
        MIC = result_data.get('MIC', np.nan)
        IC50 = result_data.get('IC50', np.nan)
        insufficient_drug = result_data.get('insufficient_drug', False)
        ic50_threshold = result_data.get('ic50_threshold', np.nan)
        mic_threshold = result_data.get('mic_threshold', np.nan)
        max_growth_scale = result_data.get('max_growth_scale', np.nan)

        # Plotting fit data if available
        x_fit, y_fit = (np.array(json.loads(result_data.get('x_fit', '[]'))), np.array(json.loads(result_data.get('y_fit', '[]'))))
        if x_fit.size and y_fit.size:
            ax.semilogx(x_fit, y_fit, linestyle='--', color='black', alpha=0.5, label='Fit')

        # Set x-ticks and x-limits dynamically based on non-zero doses from the actual data
        if all_doses.size > 0:
            unique_doses = non_zero_doses.unique()
            ax.set_xticks(unique_doses)
            # Prepare labels for the x-ticks, starting with '0' for the pseudo-zero dose
            labels = [f'{dose:.1e}' for dose in unique_doses]
            ax.set_xticklabels(labels, fontsize=8, fontweight='normal', rotation=60)
            ax.set_xlim(left=adjusted_zero_dose / 1.05, right=1.05 * unique_doses.max())
        else:
            ax.set_xticks([])
            ax.set_xticklabels([])

        ax.set_xlabel(f'{antibiotic} (μg/mL)', fontsize=10)
        ax.text(0.02, 0.96, display_text, transform=ax.transAxes, fontsize=9, verticalalignment='top')

        # Display and plot values for MIC and IC50
        ic50_text, ic50_float, plot_ic50 = display_value(IC50)
        mic_text, mic_float, plot_mic = display_value(MIC)
        # Display dagger symbol if insufficient drug flag is true
        dagger_text = '†' if insufficient_drug else ''
        ax.text(0.99, 0.9, f'IC50{dagger_text}: {ic50_text}', transform=ax.transAxes, fontsize=9, fontweight='bold', color='magenta', horizontalalignment='right')
        ax.text(0.99, 0.8, f'MIC{dagger_text}: {mic_text}', transform=ax.transAxes, fontsize=9, fontweight='bold', color='blue', horizontalalignment='right')
        # ci_upper = float(result_data['IC50_ci_upper'])

        # # Print R_squared value
        # ax.text(0.72, 0.64, f'CI95 = {ci_upper:.2f}', transform=ax.transAxes, fontsize=9, fontweight='bold',        horizontalalignment='left')
        # ax.grid(True, linestyle=':')
        
        # Finalize x-axis limits after all plot elements
        ax.figure.canvas.draw()
        final_x_limits = ax.get_xlim()

        # Calculate normalized values for MIC and IC50 for plotting horizontal lines
        normalized_MIC = 1 if not plot_mic else (np.log10(mic_float) - np.log10(final_x_limits[0])) / (np.log10(final_x_limits[1]) - np.log10(final_x_limits[0]))
        normalized_IC50 = 1 if not plot_ic50 else (np.log10(ic50_float) - np.log10(final_x_limits[0])) / (np.log10(final_x_limits[1]) - np.log10(final_x_limits[0]))

        # Drawing horizontal lines with dynamic xmax based on normalized values or full width
        ax.plot([final_x_limits[0], IC50], [ic50_threshold, ic50_threshold], color='magenta', linestyle=':', linewidth=0.85)
        ax.plot([IC50, IC50], [0, ic50_threshold], color='magenta', linestyle=':', linewidth=0.85)
        ax.plot([final_x_limits[0], MIC], [mic_threshold, mic_threshold], color='blue', linestyle=':', linewidth=0.85)
        ax.plot([MIC, MIC], [0, mic_threshold], color='blue', linestyle=':', linewidth=0.85)

        # if plot_ic50 and ic50_float:
        #     ax.axvline(x=ic50_float, color='magenta', linestyle=':', ymin=0, ymax=ic50_threshold/1.1)
        # if plot_mic and mic_float:
        #     ax.axvline(x=mic_float, color='blue', linestyle=':', ymin=0, ymax=mic_threshold)

        ax.grid(False)#, linestyle=':')

class dynamic:

    def read_biotek_data(file_path):
        try:
            # Get all sheet names
            sheet_names = pd.ExcelFile(file_path).sheet_names
            dfs = []
            # Iterate through each sheet
            for sheet_name in sheet_names:
                # Extract Plate ID from sheet name
                plate_id = sheet_name.split('-')[0].strip()
                # Read the data from the sheet
                df = pd.read_excel(file_path, sheet_name=sheet_name)
                # Identify the row number that data matrix starts from
                data_start_row = (df.iloc[:,2] == "T° 600")
                data_start_row = data_start_row[data_start_row].index.values[0]
                # Extract header
                header = df.iloc[data_start_row,:]
                # Replace header
                df.columns = header 
                # Slice out growth data matrix
                df = df.iloc[data_start_row+1:,:]
                # Rename 'T° OD:600' to 'Temperature'
                df.rename(columns={'T° OD:600': 'Temperature'}, inplace=True)
                df.rename(columns={'T° 600': 'Temperature'}, inplace=True)
                # Drop rows with missing data
                df = df.dropna(subset=['Temperature'])
                # Drop empty columns
                df = df.dropna(axis=1, how='all')
                # Filter the DataFrame to include rows with time values in the format 'HH:MM:SS'
                time_format = r'^\d{2}:\d{2}:\d{2}$'
                df = df[df['Time'].astype(str).str.match(time_format)]
                # Add Plate_ID column
                df['Plate_ID'] = plate_id
                dfs.append(df)
            # Concatenate all DataFrames
            result_df = pd.concat(dfs, ignore_index=True)
            # Convert 'Time' to timedelta since the start
            result_df['Time'] = pd.to_timedelta(result_df['Time'].astype(str))
            start_time = result_df['Time'].min()
            result_df['Hours'] = (result_df['Time'] - start_time).dt.total_seconds() / 3600
            return result_df
        except Exception as e:
            # Return the error message
            return str(e)
    
    def read_tecan_data(file_path):
        """
        To be implemented
        """
        return

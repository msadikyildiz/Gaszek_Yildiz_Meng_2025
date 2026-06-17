# Epistasis pipeline coordinator: orchestrates the sequential execution of all steps
import subprocess
from pathlib import Path
import os
import sys

# Import project configuration
from utils.config import PROJECT_PATH

def run_step(script_name):
    """Execute a Python script in the utils folder and handle any errors"""
    print(f"\nRunning {script_name}...")
    script_path = PROJECT_PATH / "src" / "utils" / script_name
    try:
        subprocess.run([sys.executable, script_path], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running {script_name}: {e}")
        raise

def main():
    """Execute the complete epistasis calculation pipeline in the correct sequence"""
    # Run each step in sequence - order is important due to dependencies between steps
    steps = [
        'generate_g_matrix.py',      # Step 1: Generate genotype mapping matrix
        'generate_hv_matrices.py',   # Step 2: Generate Hadamard and coefficient matrices
        'generate_x_matrix.py',      # Step 3: Generate genotype encoding matrix
        'generate_w_matrices.py',    # Step 4: Generate epistasis calculation matrices
        'calculate_epistasis.py'     # Step 5: Calculate epistasis for all drug conditions
    ]
    
    print("Starting TEM1 epistasis analysis pipeline...")
    print(f"Using project path: {PROJECT_PATH}")
    
    for step in steps:
        run_step(step)
    
    print("\nEpistasis analysis pipeline completed successfully!")
        
if __name__ == "__main__":
    main()

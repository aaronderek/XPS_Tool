
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import re

import sys

def generate_graphs():
    # Setup directories
    if len(sys.argv) > 1:
        target_dir = sys.argv[1]
    else:
        target_dir = os.getcwd()
        
    print(f"Target Directory: {target_dir}")
    
    if not os.path.exists(target_dir):
        print(f"Error: Directory '{target_dir}' not found.")
        return

    output_dir = os.path.join(target_dir, 'Graphs')
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")

    # style configuration for HD plots
    plt.rcParams['figure.figsize'] = (10, 8)
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['font.size'] = 14
    plt.rcParams['axes.labelsize'] = 16
    plt.rcParams['axes.titlesize'] = 18
    plt.rcParams['xtick.major.size'] = 6
    plt.rcParams['xtick.major.width'] = 1.5
    plt.rcParams['ytick.major.size'] = 6
    plt.rcParams['ytick.major.width'] = 1.5
    
    # search for csv files
    csv_files = glob.glob(os.path.join(target_dir, '*.csv'))
    
    print(f"Found {len(csv_files)} CSV files.")

    for file_path in csv_files:
        try:
            filename = os.path.basename(file_path)
            print(f"Processing: {filename}")
            
            # Read CSV
            df = pd.read_csv(file_path)
            
            # Basic validation
            if 'E' not in df.columns or 'y' not in df.columns:
                print(f"Skipping {filename}: 'E' or 'y' columns missing.")
                continue

            # Plotting
            fig, ax = plt.subplots()
            
            # Plot data
            ax.plot(df['E'], df['y'], color='#1f77b4', label='Experimental Data')
            
            # Formatting
            ax.set_xlabel('Binding Energy (eV)', fontweight='bold')
            ax.set_ylabel('Intensity (a.u.)', fontweight='bold')
            
            # Construct a readable title from filename
            # Example: In2O3 6000 shots_0106_XPS_1_20260124_C1s.csv -> In2O3 6000 shots - C1s
            base_name = os.path.splitext(filename)[0]
            
            # Attempt to extract relevant parts (This is heuristic based on user files)
            # We want the material/experiment name and the core level (e.g. C1s, O1s)
            parts = base_name.split('_')
            if len(parts) >= 5:
                # Based on: In2O3 6000 shots_0106_XPS_1_20260124_C1s
                # parts[0]: In2O3 6000 shots
                # parts[-1]: C1s
                title_text = f"{parts[0]} - {parts[-1]}"
            else:
                title_text = base_name
                
            ax.set_title(title_text, pad=15, fontweight='bold')
            
            # Clean up plot area
            ax.minorticks_on()
            ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
            
            # XPS standard: Reverse X axis (High Binding Energy -> Low Binding Energy)
            ax.set_xlim(df['E'].max(), df['E'].min())
            
            # Tight layout
            plt.tight_layout()
            
            # Save
            save_name = os.path.splitext(filename)[0] + '.png'
            save_path = os.path.join(output_dir, save_name)
            plt.savefig(save_path, bbox_inches='tight', dpi=300)
            print(f"Saved: {save_path}")
            
            plt.close(fig)

        except Exception as e:
            print(f"Error processing {file_path}: {e}")

if __name__ == "__main__":
    generate_graphs()

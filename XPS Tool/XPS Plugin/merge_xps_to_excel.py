import pandas as pd
import os
import glob
import sys

import re

def merge_csvs_to_excel():
    # Setup target directory
    if len(sys.argv) > 1:
        target_dir = sys.argv[1]
    else:
        target_dir = os.getcwd()
        
    print(f"Target Directory: {target_dir}")
    
    if not os.path.exists(target_dir):
        print(f"Error: Directory '{target_dir}' not found.")
        return

    # Find CSV files
    csv_files = glob.glob(os.path.join(target_dir, '*.csv'))
    
    if not csv_files:
        print("No CSV files found in the target directory.")
        return

    print(f"Found {len(csv_files)} CSV files. Merging...")

    output_file = os.path.join(target_dir, 'Merged_XPS_Data.xlsx')

    # Regex for user-requested format: 0103_XPS_1_Su1s
    # Filename: ..._0103_XPS_1_20260124_Su1s.csv
    # Pattern: ..._(ID)_(TYPE)_(INDEX)_DATE_(ELEMENT).csv
    pattern = re.compile(r'_(\d{4})_(XPS)_(\d+)_\d+_([A-Za-z0-9]+)\.csv$', re.IGNORECASE)

    try:
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            used_sheet_names = {}
            
            for file_path in csv_files:
                filename = os.path.basename(file_path)
                
                # Try to match the user's specific naming convention
                match = pattern.search(filename)
                
                if match:
                    # Extracted parts: 0103, XPS, 1, Su1s
                    # format: 0103_XPS_1_Su1s
                    sheet_name = f"{match.group(1)}_{match.group(2)}_{match.group(3)}_{match.group(4)}"
                else:
                    # Fallback to smart truncation if pattern doesn't match
                    name_no_ext = os.path.splitext(filename)[0]
                    if len(name_no_ext) <= 31:
                        sheet_name = name_no_ext
                    else:
                        sheet_name = name_no_ext[:14] + '...' + name_no_ext[-14:]
                
                # Ensure uniqueness
                original_sheet_name = sheet_name
                counter = 1
                while sheet_name in used_sheet_names:
                    suffix = f"_{counter}"
                    allowed_len = 31 - len(suffix)
                    sheet_name = original_sheet_name[:allowed_len] + suffix
                    counter += 1
                
                used_sheet_names[sheet_name] = True
                
                # Write to Excel
                try:
                    df = pd.read_csv(file_path)
                    df.to_excel(writer, sheet_name=sheet_name, index=False)
                    print(f"Added {filename} as sheet '{sheet_name}'")
                except Exception as e:
                    print(f"Failed to process {filename}: {e}")
                    
        print(f"Successfully created: {output_file}")
        
    except Exception as e:
        print(f"Error creating Excel file: {e}")

if __name__ == "__main__":
    merge_csvs_to_excel()

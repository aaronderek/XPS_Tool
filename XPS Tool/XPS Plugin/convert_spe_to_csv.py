from yadg.extractors import extract
import pandas as pd
import glob
import os
import sys

def convert_spe_files():
    """
    Converts .spe files (PHI Multipak) to .csv using yadg.
    
    Requirements:
        pip install "pydantic<2.0" yadg pandas xarray
        
    Note: 'pydantic<2.0' is currently required for yadg/dgbowl-schemas compatibility.
    """
    
    # Check for .spe files
    spe_files = glob.glob("*.spe")
    if not spe_files:
        print("No .spe files found in the current directory.")
        return

    print(f"Found {len(spe_files)} .spe files. Starting conversion...")

    for spe_file in spe_files:
        print(f"Processing {spe_file}...")
        try:
            # Extract data using yadg
            # 'phi.spe' is the filetype for PHI Multipak SPE files
            ds = extract(filetype='phi.spe', path=spe_file)
            
            # The result is a DataTree with children representing traces
            for trace_name in ds.keys():
                trace_node = ds[trace_name]
                
                # Convert the child node's Dataset to a DataFrame
                # We access .ds to get the xarray.Dataset
                if hasattr(trace_node, 'ds'):
                    df = trace_node.ds.to_dataframe()
                else:
                    # Fallback if structure is different (e.g. older datatree)
                    df = trace_node.to_dataframe()

                # Construct output filename
                base_name = os.path.splitext(spe_file)[0]
                output_csv = f"{base_name}_{trace_name}.csv"
                
                # Save to CSV
                df.to_csv(output_csv)
                print(f"Successfully converted {spe_file} trace {trace_name} to {output_csv}")
            
        except ImportError:
            print("Error: Missing libraries.")
            print("Run: pip install \"pydantic<2.0\" yadg pandas xarray")
            return
        except Exception as e:
            print(f"Failed to convert {spe_file}: {e}")
            if "pydantic" in str(e) or "validator" in str(e):
                print("Hint: This might be a Pydantic version issue. Try: pip install \"pydantic<2.0\"")

if __name__ == "__main__":
    convert_spe_files()

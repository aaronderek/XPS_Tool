import subprocess
import glob
import os
import h5py
import xarray as xr
import pandas as pd
import sys

def convert():
    print("Starting batch conversion...")
    
    # Locate yadg executable in current venv
    scripts_dir = os.path.dirname(sys.executable)
    yadg_exe = os.path.join(scripts_dir, "yadg.exe")
    if not os.path.exists(yadg_exe):
        # Fallback to assumming it's in PATH
        yadg_exe = "yadg"
        
    spe_files = glob.glob("*.spe")
    if not spe_files:
        print("No .spe files found.")
        return

    for spe in spe_files:
        print(f"Processing {spe}...")
        nc_file = spe + ".temp.nc"
        
        # Run yadg extraction
        try:
            # We capture output to avoid spam, but print on error
            result = subprocess.run(
                [yadg_exe, "extract", "phi.spe", spe, nc_file],
                capture_output=True, text=True
            )
            if result.returncode != 0:
                print(f"yadg failed for {spe}:")
                print(result.stderr)
                # If it failed but created the file (e.g. partial?), check existence
                if not os.path.exists(nc_file):
                    continue
            else:
                print(f"Intermediate file info: {result.stdout}")
        except Exception as e:
            print(f"Execution error: {e}")
            continue

        if not os.path.exists(nc_file):
            print(f"Failed to generate {nc_file}")
            continue
            
        # Extract traces from NeXus/HDF5 file
        try:
            with h5py.File(nc_file, 'r') as f:
                # Top level keys are usually the traces (groups)
                trace_names = list(f.keys())
            
            print(f"Found traces: {trace_names}")
            
            for trace in trace_names:
                 try:
                     # Open specific group
                     ds = xr.open_dataset(nc_file, group=trace, engine="h5netcdf")
                     df = ds.to_dataframe().reset_index()
                     
                     # Construct CSV name: Original_Trace.csv
                     out_csv = f"{os.path.splitext(spe)[0]}_{trace}.csv"
                     df.to_csv(out_csv, index=False)
                     print(f"  -> Saved {out_csv}")
                     ds.close()
                 except Exception as e:
                     print(f"  Failed to save trace {trace}: {e}")
                 
        except Exception as e:
            print(f"Error extracting traces from {nc_file}: {e}")
        finally:
            # Clean up temp file
            if os.path.exists(nc_file):
                os.remove(nc_file)

    print("Conversion complete.")

if __name__ == "__main__":
    convert()

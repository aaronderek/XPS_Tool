# XPS Data Processing Tools

This document provides instructions for using the Python tools created to process XPS CSV data.

## Prerequisites
- **Python 3.x**: Must be installed.
- **Libraries**: `pandas`, `matplotlib`, `openpyxl`.
  - Install command: `pip install pandas matplotlib openpyxl`

## Tool 1: Generate HD Graphs
**File:** `generate_xps_graphs.py`

This tool converts individual CSV files into high-definition PNG charts.

### Usage
1.  Open a terminal.
2.  Run the script with the target directory as an argument:
    ```powershell
    python generate_xps_graphs.py "path/to/data_folder"
    ```
    *Example:* `python generate_xps_graphs.py "Converted_XPS_In2O3_6000SHS"`

### Output
- Creates a `Graphs` subdirectory inside the target folder.
- Saves `.png` files for each `.csv` found.

## Tool 2: Merge to Excel
**File:** `merge_xps_to_excel.py`

This tool combines all CSV files in a folder into a single Excel file, with each CSV on a separate tab.

### Usage
1.  Open a terminal.
2.  Run the script with the target directory as an argument:
    ```powershell
    python merge_xps_to_excel.py "path/to/data_folder"
    ```
    *Example:* `python merge_xps_to_excel.py "Converted_XPS_In2O3_6000SHS"`

### Output
- Creates a file named `Merged_XPS_Data.xlsx` in the target folder.
- **Sheet Naming**: The tool attempts to extract the Experiment #, Type, Index, and Element from the filename to create concise sheet names (e.g., `0103_XPS_1_Su1s`).
    - If the filename pattern doesn't match, it falls back to truncating the filename to 31 characters.

---
**Note:** Both scripts default to the current working directory if no argument is provided.

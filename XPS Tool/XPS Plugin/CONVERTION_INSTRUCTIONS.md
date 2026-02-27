# Converting SPE Files to CSV

This guide explains how to convert PHI XPS `.spe` files to CSV format using the provided Python scripts and environment. It consolidates findings from previous troubleshooting sessions.

## 1. Quick Start

**Use the pre-configured environment.** Do not try to install dependencies in a new environment unless necessary, as the `yadg` library in `.venv311` contains critical patches for these specific files.

1.  Open a terminal in `c:\Users\robin\Files\XPS results 1`.
2.  Run the conversion script using the `.venv311` python executable:
    ```powershell
    & ".\.venv311\Scripts\python.exe" convert_spe_to_csv.py
    ```
3.  The script will:
    *   Find all `.spe` files in the current folder.
    *   Extract each spectral region (trace) (e.g., C1s, O1s).
    *   Save them as separate CSV files named `[OriginalName]_[TraceName].csv`.

---

## 2. Environment & Prerequisites

### The `.venv311` Environment
*   **Python Version**: 3.11
*   **Critical Libraries**:
    *   `yadg==6.2.2` (**Patched** to handle non-standard SPE headers)
    *   `pydantic` (Version 2.x is supported with the current script)
    *   `pandas`, `xarray`

> [!IMPORTANT]
> **Do not upgrade or reinstall `yadg`**. The version installed in `.venv311` has been manually patched to fix bugs related to Windows line endings (`\r\n`) and zero-value data offsets found in your specific `.spe` files. Reinstalling from PyPI will break the conversion.

---

## 3. Findings & Troubleshooting

### A. "ModuleNotFoundError: No module named 'yadg.extract'"
*   **Cause**: The API structure of `yadg` was misunderstood or changed.
*   **Fix**: The script now correctly imports `from yadg.extractors import extract`.

### B. "DataTree object has no attribute 'to_dataframe'"
*   **Cause**: `yadg.extract()` returns an `xarray.DataTree` object, which is a tree-like structure containing multiple datasets (one for each spectral region/trace in the file), rather than a single DataFrame.
*   **Fix**: The script iterates through the children of the `DataTree`. Each child represents a trace (e.g., 'C1s', 'O1s') and contains the actual data in a `.ds` attribute (the `xarray.Dataset`), which can then be converted to a DataFrame.

### C. File Format Deviations (Why the Patch is Needed)
Your `.spe` files differ from the standard PHI format:
1.  **Line Endings**: They use `\r\n` (CRLF) in the ASCII header, whereas standard parsers expect `\n`.
2.  **Duplicate Keys**: Headers contain duplicate `SpectralRegDef` keys which confuse standard parsers.
3.  **Sequential Data**: The binary header specifies `end_of_data = 0`, meaning data is written sequentially without absolute file offsets. Use the patched `yadg` to handle this.

### D. Pydantic Version
*   Early findings suggested `pydantic<2.0` was required.
*   **Correction**: `pydantic>=2.0` (tested with 2.12.5) works correctly with the updated script and patched `yadg` version 6.2.2.

---

## 4. Script Reference (`convert_spe_to_csv.py`)

The script has been updated to robustly handle the data structure:

```python
from yadg.extractors import extract
# ... imports ...

# ... inside loop ...
    ds = extract(filetype='phi.spe', path=spe_file)
    
    # Iterate over every trace (child node) in the DataTree
    for trace_name in ds.keys():
        trace_node = ds[trace_name]
        # Access the Dataset via .ds and convert
        df = trace_node.ds.to_dataframe()
        # Save to CSV...
```

---

## 5. Project Structure

Current layout of the software and data:

```text
.
├── .venv311/                     # [CRITICAL] Python 3.11 Env with patched yadg 6.2.2
├── convert_spe_to_csv.py         # The conversion script
├── INSTRUCTIONS.md               # This guide
├── XPS In2O3 Si111 6000shs/      # Source directory containing .spe files
└── Converted_***  # Output directory for converted .csv files
```

### Key Files
- **`.venv311/`**: Contains the specific python executable you MUST use. Do not delete or modify this folder.
- **`convert_spe_to_csv.py`**: The script that performs the work. It is configured to run in the root directory and find `.spe` files in the current or subdirectories (depending on how you run it, currently setup to look in the *current* directory but logic supports paths).


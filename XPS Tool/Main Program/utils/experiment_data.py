"""
Experimental data import and validation utilities.
"""

import numpy as np
import pandas as pd
from typing import Dict, Any, Tuple, Optional, List
import io


def validate_csv_data(df: pd.DataFrame) -> Tuple[bool, str]:
    """
    Validate imported CSV data.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe to validate

    Returns
    -------
    is_valid : bool
        True if data is valid
    message : str
        Validation message or error description
    """
    # Check if dataframe is empty
    if df.empty:
        return False, "CSV file is empty"

    # Check for NaN values
    if df.isnull().values.any():
        return False, "Data contains NaN values. Please clean your data."

    # Check data ranges
    if 'Delta_WF_eV' in df.columns:
        if np.any(np.abs(df['Delta_WF_eV']) > 3.0):
            return False, "ΔWF values exceed ±3 eV. Please check your data."

    if 'Delta_CL_eV' in df.columns:
        if np.any(np.abs(df['Delta_CL_eV']) > 2.0):
            return False, "ΔE_CL values exceed ±2 eV. Please check your data."

    return True, "Data validated successfully"


def validate_experimental_data_physics(data: Dict[str, np.ndarray]) -> List[str]:
    """
    Validate physical reasonableness of experimental data.

    Checks for common issues:
    - Negative correlation between ΔWF and ΔE_CL (should be positive)
    - Unusual η values (should be 0.5-0.95 for typical TCOs)
    - Insufficient data range

    Parameters
    ----------
    data : dict
        Data dictionary with 'Delta_WF' and 'Delta_CL' arrays

    Returns
    -------
    warnings : list of str
        List of warning messages (empty if no issues)
    """
    warnings = []

    Delta_WF = data['Delta_WF']
    Delta_CL = data['Delta_CL']

    # Need at least 3 points for meaningful analysis
    if len(Delta_WF) < 3:
        return warnings

    # Check 1: Correlation between ΔWF and ΔE_CL
    from scipy.stats import pearsonr
    try:
        corr, p_value = pearsonr(Delta_WF, Delta_CL)

        if corr < 0:
            warnings.append(
                "⚠️ **Warning: Negative correlation detected**\n"
                f"ΔWF and ΔE_CL have correlation coefficient = {corr:.3f}\n\n"
                "**Expected:** Positive correlation (both negative or both positive)\n\n"
                "**Typical annealing scenarios:**\n"
                "- Vacuum annealing (desorption): Both ΔWF and ΔE_CL < 0\n"
                "- Oxygen annealing (adsorption): Both ΔWF and ΔE_CL > 0\n\n"
                "**Possible causes:**\n"
                "- Sign convention error in data processing\n"
                "- Mixed experimental conditions\n"
                "- Incorrect reference point"
            )
        elif corr < 0.9 and p_value < 0.05:
            warnings.append(
                f"⚠️ **Notice: Low correlation** (r = {corr:.3f}, p = {p_value:.3f})\n"
                "Typical high-quality data has r > 0.95\n"
                "This might indicate experimental noise or complex behavior."
            )
    except Exception:
        pass  # Skip correlation check if it fails

    # Check 2: η value reasonableness
    try:
        # Linear fit to get slope (η)
        slope = np.polyfit(Delta_WF, Delta_CL, 1)[0]

        if slope < 0:
            warnings.append(
                f"⚠️ **Warning: Negative slope** (η = {slope:.3f})\n\n"
                "**Expected:** η ∈ (0.5, 0.95) for typical TCOs\n\n"
                "The slope η = ΔE_CL / ΔWF represents XPS sampling depth.\n"
                "A negative value is physically unreasonable."
            )
        elif slope < 0.5:
            warnings.append(
                f"⚠️ **Warning: Unusually low η** (η = {slope:.3f})\n\n"
                "**Typical range:** η ∈ (0.7, 0.9) for In₂O₃\n\n"
                "Low η might indicate:\n"
                "- Very thick depletion layer (W >> λ)\n"
                "- Surface contamination layer\n"
                "- Incorrect binding energy calibration"
            )
        elif slope > 0.95:
            warnings.append(
                f"⚠️ **Warning: Unusually high η** (η = {slope:.3f})\n\n"
                "**Typical range:** η ∈ (0.7, 0.9) for In₂O₃\n\n"
                "High η might indicate:\n"
                "- Very thin depletion layer (W << λ)\n"
                "- Long photoelectron escape depth\n"
                "- Check λ and W parameters"
            )
    except Exception:
        pass  # Skip slope check if it fails

    # Check 3: Data range
    Delta_WF_range = np.max(Delta_WF) - np.min(Delta_WF)
    if Delta_WF_range < 0.1:
        warnings.append(
            f"⚠️ **Notice: Small ΔWF range** ({Delta_WF_range:.3f} eV)\n\n"
            "For reliable fitting, aim for ΔWF range > 0.2 eV\n\n"
            "**Suggestions:**\n"
            "- Increase annealing temperature range\n"
            "- Improve vacuum conditions\n"
            "- Allow longer equilibration time"
        )

    # Check 4: Both in same direction (both positive or both negative)
    Delta_WF_mostly_positive = np.sum(Delta_WF > 0) > len(Delta_WF) / 2
    Delta_CL_mostly_positive = np.sum(Delta_CL > 0) > len(Delta_CL) / 2

    if Delta_WF_mostly_positive != Delta_CL_mostly_positive:
        warnings.append(
            "⚠️ **Warning: Mixed signs**\n\n"
            "ΔWF and ΔE_CL should have the same sign trend:\n"
            "- Both negative (typical vacuum annealing)\n"
            "- Both positive (typical oxidation)\n\n"
            "Current data shows opposite trends, suggesting:\n"
            "- Sign convention mismatch\n"
            "- Complex surface chemistry"
        )

    return warnings


def detect_format(df: pd.DataFrame) -> str:
    """
    Detect the format of the CSV file.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe

    Returns
    -------
    format_type : str
        'processed' or 'raw' or 'unknown'
    """
    # Check for processed format (Delta values already calculated)
    if 'Delta_WF_eV' in df.columns and 'Delta_CL_eV' in df.columns:
        return 'processed'

    # Check for raw format (absolute values)
    if 'WF_eV' in df.columns and ('CL_eV' in df.columns or 'In3d_eV' in df.columns):
        return 'raw'

    return 'unknown'


def process_raw_data(df: pd.DataFrame) -> Dict[str, np.ndarray]:
    """
    Process raw measurement data.

    Calculates Delta values relative to first measurement point.

    Parameters
    ----------
    df : pd.DataFrame
        Raw data with 'T_degC', 'WF_eV', and core level columns

    Returns
    -------
    data : dict
        Processed data dictionary
    """
    # Get temperature
    T = df['T_degC'].values if 'T_degC' in df.columns else np.arange(len(df))

    # Calculate Delta WF
    WF = df['WF_eV'].values
    Delta_WF = WF - WF[0]

    # Calculate Delta CL (check for different column names)
    if 'CL_eV' in df.columns:
        CL = df['CL_eV'].values
    elif 'In3d_eV' in df.columns:
        CL = df['In3d_eV'].values
    elif 'O1s_eV' in df.columns:
        CL = df['O1s_eV'].values
    else:
        raise ValueError("No core level column found (expected 'CL_eV', 'In3d_eV', or 'O1s_eV')")

    Delta_CL = CL - CL[0]

    return {
        'T_degC': T,
        'Delta_WF': Delta_WF,
        'Delta_CL': Delta_CL,
        'WF_abs': WF,
        'CL_abs': CL
    }


def process_processed_data(df: pd.DataFrame) -> Dict[str, np.ndarray]:
    """
    Process already-processed data.

    Parameters
    ----------
    df : pd.DataFrame
        Processed data with Delta values

    Returns
    -------
    data : dict
        Data dictionary
    """
    data = {
        'Delta_WF': df['Delta_WF_eV'].values,
        'Delta_CL': df['Delta_CL_eV'].values,
    }

    # Optional columns
    if 'T_degC' in df.columns:
        data['T_degC'] = df['T_degC'].values

    if 'ns_cm2' in df.columns:
        data['ns'] = df['ns_cm2'].values * 1e17  # Convert to m^-2

    return data


def import_experimental_data(file_content: bytes) -> Dict[str, Any]:
    """
    Import and process experimental data from CSV.

    Parameters
    ----------
    file_content : bytes
        Content of uploaded CSV file

    Returns
    -------
    result : dict
        Dictionary containing:
        - 'success': bool
        - 'data': dict or None
        - 'message': str
        - 'format': str
    """
    try:
        # Read CSV
        df = pd.read_csv(io.BytesIO(file_content))

        # Detect format
        format_type = detect_format(df)

        if format_type == 'unknown':
            return {
                'success': False,
                'data': None,
                'message': "Unrecognized CSV format. Expected columns: 'T_degC, WF_eV, CL_eV' or 'Delta_WF_eV, Delta_CL_eV'",
                'format': 'unknown'
            }

        # Validate
        is_valid, msg = validate_csv_data(df)
        if not is_valid:
            return {
                'success': False,
                'data': None,
                'message': msg,
                'format': format_type
            }

        # Process data
        if format_type == 'raw':
            data = process_raw_data(df)
        else:
            data = process_processed_data(df)

        return {
            'success': True,
            'data': data,
            'message': f"Successfully imported {len(df)} data points ({format_type} format)",
            'format': format_type
        }

    except Exception as e:
        return {
            'success': False,
            'data': None,
            'message': f"Error reading CSV: {str(e)}",
            'format': 'error'
        }


def get_format_example_text() -> str:
    """
    Get example CSV format text for user guidance.

    Returns
    -------
    example_text : str
        Formatted example text
    """
    example = """
**Format A: Raw Measurements**
```
T_degC,WF_eV,In3d_eV
25,4.20,444.50
100,4.28,444.42
150,4.35,444.35
200,4.42,444.27
```

**Format B: Processed Data**
```
Delta_WF_eV,Delta_CL_eV,T_degC
0.000,0.000,25
0.080,-0.080,100
0.150,-0.150,150
0.220,-0.220,200
```

**Notes:**
- First row must contain column headers
- Decimal separator: period (.)
- Temperature column optional for Format B
- Core level can be named: CL_eV, In3d_eV, O1s_eV, etc.
"""
    return example


def create_sample_data() -> pd.DataFrame:
    """
    Create sample experimental data for testing.

    Simulates In₂O₃ vacuum annealing experiment:
    - Temperature: 25-400°C
    - Effect: Water desorption → surface potential increases
    - Result: Work function decreases, core level binding energy decreases
    - Both ΔWF and ΔE_CL are negative and positively correlated
    - η = ΔE_CL / ΔWF ≈ 0.85

    Returns
    -------
    df : pd.DataFrame
        Sample data with physically reasonable values
    """
    # Temperature series
    T = np.array([25, 100, 150, 200, 250, 300, 350, 400])

    # Initial state (25°C, with H₂O adsorption)
    WF_base = 4.20  # eV, relatively high due to adsorbates
    In3d_base = 444.50  # eV

    # Simulate desorption kinetics
    # Physical picture: Water desorption → surface potential increases → WF decreases
    # The magnitude of WF change depends on both surface potential and dipole effects

    # Temperature-normalized coordinate (0 to 1)
    t_norm = (T - 25) / 375

    # Magnitude of total WF change increases with temperature
    # Realistic range for In₂O₃ vacuum annealing: 0 to -0.35 eV
    # Use smooth progression: starts slow, accelerates at 150-250°C, then saturates
    Delta_WF_magnitude = 0.35 * (1 - np.exp(-3 * t_norm))  # Exponential approach to saturation
    Delta_WF = -Delta_WF_magnitude  # Negative sign: WF decreases

    # Core level shift: ΔE_CL ≈ η × ΔWF
    # For typical TCOs with W~3nm, λ~1.8nm, η ≈ 0.85
    eta = 0.85
    Delta_CL = eta * Delta_WF

    # Add realistic noise (±0.01 eV for UPS/XPS)
    np.random.seed(42)  # for reproducibility
    noise_WF = np.random.normal(0, 0.008, len(T))
    noise_CL = np.random.normal(0, 0.008, len(T))

    Delta_WF = Delta_WF + noise_WF
    Delta_CL = Delta_CL + noise_CL

    # Calculate absolute values
    WF = WF_base + Delta_WF
    In3d = In3d_base + Delta_CL  # Both move in same direction (lower BE)

    df = pd.DataFrame({
        'T_degC': T,
        'WF_eV': WF,
        'In3d_eV': In3d
    })

    return df

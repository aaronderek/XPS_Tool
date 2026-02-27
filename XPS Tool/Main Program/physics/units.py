"""
Unit conversion utilities for 2DEG calculations.
"""

import numpy as np
from .constants import (
    EV_TO_JOULE, JOULE_TO_EV,
    NM_TO_METER, METER_TO_NM,
    CM2_TO_M2, M2_TO_CM2,
    DEBYE_TO_CM
)


def eV_to_J(energy_eV):
    """Convert energy from eV to Joules."""
    return energy_eV * EV_TO_JOULE


def J_to_eV(energy_J):
    """Convert energy from Joules to eV."""
    return energy_J * JOULE_TO_EV


def nm_to_m(length_nm):
    """Convert length from nanometers to meters."""
    return length_nm * NM_TO_METER


def m_to_nm(length_m):
    """Convert length from meters to nanometers."""
    return length_m * METER_TO_NM


def cm2_to_m2(area_cm2):
    """Convert area from cm² to m²."""
    return area_cm2 * CM2_TO_M2


def m2_to_cm2(area_m2):
    """Convert area from m² to cm²."""
    return area_m2 * M2_TO_CM2


def ns_to_display(ns_m2):
    """
    Convert sheet density from m⁻² to display units (10¹³ cm⁻²).

    Parameters
    ----------
    ns_m2 : float or array
        Sheet density in m⁻²

    Returns
    -------
    float or array
        Sheet density in units of 10¹³ cm⁻²
    """
    # For densities (counts per area), 1 m⁻² = 10^-4 cm⁻²
    # because if 1 object occupies 1 m², it occupies 10^4 cm²,
    # so the density is 1/(10^4 cm²) = 10^-4 cm⁻²
    ns_cm2 = ns_m2 * CM2_TO_M2  # Convert to cm⁻² (multiply by 1e-4)
    return ns_cm2 / 1e13  # Convert to 10¹³ cm⁻²


def ns_from_display(ns_display):
    """
    Convert sheet density from display units (10¹³ cm⁻²) to m⁻².

    Parameters
    ----------
    ns_display : float or array
        Sheet density in units of 10¹³ cm⁻²

    Returns
    -------
    float or array
        Sheet density in m⁻²
    """
    ns_cm2 = ns_display * 1e13  # Convert to cm⁻²
    return ns_cm2 * M2_TO_CM2  # Convert to m⁻² (multiply by 1e4)


def debye_to_C_m(dipole_debye):
    """Convert dipole moment from Debye to C·m."""
    return dipole_debye * DEBYE_TO_CM


def degrees_to_radians(angle_deg):
    """Convert angle from degrees to radians."""
    return np.deg2rad(angle_deg)


def validate_range(value, valid_range, param_name):
    """
    Validate that a parameter is within its allowed range.

    Parameters
    ----------
    value : float
        Parameter value to validate
    valid_range : tuple
        (min, max) allowed values
    param_name : str
        Name of parameter for error messages

    Returns
    -------
    bool
        True if valid, raises ValueError if not
    """
    if not (valid_range[0] <= value <= valid_range[1]):
        raise ValueError(
            f"{param_name} = {value} is outside valid range "
            f"[{valid_range[0]}, {valid_range[1]}]"
        )
    return True

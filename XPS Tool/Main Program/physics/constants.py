"""
Physical constants for 2DEG calculations.
All values in SI units unless otherwise specified.
"""

import numpy as np

# Fundamental constants
ELEMENTARY_CHARGE = 1.602176634e-19  # C (exact)
VACUUM_PERMITTIVITY = 8.8541878128e-12  # F/m (exact)
HBAR = 1.054571817e-34  # J·s (exact)
ELECTRON_MASS = 9.1093837015e-31  # kg (exact)
BOLTZMANN = 1.380649e-23  # J/K (exact)

# Conversion factors
EV_TO_JOULE = ELEMENTARY_CHARGE  # 1 eV = 1.602e-19 J
JOULE_TO_EV = 1.0 / EV_TO_JOULE
NM_TO_METER = 1e-9
METER_TO_NM = 1e9
CM2_TO_M2 = 1e-4
M2_TO_CM2 = 1e4
DEBYE_TO_CM = 3.33564e-30  # C·m (1 Debye)

# Commonly used combinations
Q = ELEMENTARY_CHARGE
EPS0 = VACUUM_PERMITTIVITY
M0 = ELECTRON_MASS

# Default material parameters
DEFAULT_M_STAR = 0.32 * M0  # Effective mass (typical for SrTiO3)
DEFAULT_EPSILON_R = 9  # Relative permittivity
DEFAULT_TEMPERATURE = 300  # K
DEFAULT_W = 3.0  # nm (depletion width)
DEFAULT_LAMBDA_XPS = 1.8  # nm (inelastic mean free path)

# Parameter ranges
M_STAR_RANGE = (0.30, 0.35)  # in units of m0
EPSILON_R_RANGE = (7, 12)
TEMPERATURE_RANGE = (10, 600)  # K
W_RANGE = (1.0, 10.0)  # nm
PHI_S_RANGE = (0.0, 1.0)  # eV
LAMBDA_XPS_RANGE = (0.5, 3.0)  # nm
THETA_XPS_RANGE = (0, 60)  # degrees
DELTA_PHI_DIP_RANGE = (-0.5, 0.5)  # eV

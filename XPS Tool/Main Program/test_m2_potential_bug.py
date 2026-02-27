"""
Test M2 (Fang-Howard) potential profile to identify the bug.

The issue: M2 gives η_theory = 0.5697 (slope) vs η_single = 0.2949 (point),
which suggests the ΔE_CL vs ΔWF relationship is non-linear.

This could be due to incorrect potential profile V(z).
"""

import numpy as np
from models.fang_howard import FangHowardModel
from models.triangular import TriangularModel
from physics.units import nm_to_m
from physics.constants import M0, Q

# Parameters
W_nm = 3.0
epsilon_r = 9.0
m_star_ratio = 0.32
Phi_s_eV = 0.3

# Create models
m2_model = FangHowardModel(m_star_ratio * M0, epsilon_r, W_nm)
m1_model = TriangularModel(m_star_ratio * M0, epsilon_r, W_nm)

# Create z array
z_array = np.linspace(0, nm_to_m(6), 200)

# Get potentials
V_z_m2 = m2_model.get_potential(Phi_s_eV, z_array)
V_z_m1 = m1_model.get_potential(Phi_s_eV, z_array)

# Get M2 parameters
Es, W_eff, _ = m2_model.solve_self_consistent(Phi_s_eV)

print("="*70)
print("M2 (FANG-HOWARD) POTENTIAL PROFILE ANALYSIS")
print("="*70)

print(f"\nInput: Φs = {Phi_s_eV} eV")
print(f"\nM2 self-consistent solution:")
print(f"  Es = {Es:.3e} V/m")
print(f"  W_eff = {W_eff*1e9:.3f} nm")

print(f"\nM2 potential V(z):")
print(f"  V(z=0) = {V_z_m2[0]/Q:.6f} eV")
print(f"  V(z=W_eff) = {V_z_m2[np.argmin(np.abs(z_array - W_eff))]/Q:.6f} eV")
print(f"  Expected at z=0: 0 eV (surface)")
print(f"  Expected at z=W_eff: {Phi_s_eV:.6f} eV (max potential)")

print(f"\n⚠️ PROBLEM IDENTIFIED:")
print(f"  Current code: V(z) = Es × z (increases from 0)")
print(f"  Physical V(z): Should be Φs - Es×z = Φs(1 - z/W)")
print(f"  → V(z=0) should be Φs = {Phi_s_eV} eV, not 0!")
print(f"  → V(z) should DECREASE from surface to bulk, not increase!")

print(f"\nM1 (Triangular) for comparison:")
print(f"  V(z=0) = {V_z_m1[0]/Q:.6f} eV")
print(f"  V(z=W) = {V_z_m1[np.argmin(np.abs(z_array - nm_to_m(W_nm)))]/Q:.6f} eV")

# Check sign convention
print(f"\n### SIGN CONVENTION CHECK ###")
print(f"In XPS convention:")
print(f"  - Surface potential Φs > 0 (band bending upward)")
print(f"  - Potential energy V(z) at surface should be maximum")
print(f"  - V(z) should decay from Φs (at z=0) to 0 (at bulk)")
print(f"\nFor triangular model:")
print(f"  V(z) = Φs × (1 - z/W)")
print(f"  → V(0) = Φs, V(W) = 0  ✓")
print(f"\nFor M2 current code:")
print(f"  V(z) = Es × z")
print(f"  → V(0) = 0, V(W) = Es×W = Φs  ✗ (WRONG!)")
print(f"\nCorrect M2 formula should be:")
print(f"  V(z) = Φs - Es×z = Φs × (1 - z/W_eff)")

print("\n" + "="*70)
print("FIX REQUIRED")
print("="*70)
print("""
In models/fang_howard.py, line 203:

  CURRENT (WRONG):
    V_z = Q * Es * z_array
    V_z[z_array > W_eff] = 0

  SHOULD BE:
    Phi_s_V = Phi_s_eV
    V_z = Q * (Phi_s_V - Es * z_array)
    V_z[z_array > W_eff] = 0

  OR equivalently:
    V_z = Q * Phi_s_V * (1 - z_array / W_eff)
    V_z[z_array > W_eff] = 0
""")

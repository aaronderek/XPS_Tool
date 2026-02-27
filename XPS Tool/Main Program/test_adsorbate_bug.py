#!/usr/bin/env python3
"""
Test script to verify adsorbate effect calculation in Figure 2.
This script reproduces the calculation from app.py to check if the
"with adsorbates" curve is horizontal or properly shifted.
"""

import numpy as np
import sys
sys.path.insert(0, '/home/user/2DEG-S-P-toy')

from models.triangular import TriangularModel
from physics.xps import calculate_dipole_from_coverage

# Test parameters from bug report
W_nm = 3.0
epsilon_r = 9
m_star_ratio = 0.32
Phi_s = 0.40

# Adsorbate parameters
coverage = 0.50
mu_debye = 1.50
N_site_exp = 14.00  # Slider value (10^x scale)
N_site = 10 ** N_site_exp  # Convert to actual value: 10^14.00 = 1.0e14 cm⁻²

print("=" * 70)
print("Testing Adsorbate Effect on Figure 2")
print("=" * 70)
print(f"\nParameters:")
print(f"  W = {W_nm:.2f} nm")
print(f"  ε_r = {epsilon_r}")
print(f"  Φₛ = {Phi_s:.2f} eV")
print(f"  θ = {coverage:.2f}")
print(f"  μ⊥ = {mu_debye:.2f} Debye")
print(f"  N_site_exp = {N_site_exp:.2f} (slider value)")
print(f"  N_site = 10^{N_site_exp:.2f} = {N_site:.2e} cm⁻²")

# Calculate dipole shift
Delta_Phi_dip = calculate_dipole_from_coverage(coverage, N_site, mu_debye)
print(f"\nCalculated ΔΦ_dip = {Delta_Phi_dip:+.3f} eV")
print(f"Expected ΔΦ_dip ≈ -0.283 eV")

# Initialize model
model = TriangularModel(m_star_ratio, epsilon_r, W_nm)

# Generate ns range (matching updated app.py)
# Typical 2DEG density range: 0.2 to 2.0 × 10¹³ cm⁻²
ns_range = np.linspace(2e16, 2e17, 100)  # m⁻² (0.2 to 2.0 × 10¹³ cm⁻²)

# Calculate Delta_WF for each ns
Delta_WF_array = model.calculate_Delta_WF(ns_range)

# Calculate Delta_WF with adsorbates (matching app.py line 283)
Delta_WF_with_ads = Delta_WF_array + Delta_Phi_dip

print(f"\n" + "=" * 70)
print("Results for Figure 2: ΔWF vs nₛ")
print("=" * 70)

# Convert ns to display units (10¹³ cm⁻²) using the fixed conversion
from physics.units import ns_to_display
ns_display = ns_to_display(ns_range)

# Print first, middle, and last points
indices = [0, 49, 99]
print(f"\n{'Index':<10} {'nₛ (10¹³ cm⁻²)':<20} {'ΔWF (no ads)':<20} {'ΔWF (with ads)':<20}")
print("-" * 70)
for i in indices:
    print(f"{i:<10} {ns_display[i]:<20.4f} {Delta_WF_array[i]:<20.4f} {Delta_WF_with_ads[i]:<20.4f}")

# Calculate slope for "no ads" curve
slope_no_ads = (Delta_WF_array[-1] - Delta_WF_array[0]) / (ns_display[-1] - ns_display[0])
print(f"\nSlope (no ads): {slope_no_ads:.3f} eV/(10¹³ cm⁻²)")
print(f"Expected slope: -0.603 eV/(10¹³ cm⁻²)")

# Calculate slope for "with ads" curve
slope_with_ads = (Delta_WF_with_ads[-1] - Delta_WF_with_ads[0]) / (ns_display[-1] - ns_display[0])
print(f"\nSlope (with ads): {slope_with_ads:.3f} eV/(10¹³ cm⁻²)")
print(f"Expected slope: -0.603 eV/(10¹³ cm⁻²) (same as no ads)")

# Check if slopes are equal
if abs(slope_no_ads - slope_with_ads) < 0.001:
    print(f"\n✅ PASS: Slopes are equal (curves are parallel)")
else:
    print(f"\n❌ FAIL: Slopes differ by {abs(slope_no_ads - slope_with_ads):.3f} eV/(10¹³ cm⁻²)")

# Check if curve is horizontal
if abs(slope_with_ads) < 0.01:
    print(f"❌ FAIL: 'With ads' curve is horizontal (slope ≈ 0)")
else:
    print(f"✅ PASS: 'With ads' curve has correct slope")

# Check vertical shift
vertical_shift = Delta_WF_with_ads[0] - Delta_WF_array[0]
print(f"\nVertical shift: {vertical_shift:+.3f} eV")
print(f"Expected shift: {Delta_Phi_dip:+.3f} eV (= ΔΦ_dip)")
if abs(vertical_shift - Delta_Phi_dip) < 0.001:
    print(f"✅ PASS: Vertical shift matches ΔΦ_dip")
else:
    print(f"❌ FAIL: Vertical shift does not match ΔΦ_dip")

print("\n" + "=" * 70)
print("Summary")
print("=" * 70)
print("Expected behavior:")
print("  1. ΔΦ_dip should be approximately -0.283 eV (constant)")
print("  2. Both curves should have slope -0.603 eV/(10¹³ cm⁻²)")
print("  3. 'With ads' curve should be shifted down by ΔΦ_dip")
print("  4. 'With ads' curve should NOT be horizontal")

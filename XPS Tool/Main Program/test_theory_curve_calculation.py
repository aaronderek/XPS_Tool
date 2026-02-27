"""
Test the actual theory curve calculation as done in create_publication_comparison_figure.

This reproduces the exact calculation flow to diagnose why η_theory = 0.571
instead of expected ~0.85.
"""

import numpy as np
from models.parabolic import ParabolicModel
from physics.xps import calculate_core_level_shift
from physics.units import nm_to_m
from physics.constants import M0
from scipy.stats import linregress

# Sample data parameters
W_nm = 3.0
lambda_nm = 1.8
theta_deg = 0.0
epsilon_r = 9.0
m_star_ratio = 0.32

# Create model
model = ParabolicModel(m_star_ratio * M0, epsilon_r, W_nm)

# Simulate sample data range (from experiment_data.py)
# Sample data has Delta_WF ranging from approximately -0.35 to 0
Delta_WF_min = -0.35 - 0.05
Delta_WF_max = 0.0 + 0.05

print("="*70)
print("THEORY CURVE CALCULATION TEST")
print("="*70)

print(f"\nParameters:")
print(f"  W = {W_nm} nm")
print(f"  λ = {lambda_nm} nm")
print(f"  θ = {theta_deg}°")
print(f"  Model: M3-Parabolic")
print(f"\nΔWF range: [{Delta_WF_min:.3f}, {Delta_WF_max:.3f}] eV")

# Generate theoretical curve (EXACT CODE FROM create_publication_comparison_figure)
Delta_WF_theory = np.linspace(Delta_WF_min, Delta_WF_max, 100)
Delta_CL_theory = []

lambda_m = nm_to_m(lambda_nm)

for i, Delta_WF in enumerate(Delta_WF_theory):
    # From ΔWF to Phi_s: ΔWF = -Phi_s (ignoring adsorbates)
    Phi_s_eV = -Delta_WF

    # Skip if Phi_s is too small or negative
    if Phi_s_eV < 0.01:
        Delta_CL_theory.append(0)
        if i < 5 or i >= len(Delta_WF_theory) - 5:
            print(f"  Point {i}: ΔWF={Delta_WF:.4f}, Φs={Phi_s_eV:.4f} → SKIPPED")
        continue

    # Create z array for integration
    z_max = 3 * lambda_m
    z_array = np.linspace(0, z_max, 500)

    # Get potential profile from model
    V_z = model.get_potential(Phi_s_eV, z_array)

    # Calculate XPS shift
    Delta_E_CL, _ = calculate_core_level_shift(
        V_z, z_array, lambda_nm, theta_deg
    )

    Delta_CL_theory.append(Delta_E_CL)

    # Print first, middle, and last few points
    if i < 3 or i == len(Delta_WF_theory)//2 or i >= len(Delta_WF_theory) - 3:
        print(f"  Point {i}: ΔWF={Delta_WF:.4f}, Φs={Phi_s_eV:.4f}, ΔE_CL={Delta_E_CL:.4f}")

Delta_CL_theory = np.array(Delta_CL_theory)

# Calculate η from slope (EXACT CODE FROM create_publication_comparison_figure)
if len(Delta_WF_theory) > 2:
    eta_theory = np.polyfit(Delta_WF_theory, Delta_CL_theory, 1)[0]

print(f"\n### RESULT ###")
print(f"η_theory (from slope) = {eta_theory:.4f}")

# Also calculate using scipy linregress for more info
slope, intercept, r_value, p_value, std_err = linregress(Delta_WF_theory, Delta_CL_theory)
print(f"Linear fit: ΔE_CL = {slope:.4f} × ΔWF + {intercept:.4f}")
print(f"R² = {r_value**2:.6f}")
print(f"Intercept = {intercept:.6f} eV (should be near 0)")

# Calculate η at a single point for comparison
Phi_s_test = 0.3
z_array_test = np.linspace(0, 3 * lambda_m, 500)
V_z_test = model.get_potential(Phi_s_test, z_array_test)
Delta_E_CL_test, eta_single = calculate_core_level_shift(
    V_z_test, z_array_test, lambda_nm, theta_deg
)

print(f"\n### SINGLE POINT COMPARISON ###")
print(f"At Φs = {Phi_s_test} eV:")
print(f"  ΔE_CL = {Delta_E_CL_test:.4f} eV")
print(f"  η_single = |ΔE_CL| / Φs = {eta_single:.4f}")
print(f"  η_theory (from slope) = {eta_theory:.4f}")

# Check if the issue is the skipped points near zero
num_skipped = np.sum(np.array(Delta_CL_theory) == 0)
print(f"\n### DIAGNOSTICS ###")
print(f"Number of skipped points (Φs < 0.01): {num_skipped}")
print(f"Skipped percentage: {num_skipped/len(Delta_WF_theory)*100:.1f}%")

if num_skipped > len(Delta_WF_theory) * 0.1:
    print("⚠️ WARNING: Many points skipped due to Phi_s < 0.01 cutoff")
    print("   This could affect the slope calculation")

# Calculate simplified η
lambda_eff = lambda_m * np.cos(0)  # theta = 0
eta_simple = 1 - np.exp(-nm_to_m(W_nm) / lambda_eff)

print(f"\n### COMPARISON ###")
print(f"η_simple (uniform potential) = {eta_simple:.4f}")
print(f"η_single (parabolic, Φs=0.3) = {eta_single:.4f}")
print(f"η_theory (from slope)        = {eta_theory:.4f}")
print(f"η_expected (sample data)     = 0.8500")

print(f"\n### DIAGNOSIS ###")
if abs(eta_theory - 0.85) / 0.85 < 0.1:
    print("✓ η_theory is correct (within 10% of expected)")
elif abs(eta_theory - eta_single) / eta_single < 0.05:
    print("⚠️ η_theory matches single-point calculation")
    print("   BUT both are lower than expected for sample data")
    print("   → This is NORMAL for parabolic potential!")
    print("   → Sample data was generated with η=0.85, not physical η")
else:
    print("✗ η_theory has unexpected value")
    print("   → Check for bugs in the calculation")

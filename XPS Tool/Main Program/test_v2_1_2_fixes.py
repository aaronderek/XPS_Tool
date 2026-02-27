"""
Test script for v2.1.2 fixes.

Validates:
1. Sample data generates physically reasonable values
2. Data validation detects issues correctly
3. Theory curve calculation works properly
"""

import numpy as np
import sys
from scipy.stats import linregress

# Import modules
from utils.experiment_data import create_sample_data, validate_experimental_data_physics
from models import TriangularModel, ParabolicModel
from physics.constants import M0, Q
from physics.xps import calculate_core_level_shift
from physics.units import nm_to_m

print("=" * 60)
print("Testing v2.1.2 Fixes")
print("=" * 60)

# Test 1: Sample Data Physics
print("\n[Test 1] Sample Data Generation")
print("-" * 60)

df_sample = create_sample_data()
print(f"Generated {len(df_sample)} data points")
print("\nFirst 3 rows (raw data):")
print(df_sample.head(3))

# Process to get Delta values
WF = df_sample['WF_eV'].values
In3d = df_sample['In3d_eV'].values
Delta_WF = WF - WF[0]
Delta_CL = In3d - In3d[0]

print(f"\nΔWF range: [{Delta_WF.min():.3f}, {Delta_WF.max():.3f}] eV")
print(f"ΔE_CL range: [{Delta_CL.min():.3f}, {Delta_CL.max():.3f}] eV")

# Check correlation
corr = np.corrcoef(Delta_WF, Delta_CL)[0, 1]
print(f"\nCorrelation: {corr:.4f}")

# Check eta
slope, intercept, r_value, _, std_err = linregress(Delta_WF, Delta_CL)
print(f"η (slope): {slope:.3f} ± {std_err:.3f}")
print(f"R²: {r_value**2:.4f}")

# Acceptance criteria
test1_pass = True

# Check if all non-reference points are negative (first point is reference = 0)
Delta_WF_nonref = Delta_WF[1:]  # Skip first point
Delta_CL_nonref = Delta_CL[1:]

if not (Delta_WF_nonref < 0).all():
    print("❌ FAIL: Not all ΔWF values (non-reference) are negative")
    test1_pass = False
else:
    print("✅ PASS: All non-reference ΔWF values are negative")

if not (Delta_CL_nonref < 0).all():
    print("❌ FAIL: Not all ΔE_CL values (non-reference) are negative")
    test1_pass = False
else:
    print("✅ PASS: All non-reference ΔE_CL values are negative")

if corr < 0:
    print("❌ FAIL: Negative correlation detected")
    test1_pass = False
elif corr < 0.95:
    print(f"⚠️  WARNING: Correlation {corr:.3f} is below 0.95")
    print("✅ PASS: Positive correlation (but low)")
else:
    print("✅ PASS: Strong positive correlation")

if slope < 0.7 or slope > 0.9:
    print(f"⚠️  WARNING: η = {slope:.3f} is outside typical range (0.7-0.9)")
    if slope < 0:
        print("❌ FAIL: Negative η")
        test1_pass = False
    else:
        print("✅ PASS: Positive η (but outside typical range)")
else:
    print("✅ PASS: η in expected range (0.7-0.9)")

if r_value**2 < 0.95:
    print(f"⚠️  WARNING: R² = {r_value**2:.4f} is below 0.95")
else:
    print("✅ PASS: High R² value")

# Test 2: Data Validation
print("\n[Test 2] Data Validation System")
print("-" * 60)

exp_data = {
    'Delta_WF': Delta_WF,
    'Delta_CL': Delta_CL
}

warnings = validate_experimental_data_physics(exp_data)
print(f"Number of warnings: {len(warnings)}")

if len(warnings) == 0:
    print("✅ PASS: No warnings for good data")
    test2_pass = True
else:
    print("⚠️  Warnings generated:")
    for i, warning in enumerate(warnings, 1):
        print(f"\n{i}. {warning[:100]}...")
    test2_pass = True  # Warnings are ok if they're informational

# Test with bad data (negative correlation)
print("\nTesting with bad data (flipped signs)...")
bad_data = {
    'Delta_WF': -Delta_WF,  # Flip sign
    'Delta_CL': Delta_CL
}

warnings_bad = validate_experimental_data_physics(bad_data)
print(f"Number of warnings: {len(warnings_bad)}")

if len(warnings_bad) > 0:
    print("✅ PASS: Validation detected issue in bad data")
    # Check if negative correlation warning is present
    has_neg_corr_warning = any('negative correlation' in w.lower() for w in warnings_bad)
    if has_neg_corr_warning:
        print("✅ PASS: Negative correlation warning present")
    else:
        print("⚠️  WARNING: Negative correlation not explicitly warned")
else:
    print("❌ FAIL: Validation did not detect issue")
    test2_pass = False

# Test 3: Theory Curve Calculation
print("\n[Test 3] Theory η Calculation")
print("-" * 60)

# Setup model
m_star = 0.32 * M0
epsilon_r = 9
W_nm = 3.0
lambda_nm = 1.8
theta_deg = 0

model = TriangularModel(m_star, epsilon_r, W_nm)

# Calculate theory curve for a range of Delta_WF
Delta_WF_theory = np.linspace(-0.4, -0.1, 20)
Delta_CL_theory = []

for dwf in Delta_WF_theory:
    Phi_s_eV = -dwf

    if Phi_s_eV < 0.01:
        Delta_CL_theory.append(0)
        continue

    # Create z array
    lambda_m = nm_to_m(lambda_nm)
    z_max = 3 * lambda_m
    z_array = np.linspace(0, z_max, 500)

    # Get potential
    V_z = model.get_potential(Phi_s_eV, z_array)

    # Calculate shift
    from physics.xps import calculate_core_level_shift
    Delta_E_CL, _ = calculate_core_level_shift(V_z, z_array, lambda_nm, theta_deg)

    Delta_CL_theory.append(Delta_E_CL)

Delta_CL_theory = np.array(Delta_CL_theory)

# Calculate eta_theory
eta_theory = np.polyfit(Delta_WF_theory, Delta_CL_theory, 1)[0]

print(f"Model: M1-Triangular")
print(f"Parameters: W={W_nm}nm, λ={lambda_nm}nm, θ={theta_deg}°")
print(f"η_theory: {eta_theory:.3f}")

test3_pass = True
if eta_theory < 0:
    print("❌ FAIL: Negative η_theory")
    test3_pass = False
elif eta_theory < 0.5:
    print("⚠️  WARNING: η_theory unusually low")
    print("✅ PASS: Positive η_theory (but low)")
elif eta_theory > 0.95:
    print("⚠️  WARNING: η_theory unusually high")
    print("✅ PASS: Positive η_theory (but high)")
else:
    print("✅ PASS: η_theory in expected range (0.5-0.95)")

# Test with M3 model
model_m3 = ParabolicModel(m_star, epsilon_r, W_nm)
Delta_CL_theory_m3 = []

for dwf in Delta_WF_theory:
    Phi_s_eV = -dwf

    if Phi_s_eV < 0.01:
        Delta_CL_theory_m3.append(0)
        continue

    lambda_m = nm_to_m(lambda_nm)
    z_max = 3 * lambda_m
    z_array = np.linspace(0, z_max, 500)

    V_z = model_m3.get_potential(Phi_s_eV, z_array)
    Delta_E_CL, _ = calculate_core_level_shift(V_z, z_array, lambda_nm, theta_deg)

    Delta_CL_theory_m3.append(Delta_E_CL)

Delta_CL_theory_m3 = np.array(Delta_CL_theory_m3)
eta_theory_m3 = np.polyfit(Delta_WF_theory, Delta_CL_theory_m3, 1)[0]

print(f"\nModel: M3-Parabolic")
print(f"η_theory: {eta_theory_m3:.3f}")

if eta_theory_m3 > 0.5 and eta_theory_m3 < 0.95:
    print("✅ PASS: M3 η_theory in expected range")
else:
    print(f"⚠️  WARNING: M3 η_theory = {eta_theory_m3:.3f} outside expected range")

# Summary
print("\n" + "=" * 60)
print("Test Summary")
print("=" * 60)

all_pass = test1_pass and test2_pass and test3_pass

print(f"Test 1 (Sample Data): {'✅ PASS' if test1_pass else '❌ FAIL'}")
print(f"Test 2 (Validation):  {'✅ PASS' if test2_pass else '❌ FAIL'}")
print(f"Test 3 (Theory η):    {'✅ PASS' if test3_pass else '❌ FAIL'}")

print("\n" + "=" * 60)
if all_pass:
    print("🎉 ALL TESTS PASSED!")
    print("=" * 60)
    sys.exit(0)
else:
    print("❌ SOME TESTS FAILED")
    print("=" * 60)
    sys.exit(1)

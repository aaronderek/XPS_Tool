"""
Debug script to test theoretical η calculation.

Based on the issue description:
- W = 3 nm
- λ = 1.8 nm
- θ = 0°
- Expected η_simple = 1 - exp(-W/λ) ≈ 0.811
- Current η_theory = 0.571 (WRONG)

This script will:
1. Calculate simplified η
2. Calculate full XPS integration η
3. Find where the bug is
"""

import numpy as np
from models.parabolic import ParabolicModel
from models.triangular import TriangularModel
from models.fang_howard import FangHowardModel
from physics.xps import calculate_core_level_shift, calculate_xps_weight
from physics.units import nm_to_m, degrees_to_radians
from physics.constants import M0

# Test parameters (same as sample data)
W_nm = 3.0
lambda_nm = 1.8
theta_deg = 0.0
epsilon_r = 9.0
m_star_ratio = 0.32
Phi_s_eV = 0.3  # Test surface potential

print("="*70)
print("THEORETICAL η CALCULATION DEBUG")
print("="*70)

print("\n### INPUT PARAMETERS ###")
print(f"W = {W_nm} nm")
print(f"λ = {lambda_nm} nm")
print(f"θ = {theta_deg}°")
print(f"Φs = {Phi_s_eV} eV (test value)")
print(f"εᵣ = {epsilon_r}")
print(f"m*/m₀ = {m_star_ratio}")

# Calculate simplified η
print("\n### SIMPLIFIED FORMULA ###")
lambda_m = nm_to_m(lambda_nm)
W_m = nm_to_m(W_nm)
theta_rad = degrees_to_radians(theta_deg)
lambda_eff = lambda_m * np.cos(theta_rad)

print(f"λ = {lambda_m} m = {lambda_m*1e9:.3f} nm")
print(f"W = {W_m} m = {W_m*1e9:.3f} nm")
print(f"θ = {theta_rad} rad = {theta_rad*180/np.pi:.1f}°")
print(f"λ_eff = λ × cos(θ) = {lambda_eff*1e9:.3f} nm")

eta_simple = 1 - np.exp(-W_m / lambda_eff)
print(f"\nη_simple = 1 - exp(-W/λ_eff)")
print(f"         = 1 - exp(-{W_nm}/{lambda_eff*1e9:.3f})")
print(f"         = 1 - exp({-W_m/lambda_eff:.4f})")
print(f"         = 1 - {np.exp(-W_m/lambda_eff):.4f}")
print(f"         = {eta_simple:.4f}")

# Test all three models
models = {
    'M1-Triangular': TriangularModel(m_star_ratio * M0, epsilon_r, W_nm),
    'M2-Fang-Howard': FangHowardModel(m_star_ratio * M0, epsilon_r, W_nm),
    'M3-Parabolic': ParabolicModel(m_star_ratio * M0, epsilon_r, W_nm)
}

print("\n### XPS INTEGRATION CALCULATION ###")
for model_name, model in models.items():
    print(f"\n--- {model_name} ---")

    # Create z array for integration
    z_max = 3 * lambda_m
    z_array = np.linspace(0, z_max, 500)

    print(f"z_max = {z_max*1e9:.3f} nm (3 × λ)")
    print(f"Integration points: {len(z_array)}")

    # Get potential profile
    V_z = model.get_potential(Phi_s_eV, z_array)

    # Check units
    print(f"\nPotential profile V(z):")
    print(f"  V(z=0) = {V_z[0]:.6e} J")
    print(f"  V(z=W) = {V_z[np.argmin(np.abs(z_array - W_m))]:.6e} J")
    print(f"  V(z=3λ) = {V_z[-1]:.6e} J")

    # Get weight function
    w_z = calculate_xps_weight(z_array, lambda_nm, theta_deg)

    print(f"\nWeight function w(z):")
    print(f"  w(z=0) = {w_z[0]:.6e} m⁻¹")
    print(f"  ∫w(z)dz = {np.trapezoid(w_z, z_array):.6f} (should ≈ 1)")

    # Calculate XPS shift
    Delta_E_CL, eta_xps = calculate_core_level_shift(
        V_z, z_array, lambda_nm, theta_deg
    )

    print(f"\nXPS core level shift:")
    print(f"  ΔE_CL = {Delta_E_CL:.6f} eV")
    print(f"  η_xps = ΔE_CL / Φs = {Delta_E_CL:.6f} / {Phi_s_eV:.6f} = {eta_xps:.4f}")

    # Compare with simplified
    print(f"\nComparison:")
    print(f"  η_simple = {eta_simple:.4f}")
    print(f"  η_xps    = {eta_xps:.4f}")
    print(f"  Difference = {abs(eta_simple - eta_xps)/eta_simple*100:.1f}%")

    # Expected: difference should be < 10% for reasonable models
    if abs(eta_simple - eta_xps) / eta_simple > 0.1:
        print(f"  ⚠️ WARNING: Large difference detected!")
        print(f"     This suggests a bug in the XPS integration.")
    else:
        print(f"  ✓ Good agreement")

print("\n" + "="*70)
print("DIAGNOSIS")
print("="*70)

print(f"""
If η_xps ≈ {eta_simple:.4f} for all models:
  ✓ XPS integration is working correctly
  → Bug must be in how create_publication_comparison_figure calls this code

If η_xps ≈ 0.571 (the reported wrong value):
  ✗ Bug is in calculate_core_level_shift or weight function
  → Check sign conventions, unit conversions, integration formula

If η_xps varies wildly between models:
  ✗ Bug might be in model.get_potential
  → Check potential energy sign and magnitude
""")

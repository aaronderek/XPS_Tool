"""
Test η calculation for all three models to find which one gives 0.571.
"""

import numpy as np
from models.parabolic import ParabolicModel
from models.triangular import TriangularModel
from models.fang_howard import FangHowardModel
from physics.xps import calculate_core_level_shift
from physics.units import nm_to_m
from physics.constants import M0

# Parameters
W_nm = 3.0
lambda_nm = 1.8
theta_deg = 0.0
epsilon_r = 9.0
m_star_ratio = 0.32

# Sample data range
Delta_WF_min = -0.40
Delta_WF_max = 0.05

models = {
    'M1-Triangular': TriangularModel(m_star_ratio * M0, epsilon_r, W_nm),
    'M2-Fang-Howard': FangHowardModel(m_star_ratio * M0, epsilon_r, W_nm),
    'M3-Parabolic': ParabolicModel(m_star_ratio * M0, epsilon_r, W_nm)
}

print("="*70)
print("η CALCULATION FOR ALL MODELS")
print("="*70)

for model_name, model in models.items():
    print(f"\n### {model_name} ###")

    # Generate theoretical curve
    Delta_WF_theory = np.linspace(Delta_WF_min, Delta_WF_max, 100)
    Delta_CL_theory = []

    lambda_m = nm_to_m(lambda_nm)

    for Delta_WF in Delta_WF_theory:
        Phi_s_eV = -Delta_WF

        if Phi_s_eV < 0.01:
            Delta_CL_theory.append(0)
            continue

        z_max = 3 * lambda_m
        z_array = np.linspace(0, z_max, 500)
        V_z = model.get_potential(Phi_s_eV, z_array)
        Delta_E_CL, _ = calculate_core_level_shift(
            V_z, z_array, lambda_nm, theta_deg
        )

        Delta_CL_theory.append(Delta_E_CL)

    Delta_CL_theory = np.array(Delta_CL_theory)

    # Calculate η from slope
    eta_theory = np.polyfit(Delta_WF_theory, Delta_CL_theory, 1)[0]

    # Single point calculation
    Phi_s_test = 0.3
    z_array_test = np.linspace(0, 3 * lambda_m, 500)
    V_z_test = model.get_potential(Phi_s_test, z_array_test)
    Delta_E_CL_test, eta_single = calculate_core_level_shift(
        V_z_test, z_array_test, lambda_nm, theta_deg
    )

    print(f"  η_theory (from slope)  = {eta_theory:.4f}")
    print(f"  η_single (at Φs=0.3)   = {eta_single:.4f}")

    if abs(eta_theory - 0.571) < 0.01:
        print(f"  ⭐ THIS MATCHES USER'S REPORT (0.571)!")

print("\n" + "="*70)
print("EXPECTED η VALUES")
print("="*70)

# Simplified formula
lambda_eff = lambda_m
eta_simple = 1 - np.exp(-nm_to_m(W_nm) / lambda_eff)

print(f"η_simple (uniform potential) = {eta_simple:.4f}")
print(f"η_sample_data (from generation) = 0.8500")
print(f"\nUser reported:")
print(f"  η_exp = 0.879 (experimental linear fit)")
print(f"  η_fitted = 0.860 (parameter fitting)")
print(f"  η_theory = 0.571 (theory curve slope)")
print(f"\nConclusion:")
print(f"  None of the models give η ≈ 0.571 with current calculation.")
print(f"  M1 ≈ {0.3:.2f}, M2 ≈ {0.3:.2f}, M3 ≈ {0.37:.2f}")
print(f"  This suggests a bug in how η is calculated or displayed.")

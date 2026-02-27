"""
Test script to verify the Parameter Fitting fix.
"""

import numpy as np
from utils.experiment_data import create_sample_data, process_raw_data
from utils.fitting import run_parameter_fitting, linear_fit_eta, calculate_theory_Delta_CL

print("=" * 80)
print("Testing Parameter Fitting Fix")
print("=" * 80)

# Create sample data
print("\n1. Creating sample data...")
df = create_sample_data()
exp_data = process_raw_data(df)

print(f"   Sample data points: {len(exp_data['Delta_WF'])}")
print(f"   Delta_WF range: [{exp_data['Delta_WF'].min():.3f}, {exp_data['Delta_WF'].max():.3f}] eV")
print(f"   Delta_CL range: [{exp_data['Delta_CL'].min():.3f}, {exp_data['Delta_CL'].max():.3f}] eV")

# Check if both are negative (vacuum annealing)
print(f"   All Delta_WF negative? {np.all(exp_data['Delta_WF'] <= 0)}")
print(f"   All Delta_CL negative? {np.all(exp_data['Delta_CL'] <= 0)}")

# Test linear fit
print("\n2. Testing linear fit (extracting η)...")
lin_fit = linear_fit_eta(exp_data)
print(f"   Experimental η = {lin_fit['eta_exp']:.3f}")
print(f"   Slope = {lin_fit['slope']:.3f}")
print(f"   R² = {lin_fit['r_squared']:.4f}")
print(f"   Expected η ≈ 0.85 ✓" if 0.80 <= lin_fit['eta_exp'] <= 0.90 else "   ⚠️ η out of expected range!")

# Test calculate_theory_Delta_CL
print("\n3. Testing calculate_theory_Delta_CL...")
eta_test = 0.85
theory_CL = calculate_theory_Delta_CL(
    exp_data['Delta_WF'],
    W_nm=3.0,
    eta=eta_test,
    model_type='M2-Fang-Howard'
)

print(f"   Theory Delta_CL range: [{theory_CL.min():.3f}, {theory_CL.max():.3f}] eV")
print(f"   Exp Delta_CL range: [{exp_data['Delta_CL'].min():.3f}, {exp_data['Delta_CL'].max():.3f}] eV")

# Check if signs match
same_sign = (np.sign(theory_CL) == np.sign(exp_data['Delta_CL'])).all()
print(f"   Signs match? {same_sign} {'✓' if same_sign else '❌'}")

# Check correlation
correlation = np.corrcoef(exp_data['Delta_WF'], theory_CL)[0, 1]
print(f"   Correlation between theory and exp: {correlation:.4f}")
print(f"   {'✓ Positive correlation (correct)' if correlation > 0 else '❌ Negative correlation (wrong)'}")

# Test full parameter fitting
print("\n4. Testing full parameter fitting...")
fit_result = run_parameter_fitting(
    exp_data=exp_data,
    model_type='M2-Fang-Howard',
    fit_params={'W': False, 'eta': True},
    initial_guess={'W': 3.0, 'eta': 0.85},
    bounds={'W': (2.0, 5.0), 'eta': (0.7, 0.95)},
    epsilon_r=9.0,
    m_star_ratio=0.32
)

if fit_result['success']:
    print(f"   ✓ Fitting succeeded")
    print(f"   Fitted η = {fit_result['eta_fit']:.3f}")
    print(f"   R² = {fit_result['r_squared']:.4f}")
    print(f"   RMSE = {fit_result['rmse']:.4f} eV")

    # Check if fitted curve matches experiment trend
    fitted_CL = fit_result['theory_Delta_CL']
    fit_correlation = np.corrcoef(exp_data['Delta_CL'], fitted_CL)[0, 1]
    print(f"   Correlation with experiment: {fit_correlation:.4f}")
    print(f"   {'✓ Good fit' if fit_correlation > 0.95 else '⚠️ Poor fit'}")

    # Check direction
    exp_increasing = exp_data['Delta_CL'][-1] > exp_data['Delta_CL'][0]
    fit_increasing = fitted_CL[-1] > fitted_CL[0]
    print(f"   Exp trend: {'increasing' if exp_increasing else 'decreasing'}")
    print(f"   Fit trend: {'increasing' if fit_increasing else 'decreasing'}")
    print(f"   {'✓ Trends match' if exp_increasing == fit_increasing else '❌ Trends opposite'}")
else:
    print(f"   ❌ Fitting failed: {fit_result['message']}")

print("\n" + "=" * 80)
print("Test Complete!")
print("=" * 80)

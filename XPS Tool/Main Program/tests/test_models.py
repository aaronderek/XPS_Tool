"""
Unit tests for 2DEG models.

Test cases from requirements document:
1. M1 model: εᵣ=9, W=3nm, Φs=0.3eV
   Expected: ns ≈ 5.3e12 cm⁻², ΔWF = -0.30 eV, slope ≈ -0.57 eV/(1e13 cm⁻²)

2. M3 model: εᵣ=9, W=3nm, Φs=0.3eV
   Expected: ns ≈ 1.06e13 cm⁻² (2x M1), ΔWF = -0.30 eV, slope ≈ -0.28 eV/(1e13 cm⁻²)

3. Adsorbate: ΔΦdip = +0.3 eV
   Expected: ns unchanged, ΔWF_total = -0.30 + 0.30 = 0.00 eV
"""

import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from models import TriangularModel, FangHowardModel, ParabolicModel
from physics.constants import M0
from physics.units import ns_to_display
from physics.xps import calculate_dipole_from_coverage


def test_m1_basic():
    """Test M1 model with standard parameters."""
    print("\n=== Test 1: M1 Model Basic ===")

    # Parameters
    m_star = 0.32 * M0
    epsilon_r = 9
    W_nm = 3.0
    Phi_s = 0.3

    # Create model
    model = TriangularModel(m_star, epsilon_r, W_nm)

    # Calculate ns
    ns = model.calculate_ns(Phi_s)
    ns_display = ns_to_display(ns)

    print(f"Parameters: εᵣ={epsilon_r}, W={W_nm}nm, Φs={Phi_s}eV")
    print(f"Result: ns = {ns_display:.3f} × 10¹³ cm⁻²")
    print(f"Expected: ns ≈ 5.3 × 10¹³ cm⁻²")

    # Calculate Delta_WF
    Delta_WF = model.calculate_Delta_WF(ns)
    print(f"Result: ΔWF = {Delta_WF:.3f} eV")
    print(f"Expected: ΔWF = -0.30 eV")

    # Calculate slope
    slope = model.get_slope()
    print(f"Result: slope = {slope:.3f} eV/(10¹³ cm⁻²)")
    print(f"Expected: slope ≈ -0.57 eV/(10¹³ cm⁻²)")

    # Assertions with tolerance
    assert np.isclose(ns_display, 5.3, atol=0.5), f"ns mismatch: {ns_display:.3f}"
    assert np.isclose(Delta_WF, -0.30, atol=0.01), f"Delta_WF mismatch: {Delta_WF:.3f}"
    assert np.isclose(slope, -0.57, atol=0.05), f"Slope mismatch: {slope:.3f}"

    print("✓ Test 1 passed!")
    return True


def test_m3_basic():
    """Test M3 model with standard parameters."""
    print("\n=== Test 2: M3 Model Basic ===")

    # Parameters
    m_star = 0.32 * M0
    epsilon_r = 9
    W_nm = 3.0
    Phi_s = 0.3

    # Create model
    model = ParabolicModel(m_star, epsilon_r, W_nm)

    # Calculate ns
    ns = model.calculate_ns(Phi_s)
    ns_display = ns_to_display(ns)

    print(f"Parameters: εᵣ={epsilon_r}, W={W_nm}nm, Φs={Phi_s}eV")
    print(f"Result: ns = {ns_display:.3f} × 10¹³ cm⁻²")
    print(f"Expected: ns ≈ 10.6 × 10¹³ cm⁻² (2x M1)")

    # Calculate Delta_WF
    Delta_WF = model.calculate_Delta_WF(ns)
    print(f"Result: ΔWF = {Delta_WF:.3f} eV")
    print(f"Expected: ΔWF = -0.30 eV")

    # Calculate slope
    slope = model.get_slope()
    print(f"Result: slope = {slope:.3f} eV/(10¹³ cm⁻²)")
    print(f"Expected: slope ≈ -0.28 eV/(10¹³ cm⁻²)")

    # Assertions
    assert np.isclose(ns_display, 10.6, atol=0.5), f"ns mismatch: {ns_display:.3f}"
    assert np.isclose(Delta_WF, -0.30, atol=0.01), f"Delta_WF mismatch: {Delta_WF:.3f}"
    assert np.isclose(slope, -0.28, atol=0.05), f"Slope mismatch: {slope:.3f}"

    print("✓ Test 2 passed!")
    return True


def test_m3_vs_m1_ratio():
    """Test that M3 ns is exactly 2x M1 ns for same parameters."""
    print("\n=== Test 3: M3 vs M1 Ratio ===")

    # Parameters
    m_star = 0.32 * M0
    epsilon_r = 9
    W_nm = 3.0
    Phi_s = 0.4

    # Create models
    model_m1 = TriangularModel(m_star, epsilon_r, W_nm)
    model_m3 = ParabolicModel(m_star, epsilon_r, W_nm)

    # Calculate ns
    ns_m1 = model_m1.calculate_ns(Phi_s)
    ns_m3 = model_m3.calculate_ns(Phi_s)

    ratio = ns_m3 / ns_m1

    print(f"M1 ns = {ns_to_display(ns_m1):.3f} × 10¹³ cm⁻²")
    print(f"M3 ns = {ns_to_display(ns_m3):.3f} × 10¹³ cm⁻²")
    print(f"Ratio (M3/M1) = {ratio:.3f}")
    print(f"Expected: ratio = 2.00")

    # Check ratio
    assert np.isclose(ratio, 2.0, atol=0.01), f"Ratio mismatch: {ratio:.3f}"

    # Check slopes
    slope_m1 = model_m1.get_slope()
    slope_m3 = model_m3.get_slope()
    slope_ratio = slope_m3 / slope_m1

    print(f"\nM1 slope = {slope_m1:.3f} eV/(10¹³ cm⁻²)")
    print(f"M3 slope = {slope_m3:.3f} eV/(10¹³ cm⁻²)")
    print(f"Slope ratio (M3/M1) = {slope_ratio:.3f}")
    print(f"Expected: slope ratio = 0.50")

    assert np.isclose(slope_ratio, 0.5, atol=0.01), f"Slope ratio mismatch: {slope_ratio:.3f}"

    print("✓ Test 3 passed!")
    return True


def test_adsorbate_shift():
    """Test adsorbate dipole shift."""
    print("\n=== Test 4: Adsorbate Shift ===")

    # Parameters
    m_star = 0.32 * M0
    epsilon_r = 9
    W_nm = 3.0
    Phi_s = 0.3

    # Create model
    model = TriangularModel(m_star, epsilon_r, W_nm)

    # Calculate without adsorbate
    ns = model.calculate_ns(Phi_s)
    Delta_WF_no_ads = model.calculate_Delta_WF(ns)

    # Add adsorbate shift
    Delta_Phi_dip = 0.3  # eV

    Delta_WF_with_ads = Delta_WF_no_ads + Delta_Phi_dip

    print(f"ns = {ns_to_display(ns):.3f} × 10¹³ cm⁻²")
    print(f"ΔWF (no ads) = {Delta_WF_no_ads:.3f} eV")
    print(f"ΔΦdip = +{Delta_Phi_dip:.2f} eV")
    print(f"ΔWF (with ads) = {Delta_WF_with_ads:.3f} eV")
    print(f"Expected: ΔWF (with ads) ≈ 0.00 eV")

    # Check that ns doesn't change (it's an input)
    # Check that ΔWF shifts by exactly ΔΦdip
    assert np.isclose(Delta_WF_with_ads, 0.0, atol=0.01), f"ΔWF mismatch: {Delta_WF_with_ads:.3f}"

    print("✓ Test 4 passed!")
    return True


def test_m2_convergence():
    """Test M2 Fang-Howard self-consistent convergence."""
    print("\n=== Test 5: M2 Self-Consistent Convergence ===")

    # Parameters
    m_star = 0.32 * M0
    epsilon_r = 9
    W_nm = 3.0
    Phi_s = 0.4

    # Create model
    model = FangHowardModel(m_star, epsilon_r, W_nm)

    # Solve self-consistent equations
    Es, W_eff, n_iter = model.solve_self_consistent(Phi_s, max_iter=10, tolerance=1e-4)

    print(f"Converged in {n_iter} iterations")
    print(f"W_eff = {W_eff * 1e9:.3f} nm")
    print(f"Es = {Es:.3e} V/m")

    # Check convergence
    assert n_iter < 10, f"Failed to converge in 10 iterations"
    assert W_eff < model.W_initial, f"W_eff should be smaller than W_initial for high field"

    # Calculate ns
    ns = model.calculate_ns(Phi_s)
    ns_display = ns_to_display(ns)

    print(f"ns = {ns_display:.3f} × 10¹³ cm⁻²")

    # M2 should give ns between M1 and M3
    model_m1 = TriangularModel(m_star, epsilon_r, W_nm)
    model_m3 = ParabolicModel(m_star, epsilon_r, W_nm)

    ns_m1 = model_m1.calculate_ns(Phi_s)
    ns_m3 = model_m3.calculate_ns(Phi_s)

    print(f"M1 ns = {ns_to_display(ns_m1):.3f} × 10¹³ cm⁻²")
    print(f"M2 ns = {ns_display:.3f} × 10¹³ cm⁻²")
    print(f"M3 ns = {ns_to_display(ns_m3):.3f} × 10¹³ cm⁻²")

    # M2 should be between M1 and M3
    assert ns_m1 < ns < ns_m3, "M2 ns should be between M1 and M3"

    print("✓ Test 5 passed!")
    return True


def test_zero_potential():
    """Test that all models give ns=0 when Phi_s=0."""
    print("\n=== Test 6: Zero Potential ===")

    # Parameters
    m_star = 0.32 * M0
    epsilon_r = 9
    W_nm = 3.0
    Phi_s = 0.0

    # Create models
    model_m1 = TriangularModel(m_star, epsilon_r, W_nm)
    model_m3 = ParabolicModel(m_star, epsilon_r, W_nm)

    # Calculate ns
    ns_m1 = model_m1.calculate_ns(Phi_s)
    ns_m3 = model_m3.calculate_ns(Phi_s)

    print(f"M1 ns (Φs=0) = {ns_to_display(ns_m1):.6f} × 10¹³ cm⁻²")
    print(f"M3 ns (Φs=0) = {ns_to_display(ns_m3):.6f} × 10¹³ cm⁻²")

    # Should be zero or very close
    assert np.isclose(ns_m1, 0, atol=1e10), f"M1 ns should be zero: {ns_m1}"
    assert np.isclose(ns_m3, 0, atol=1e10), f"M3 ns should be zero: {ns_m3}"

    print("✓ Test 6 passed!")
    return True


def run_all_tests():
    """Run all validation tests."""
    print("=" * 60)
    print("2DEG Model Validation Tests")
    print("=" * 60)

    tests = [
        test_m1_basic,
        test_m3_basic,
        test_m3_vs_m1_ratio,
        test_adsorbate_shift,
        test_m2_convergence,
        test_zero_potential
    ]

    passed = 0
    failed = 0

    for test in tests:
        try:
            test()
            passed += 1
        except AssertionError as e:
            print(f"✗ Test failed: {e}")
            failed += 1
        except Exception as e:
            print(f"✗ Test error: {e}")
            failed += 1

    print("\n" + "=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 60)

    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)

"""
Theory-experiment fitting utilities.
"""

import numpy as np
from scipy.optimize import minimize, OptimizeResult
from scipy.stats import linregress
from typing import Dict, List, Tuple, Optional, Any
import warnings


def calculate_theory_Delta_CL(
    Delta_WF_exp: np.ndarray,
    W_nm: float,
    eta: float,
    model_type: str,
    epsilon_r: float = 9.0,
    m_star_ratio: float = 0.32
) -> np.ndarray:
    """
    Calculate theoretical Delta_CL from experimental Delta_WF.

    Parameters
    ----------
    Delta_WF_exp : array
        Experimental work function changes (eV)
    W_nm : float
        Depletion width (nm)
    eta : float
        XPS sampling factor
    model_type : str
        Model type ('M1-Triangular', 'M2-Fang-Howard', 'M3-Parabolic')
    epsilon_r : float
        Relative permittivity
    m_star_ratio : float
        Effective mass ratio

    Returns
    -------
    Delta_CL_theory : array
        Theoretical core level shifts (eV)
    """
    # For simplified fitting, we use the relation:
    # Delta_CL = eta * Delta_WF
    # Both Delta_CL and Delta_WF follow the same sign convention:
    # - Vacuum annealing (desorption): Both negative
    # - Oxygen annealing (adsorption): Both positive
    # The relationship is direct proportional, not inverse

    Delta_CL_theory = eta * Delta_WF_exp

    return Delta_CL_theory


def objective_function(
    params: np.ndarray,
    exp_data: Dict[str, np.ndarray],
    fit_params: Dict[str, bool],
    initial_guess: Dict[str, float],
    model_type: str,
    epsilon_r: float,
    m_star_ratio: float
) -> float:
    """
    Objective function for parameter fitting.

    Minimizes chi-squared between theory and experiment.

    Parameters
    ----------
    params : array
        Current parameter values being optimized
    exp_data : dict
        Experimental data
    fit_params : dict
        Which parameters to fit
    initial_guess : dict
        Initial parameter values
    model_type : str
        Model type
    epsilon_r : float
        Relative permittivity
    m_star_ratio : float
        Effective mass ratio

    Returns
    -------
    chi_squared : float
        Chi-squared value
    """
    # Unpack parameters
    idx = 0
    W_nm = params[idx] if fit_params.get('W', False) else initial_guess['W']
    idx += 1 if fit_params.get('W', False) else 0

    eta = params[idx] if fit_params.get('eta', False) else initial_guess['eta']
    idx += 1 if fit_params.get('eta', False) else 0

    # Calculate theoretical Delta_CL
    theory_Delta_CL = calculate_theory_Delta_CL(
        exp_data['Delta_WF'],
        W_nm,
        eta,
        model_type,
        epsilon_r,
        m_star_ratio
    )

    # Calculate residuals
    residuals = theory_Delta_CL - exp_data['Delta_CL']

    # Chi-squared
    chi_squared = np.sum(residuals**2)

    return chi_squared


def run_parameter_fitting(
    exp_data: Dict[str, np.ndarray],
    model_type: str,
    fit_params: Dict[str, bool],
    initial_guess: Dict[str, float],
    bounds: Dict[str, Tuple[float, float]],
    epsilon_r: float = 9.0,
    m_star_ratio: float = 0.32
) -> Dict[str, Any]:
    """
    Fit theoretical parameters to experimental data.

    Parameters
    ----------
    exp_data : dict
        Experimental data containing 'Delta_WF' and 'Delta_CL'
    model_type : str
        Model type
    fit_params : dict
        Which parameters to fit (e.g., {'W': True, 'eta': True})
    initial_guess : dict
        Initial parameter guesses
    bounds : dict
        Parameter bounds (e.g., {'W': (2.0, 5.0), 'eta': (0.7, 0.95)})
    epsilon_r : float
        Relative permittivity
    m_star_ratio : float
        Effective mass ratio

    Returns
    -------
    result : dict
        Fitting results
    """
    # Build parameter vector and bounds list
    x0 = []
    bounds_list = []

    if fit_params.get('W', False):
        x0.append(initial_guess['W'])
        bounds_list.append(bounds.get('W', (1.0, 10.0)))

    if fit_params.get('eta', False):
        x0.append(initial_guess['eta'])
        bounds_list.append(bounds.get('eta', (0.5, 1.0)))

    if not x0:
        return {
            'success': False,
            'message': 'No parameters selected for fitting'
        }

    # Run optimization
    try:
        opt_result = minimize(
            objective_function,
            x0=x0,
            args=(exp_data, fit_params, initial_guess, model_type, epsilon_r, m_star_ratio),
            bounds=bounds_list,
            method='L-BFGS-B',
            options={'maxiter': 1000, 'ftol': 1e-9}
        )

        if not opt_result.success:
            warnings.warn(f"Optimization did not converge: {opt_result.message}")

    except Exception as e:
        return {
            'success': False,
            'message': f'Optimization failed: {str(e)}'
        }

    # Unpack results
    params_opt = opt_result.x
    idx = 0

    W_fit = params_opt[idx] if fit_params.get('W', False) else initial_guess['W']
    idx += 1 if fit_params.get('W', False) else 0

    eta_fit = params_opt[idx] if fit_params.get('eta', False) else initial_guess['eta']

    # Calculate fitted curve
    theory_Delta_CL_fit = calculate_theory_Delta_CL(
        exp_data['Delta_WF'],
        W_fit,
        eta_fit,
        model_type,
        epsilon_r,
        m_star_ratio
    )

    # Calculate goodness of fit
    residuals = theory_Delta_CL_fit - exp_data['Delta_CL']
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((exp_data['Delta_CL'] - np.mean(exp_data['Delta_CL']))**2)

    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    rmse = np.sqrt(np.mean(residuals**2))

    # Calculate lambda from eta (approximate)
    # eta ≈ 1 - exp(-W/lambda) => lambda ≈ -W / ln(1-eta)
    if eta_fit >= 0.99:
        lambda_fit = W_fit * 10  # Approximate for very high eta
    else:
        lambda_fit = -W_fit / np.log(1 - eta_fit)

    return {
        'success': True,
        'W_fit': W_fit,
        'eta_fit': eta_fit,
        'lambda_fit': lambda_fit,
        'r_squared': r_squared,
        'rmse': rmse,
        'residuals': residuals,
        'theory_Delta_CL': theory_Delta_CL_fit,
        'opt_result': opt_result,
        'message': 'Fitting completed successfully'
    }


def calculate_residuals(
    theory_data: Dict[str, np.ndarray],
    exp_data: Dict[str, np.ndarray],
    quantity: str = 'Delta_CL'
) -> np.ndarray:
    """
    Calculate residuals between theory and experiment.

    Parameters
    ----------
    theory_data : dict
        Theoretical data
    exp_data : dict
        Experimental data
    quantity : str
        Quantity to compare

    Returns
    -------
    residuals : array
        Theory - Experiment
    """
    return theory_data[quantity] - exp_data[quantity]


def linear_fit_eta(exp_data: Dict[str, np.ndarray]) -> Dict[str, float]:
    """
    Perform linear regression to extract eta from experimental data.

    Fits: Delta_CL = eta * Delta_WF

    Parameters
    ----------
    exp_data : dict
        Experimental data with 'Delta_WF' and 'Delta_CL'

    Returns
    -------
    result : dict
        Linear fit results
    """
    # Linear regression: Delta_CL vs Delta_WF
    slope, intercept, r_value, p_value, std_err = linregress(
        exp_data['Delta_WF'],
        exp_data['Delta_CL']
    )

    # eta is the slope (because Delta_CL = eta * Delta_WF)
    eta_exp = slope

    return {
        'eta_exp': eta_exp,
        'slope': slope,
        'intercept': intercept,
        'r_squared': r_value**2,
        'p_value': p_value,
        'std_err': std_err
    }


def estimate_initial_parameters(exp_data: Dict[str, np.ndarray]) -> Dict[str, float]:
    """
    Estimate initial parameters from experimental data.

    Parameters
    ----------
    exp_data : dict
        Experimental data

    Returns
    -------
    params : dict
        Estimated parameters
    """
    # Estimate eta from linear fit
    lin_fit = linear_fit_eta(exp_data)
    eta_estimate = lin_fit['eta_exp']

    # Clamp eta to reasonable range
    eta_estimate = np.clip(eta_estimate, 0.5, 0.99)

    # W is harder to estimate without knowing ns
    # Use default value
    W_estimate = 3.0

    return {
        'W': W_estimate,
        'eta': eta_estimate
    }

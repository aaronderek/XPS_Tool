"""
Plotting functions for 2DEG visualization.
"""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from physics.units import ns_to_display


def create_ns_vs_Phi_s_plot(curves_data, title="Sheet Density vs Surface Potential",
                            show_uncertainty=False, model_obj=None, m_star_range=(0.30, 0.35)):
    """
    Create Figure 1: ns vs Φs plot.

    Parameters
    ----------
    curves_data : list of dict
        List of curve data, each dict containing:
        - 'name': curve label
        - 'Phi_s': array of surface potential values (eV)
        - 'ns': array of sheet density values (m⁻²)
        - 'color': optional color
    title : str
        Plot title
    show_uncertainty : bool
        If True, show m* uncertainty band
    model_obj : object
        Model object to calculate uncertainty (needed if show_uncertainty=True)
    m_star_range : tuple
        Range of m*/m0 ratios for uncertainty (default: (0.30, 0.35))

    Returns
    -------
    fig : plotly.graph_objects.Figure
        Plotly figure object
    """
    from physics.constants import M0

    fig = go.Figure()

    for i, curve in enumerate(curves_data):
        Phi_s = curve['Phi_s']
        ns = curve['ns']
        name = curve['name']
        color = curve.get('color', None)

        # Convert ns to display units (10¹³ cm⁻²)
        ns_display = ns_to_display(ns)

        # Add uncertainty band if requested (only for first curve)
        if show_uncertainty and i == 0 and model_obj is not None:
            from models import TriangularModel, FangHowardModel, ParabolicModel

            # Get model class and parameters
            model_class = type(model_obj)
            epsilon_r = model_obj.epsilon_r
            W_nm = getattr(model_obj, 'W_nm', None)
            if W_nm is None:
                W_nm = getattr(model_obj, 'W_initial', 3.0) * 1e9  # Convert to nm

            # Calculate bounds
            model_low = model_class(m_star_range[0] * M0, epsilon_r, W_nm)
            model_high = model_class(m_star_range[1] * M0, epsilon_r, W_nm)

            ns_low = model_low.calculate_ns(Phi_s)
            ns_high = model_high.calculate_ns(Phi_s)

            ns_low_display = ns_to_display(ns_low)
            ns_high_display = ns_to_display(ns_high)

            # Add filled area
            fig.add_trace(go.Scatter(
                x=np.concatenate([Phi_s, Phi_s[::-1]]),
                y=np.concatenate([ns_low_display, ns_high_display[::-1]]),
                fill='toself',
                fillcolor='rgba(100, 150, 250, 0.2)',
                line=dict(width=0),
                showlegend=True,
                name=f'm* uncertainty ({m_star_range[0]:.2f}-{m_star_range[1]:.2f} m₀)',
                hoverinfo='skip'
            ))

        fig.add_trace(go.Scatter(
            x=Phi_s,
            y=ns_display,
            mode='lines',
            name=name,
            line=dict(color=color, width=2),
            hovertemplate='Φₛ: %{x:.3f} eV<br>nₛ: %{y:.2f} × 10¹³ cm⁻²<extra></extra>'
        ))

    fig.update_layout(
        title=title,
        xaxis_title="Surface Potential Φₛ (eV)",
        yaxis_title="Sheet Density nₛ (10¹³ cm⁻²)",
        hovermode='closest',
        template='plotly_white',
        font=dict(size=12),
        width=800,
        height=500,
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01
        )
    )

    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')

    return fig


def create_Delta_WF_vs_ns_plot(curves_data, title="Work Function Change vs Sheet Density",
                               show_uncertainty=False, model_obj=None, m_star_range=(0.30, 0.35)):
    """
    Create Figure 2: ΔWF vs ns plot.

    Parameters
    ----------
    curves_data : list of dict
        List of curve data, each dict containing:
        - 'name': curve label
        - 'ns': array of sheet density values (m⁻²)
        - 'Delta_WF': array of work function change (eV)
        - 'with_adsorbate': bool, whether adsorbates are included
        - 'color': optional color
    title : str
        Plot title
    show_uncertainty : bool
        If True, show m* uncertainty band
    model_obj : object
        Model object to calculate uncertainty (needed if show_uncertainty=True)
    m_star_range : tuple
        Range of m*/m0 ratios for uncertainty (default: (0.30, 0.35))

    Returns
    -------
    fig : plotly.graph_objects.Figure
        Plotly figure object
    """
    from physics.constants import M0

    fig = go.Figure()

    for i, curve in enumerate(curves_data):
        ns = curve['ns']
        Delta_WF = curve['Delta_WF']
        name = curve['name']
        with_ads = curve.get('with_adsorbate', False)
        color = curve.get('color', None)

        # Convert ns to display units
        ns_display = ns_to_display(ns)

        # Sanity checks to catch unit conversion errors
        if np.any(np.abs(Delta_WF) > 5.0):
            import warnings
            warnings.warn(f"ΔWF values exceed ±5 eV, possible unit error in curve '{name}'")
        if np.any(ns_display < 1e-3) or np.any(ns_display > 1e3):
            import warnings
            warnings.warn(f"nₛ display values outside plausible range (0.001-1000 × 10¹³ cm⁻²) in curve '{name}'")

        # Add uncertainty band if requested (only for first curve without adsorbates)
        if show_uncertainty and i == 0 and not with_ads and model_obj is not None:
            from models import TriangularModel, FangHowardModel, ParabolicModel

            # Get model class and parameters
            model_class = type(model_obj)
            epsilon_r = model_obj.epsilon_r
            W_nm = getattr(model_obj, 'W_nm', None)
            if W_nm is None:
                W_nm = getattr(model_obj, 'W_initial', 3.0) * 1e9  # Convert to nm

            # Calculate bounds
            model_low = model_class(m_star_range[0] * M0, epsilon_r, W_nm)
            model_high = model_class(m_star_range[1] * M0, epsilon_r, W_nm)

            Delta_WF_low = model_low.calculate_Delta_WF(ns)
            Delta_WF_high = model_high.calculate_Delta_WF(ns)

            # Add filled area
            fig.add_trace(go.Scatter(
                x=np.concatenate([ns_display, ns_display[::-1]]),
                y=np.concatenate([Delta_WF_low, Delta_WF_high[::-1]]),
                fill='toself',
                fillcolor='rgba(100, 150, 250, 0.2)',
                line=dict(width=0),
                showlegend=True,
                name=f'm* uncertainty ({m_star_range[0]:.2f}-{m_star_range[1]:.2f} m₀)',
                hoverinfo='skip'
            ))

        # Use different line styles for with/without adsorbates
        line_style = 'dash' if with_ads else 'solid'

        fig.add_trace(go.Scatter(
            x=ns_display,
            y=Delta_WF,
            mode='lines',
            name=name,
            line=dict(color=color, width=2, dash=line_style),
            hovertemplate='nₛ: %{x:.2f} × 10¹³ cm⁻²<br>ΔWF: %{y:.3f} eV<extra></extra>'
        ))

    fig.update_layout(
        title=title,
        xaxis_title="Sheet Density nₛ (10¹³ cm⁻²)",
        yaxis_title="Work Function Change ΔWF (eV)",
        hovermode='closest',
        template='plotly_white',
        font=dict(size=12),
        width=800,
        height=500,
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99
        )
    )

    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')

    return fig


def create_potential_profile_plot(z_nm, V_eV, title="Potential Profile"):
    """
    Create potential energy profile V(z) plot.

    Parameters
    ----------
    z_nm : array
        Depth positions in nm
    V_eV : array or dict of arrays
        Potential energy in eV. Can be single array or dict with multiple curves
    title : str
        Plot title

    Returns
    -------
    fig : plotly.graph_objects.Figure
        Plotly figure object
    """
    fig = go.Figure()

    if isinstance(V_eV, dict):
        for name, V in V_eV.items():
            fig.add_trace(go.Scatter(
                x=z_nm,
                y=V,
                mode='lines',
                name=name,
                line=dict(width=2),
                hovertemplate='z: %{x:.2f} nm<br>V: %{y:.3f} eV<extra></extra>'
            ))
    else:
        fig.add_trace(go.Scatter(
            x=z_nm,
            y=V_eV,
            mode='lines',
            line=dict(width=2, color='blue'),
            hovertemplate='z: %{x:.2f} nm<br>V: %{y:.3f} eV<extra></extra>',
            showlegend=False
        ))

    fig.update_layout(
        title=title,
        xaxis_title="Depth z (nm)",
        yaxis_title="Potential Energy V (eV)",
        hovermode='closest',
        template='plotly_white',
        font=dict(size=12),
        width=700,
        height=400
    )

    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')

    return fig


def create_electron_density_plot(z_nm, n_z, title="Electron Density Distribution"):
    """
    Create electron density n(z) profile plot.

    Parameters
    ----------
    z_nm : array
        Depth positions in nm
    n_z : array
        Electron density (arbitrary units)
    title : str
        Plot title

    Returns
    -------
    fig : plotly.graph_objects.Figure
        Plotly figure object
    """
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=z_nm,
        y=n_z,
        mode='lines',
        line=dict(width=2, color='red'),
        fill='tozeroy',
        fillcolor='rgba(255, 0, 0, 0.1)',
        hovertemplate='z: %{x:.2f} nm<br>n(z): %{y:.2e}<extra></extra>',
        showlegend=False
    ))

    fig.update_layout(
        title=title,
        xaxis_title="Depth z (nm)",
        yaxis_title="Electron Density n(z) (a.u.)",
        hovermode='closest',
        template='plotly_white',
        font=dict(size=12),
        width=700,
        height=400
    )

    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')

    return fig


def create_xps_weight_plot(z_nm, w_z, title="XPS Sampling Weight"):
    """
    Create XPS weight function w(z) plot.

    Parameters
    ----------
    z_nm : array
        Depth positions in nm
    w_z : array
        Weight function values
    title : str
        Plot title

    Returns
    -------
    fig : plotly.graph_objects.Figure
        Plotly figure object
    """
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=z_nm,
        y=w_z,
        mode='lines',
        line=dict(width=2, color='green'),
        fill='tozeroy',
        fillcolor='rgba(0, 255, 0, 0.1)',
        hovertemplate='z: %{x:.2f} nm<br>w(z): %{y:.2e}<extra></extra>',
        showlegend=False
    ))

    fig.update_layout(
        title=title,
        xaxis_title="Depth z (nm)",
        yaxis_title="XPS Weight w(z) (nm⁻¹)",
        hovermode='closest',
        template='plotly_white',
        font=dict(size=12),
        width=700,
        height=400
    )

    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')

    return fig


def create_combined_profile_plot(z_nm, V_eV, n_z=None, w_z=None):
    """
    Create combined plot with V(z), n(z), and w(z).

    Parameters
    ----------
    z_nm : array
        Depth positions in nm
    V_eV : array
        Potential energy in eV
    n_z : array, optional
        Electron density
    w_z : array, optional
        XPS weight function

    Returns
    -------
    fig : plotly.graph_objects.Figure
        Plotly figure with subplots
    """
    n_plots = 1 + (n_z is not None) + (w_z is not None)

    fig = make_subplots(
        rows=n_plots, cols=1,
        subplot_titles=("Potential V(z)",
                       "Electron Density n(z)" if n_z is not None else None,
                       "XPS Weight w(z)" if w_z is not None else None),
        vertical_spacing=0.12
    )

    # Plot V(z)
    fig.add_trace(
        go.Scatter(x=z_nm, y=V_eV, mode='lines',
                  line=dict(color='blue', width=2), name='V(z)'),
        row=1, col=1
    )

    row = 2
    if n_z is not None:
        fig.add_trace(
            go.Scatter(x=z_nm, y=n_z, mode='lines',
                      line=dict(color='red', width=2), name='n(z)',
                      fill='tozeroy', fillcolor='rgba(255, 0, 0, 0.1)'),
            row=row, col=1
        )
        row += 1

    if w_z is not None:
        fig.add_trace(
            go.Scatter(x=z_nm, y=w_z, mode='lines',
                      line=dict(color='green', width=2), name='w(z)',
                      fill='tozeroy', fillcolor='rgba(0, 255, 0, 0.1)'),
            row=row, col=1
        )

    fig.update_xaxes(title_text="Depth z (nm)", row=n_plots, col=1)
    fig.update_yaxes(title_text="V (eV)", row=1, col=1)

    if n_z is not None:
        fig.update_yaxes(title_text="n(z) (a.u.)", row=2, col=1)

    if w_z is not None:
        fig.update_yaxes(title_text="w(z) (nm⁻¹)", row=n_plots, col=1)

    fig.update_layout(
        height=300 * n_plots,
        width=700,
        template='plotly_white',
        showlegend=False
    )

    return fig


def create_comparison_CL_vs_WF_plot(
    exp_data,
    theory_data=None,
    fit_result=None,
    model=None,
    lambda_nm=None,
    theta_deg=None,
    title="Core Level Shift vs Work Function Change"
):
    """
    Create comparison plot: ΔE_CL vs ΔWF.

    Parameters
    ----------
    exp_data : dict
        Experimental data with 'Delta_WF' and 'Delta_CL'
    theory_data : dict, optional
        Theory data with 'Delta_WF' and 'Delta_CL'
    fit_result : dict, optional
        Fitting results
    model : object, optional
        Model object (TriangularModel, FangHowardModel, or ParabolicModel)
        to calculate theory curve
    lambda_nm : float, optional
        XPS mean free path in nm (needed if model is provided)
    theta_deg : float, optional
        XPS detection angle in degrees (needed if model is provided)
    title : str
        Plot title

    Returns
    -------
    fig : plotly.graph_objects.Figure
        Plotly figure object
    """
    from physics.xps import calculate_core_level_shift
    from physics.units import nm_to_m

    fig = go.Figure()

    # ========== 1. Calculate Theory Curve (if model provided) ==========
    eta_theory = None
    if model is not None and lambda_nm is not None and theta_deg is not None:
        # DEBUG: Print model parameters
        print("=" * 60)
        print("=== DEBUG: Theory Curve Calculation ===")
        print(f"Model: {model}")
        print(f"Model type: {type(model).__name__}")

        # Get model parameters
        W_nm_model = getattr(model, 'W_nm', None)
        W_m_model = getattr(model, 'W', None)
        print(f"W_nm (from model) = {W_nm_model} nm")
        W_m_str = f"{W_m_model*1e9:.3f}" if W_m_model else "N/A"
        print(f"W_m (from model) = {W_m_str} nm")
        print(f"λ = {lambda_nm} nm")
        print(f"θ = {theta_deg}°")

        # Calculate λeff
        theta_rad = theta_deg * np.pi / 180
        lambda_eff = lambda_nm * np.cos(theta_rad)
        print(f"λeff = λ·cos(θ) = {lambda_eff:.3f} nm")

        # Expected η from simple formula
        if W_nm_model:
            eta_expected = 1 - np.exp(-W_nm_model / lambda_eff)
            print(f"Expected η ≈ 1 - exp(-W/λeff) = 1 - exp(-{W_nm_model}/{lambda_eff:.3f}) = {eta_expected:.3f}")
        print("=" * 60)

        # Determine Delta_WF range from experimental data
        if exp_data is not None and len(exp_data['Delta_WF']) > 0:
            Delta_WF_min = min(exp_data['Delta_WF']) - 0.05
            Delta_WF_max = max(exp_data['Delta_WF']) + 0.05
        else:
            Delta_WF_min = -0.6
            Delta_WF_max = 0.0

        # Generate theoretical curve
        Delta_WF_theory = np.linspace(Delta_WF_min, Delta_WF_max, 100)
        Delta_CL_theory = []

        for i, Delta_WF in enumerate(Delta_WF_theory):
            # From ΔWF to Phi_s: ΔWF = -Phi_s (ignoring adsorbates)
            Phi_s_eV = -Delta_WF

            # Skip if Phi_s is too small or negative
            if Phi_s_eV < 0.01:
                Delta_CL_theory.append(0)
                continue

            # Create z array for integration
            lambda_m = nm_to_m(lambda_nm)
            z_max = 3 * lambda_m
            z_array = np.linspace(0, z_max, 500)

            # Get potential profile from model
            V_z = model.get_potential(Phi_s_eV, z_array)

            # DEBUG: First point details
            if i == 0:
                print(f"\n=== First Point Calculation ===")
                print(f"ΔWF = {Delta_WF:.3f} eV")
                print(f"Φs = -ΔWF = {Phi_s_eV:.3f} eV")
                print(f"z_max = 3λ = {z_max*1e9:.3f} nm")
                print(f"z_array: {len(z_array)} points from 0 to {z_max*1e9:.3f} nm")
                print(f"V_z shape: {V_z.shape}")
                print(f"V(0) = {V_z[0]:.6e} J = {V_z[0]/1.602e-19:.3f} eV")
                if W_m_model and W_m_model < z_max:
                    idx_W = np.argmin(np.abs(z_array - W_m_model))
                    print(f"V(W) at z={z_array[idx_W]*1e9:.3f}nm = {V_z[idx_W]:.6e} J = {V_z[idx_W]/1.602e-19:.3f} eV")
                print(f"V_max = {np.max(V_z):.6e} J = {np.max(V_z)/1.602e-19:.3f} eV")
                print(f"V_min = {np.min(V_z):.6e} J = {np.min(V_z)/1.602e-19:.3f} eV")

            # Calculate XPS shift
            Delta_E_CL, _ = calculate_core_level_shift(
                V_z, z_array, lambda_nm, theta_deg
            )

            # DEBUG: First point result
            if i == 0:
                print(f"Result: ΔE_CL = {Delta_E_CL:.3f} eV")
                print(f"Ratio: |ΔE_CL/Φs| = {abs(Delta_E_CL/Phi_s_eV):.3f}")
                print("=" * 60)

            Delta_CL_theory.append(Delta_E_CL)

        Delta_CL_theory = np.array(Delta_CL_theory)

        # Calculate theoretical η from slope
        if len(Delta_WF_theory) > 2:
            eta_theory = np.polyfit(Delta_WF_theory, Delta_CL_theory, 1)[0]

            print(f"\n=== Final Result ===")
            print(f"η_theory (from linear fit slope) = {eta_theory:.3f}")
            print("=" * 60)

        # Plot theory curve
        fig.add_trace(go.Scatter(
            x=Delta_WF_theory,
            y=Delta_CL_theory,
            mode='lines',
            name='Theory',
            line=dict(color='blue', width=3),
            hovertemplate='ΔWF: %{x:.3f} eV<br>ΔE_CL: %{y:.3f} eV (theory)<extra></extra>'
        ))

    # ========== 2. Experimental Data ==========
    if exp_data is not None and len(exp_data['Delta_WF']) > 0:
        fig.add_trace(go.Scatter(
            x=exp_data['Delta_WF'],
            y=exp_data['Delta_CL'],
            mode='markers',
            name='Experiment',
            marker=dict(size=12, color='red', symbol='circle',
                       line=dict(width=1.5, color='black')),
            hovertemplate='ΔWF: %{x:.3f} eV<br>ΔE_CL: %{y:.3f} eV (exp)<extra></extra>'
        ))

        # ========== 3. Experimental Linear Fit ==========
        if len(exp_data['Delta_WF']) >= 3:
            from scipy.stats import linregress

            slope, intercept, r_value, _, std_err = linregress(
                exp_data['Delta_WF'],
                exp_data['Delta_CL']
            )

            # Fit line
            x_fit = np.array([exp_data['Delta_WF'].min(), exp_data['Delta_WF'].max()])
            y_fit = slope * x_fit + intercept

            fig.add_trace(go.Scatter(
                x=x_fit,
                y=y_fit,
                mode='lines',
                name=f'Exp. fit (η={slope:.3f})',
                line=dict(color='red', width=2, dash='dash'),
                hovertemplate=f'Linear fit<br>η={slope:.3f}<br>R²={r_value**2:.4f}<extra></extra>'
            ))

            # Add annotation with fit statistics (positioned at top-right to avoid overlapping theory curve)
            annotation_text = f"<b>Experimental Fit:</b><br>η_exp = {slope:.3f} ± {std_err:.3f}<br>R² = {r_value**2:.4f}"
            if eta_theory is not None:
                annotation_text += f"<br><br><b>Theory:</b><br>η_theory = {eta_theory:.3f}<br>Difference: {abs(slope-eta_theory)/eta_theory*100:.1f}%"

            fig.add_annotation(
                text=annotation_text,
                xref="paper", yref="paper",
                x=0.98, y=0.98,  # Top-right corner
                xanchor='right', yanchor='top',  # Right-aligned
                showarrow=False,
                font=dict(size=11),
                bgcolor="rgba(255,255,255,0.95)",  # Slightly more opaque
                bordercolor="black",
                borderwidth=1.5,
                borderpad=8
            )

    # ========== 4. Legacy theory_data support ==========
    elif theory_data is not None:
        fig.add_trace(go.Scatter(
            x=theory_data['Delta_WF'],
            y=theory_data['Delta_CL'],
            mode='lines',
            name='Theory',
            line=dict(color='blue', width=3),
            hovertemplate='ΔWF: %{x:.3f} eV<br>ΔE_CL: %{y:.3f} eV<extra></extra>'
        ))

    # ========== 5. Fitted Curve (from parameter fitting) ==========
    if fit_result is not None and 'theory_Delta_CL' in fit_result and exp_data is not None:
        # Sort for smooth line
        sort_idx = np.argsort(exp_data['Delta_WF'])
        fig.add_trace(go.Scatter(
            x=exp_data['Delta_WF'][sort_idx],
            y=fit_result['theory_Delta_CL'][sort_idx],
            mode='lines',
            name=f"Fit (η={fit_result['eta_fit']:.3f}, R²={fit_result['r_squared']:.3f})",
            line=dict(color='green', width=3, dash='dash'),
            hovertemplate='ΔWF: %{x:.3f} eV<br>ΔE_CL (fit): %{y:.3f} eV<extra></extra>'
        ))

    fig.update_layout(
        title=title,
        xaxis_title="Work Function Change ΔWF (eV)",
        yaxis_title="Core Level Shift ΔE_CL (eV)",
        hovermode='closest',
        template='plotly_white',
        font=dict(size=12),
        width=800,
        height=600,
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01
        )
    )

    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')

    return fig


def create_residual_analysis_plot(exp_data, fit_result):
    """
    Create comprehensive residual analysis plots.

    Parameters
    ----------
    exp_data : dict
        Experimental data
    fit_result : dict
        Fitting results with residuals

    Returns
    -------
    fig : plotly.graph_objects.Figure
        Plotly figure with subplots
    """
    from plotly.subplots import make_subplots
    from scipy import stats

    residuals = fit_result['residuals']

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            'Residuals vs ΔWF',
            'Residuals Distribution',
            'Predicted vs Actual',
            'Residual Statistics'
        ),
        specs=[[{"type": "scatter"}, {"type": "histogram"}],
               [{"type": "scatter"}, {"type": "table"}]]
    )

    # 1. Residuals vs ΔWF
    fig.add_trace(
        go.Scatter(
            x=exp_data['Delta_WF'],
            y=residuals,
            mode='markers',
            marker=dict(size=10, color='blue'),
            name='Residuals',
            showlegend=False
        ),
        row=1, col=1
    )
    fig.add_hline(y=0, line_dash="dash", line_color="red", row=1, col=1)

    # 2. Histogram
    fig.add_trace(
        go.Histogram(
            x=residuals,
            nbinsx=15,
            marker=dict(color='lightblue', line=dict(color='black', width=1)),
            name='Distribution',
            showlegend=False
        ),
        row=1, col=2
    )

    # 3. Predicted vs Actual
    predicted = exp_data['Delta_CL'] - residuals
    fig.add_trace(
        go.Scatter(
            x=exp_data['Delta_CL'],
            y=predicted,
            mode='markers',
            marker=dict(size=10, color='green'),
            name='Data',
            showlegend=False
        ),
        row=2, col=1
    )
    # Perfect fit line
    min_val = min(exp_data['Delta_CL'].min(), predicted.min())
    max_val = max(exp_data['Delta_CL'].max(), predicted.max())
    fig.add_trace(
        go.Scatter(
            x=[min_val, max_val],
            y=[min_val, max_val],
            mode='lines',
            line=dict(color='red', dash='dash'),
            name='Perfect fit',
            showlegend=False
        ),
        row=2, col=1
    )

    # 4. Statistics table
    mean_res = np.mean(residuals)
    std_res = np.std(residuals)
    max_res = np.max(np.abs(residuals))

    fig.add_trace(
        go.Table(
            header=dict(values=['Statistic', 'Value'],
                       fill_color='lightgray',
                       align='left'),
            cells=dict(values=[
                ['Mean residual', 'Std residual', 'Max |residual|', 'RMSE', 'R²'],
                [f'{mean_res:.4f} eV', f'{std_res:.4f} eV', f'{max_res:.4f} eV',
                 f"{fit_result['rmse']:.4f} eV", f"{fit_result['r_squared']:.4f}"]
            ],
            fill_color='white',
            align='left')
        ),
        row=2, col=2
    )

    # Update axes labels
    fig.update_xaxes(title_text="ΔWF (eV)", row=1, col=1)
    fig.update_yaxes(title_text="Residual (eV)", row=1, col=1)

    fig.update_xaxes(title_text="Residual (eV)", row=1, col=2)
    fig.update_yaxes(title_text="Count", row=1, col=2)

    fig.update_xaxes(title_text="Actual ΔE_CL (eV)", row=2, col=1)
    fig.update_yaxes(title_text="Predicted ΔE_CL (eV)", row=2, col=1)

    fig.update_layout(
        height=800,
        showlegend=False,
        template='plotly_white'
    )

    return fig


def create_annealing_trajectory_plot(exp_data, title="Annealing Trajectory"):
    """
    Create annealing trajectory plot showing evolution with temperature.

    Parameters
    ----------
    exp_data : dict
        Experimental data with 'T_degC', 'Delta_WF', 'Delta_CL'
    title : str
        Plot title

    Returns
    -------
    fig : plotly.graph_objects.Figure
        Plotly figure object
    """
    if 'T_degC' not in exp_data:
        return None

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=('ΔWF vs Temperature', 'ΔE_CL vs Temperature')
    )

    # ΔWF vs T
    fig.add_trace(
        go.Scatter(
            x=exp_data['T_degC'],
            y=exp_data['Delta_WF'],
            mode='lines+markers',
            marker=dict(size=10, color='red'),
            line=dict(color='red', width=2),
            name='ΔWF',
            hovertemplate='T: %{x}°C<br>ΔWF: %{y:.3f} eV<extra></extra>'
        ),
        row=1, col=1
    )

    # ΔE_CL vs T
    fig.add_trace(
        go.Scatter(
            x=exp_data['T_degC'],
            y=exp_data['Delta_CL'],
            mode='lines+markers',
            marker=dict(size=10, color='blue'),
            line=dict(color='blue', width=2),
            name='ΔE_CL',
            hovertemplate='T: %{x}°C<br>ΔE_CL: %{y:.3f} eV<extra></extra>'
        ),
        row=1, col=2
    )

    # Update axes
    fig.update_xaxes(title_text="Temperature (°C)", row=1, col=1)
    fig.update_yaxes(title_text="ΔWF (eV)", row=1, col=1)

    fig.update_xaxes(title_text="Temperature (°C)", row=1, col=2)
    fig.update_yaxes(title_text="ΔE_CL (eV)", row=1, col=2)

    fig.update_layout(
        title=title,
        height=400,
        showlegend=False,
        template='plotly_white'
    )

    return fig

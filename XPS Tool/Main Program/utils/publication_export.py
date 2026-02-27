"""
Publication-quality figure export with journal styles.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.figure import Figure
from typing import Dict, Tuple, Optional, List, Any
import io


def get_journal_style(style_name: str) -> Dict[str, Any]:
    """
    Get journal-specific style configuration.

    Parameters
    ----------
    style_name : str
        Journal name ('Nature', 'Science', 'ACS', 'Grayscale', 'Default')

    Returns
    -------
    style : dict
        Style configuration dictionary
    """
    styles = {
        'Nature': {
            'font_family': 'Arial',
            'colors': ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F'],
            'marker_symbols': ['o', 's', '^', 'D', 'v'],
            'grid': True,
            'grid_alpha': 0.3,
            'grid_style': '--',
            'background': 'white',
            'legend_frame': True,
            'spine_width': 1.0
        },
        'Science': {
            'font_family': 'Helvetica',
            'colors': ['#BC3C29', '#0072B5', '#E18727', '#20854E', '#7876B1'],
            'marker_symbols': ['o', 's', '^', 'D', 'v'],
            'grid': False,
            'background': 'white',
            'legend_frame': False,
            'spine_width': 1.0
        },
        'ACS': {
            'font_family': 'Arial',
            'colors': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'],
            'marker_symbols': ['o', 's', '^', 'D', 'v'],
            'grid': True,
            'grid_alpha': 0.2,
            'grid_style': ':',
            'background': 'white',
            'legend_frame': True,
            'spine_width': 1.2
        },
        'Grayscale': {
            'font_family': 'Arial',
            'colors': ['#000000', '#404040', '#808080', '#A0A0A0', '#C0C0C0'],
            'marker_symbols': ['o', 's', '^', 'D', 'v', 'p', '*'],
            'grid': True,
            'grid_alpha': 0.3,
            'grid_style': '--',
            'background': 'white',
            'legend_frame': True,
            'spine_width': 1.0,
            'use_different_markers': True
        },
        'Default': {
            'font_family': 'Arial',
            'colors': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'],
            'marker_symbols': ['o', 's', '^', 'D', 'v'],
            'grid': True,
            'grid_alpha': 0.3,
            'grid_style': '--',
            'background': 'white',
            'legend_frame': True,
            'spine_width': 1.0
        }
    }

    return styles.get(style_name, styles['Default'])


def setup_matplotlib_style(
    font_size: int = 12,
    line_width: float = 2.0,
    journal_style: Optional[Dict] = None
):
    """
    Configure matplotlib global style settings.

    Parameters
    ----------
    font_size : int
        Base font size in points
    line_width : float
        Line width in points
    journal_style : dict, optional
        Journal-specific style settings
    """
    if journal_style is None:
        journal_style = get_journal_style('Default')

    # Reset to defaults first
    mpl.rcParams.update(mpl.rcParamsDefault)

    # Set font
    mpl.rcParams['font.family'] = journal_style['font_family']
    mpl.rcParams['font.size'] = font_size

    # Set line widths
    mpl.rcParams['axes.linewidth'] = journal_style.get('spine_width', 1.0)
    mpl.rcParams['lines.linewidth'] = line_width
    mpl.rcParams['lines.markersize'] = 8

    # Set tick widths
    mpl.rcParams['xtick.major.width'] = journal_style.get('spine_width', 1.0)
    mpl.rcParams['ytick.major.width'] = journal_style.get('spine_width', 1.0)
    mpl.rcParams['xtick.major.size'] = 6
    mpl.rcParams['ytick.major.size'] = 6

    # Grid
    mpl.rcParams['grid.linewidth'] = 0.5
    mpl.rcParams['grid.alpha'] = journal_style.get('grid_alpha', 0.3)
    mpl.rcParams['grid.linestyle'] = journal_style.get('grid_style', '--')

    # Legend
    mpl.rcParams['legend.frameon'] = journal_style.get('legend_frame', True)
    mpl.rcParams['legend.fontsize'] = font_size - 1

    # Figure
    mpl.rcParams['figure.facecolor'] = journal_style.get('background', 'white')
    mpl.rcParams['axes.facecolor'] = journal_style.get('background', 'white')


def create_publication_figure_1(
    curves_data: List[Dict],
    exp_data: Optional[Dict] = None,
    style_config: Optional[Dict] = None,
    size: Tuple[float, float] = (7.0, 5.0)
) -> Figure:
    """
    Create publication-quality Figure 1: ns vs Φs.

    Parameters
    ----------
    curves_data : list of dict
        Theory curve data
    exp_data : dict, optional
        Experimental data
    style_config : dict, optional
        Style configuration
    size : tuple
        Figure size in inches (width, height)

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object
    """
    if style_config is None:
        style_config = get_journal_style('Default')

    fig, ax = plt.subplots(figsize=size)

    # Plot theory curves
    for i, curve in enumerate(curves_data):
        color_idx = i % len(style_config['colors'])
        ax.plot(
            curve['Phi_s'],
            curve['ns'] / 1e17,  # Convert to 10^13 cm^-2
            color=style_config['colors'][color_idx],
            linewidth=2.5,
            label=curve['name'],
            zorder=2
        )

    # Plot experimental data if available
    if exp_data is not None and 'Phi_s_exp' in exp_data and 'ns_exp' in exp_data:
        ax.scatter(
            exp_data['Phi_s_exp'],
            exp_data['ns_exp'] / 1e17,
            color='red',
            s=100,
            marker='o',
            edgecolors='black',
            linewidths=1.5,
            label='Experiment',
            zorder=10
        )

    # Labels and formatting
    ax.set_xlabel(r'Surface Potential $\Phi_s$ (eV)', fontsize=14, fontweight='bold')
    ax.set_ylabel(r'Sheet Density $n_s$ ($10^{13}$ cm$^{-2}$)', fontsize=14, fontweight='bold')

    # Grid
    if style_config.get('grid', True):
        ax.grid(True, alpha=style_config.get('grid_alpha', 0.3),
                linestyle=style_config.get('grid_style', '--'))

    # Legend
    ax.legend(frameon=style_config.get('legend_frame', True), loc='best')

    # Tight layout
    plt.tight_layout()

    return fig


def create_publication_figure_2(
    curves_data: List[Dict],
    exp_data: Optional[Dict] = None,
    style_config: Optional[Dict] = None,
    size: Tuple[float, float] = (7.0, 5.0)
) -> Figure:
    """
    Create publication-quality Figure 2: ΔWF vs ns.

    Parameters
    ----------
    curves_data : list of dict
        Theory curve data
    exp_data : dict, optional
        Experimental data
    style_config : dict, optional
        Style configuration
    size : tuple
        Figure size in inches

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object
    """
    if style_config is None:
        style_config = get_journal_style('Default')

    fig, ax = plt.subplots(figsize=size)

    # Plot theory curves
    for i, curve in enumerate(curves_data):
        color_idx = i % len(style_config['colors'])
        line_style = '--' if curve.get('with_adsorbate', False) else '-'

        ax.plot(
            curve['ns'] / 1e17,
            curve['Delta_WF'],
            color=style_config['colors'][color_idx],
            linewidth=2.5,
            linestyle=line_style,
            label=curve['name'],
            zorder=2
        )

    # Plot experimental data if available
    if exp_data is not None and 'ns_exp' in exp_data and 'Delta_WF' in exp_data:
        ax.scatter(
            exp_data['ns_exp'] / 1e17,
            exp_data['Delta_WF'],
            color='red',
            s=100,
            marker='o',
            edgecolors='black',
            linewidths=1.5,
            label='Experiment',
            zorder=10
        )

        # Error bars if available
        if 'Delta_WF_err' in exp_data:
            ax.errorbar(
                exp_data['ns_exp'] / 1e17,
                exp_data['Delta_WF'],
                yerr=exp_data['Delta_WF_err'],
                fmt='none',
                ecolor='red',
                capsize=4,
                capthick=1.5,
                zorder=9
            )

    # Labels
    ax.set_xlabel(r'Sheet Density $n_s$ ($10^{13}$ cm$^{-2}$)', fontsize=14, fontweight='bold')
    ax.set_ylabel(r'Work Function Change $\Delta$WF (eV)', fontsize=14, fontweight='bold')

    # Grid
    if style_config.get('grid', True):
        ax.grid(True, alpha=style_config.get('grid_alpha', 0.3),
                linestyle=style_config.get('grid_style', '--'))

    # Legend
    ax.legend(frameon=style_config.get('legend_frame', True), loc='best')

    plt.tight_layout()

    return fig


def create_publication_comparison_figure(
    exp_data: Dict[str, np.ndarray],
    theory_data: Optional[Dict[str, np.ndarray]] = None,
    fit_result: Optional[Dict] = None,
    style_config: Optional[Dict] = None,
    size: Tuple[float, float] = (7.0, 5.0),
    model: Optional[Any] = None,
    lambda_nm: Optional[float] = None,
    theta_deg: Optional[float] = None
) -> Figure:
    """
    Create publication-quality comparison figure: ΔE_CL vs ΔWF.

    Parameters
    ----------
    exp_data : dict
        Experimental data with 'Delta_WF' and 'Delta_CL'
    theory_data : dict, optional
        Theory curve data (legacy support)
    fit_result : dict, optional
        Fitting results
    style_config : dict, optional
        Style configuration
    size : tuple
        Figure size in inches
    model : object, optional
        Model object (TriangularModel, etc.) to calculate theory curve
    lambda_nm : float, optional
        XPS mean free path in nm (needed if model provided)
    theta_deg : float, optional
        XPS detection angle in degrees (needed if model provided)

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object
    """
    from physics.xps import calculate_core_level_shift
    from physics.units import nm_to_m
    from scipy.stats import linregress

    if style_config is None:
        style_config = get_journal_style('Default')

    fig, ax = plt.subplots(figsize=size)

    # ========== 1. Calculate Theory Curve (if model provided) ==========
    eta_theory = None
    if model is not None and lambda_nm is not None and theta_deg is not None:
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

        for Delta_WF in Delta_WF_theory:
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

            # Calculate XPS shift
            Delta_E_CL, _ = calculate_core_level_shift(
                V_z, z_array, lambda_nm, theta_deg
            )

            Delta_CL_theory.append(Delta_E_CL)

        Delta_CL_theory = np.array(Delta_CL_theory)

        # Calculate theoretical η from slope
        if len(Delta_WF_theory) > 2:
            eta_theory = np.polyfit(Delta_WF_theory, Delta_CL_theory, 1)[0]

        # Plot theory curve
        ax.plot(
            Delta_WF_theory,
            Delta_CL_theory,
            color=style_config['colors'][0],
            linewidth=2.8,
            linestyle='-',
            label='Theory (initial params)',
            zorder=2
        )

    # ========== 2. Plot experimental data ==========
    ax.scatter(
        exp_data['Delta_WF'],
        exp_data['Delta_CL'],
        color='red',
        s=120,  # Larger markers for better visibility
        marker='o',
        edgecolors='black',
        linewidths=1.8,
        label='Experimental data',
        zorder=10
    )

    # Error bars if available
    if 'Delta_CL_err' in exp_data:
        ax.errorbar(
            exp_data['Delta_WF'],
            exp_data['Delta_CL'],
            yerr=exp_data['Delta_CL_err'],
            fmt='none',
            ecolor='red',
            capsize=4,
            capthick=1.5,
            zorder=9
        )

    # ========== 3. Experimental Linear Fit ==========
    # Calculate experimental fit statistics (but don't plot the line to avoid clutter)
    eta_exp = None
    r_squared_exp = None
    std_err_exp = None

    if len(exp_data['Delta_WF']) >= 3:
        slope, intercept, r_value, _, std_err = linregress(
            exp_data['Delta_WF'],
            exp_data['Delta_CL']
        )

        eta_exp = slope
        r_squared_exp = r_value**2
        std_err_exp = std_err

        # Add annotation with simplified statistics
        annotation_text = f"Experimental Fit:\n"
        annotation_text += f"η_exp = {slope:.3f} ± {std_err:.3f}\n"
        annotation_text += f"R² = {r_value**2:.4f}"

        if eta_theory is not None:
            diff_pct = abs(slope-eta_theory)/max(abs(eta_theory), 0.001)*100
            annotation_text += f"\n\nInitial Theory:\n"
            annotation_text += f"η_init = {eta_theory:.3f}\n"
            annotation_text += f"Δ = {diff_pct:.1f}%"

        # Text box positioned at left side to avoid overlapping curves
        props = dict(boxstyle='round', facecolor='white',
                    edgecolor='black', alpha=0.95, linewidth=1.5, pad=0.8)
        ax.text(0.02, 0.98, annotation_text, transform=ax.transAxes,
               fontsize=11, verticalalignment='top', horizontalalignment='left', bbox=props)

    # ========== 4. Legacy theory_data support ==========
    elif theory_data is not None:
        ax.plot(
            theory_data['Delta_WF'],
            theory_data['Delta_CL'],
            color=style_config['colors'][0],
            linewidth=2.5,
            linestyle='-',
            label='Theory',
            zorder=2
        )

    # ========== 5. Fitted Curve (from parameter fitting) ==========
    if fit_result is not None and 'theory_Delta_CL' in fit_result:
        # Sort by Delta_WF for smooth line
        sort_idx = np.argsort(exp_data['Delta_WF'])
        ax.plot(
            exp_data['Delta_WF'][sort_idx],
            fit_result['theory_Delta_CL'][sort_idx],
            color=style_config['colors'][2] if len(style_config['colors']) > 2 else 'green',
            linewidth=3.0,
            linestyle='-',  # Solid line for fitted curve (more important than initial theory)
            label=f"Best fit (η={fit_result['eta_fit']:.3f}, R²={fit_result['r_squared']:.3f})",
            zorder=4
        )

    # Labels with larger font
    ax.set_xlabel(r'Work Function Change $\Delta$WF (eV)', fontsize=15, fontweight='bold')
    ax.set_ylabel(r'Core Level Shift $\Delta E_{CL}$ (eV)', fontsize=15, fontweight='bold')

    # Adjust X-axis range to minimize empty space
    if exp_data is not None and len(exp_data['Delta_WF']) > 0:
        x_min = np.min(exp_data['Delta_WF']) - 0.02
        x_max = np.max(exp_data['Delta_WF']) + 0.02
        ax.set_xlim(x_min, x_max)

    # Grid
    if style_config.get('grid', True):
        ax.grid(True, alpha=style_config.get('grid_alpha', 0.3),
                linestyle=style_config.get('grid_style', '--'))

    # Legend - positioned to avoid annotation box
    ax.legend(frameon=style_config.get('legend_frame', True), loc='lower right',
             fontsize=11, framealpha=0.95, edgecolor='black', fancybox=True, shadow=True)

    # Add zero lines
    ax.axhline(y=0, color='k', linestyle='-', linewidth=0.8, alpha=0.5)
    ax.axvline(x=0, color='k', linestyle='-', linewidth=0.8, alpha=0.5)

    plt.tight_layout()

    return fig


def save_figure(
    fig: Figure,
    format_type: str = 'png',
    dpi: int = 300
) -> bytes:
    """
    Save figure to bytes buffer.

    Parameters
    ----------
    fig : Figure
        Matplotlib figure
    format_type : str
        Format ('png', 'svg', 'pdf')
    dpi : int
        DPI for raster formats

    Returns
    -------
    img_bytes : bytes
        Image data
    """
    buf = io.BytesIO()

    if format_type.lower() == 'svg':
        fig.savefig(buf, format='svg', bbox_inches='tight')
    elif format_type.lower() == 'pdf':
        fig.savefig(buf, format='pdf', bbox_inches='tight')
    else:  # png
        fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight')

    buf.seek(0)
    return buf.getvalue()


def get_figure_size_presets() -> Dict[str, Tuple[float, float]]:
    """
    Get common figure size presets for journals.

    Returns
    -------
    presets : dict
        Dictionary of preset names to (width, height) in inches
    """
    return {
        'Single column (3.5")': (3.5, 2.6),
        'Double column (7.0")': (7.0, 5.0),
        'Full page (7.0" x 9.0")': (7.0, 9.0),
        'Square (5.0" x 5.0")': (5.0, 5.0)
    }

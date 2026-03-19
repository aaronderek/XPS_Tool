"""
2DEG Surface Electron Gas Visualization Tool

Interactive tool connecting surface band bending (Φs) → 2DEG sheet density (ns)
→ XPS/UPS observables (ΔWF, ΔE_CL)
"""

import streamlit as st
import numpy as np
import pandas as pd
from datetime import datetime
import io
import os
import hashlib
import matplotlib.pyplot as plt

# Import models
from models import TriangularModel, FangHowardModel, ParabolicModel
from physics.constants import (
    M0, DEFAULT_M_STAR, DEFAULT_EPSILON_R, DEFAULT_TEMPERATURE,
    DEFAULT_W, DEFAULT_LAMBDA_XPS, M_STAR_RANGE, EPSILON_R_RANGE,
    TEMPERATURE_RANGE, W_RANGE, PHI_S_RANGE, LAMBDA_XPS_RANGE,
    THETA_XPS_RANGE, DELTA_PHI_DIP_RANGE, EPS0, Q, HBAR
)
from physics.xps import XPSModel, calculate_xps_weight
from physics.units import ns_to_display, ns_from_display, m_to_nm, nm_to_m, J_to_eV
from ui.plots import (
    create_ns_vs_Phi_s_plot, create_Delta_WF_vs_ns_plot,
    create_potential_profile_plot, create_electron_density_plot,
    create_xps_weight_plot, create_combined_profile_plot,
    create_comparison_CL_vs_WF_plot, create_residual_analysis_plot,
    create_annealing_trajectory_plot
)
from utils.experiment_data import (
    import_experimental_data, get_format_example_text, create_sample_data,
    validate_experimental_data_physics
)
from utils.fitting import (
    run_parameter_fitting, linear_fit_eta, estimate_initial_parameters
)
from utils.publication_export import (
    get_journal_style, setup_matplotlib_style, create_publication_figure_1,
    create_publication_figure_2, create_publication_comparison_figure,
    save_figure, get_figure_size_presets
)
from utils.xps_spectrum import (
    load_uploaded_xps_file,
    detect_spectrum_peaks,
    fit_spectrum_peaks,
    create_xps_spectrum_figure,
    export_xps_figure_bytes,
    get_element_reference_library,
    suggest_reference_labels,
    match_fitted_peaks,
    build_spectra_zip_bytes
)
from utils.ups_spectrum import (
    load_uploaded_ups_file,
    detect_secondary_cutoff,
    detect_valence_band_edge,
    calculate_work_function,
    create_ups_spectrum_figure,
    export_ups_figure_bytes,
    PHOTON_SOURCES,
)

# Page configuration
st.set_page_config(
    page_title="2DEG Visualization Tool",
    page_icon="⚛️",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Title
st.title("2DEG Surface Electron Gas Visualization Tool")
st.markdown("**Interactive tool for exploring 2DEG physics and XPS measurements**")

# Initialize session state for comparison curves
if 'comparison_curves' not in st.session_state:
    st.session_state.comparison_curves = []

# Initialize session state for experimental data
if 'exp_data' not in st.session_state:
    st.session_state.exp_data = None
if 'exp_data_loaded' not in st.session_state:
    st.session_state.exp_data_loaded = False
if 'fit_result' not in st.session_state:
    st.session_state.fit_result = None
if 'xps_fit_store' not in st.session_state:
    st.session_state.xps_fit_store = {}
if 'xps_batch_fit_file' not in st.session_state:
    st.session_state.xps_batch_fit_file = ""
if 'xps_batch_fit_result' not in st.session_state:
    st.session_state.xps_batch_fit_result = None
if 'ups_wf_result' not in st.session_state:
    st.session_state.ups_wf_result = None
if 'ups_cutoff_result' not in st.session_state:
    st.session_state.ups_cutoff_result = None
if 'ups_vbm_result' not in st.session_state:
    st.session_state.ups_vbm_result = None

# Sidebar - Parameter Controls
st.sidebar.header("Parameters")

# Model selection
model_type = st.sidebar.selectbox(
    "Model Selection",
    options=["M1-Triangular", "M2-Fang-Howard", "M3-Parabolic"],
    index=0,
    help="Select the potential well model"
)

st.sidebar.markdown("---")
st.sidebar.subheader("Material Parameters")

# Effective mass
m_star_ratio = st.sidebar.slider(
    "m*/m₀",
    min_value=M_STAR_RANGE[0],
    max_value=M_STAR_RANGE[1],
    value=DEFAULT_M_STAR / M0,
    step=0.01,
    help="Effective mass ratio: m*/m₀ represents the ratio of the electron's effective mass in the oxide to the free electron mass. In SrTiO₃, electrons move in the conduction band with reduced inertia (m* ≈ 0.32m₀) due to the crystal lattice periodic potential. This parameter directly affects quantum confinement energies E_n ∝ (m*)^(-1/3) and the 2DEG formation. Lower m* leads to higher mobility and stronger quantum effects."
)
m_star = m_star_ratio * M0

# Relative permittivity
epsilon_r = st.sidebar.slider(
    "εᵣ (Relative Permittivity)",
    min_value=EPSILON_R_RANGE[0],
    max_value=EPSILON_R_RANGE[1],
    value=DEFAULT_EPSILON_R,
    step=1,
    help="Relative dielectric constant: εᵣ describes how effectively the material screens electric fields compared to vacuum (ε₀). For SrTiO₃, εᵣ ≈ 9 at room temperature (static value ~300 at low T). This parameter controls the relationship between surface field and carrier density: n_s = (ε₀·εᵣ·E_s)/q. Higher εᵣ means stronger screening and higher carrier density for the same band bending. It also affects the confinement width W_eff and quantum subband energies."
)

# Temperature
temperature = st.sidebar.slider(
    "T (K)",
    min_value=TEMPERATURE_RANGE[0],
    max_value=TEMPERATURE_RANGE[1],
    value=DEFAULT_TEMPERATURE,
    step=10,
    help="Temperature in Kelvin: Sets the thermal energy k_B·T which controls the Fermi-Dirac distribution tail and thermal broadening of quantum states. At room temperature (300K), k_B·T ≈ 26 meV. For typical 2DEG with E_F ~ 100-300 meV, the system is partially degenerate. Low T (< 100K) sharpens the Fermi edge and enables observation of quantum oscillations. High T (> 400K) increases thermal smearing and can affect XPS peak widths. Currently used for future thermal effects modeling."
)

st.sidebar.markdown("---")
st.sidebar.subheader("Surface Electrostatics")

# Surface potential
Phi_s = st.sidebar.slider(
    "Φₛ (eV)",
    min_value=PHI_S_RANGE[0],
    max_value=PHI_S_RANGE[1],
    value=0.4,
    step=0.01,
    help="Surface potential (band bending): Φ_s represents the electrostatic potential drop from the bulk to the surface, driving 2DEG formation through downward band bending. Positive Φ_s creates an electron accumulation layer at the surface. This parameter directly determines: (1) the surface electric field E_s, (2) the 2D carrier density n_s, (3) the quantum confinement strength, and (4) the work function change ΔWF. Typical values: 0.2-0.6 eV for oxide 2DEG systems. Measured via work function shifts in UPS/KPFM or core level shifts in XPS."
)

# Depletion width
W_nm = st.sidebar.slider(
    "W (nm)",
    min_value=W_RANGE[0],
    max_value=W_RANGE[1],
    value=DEFAULT_W,
    step=0.1,
    help="Depletion/confinement width: W defines the characteristic length scale over which the surface band bending potential extends into the bulk. Physical interpretation depends on the model: (M1-Triangular) W is the depletion width where V(W)=0; (M2-Fang-Howard) W sets initial guess for self-consistent W_eff=6/b; (M3-Parabolic) W is the full extent of the parabolic well. Typical values: 2-5 nm for strong confinement. Affects: n_s ∝ 1/W, confinement energies, and the slope dΔWF/dn_s. Can be fitted from experimental ΔWF vs n_s data."
)

st.sidebar.markdown("---")
st.sidebar.subheader("Surface Adsorbates")

# Adsorbate checkbox
show_adsorbates = st.sidebar.checkbox(
    "Show with adsorbates",
    value=False,
    help="Include adsorbate dipole layer effects: Enables modeling of surface dipole layers from adsorbed species (O₂, H₂O, organic molecules, etc.). Adsorbates with perpendicular dipole moments create a Helmholtz-type potential step ΔΦ_dip = -(N_ads·μ⊥)/ε₀ that shifts the work function without changing the underlying 2DEG density. This decouples ΔWF from n_s, which is crucial for interpreting UPS/XPS data during annealing experiments. Oxygen typically has μ⊥ ~ 1-2 Debye pointing outward (increases WF), while hydrogen points inward (decreases WF)."
)

if show_adsorbates:
    # Coverage
    coverage = st.sidebar.slider(
        "θ (Coverage)",
        min_value=0.0,
        max_value=1.0,
        value=0.5,
        step=0.05,
        help="Adsorbate coverage fraction: θ = N_ads/N_site represents the fraction of available surface sites occupied by adsorbates. θ=0 means clean surface (vacuum annealing), θ=1 means saturated monolayer (full oxygen exposure). The actual adsorbate density is N_ads = θ·N_site. Coverage controls the magnitude of the dipole shift: ΔΦ_dip ∝ θ. In temperature-programmed experiments, θ typically decreases from ~1 at low T to ~0 at high T as adsorbates desorb. Measured via XPS intensity ratios (e.g., O1s/Ti2p)."
    )

    # Dipole moment
    mu_debye = st.sidebar.slider(
        "μ⊥ (Debye)",
        min_value=0.0,
        max_value=3.0,
        value=1.5,
        step=0.1,
        help="Perpendicular dipole moment: μ⊥ is the component of the molecular dipole moment normal to the surface (in Debye units, 1 D ≈ 3.34×10⁻³⁰ C·m). Determines the magnitude of work function shift per adsorbate: ΔΦ_dip = -(N_ads·μ⊥)/ε₀. Positive μ⊥ (electron-rich end outward, e.g., O⁻-Ti) increases work function. Typical values: O₂ on oxides ~1-2 D, H₂O ~1.85 D, organic molecules 0.5-3 D. Can be estimated from DFT calculations or fitted from experimental WF vs coverage data. The sign convention here: positive μ⊥ → negative ΔΦ_dip in the Helmholtz equation used in this code."
    )

    # Surface site density
    N_site_exp = st.sidebar.slider(
        "N_site (10^x cm⁻²)",
        min_value=13.0,
        max_value=15.0,
        value=14.0,
        step=0.1,
        help="Surface site density: N_site is the total number of available adsorption sites per unit area (displayed as 10^x cm⁻²). For crystalline surfaces, this is related to the surface atomic density. For SrTiO₃ (001), the TiO₂-terminated surface has N_site ~ 6.5×10¹⁴ cm⁻² (one Ti per unit cell). For SrO termination, N_site ~ 6.5×10¹⁴ cm⁻² (one Sr per unit cell). Polycrystalline or amorphous surfaces may have lower effective site densities (10¹³-10¹⁴ cm⁻²). Together with coverage θ and dipole μ⊥, determines the total work function shift from adsorbates."
    )
    N_site = 10 ** N_site_exp  # Convert from log scale to actual value
    from physics.xps import calculate_dipole_from_coverage
    Delta_Phi_dip = calculate_dipole_from_coverage(coverage, N_site, mu_debye)

    # Warning for unusually large dipole shifts
    if abs(Delta_Phi_dip) > 2.0:
        st.sidebar.warning(f"⚠️ Dipole shift unusually large: {Delta_Phi_dip:+.3f} eV. Check Debye or N_site units.")
    else:
        st.sidebar.info(f"ΔΦ_dip = {Delta_Phi_dip:+.3f} eV")
else:
    Delta_Phi_dip = 0.0
    coverage = 0.0
    mu_debye = 0.0

st.sidebar.markdown("---")
st.sidebar.subheader("XPS Parameters")

# XPS mean free path
lambda_xps = st.sidebar.slider(
    "λ (nm)",
    min_value=LAMBDA_XPS_RANGE[0],
    max_value=LAMBDA_XPS_RANGE[1],
    value=DEFAULT_LAMBDA_XPS,
    step=0.1,
    help="Inelastic mean free path (IMFP): λ is the average distance a photoelectron travels in the solid before undergoing inelastic scattering, defining the XPS probing depth. For typical XPS kinetic energies (200-1500 eV), λ ≈ 0.5-3 nm in oxides. The effective sampling depth is λ_eff = λ·cos(θ) where θ is the emission angle. XPS is surface-sensitive due to small λ: 95% of signal comes from 3λ_eff. Controls the weighting function w(z) = (1/λ_eff)·exp(-z/λ_eff) and thus the measured core level shift ΔE_CL. Shorter λ → more surface-sensitive, larger η = |ΔE_CL|/Φ_s."
)

# XPS detection angle
theta_xps = st.sidebar.slider(
    "θ (degrees)",
    min_value=THETA_XPS_RANGE[0],
    max_value=THETA_XPS_RANGE[1],
    value=0,
    step=5,
    help="Detection angle (XPS emission angle): θ is the angle between the surface normal and the photoelectron detection direction. θ=0° (normal emission) samples deeper into the bulk. θ=60° (grazing emission) enhances surface sensitivity. The effective probing depth scales as λ_eff = λ·cos(θ), so increasing θ reduces the sampling depth and increases surface weighting. Angle-resolved XPS (ARXPS) at different θ can probe the depth profile of band bending. In this model, affects the XPS weight function w(z) and thus the core level shift η parameter."
)

st.sidebar.markdown("---")
st.sidebar.subheader("Visualization Options")

# Visualization checkboxes
show_subbands = st.sidebar.checkbox("Show subband levels (Eₙ)", value=False)
show_n_z = st.sidebar.checkbox("Show n(z) distribution", value=False)
show_w_z = st.sidebar.checkbox("Show XPS weight w(z)", value=False)
show_uncertainty = st.sidebar.checkbox("Show m* uncertainty band", value=False)

# Buttons
col1, col2 = st.sidebar.columns(2)
with col1:
    add_comparison = st.button("Add to Compare", use_container_width=True)
with col2:
    clear_comparison = st.button("Clear Compare", use_container_width=True)

# ============================================================================
# INITIALIZE MODEL AND DATA (before tabs, so all tabs can access)
# ============================================================================

# Initialize the selected model
if model_type == "M1-Triangular":
    model = TriangularModel(m_star, epsilon_r, W_nm)
elif model_type == "M2-Fang-Howard":
    model = FangHowardModel(m_star, epsilon_r, W_nm)
else:  # M3-Parabolic
    model = ParabolicModel(m_star, epsilon_r, W_nm)

# Initialize XPS model
xps_model = XPSModel(lambda_xps, theta_xps)

# Generate curve data for Figure 1: ns vs Phi_s
Phi_s_range = np.linspace(0.1, 0.6, 100)
ns_array = model.calculate_ns(Phi_s_range)

curves_fig1 = [{
    'name': f"{model_type}, W={W_nm:.1f}nm",
    'Phi_s': Phi_s_range,
    'ns': ns_array,
    'color': 'blue'
}]

# Generate curve data for Figure 2: Delta_WF vs ns
# Typical 2DEG density range: 0.2 to 2.0 × 10¹³ cm⁻²
# Converting to m⁻²: multiply by 1e13 * 1e4 = 1e17
ns_range = np.linspace(2e16, 2e17, 100)  # m⁻² (0.2 to 2.0 × 10¹³ cm⁻²)

# Calculate Delta_WF for each ns
Delta_WF_array = model.calculate_Delta_WF(ns_range)

# Without adsorbates
curves_fig2 = [{
    'name': f"{model_type}, W={W_nm:.1f}nm (no ads)",
    'ns': ns_range,
    'Delta_WF': Delta_WF_array,
    'with_adsorbate': False,
    'color': 'blue'
}]

# With adsorbates
if show_adsorbates:
    Delta_WF_with_ads = Delta_WF_array + Delta_Phi_dip
    curves_fig2.append({
        'name': f"{model_type}, W={W_nm:.1f}nm (with ads)",
        'ns': ns_range,
        'Delta_WF': Delta_WF_with_ads,
        'with_adsorbate': True,
        'color': 'blue'
    })

# Add comparison curves if any
for comp_curve in st.session_state.comparison_curves:
    # For Figure 1
    curves_fig1.append({
        'name': comp_curve['name'],
        'Phi_s': comp_curve['Phi_s_range'],
        'ns': comp_curve['ns_array'],
        'color': comp_curve.get('color', 'gray')
    })

    # For Figure 2
    curves_fig2.append({
        'name': comp_curve['name'] + " (no ads)",
        'ns': comp_curve['ns_range'],
        'Delta_WF': comp_curve['Delta_WF_array'],
        'with_adsorbate': False,
        'color': comp_curve.get('color', 'gray')
    })

    if comp_curve.get('with_adsorbates', False):
        curves_fig2.append({
            'name': comp_curve['name'] + " (with ads)",
            'ns': comp_curve['ns_range'],
            'Delta_WF': comp_curve['Delta_WF_with_ads'],
            'with_adsorbate': True,
            'color': comp_curve.get('color', 'gray')
        })

# Main area - Create tabs
tab1, tab2, tab3, tab4, tab6, tab7, tab8, tab5 = st.tabs([
    "📊 Core Figures",
    "📈 Additional Plots",
    "🔬 Experiment Comparison",
    "📤 Publication Export",
    "🧲 XPS Spectrum",
    "🔮 UPS Spectrum",
    "ℹ️ About",
    "🧪 Beta Features"
])

# ============================================================================
# TAB 1: CORE FIGURES
# ============================================================================

with tab1:

    # Display physics formulas for the selected model
    st.markdown("### 🔬 Physics Formulas for Selected Model")

    if model_type == "M1-Triangular":
        with st.expander("📐 **M1: Triangular Potential (Constant Electric Field)**", expanded=True):
            st.markdown("""
            **Physical Picture:** Constant electric field at the surface, creating a triangular potential well.

            **Key Equations:**
            """)

            col_f1, col_f2 = st.columns(2)
            with col_f1:
                st.latex(r"V(z) = q \cdot E_s \cdot z \quad (0 \leq z \leq W)")
                st.caption("Potential energy profile")

                st.latex(r"E_s = \frac{\Phi_s}{W}")
                st.caption("Surface electric field")

                st.latex(r"n_s = \frac{\varepsilon_0 \varepsilon_r E_s}{q} = \frac{\varepsilon_0 \varepsilon_r \Phi_s}{q W}")
                st.caption("2D electron sheet density")

            with col_f2:
                st.latex(r"\Delta WF = -\Phi_s")
                st.caption("Work function change")

                st.latex(r"\frac{d\Delta WF}{dn_s} = -\frac{q W}{\varepsilon_0 \varepsilon_r}")
                st.caption("Slope of ΔWF vs n_s")

                st.latex(r"E_n = a_n \left(\frac{\hbar^2}{2m^*}\right)^{1/3} (q E_s)^{2/3}")
                st.caption("Quantum subband energies (Airy zeros: a₀=2.338, a₁=4.088, ...)")

    elif model_type == "M2-Fang-Howard":
        with st.expander("🌊 **M2: Fang-Howard (Self-Consistent Variational)**", expanded=True):
            st.markdown("""
            **Physical Picture:** Self-consistent variational approach with exponentially decaying electron wavefunction.

            **Key Equations:**
            """)

            col_f1, col_f2 = st.columns(2)
            with col_f1:
                st.latex(r"\psi(z) = \sqrt{\frac{b^3}{2}} \cdot z \cdot e^{-bz/2}")
                st.caption("Variational wavefunction (normalized)")

                st.latex(r"b = \left(\frac{12 m^* q E_s}{\hbar^2}\right)^{1/3}")
                st.caption("Variational parameter [m⁻¹]")

                st.latex(r"W_{eff} = \frac{6}{b}")
                st.caption("Effective confinement width")

            with col_f2:
                st.latex(r"n(z) = n_s |\psi(z)|^2 = n_s \frac{b^3}{2} z^2 e^{-bz}")
                st.caption("Electron density distribution")

                st.latex(r"E_s = \frac{2\Phi_s}{W_{eff}}")
                st.caption("Self-consistent surface field")

                st.latex(r"n_s = \frac{\varepsilon_0 \varepsilon_r E_s}{q}")
                st.caption("2D electron sheet density")

            st.info("⚙️ **Self-consistency:** Iteratively solve for E_s until W_eff converges (typically 3-5 iterations)")

    else:  # M3-Parabolic
        with st.expander("📊 **M3: Parabolic Potential (Linear Field Decay)**", expanded=True):
            st.markdown("""
            **Physical Picture:** Linearly decaying electric field, creating a parabolic (harmonic) potential well.

            **Key Equations:**
            """)

            col_f1, col_f2 = st.columns(2)
            with col_f1:
                st.latex(r"V(z) = -\Phi_s \left(1 - \frac{z}{W}\right)^2 \quad (0 \leq z \leq W)")
                st.caption("Parabolic potential profile")

                st.latex(r"E(z) = \frac{2\Phi_s}{W}\left(1 - \frac{z}{W}\right)")
                st.caption("Electric field (linear decay)")

                st.latex(r"E_s = \frac{2\Phi_s}{W}")
                st.caption("Surface field (2× M1!)")

            with col_f2:
                st.latex(r"n_s = \frac{\varepsilon_0 \varepsilon_r E_s}{q} = \frac{2\varepsilon_0 \varepsilon_r \Phi_s}{q W}")
                st.caption("2D electron sheet density (2× M1!)")

                st.latex(r"\frac{d\Delta WF}{dn_s} = -\frac{q W}{2\varepsilon_0 \varepsilon_r}")
                st.caption("Slope of ΔWF vs n_s (half of M1!)")

                st.latex(r"E_n = \hbar\omega\left(n + \frac{1}{2}\right), \quad \omega = \sqrt{\frac{2qE_s}{m^* W}}")
                st.caption("Harmonic oscillator energy levels")

    st.markdown("---")

    # Create and display plots
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Figure 1: Sheet Density vs Surface Potential")
        with st.expander("ℹ️ **Chart Explanation**", expanded=False):
            st.markdown("""
            **What this chart shows:** The relationship between surface band bending (Φₛ) and 2D electron sheet density (nₛ).

            **X-axis (Φₛ):** Surface potential in eV, representing the electrostatic potential drop from bulk to surface. Positive values indicate downward band bending that creates electron accumulation.

            **Y-axis (nₛ):** Two-dimensional electron sheet density in units of 10¹³ cm⁻². This is the total number of electrons per unit area confined at the surface.

            **Physical interpretation:** As band bending increases (larger Φₛ), more electrons accumulate at the surface. The slope depends on the model and depletion width W. Linear relationship for M1/M3, slightly nonlinear for M2 due to self-consistency.

            **Typical values:** For oxide 2DEGs, Φₛ ~ 0.2-0.6 eV corresponds to nₛ ~ 0.5-2.0 × 10¹³ cm⁻².
            """)

        fig1 = create_ns_vs_Phi_s_plot(curves_fig1, show_uncertainty=show_uncertainty, model_obj=model)
        st.plotly_chart(fig1, use_container_width=True)

        # Display current value
        current_ns = model.calculate_ns(Phi_s)
        st.info(f"**Current:** Φₛ = {Phi_s:.3f} eV → nₛ = {ns_to_display(current_ns):.3f} × 10¹³ cm⁻²")

    with col2:
        st.subheader("Figure 2: Work Function Change vs Sheet Density")
        with st.expander("ℹ️ **Chart Explanation**", expanded=False):
            st.markdown("""
            **What this chart shows:** The relationship between 2D electron density (nₛ) and work function change (ΔWF), which is measurable by UPS/KPFM.

            **X-axis (nₛ):** Two-dimensional electron sheet density in units of 10¹³ cm⁻². Represents the carrier density confined at the surface.

            **Y-axis (ΔWF):** Work function change in eV, relative to the clean surface. Negative values indicate work function decrease (easier to extract electrons), which occurs when electrons accumulate at the surface.

            **Physical interpretation:** This is the key experimental observable! The slope dΔWF/dnₛ is model-dependent and directly relates to the confinement width W. Solid lines show intrinsic band bending effect. Dashed lines (when adsorbates enabled) include additional Helmholtz dipole shift that decouples ΔWF from nₛ.

            **Typical values:** ΔWF ~ -0.2 to -0.6 eV for oxide 2DEGs. Slope ranges from -0.2 to -0.6 eV/(10¹³ cm⁻²) depending on W and model.
            """)

        fig2 = create_Delta_WF_vs_ns_plot(curves_fig2, show_uncertainty=show_uncertainty, model_obj=model)
        st.plotly_chart(fig2, use_container_width=True)

        # Display current value
        current_ns = model.calculate_ns(Phi_s)
        current_Delta_WF = model.calculate_Delta_WF(current_ns)
        if show_adsorbates:
            current_Delta_WF_total = current_Delta_WF + Delta_Phi_dip
            st.info(f"**Current:** nₛ = {ns_to_display(current_ns):.3f} × 10¹³ cm⁻² → ΔWF = {current_Delta_WF:.3f} eV (no ads) / {current_Delta_WF_total:.3f} eV (with ads)")
        else:
            st.info(f"**Current:** nₛ = {ns_to_display(current_ns):.3f} × 10¹³ cm⁻² → ΔWF = {current_Delta_WF:.3f} eV")

    # Calculate and display slope
    if hasattr(model, 'get_slope'):
        slope = model.get_slope()
        st.success(f"**Slope of ΔWF vs nₛ:** {slope:.3f} eV/(10¹³ cm⁻²)")

    # XPS core level shift
    if model_type != "M2-Fang-Howard":  # Simpler calculation for M1 and M3
        Delta_E_CL, eta = xps_model.calculate_shift(model, Phi_s)
        st.success(f"**XPS Core Level Shift:** ΔE_CL = {Delta_E_CL:.3f} eV (η = {eta:.2f})")

# ============================================================================
# TAB 2: ADDITIONAL PLOTS
# ============================================================================

with tab2:
    st.subheader("Depth Profiles and Distributions")

    # Display XPS formulas
    with st.expander("📊 **XPS Core Level Shift Formulas**", expanded=False):
        st.markdown("""
        **Physical Picture:** X-ray photoelectron spectroscopy (XPS) probes the surface band bending through depth-weighted sampling of the potential profile.

        **Key Equations:**
        """)

        col_x1, col_x2 = st.columns(2)
        with col_x1:
            st.latex(r"w(z) = \frac{1}{\lambda_{eff}} \exp\left(-\frac{z}{\lambda_{eff}}\right)")
            st.caption("XPS sampling weight (exponential decay)")

            st.latex(r"\lambda_{eff} = \lambda \cdot \cos(\theta)")
            st.caption("Effective probing depth (angle-dependent)")

            st.latex(r"\int_0^\infty w(z) \, dz = 1")
            st.caption("Normalization condition")

        with col_x2:
            st.latex(r"\Delta E_{CL} = -\int_0^\infty w(z) \cdot V(z) \, dz")
            st.caption("Core level shift (depth-weighted potential)")

            st.latex(r"\eta = \frac{|\Delta E_{CL}|}{\Phi_s}")
            st.caption("Effective sampling factor (0 < η < 1)")

            st.latex(r"\Delta \Phi_{dip} = -\frac{N_{ads} \cdot \mu_\perp}{\varepsilon_0}")
            st.caption("Adsorbate dipole shift (Helmholtz equation)")

        st.info("💡 **Interpretation:** η ≈ 0.85 typical for λ ~ 2 nm. Smaller λ or larger θ → more surface-sensitive → larger η.")

    st.markdown("---")

    # Create z array for profiles
    z_max = 3 * W_nm
    z_nm_array = np.linspace(0, z_max, 500)
    z_m_array = nm_to_m(z_nm_array)

    # Get potential profile
    V_z_J = model.get_potential(Phi_s, z_m_array)
    V_z_eV = J_to_eV(V_z_J)

    # Potential profile
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("### Potential Profile V(z)")
        with st.expander("ℹ️ **Chart Explanation**", expanded=False):
            st.markdown("""
            **What this chart shows:** The electrostatic potential energy profile as a function of depth from the surface.

            **X-axis (z):** Depth from the surface in nanometers. z=0 is the surface, increasing z goes into the bulk.

            **Y-axis (V):** Potential energy in eV. The shape depends on the model: triangular (M1) shows linear increase, Fang-Howard (M2) shows soft triangular with exponential decay, parabolic (M3) shows quadratic profile.

            **Physical interpretation:** This potential well confines electrons near the surface. The steeper the gradient near z=0, the stronger the surface electric field E_s = -dV/dz. The depth where V(z)→0 defines the depletion width W. Electrons are confined within the first few nanometers.

            **Model comparison:** M1 (constant field) vs M2 (self-consistent decay) vs M3 (linear field decay). The choice affects n_s(Φ_s) relationship and XPS core level shifts.
            """)

        fig_V = create_potential_profile_plot(z_nm_array, V_z_eV)
        st.plotly_chart(fig_V, use_container_width=True)

    with col2:
        if show_w_z:
            st.markdown("### XPS Sampling Weight w(z)")
            with st.expander("ℹ️ **Chart Explanation**", expanded=False):
                st.markdown("""
                **What this chart shows:** The depth-dependent weighting function for XPS measurements, showing how much each depth contributes to the measured signal.

                **X-axis (z):** Depth from the surface in nanometers. z=0 is the surface.

                **Y-axis (w):** Sampling weight in units of nm⁻¹. This represents the probability density that a photoelectron detected in XPS originated from depth z.

                **Physical interpretation:** w(z) = (1/λ_eff)·exp(-z/λ_eff) describes exponential attenuation of photoelectrons due to inelastic scattering. Most signal (63%) comes from 0 to λ_eff, and 95% from 0 to 3λ_eff. The effective depth λ_eff = λ·cos(θ) decreases with increasing emission angle θ, making grazing angle more surface-sensitive.

                **Use in XPS:** The measured core level shift ΔE_CL is the integral of w(z)·V(z), not just the surface value. This is why η = |ΔE_CL|/Φ_s < 1.
                """)

            w_z = calculate_xps_weight(z_m_array, lambda_xps, theta_xps)
            fig_w = create_xps_weight_plot(z_nm_array, w_z)
            st.plotly_chart(fig_w, use_container_width=True)

    # Electron density for M2
    if show_n_z and model_type == "M2-Fang-Howard":
        st.markdown("### Electron Density Distribution n(z)")
        with st.expander("ℹ️ **Chart Explanation**", expanded=False):
            st.markdown("""
            **What this chart shows:** The spatial distribution of electron density as a function of depth for the M2 Fang-Howard model.

            **X-axis (z):** Depth from the surface in nanometers. z=0 is the surface.

            **Y-axis (n):** Electron density in arbitrary units (normalized). Shows the probability density of finding electrons at each depth.

            **Physical interpretation:** For M2 Fang-Howard, n(z) = n_s·|ψ(z)|² = n_s·(b³/2)·z²·exp(-bz). The distribution peaks at depth z_peak = 2/b, then decays exponentially. Average depth is ⟨z⟩ = 3/b. This is more realistic than M1/M3 which assume sharp or infinite wells.

            **Parameters:** The variational parameter b determines the confinement: larger b → tighter confinement near surface. b is calculated self-consistently from the electric field and material properties.

            **Note:** Only available for M2 Fang-Howard model which explicitly calculates the wavefunction.
            """)

        st.info("""
        **M2 Fang-Howard Model:**
        Variational approximation using single effective wavefunction.
        n(z) is calculated from |ψ(z)|² = (b³/2) z² exp(-bz).
        """)
        n_z = model.get_electron_density(Phi_s, z_m_array)
        fig_n = create_electron_density_plot(z_nm_array, n_z)
        st.plotly_chart(fig_n, use_container_width=True)

        # Display key parameters for M2
        Es, W_eff, _ = model.solve_self_consistent(Phi_s)
        b_M2 = model.calculate_b_parameter(Es)
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Variational parameter b", f"{b_M2*1e-9:.2f} nm⁻¹")
        with col2:
            st.metric("Effective width W_eff", f"{W_eff*1e9:.2f} nm")
        with col3:
            st.metric("Average depth ⟨z⟩", f"{3/b_M2*1e9:.2f} nm")

    # Subband energies (only for M1 and M3)
    if show_subbands:
        # Only show if model has subband methods (M1 and M3)
        if hasattr(model, 'get_subband_energies'):
            st.markdown("### Quantum Subband Energy Levels")
            energies = model.get_subband_energies(Phi_s, n_levels=3)
            df_energies = pd.DataFrame({
                'Level': ['E₁', 'E₂', 'E₃'],
                'Energy (eV)': energies
            })
            st.table(df_energies)
        elif hasattr(model, 'get_harmonic_levels'):
            st.markdown("### Quantum Subband Energy Levels")
            energies = model.get_harmonic_levels(Phi_s, n_levels=3)
            df_energies = pd.DataFrame({
                'Level': ['E₀', 'E₁', 'E₂'],
                'Energy (eV)': energies
            })
            st.table(df_energies)
        elif model_type == "M2-Fang-Howard":
            # Show explanation for M2
            st.info("""
            **ℹ️ Why no subband levels for M2?**

            The Fang-Howard model is a **variational approximation** that provides:
            - ✅ An effective single wavefunction ψ(z)
            - ✅ Electron density distribution n(z) = ns·|ψ(z)|²
            - ✅ Ground state energy (first subband approximation)

            It does **not** calculate discrete subband levels (E₀, E₁, E₂...)
            because it's not solving the full Schrödinger eigenvalue problem.

            👉 For multi-subband analysis, use **M1 (Triangular)** or **M3 (Parabolic)**.
            """)

    # Self-Consistent S-P Diagnostic (M2 only) - TEMPORARILY HIDDEN
    # TODO: Re-enable after fixing the visualization
    pass

# ============================================================================
# TAB 3: EXPERIMENT COMPARISON
# ============================================================================

with tab3:
    st.header("🔬 Theory-Experiment Comparison")

    # Data import section
    st.subheader("1. Import Experimental Data")

    col1, col2 = st.columns([3, 1])

    with col1:
        uploaded_file = st.file_uploader(
            "Upload CSV file with experimental data",
            type=['csv'],
            help="See format requirements below"
        )

    with col2:
        if st.button("📋 Show Format"):
            st.session_state['show_format'] = not st.session_state.get('show_format', False)

        if st.button("📥 Load Sample Data"):
            # Create sample data
            sample_df = create_sample_data()
            sample_csv = sample_df.to_csv(index=False).encode('utf-8')

            # Import sample data
            import_result = import_experimental_data(sample_csv)
            if import_result['success']:
                st.session_state.exp_data = import_result['data']
                st.session_state.exp_data_loaded = True
                st.success(import_result['message'])
                st.rerun()

    # Show format example
    if st.session_state.get('show_format', False):
        st.markdown(get_format_example_text())

    # Process uploaded file
    if uploaded_file is not None:
        file_content = uploaded_file.read()
        import_result = import_experimental_data(file_content)

        if import_result['success']:
            st.session_state.exp_data = import_result['data']
            st.session_state.exp_data_loaded = True
            st.success(import_result['message'])
        else:
            st.error(import_result['message'])

    # Clear data button
    if st.session_state.exp_data_loaded:
        col1, col2 = st.columns([5, 1])
        with col2:
            if st.button("🗑️ Clear Data"):
                st.session_state.exp_data = None
                st.session_state.exp_data_loaded = False
                st.session_state.fit_result = None
                st.rerun()

    # Display data and analysis
    if st.session_state.exp_data_loaded:
        exp_data = st.session_state.exp_data

        st.markdown("---")
        st.subheader("2. Data Preview")

        # Create preview dataframe
        preview_df = pd.DataFrame({
            'T (°C)': exp_data.get('T_degC', np.arange(len(exp_data['Delta_WF']))),
            'ΔWF (eV)': exp_data['Delta_WF'],
            'ΔE_CL (eV)': exp_data['Delta_CL']
        })
        st.dataframe(preview_df, use_container_width=True)

        # Quick statistics
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Data Points", len(exp_data['Delta_WF']))
        with col2:
            st.metric("ΔWF Range", f"{exp_data['Delta_WF'].min():.3f} to {exp_data['Delta_WF'].max():.3f} eV")
        with col3:
            # Estimate eta from slope
            lin_fit = linear_fit_eta(exp_data)
            st.metric("Estimated η", f"{lin_fit['eta_exp']:.3f}")

        # Physics validation
        validation_warnings = validate_experimental_data_physics(exp_data)
        if validation_warnings:
            st.markdown("---")
            st.subheader("⚠️ Data Quality Checks")
            for warning in validation_warnings:
                st.warning(warning)

        st.markdown("---")
        st.subheader("3. Parameter Fitting")

        with st.expander("🔧 Automatic Fitting Options", expanded=True):
            col1, col2, col3 = st.columns(3)

            with col1:
                fit_W = st.checkbox("Fit W (depletion width)", value=False,
                                   help="Optimize depletion width parameter")
                fit_eta = st.checkbox("Fit η (XPS sampling factor)", value=True,
                                     help="Optimize XPS sampling depth factor")

            with col2:
                W_guess = st.number_input("W initial (nm)", 1.0, 10.0, W_nm, 0.1)
                W_bounds_range = st.slider("W bounds (nm)", 1.0, 10.0, (2.0, 5.0), 0.1)

            with col3:
                eta_guess = st.number_input("η initial", 0.5, 0.99, 0.85, 0.01)
                eta_bounds_range = st.slider("η bounds", 0.5, 0.99, (0.7, 0.95), 0.01)

            if st.button("🎯 Run Fitting", type="primary", use_container_width=True):
                with st.spinner("Optimizing parameters..."):
                    fit_result = run_parameter_fitting(
                        exp_data=exp_data,
                        model_type=model_type,
                        fit_params={'W': fit_W, 'eta': fit_eta},
                        initial_guess={'W': W_guess, 'eta': eta_guess},
                        bounds={'W': W_bounds_range, 'eta': eta_bounds_range},
                        epsilon_r=epsilon_r,
                        m_star_ratio=m_star_ratio
                    )

                if fit_result['success']:
                    st.session_state.fit_result = fit_result
                    st.success("Fitting completed successfully!")
                    st.rerun()
                else:
                    st.error(fit_result['message'])

        # Display fitting results
        if st.session_state.fit_result is not None:
            fit_result = st.session_state.fit_result

            st.subheader("Fitting Results")

            col1, col2, col3, col4 = st.columns(4)
            with col1:
                delta_W = fit_result['W_fit'] - W_guess
                st.metric("Fitted W", f"{fit_result['W_fit']:.2f} nm",
                         delta=f"{delta_W:+.2f} nm")
            with col2:
                delta_eta = fit_result['eta_fit'] - eta_guess
                st.metric("Fitted η", f"{fit_result['eta_fit']:.3f}",
                         delta=f"{delta_eta:+.3f}")
            with col3:
                st.metric("R² (goodness)", f"{fit_result['r_squared']:.4f}")
            with col4:
                st.metric("RMSE", f"{fit_result['rmse']:.4f} eV")

            # Apply fitted parameters
            if st.button("✅ Apply Fitted Parameters to Model"):
                # This would require updating sidebar values - show info instead
                st.info(f"**Fitted parameters:**\n- W = {fit_result['W_fit']:.2f} nm\n- η ≈ λ = {fit_result['lambda_fit']:.2f} nm")

        st.markdown("---")
        st.subheader("4. Comparison Plots")

        # Chart explanation for comparison plot
        with st.expander("ℹ️ **Comparison Chart Explanation**", expanded=False):
            st.markdown("""
            **What this chart shows:** Direct comparison between experimental measurements (ΔE_CL vs ΔWF) and theoretical predictions.

            **X-axis (ΔWF):** Work function change in eV, measured by UPS or KPFM. Represents the change in surface electronic structure.

            **Y-axis (ΔE_CL):** Core level shift in eV, measured by XPS. Represents the depth-weighted average of the surface potential.

            **Data points (red):** Your experimental measurements. Each point typically represents a different annealing temperature or treatment condition.

            **Theory line (blue):** Model prediction using the relationship ΔE_CL = η × ΔWF, where η is the XPS sampling factor (typically 0.7-0.95).

            **Fitted line (green):** Best-fit line after parameter optimization. The closer to the theory line, the better the model describes your data.

            **Physical interpretation:** Both ΔE_CL and ΔWF should follow the same band bending, but XPS (ΔE_CL) samples deeper than UPS (ΔWF), hence |ΔE_CL| < |ΔWF|. The slope η depends on λ (mean free path) and θ (emission angle).
            """)

        # Main comparison plot: ΔE_CL vs ΔWF
        fig_comparison = create_comparison_CL_vs_WF_plot(
            exp_data=exp_data,
            theory_data=None,
            fit_result=st.session_state.fit_result,
            model=model,
            lambda_nm=lambda_xps,
            theta_deg=theta_xps
        )
        st.plotly_chart(fig_comparison, use_container_width=True)

        # Annealing trajectory if temperature data available
        if 'T_degC' in exp_data:
            with st.expander("ℹ️ **Annealing Trajectory Explanation**", expanded=False):
                st.markdown("""
                **What this chart shows:** Evolution of surface electronic structure during temperature-programmed annealing experiments.

                **X-axis (ΔWF):** Work function change in eV. Typically starts negative (electron-rich) and moves toward zero as temperature increases.

                **Y-axis (ΔE_CL):** Core level shift in eV. Follows ΔWF but with reduced magnitude due to depth averaging.

                **Color scale (T):** Temperature in °C. Shows the thermal treatment sequence. Typically progresses from low T (dark blue) to high T (yellow/red).

                **Physical interpretation:** During vacuum annealing, adsorbates desorb with increasing temperature, causing work function to increase (ΔWF → 0) and 2DEG density to decrease. The trajectory should follow a straight line with slope η if only band bending changes. Deviations indicate additional effects (e.g., adsorbate dipoles, surface reconstruction).

                **Use:** Helps visualize whether your experimental trajectory is consistent with the simple band bending + adsorbate model.
                """)

            fig_trajectory = create_annealing_trajectory_plot(exp_data)
            if fig_trajectory is not None:
                st.plotly_chart(fig_trajectory, use_container_width=True)

        # Residual analysis
        if st.session_state.fit_result is not None:
            st.markdown("---")
            st.subheader("5. Residual Analysis")

            with st.expander("ℹ️ **Residual Analysis Explanation**", expanded=False):
                st.markdown("""
                **What this chart shows:** Statistical diagnostic plots to assess the quality of the fit between theory and experiment.

                **Left panel - Residuals vs Fitted:**
                - **X-axis:** Fitted ΔE_CL values (eV) from the theoretical model.
                - **Y-axis:** Residuals (eV) = Experimental - Fitted values.
                - **Interpretation:** Points should scatter randomly around zero. Patterns indicate systematic deviations (e.g., model inadequacy, nonlinear effects).

                **Middle panel - Residual Histogram:**
                - **X-axis:** Residual values (eV).
                - **Y-axis:** Frequency count.
                - **Interpretation:** Should be approximately symmetric and centered at zero. Width indicates scatter magnitude (experimental noise + model error).

                **Right panel - Q-Q Plot (Normal Probability):**
                - **X-axis:** Theoretical quantiles (normal distribution).
                - **Y-axis:** Sample quantiles (your residuals).
                - **Interpretation:** Points should follow the diagonal red line if residuals are normally distributed. Deviations suggest outliers or non-Gaussian errors.

                **Good fit indicators:** Random scatter, centered histogram, linear Q-Q plot, small RMSE, R² close to 1.0.
                """)

            fig_residuals = create_residual_analysis_plot(
                exp_data=exp_data,
                fit_result=st.session_state.fit_result
            )
            st.plotly_chart(fig_residuals, use_container_width=True)

# ============================================================================
# TAB 4: PUBLICATION EXPORT
# ============================================================================

with tab4:
    st.header("📤 Publication-Quality Export")

    st.markdown("Export high-quality figures optimized for journal submission")

    # Export settings
    st.subheader("Export Settings")

    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("**Format**")
        format_type = st.radio(
            "Format type",
            ["PNG (raster)", "SVG (vector)", "PDF (document)"],
            label_visibility="collapsed"
        )

        if "PNG" in format_type:
            dpi = st.selectbox("Resolution (DPI)", [300, 600, 1200], index=0)
        else:
            dpi = 300

    with col2:
        st.markdown("**Style**")
        journal_style = st.selectbox(
            "Journal style",
            ["Default", "Nature", "Science", "ACS", "Grayscale"]
        )

        font_size = st.slider("Font size (pt)", 8, 16, 12)

    with col3:
        st.markdown("**Size**")
        size_presets = get_figure_size_presets()
        size_preset = st.selectbox(
            "Size preset",
            list(size_presets.keys())
        )

        fig_size = size_presets[size_preset]
        st.caption(f"Size: {fig_size[0]}\" × {fig_size[1]}\"")

    st.markdown("---")

    # Configure matplotlib style
    style_config = get_journal_style(journal_style)
    setup_matplotlib_style(font_size=font_size, line_width=2.0, journal_style=style_config)

    # Figure selection
    st.subheader("Select Figures to Export")

    col1, col2, col3 = st.columns(3)

    with col1:
        export_fig1 = st.checkbox("Figure 1: ns vs Φs", value=True)
    with col2:
        export_fig2 = st.checkbox("Figure 2: ΔWF vs ns", value=True)
    with col3:
        export_comparison = st.checkbox(
            "Figure 3: Comparison",
            value=st.session_state.exp_data_loaded,
            disabled=not st.session_state.exp_data_loaded
        )

    # Generate and export button
    if st.button("🎨 Generate Publication Figures", type="primary", use_container_width=True):
        with st.spinner("Generating high-quality figures..."):
            exported_files = []

            # Figure 1
            if export_fig1:
                try:
                    mpl_fig1 = create_publication_figure_1(
                        curves_data=curves_fig1,
                        exp_data=None,
                        style_config=style_config,
                        size=fig_size
                    )

                    fmt = 'png' if 'PNG' in format_type else ('svg' if 'SVG' in format_type else 'pdf')
                    img_bytes = save_figure(mpl_fig1, format_type=fmt, dpi=dpi)

                    exported_files.append(('figure1_ns_vs_Phis.' + fmt, img_bytes))
                    plt.close(mpl_fig1)
                except Exception as e:
                    st.error(f"Error generating Figure 1: {str(e)}")

            # Figure 2
            if export_fig2:
                try:
                    mpl_fig2 = create_publication_figure_2(
                        curves_data=curves_fig2,
                        exp_data=None,
                        style_config=style_config,
                        size=fig_size
                    )

                    fmt = 'png' if 'PNG' in format_type else ('svg' if 'SVG' in format_type else 'pdf')
                    img_bytes = save_figure(mpl_fig2, format_type=fmt, dpi=dpi)

                    exported_files.append(('figure2_Delta_WF_vs_ns.' + fmt, img_bytes))
                    plt.close(mpl_fig2)
                except Exception as e:
                    st.error(f"Error generating Figure 2: {str(e)}")

            # Comparison figure
            if export_comparison and st.session_state.exp_data_loaded:
                try:
                    mpl_fig3 = create_publication_comparison_figure(
                        exp_data=st.session_state.exp_data,
                        theory_data=None,
                        fit_result=st.session_state.fit_result,
                        style_config=style_config,
                        size=fig_size,
                        model=model,
                        lambda_nm=lambda_xps,
                        theta_deg=theta_xps
                    )

                    fmt = 'png' if 'PNG' in format_type else ('svg' if 'SVG' in format_type else 'pdf')
                    img_bytes = save_figure(mpl_fig3, format_type=fmt, dpi=dpi)

                    exported_files.append(('figure3_comparison.' + fmt, img_bytes))
                    plt.close(mpl_fig3)
                except Exception as e:
                    st.error(f"Error generating Figure 3: {str(e)}")

        # Display download buttons
        if exported_files:
            st.success(f"✅ Generated {len(exported_files)} figure(s)")

            st.markdown("---")
            st.subheader("Download Figures")

            for filename, img_bytes in exported_files:
                mime_type = 'image/png' if '.png' in filename else ('image/svg+xml' if '.svg' in filename else 'application/pdf')

                st.download_button(
                    label=f"⬇️ Download {filename}",
                    data=img_bytes,
                    file_name=filename,
                    mime=mime_type,
                    use_container_width=True
                )

    st.markdown("---")
    st.subheader("Data Export")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("### Export Data")

        # Prepare CSV data
        csv_data = []
        for i, (phi, ns_val) in enumerate(zip(Phi_s_range, ns_array)):
            csv_data.append({
                'model': model_type,
                'W_nm': W_nm,
                'Phi_s_eV': phi,
                'm_star_m0': m_star_ratio,
                'epsilon_r': epsilon_r,
                'T_K': temperature,
                'ns_cm2': ns_to_display(ns_val) * 1e13,
                'Delta_WF_eV': model.calculate_Delta_WF(ns_val),
            })

        df_export = pd.DataFrame(csv_data)

        # CSV download
        csv = df_export.to_csv(index=False)
        st.download_button(
            label="Download CSV",
            data=csv,
            file_name=f"2deg_data_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
            mime="text/csv",
            use_container_width=True
        )

        # JSON export with parameters
        json_data = {
            'timestamp': datetime.now().isoformat(),
            'parameters': {
                'model': model_type,
                'm_star': m_star_ratio,
                'epsilon_r': epsilon_r,
                'T_K': temperature,
                'W_nm': W_nm,
                'Phi_s': Phi_s,
                'lambda_xps_nm': lambda_xps,
                'theta_xps_deg': theta_xps,
                'Delta_Phi_dip_eV': Delta_Phi_dip
            }
        }

        import json
        json_str = json.dumps(json_data, indent=2)
        st.download_button(
            label="Download JSON Config",
            data=json_str,
            file_name=f"2deg_config_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
            mime="application/json",
            use_container_width=True
        )

    with col2:
        st.markdown("### Quick Export (Plotly)")

        # Quick export buttons for interactive figures
        st.caption("Export interactive Plotly figures (lower quality)")

        if st.button("Download Fig 1 (HTML)", use_container_width=True):
            fig1_html = fig1.to_html()
            st.download_button(
                label="⬇️ Save HTML",
                data=fig1_html,
                file_name=f"figure1_{datetime.now().strftime('%Y%m%d_%H%M%S')}.html",
                mime="text/html",
                use_container_width=True
            )

        if st.button("Download Fig 2 (HTML)", use_container_width=True):
            fig2_html = fig2.to_html()
            st.download_button(
                label="⬇️ Save HTML",
                data=fig2_html,
                file_name=f"figure2_{datetime.now().strftime('%Y%m%d_%H%M%S')}.html",
                mime="text/html",
                use_container_width=True
            )

# ============================================================================
# TAB 5: BETA FEATURES
# ============================================================================

with tab5:
    st.header("🧪 Beta Features")
    st.markdown("*Experimental features under development*")

    st.info("These features are currently under development and will be added in future versions.")

    # Feature 1: Self-Consistent S-P Diagnostic
    with st.expander("🔍 Self-Consistent S-P Diagnostic (Coming Soon)", expanded=False):
        st.markdown("""
        **Feature:** Visualize Schrödinger-Poisson self-consistent solution

        This will show:
        - Gauss law: nₛ(F) = εF/q
        - Quantum DOS: nₛ(F) from subband occupation
        - Self-consistent intersection point

        **Status:** Planned for M2 (Fang-Howard) model
        """)

    # Feature 2: Experimental Guidance Mode
    with st.expander("🎯 Experimental Guidance Mode (Coming Soon)", expanded=False):
        st.markdown("""
        **Feature:** Design experiments and get theoretical predictions

        This will provide:
        - Sample specification input (material, thickness, growth method)
        - Annealing plan (temperature range, atmosphere)
        - Predicted measurement trajectories
        - Experimental recommendations and warnings

        **Status:** Planned
        """)

    # Feature 3: Uncertainty Analysis
    with st.expander("📊 Uncertainty Propagation Analysis (Coming Soon)", expanded=False):
        st.markdown("""
        **Feature:** Propagate parameter uncertainties to predictions

        This will show:
        - Error bands on theory curves
        - Monte Carlo confidence intervals
        - Sensitivity analysis

        **Status:** Planned
        """)

    # Feature 4: Material Database
    with st.expander("📚 Material Database (Coming Soon)", expanded=False):
        st.markdown("""
        **Feature:** Pre-configured parameters for common materials

        Will include:
        - In₂O₃, ZnO, SnO₂, ITO, etc.
        - Literature values for m*, εᵣ
        - One-click parameter loading

        **Status:** Planned
        """)

    # Feature 5: Annealing Animation
    with st.expander("🎬 Annealing Trajectory Animation (Coming Soon)", expanded=False):
        st.markdown("""
        **Feature:** Animated visualization of annealing process

        This will show:
        - Dynamic evolution of ΔWF and ns with temperature
        - Trajectory path in parameter space
        - Time-lapse animation

        **Status:** Planned
        """)

# ============================================================================
# TAB 6: XPS SPECTRUM
# ============================================================================

with tab6:
    XPS_UI_LANGUAGE = "en"
    XPS_TEXT = {
        "intro": {
            "en": "Upload XPS files (`.spe` / `.csv`) to auto-convert, create spectra, annotate element peaks, and run multi-peak fitting.",
            "zh": "拖拽或上传 XPS 文件（`.spe` / `.csv`），自动转换、出图、峰位标注并执行多峰拟合。",
        },
        "workflow_title": {
            "en": "ℹ️ Workflow Guide",
            "zh": "ℹ️ 控件说明 / Workflow Guide",
        },
        "workflow_body": {
            "en": (
                "1. Upload `.spe` or `.csv`. SPE files are split into trace-level spectra automatically.  \n"
                "2. Choose analysis scope (full trace / element region / manual window).  \n"
                "3. Tune detection settings, then set the fitting range until `Peaks inside fit range` > 0.  \n"
                "4. Run fitting and export a publication-ready figure with optional fitting curves and peak labels."
            ),
            "zh": (
                "1. 上传 `.spe` 或 `.csv`，`SPE` 会自动拆分为多个 trace。  \n"
                "2. 选择分析范围（全谱 / 元素窗口 / 自定义窗口）。  \n"
                "3. 先调峰检测参数，再设定拟合区间，确保 `Peaks inside fit range` > 0。  \n"
                "4. 导出面板支持 `PNG/SVG/PDF`，可选择是否携带拟合和峰名。"
            ),
        },
        "upload_help": {
            "en": "Supports raw `.spe` source files and converted `.csv` files.",
            "zh": "支持原始 `.spe` 文件与已转换 `.csv` 文件。",
        },
        "upload_required": {
            "en": "Upload one XPS file (`.spe` or `.csv`) to start.",
            "zh": "请先上传一个 XPS 文件（.spe 或 .csv）。",
        },
        "phi_parser_note": {
            "en": "Automatically switched to the built-in PHI parser for non-standard `.spe` headers.",
            "zh": "已自动切换到内置 PHI 解析器（兼容非标准 .spe 头格式）。",
        },
        "trace_help": {
            "en": "A single `.spe` usually contains multiple element regions (for example C1s / O1s / In3d). Choose one region to analyze.",
            "zh": "同一个 .spe 里通常包含多个元素区间（如 C1s/O1s/In3d）。在这里切换当前分析谱线。",
        },
        "reference_help": {
            "en": "Reference peak library used for line annotation and fitted-peak matching.",
            "zh": "用于标注和匹配的参考元素峰位库。",
        },
        "scope_help": {
            "en": "Choose whether to analyze the full trace, one element-centered region, or a custom energy window.",
            "zh": "选择分析全谱、元素窗口或自定义能量窗口。",
        },
        "scope_points_warning": {
            "en": "Selected scope has too few points. Widen the window for stable detection/fitting.",
            "zh": "当前分析窗口点数过少，请扩大范围以保证检测与拟合稳定。",
        },
    }

    def xps_t(key: str) -> str:
        entry = XPS_TEXT.get(key, {})
        return str(entry.get(XPS_UI_LANGUAGE) or entry.get("en") or key)

    def stable_token(raw: str, max_len: int = 80) -> str:
        out = ''.join(c if c.isalnum() else '_' for c in str(raw))
        out = out.strip('_')
        if not out:
            out = "xps"
        return out[:max_len]

    def short_hash(raw: str, n: int = 12) -> str:
        return hashlib.md5(raw.encode("utf-8")).hexdigest()[:n]

    st.header("🧲 XPS Spectrum Analyzer")
    st.markdown(xps_t("intro"))

    col_left, col_right = st.columns([1.2, 2.0], gap="large")
    left_panel = col_left.container(height=860, border=False, key="xps_left_scroll_panel")
    right_panel = col_right.container(height=860, border=False, key="xps_right_scroll_panel")


    with left_panel:
        with st.expander(xps_t("workflow_title"), expanded=False):
            st.markdown(xps_t("workflow_body"))
        with st.container(border=True):
            st.subheader("Step 1: Data Import 📁")
            
            with st.expander("🛠️ Developer / AI Agent Auto-Loader", expanded=False):
                dev_path = st.text_input("Absolute path to .spe or .csv file:", key="xps_dev_path_input")
                c1, c2 = st.columns(2)
                with c1:
                    if st.button("Auto-Load File", key="xps_btn_autoload", use_container_width=True):
                        if os.path.exists(dev_path):
                            st.session_state['xps_dev_autoload_path'] = dev_path
                        else:
                            st.error("File not found.")
                with c2:
                    if st.button("Clear loaded file", key="xps_btn_clear_autoload", use_container_width=True):
                        if 'xps_dev_autoload_path' in st.session_state:
                            del st.session_state['xps_dev_autoload_path']

            dev_file_path = st.session_state.get('xps_dev_autoload_path', '')
            
            uploaded_xps = st.file_uploader(
                "Upload XPS spectrum file",
                type=['csv', 'spe'],
                key='xps_upload_any',
                help=xps_t("upload_help")
            )

    has_file = (uploaded_xps is not None) or (dev_file_path and os.path.exists(dev_file_path))

    if not has_file:
        with left_panel:
            st.info(xps_t("upload_required"))
            
    else:
        try:
            if uploaded_xps is not None:
                file_content = uploaded_xps.getvalue()
                filename = uploaded_xps.name
            else:
                with open(dev_file_path, 'rb') as f:
                    file_content = f.read()
                filename = os.path.basename(dev_file_path)

            loaded_payload = load_uploaded_xps_file(
                file_content=file_content,
                filename=filename
            )
        except Exception as exc:
            with left_panel:
                st.error(f"Failed to parse XPS file: {exc}")
        else:
            all_spectra = loaded_payload['spectra']
            file_type = loaded_payload['file_type']

            with left_panel:
                if file_type == 'spe':
                    st.success(
                        f"SPE converted successfully: {len(all_spectra)} trace(s) "
                        f"via {loaded_payload.get('method', 'auto')}."
                    )

                    log_text = (loaded_payload.get('logs') or '').strip()
                    if loaded_payload.get('method') == "local_phi_parser":
                        st.info(xps_t("phi_parser_note"))
                    if log_text:
                        with st.expander("Conversion Log", expanded=False):
                            st.code(log_text)

                    trace_labels = [
                        f"{idx + 1}. {spec.get('trace_name', f'trace_{idx + 1}')}"
                        for idx, spec in enumerate(all_spectra)
                    ]
                    selected_trace_label = st.selectbox(
                        "Select trace from converted SPE",
                        options=trace_labels,
                        key='xps_trace_select',
                        help=xps_t("trace_help")
                    )
                    selected_idx = trace_labels.index(selected_trace_label)

                    zip_bytes = build_spectra_zip_bytes(
                        spectra=all_spectra,
                        base_name=filename.rsplit('.', 1)[0]
                    )
                    st.download_button(
                        label="⬇️ Download Converted Traces (.zip)",
                        data=zip_bytes,
                        file_name=f"{filename.rsplit('.', 1)[0]}_converted_traces.zip",
                        mime="application/zip",
                        use_container_width=True,
                        key='xps_download_converted_zip'
                    )
                else:
                    selected_idx = 0

                xps_data = all_spectra[selected_idx]
                energy = xps_data['energy']
                intensity = xps_data['intensity']
                trace_name = str(xps_data.get('trace_name', filename))
                trace_key = stable_token(f"{filename}_{trace_name}", max_len=100)

                st.caption(
                    f"File: {filename} | Trace: {trace_name} | Points: {xps_data['n_points']} | "
                    f"Range: {xps_data['energy_min']:.2f} - {xps_data['energy_max']:.2f} eV"
                )

            with left_panel:
                with st.container(border=True):
                    st.subheader("Step 2: Region & Elements 🔍")
                    reference_library = get_element_reference_library()
                    ref_options = [item['label'] for item in reference_library]
                    suggested_refs = [
                        x for x in suggest_reference_labels(f"{filename}_{trace_name}") if x in ref_options
                    ]
                    selected_ref_labels = st.multiselect(
                        "Reference element peaks",
                        options=ref_options,
                        default=suggested_refs,
                        key=f'xps_reference_labels_{trace_key}',
                        help=xps_t("reference_help")
                    )
                    selected_reference_peaks = [item for item in reference_library if item['label'] in selected_ref_labels]

                    energy_min = float(np.min(energy))
                    energy_max = float(np.max(energy))

                    scope_mode = st.selectbox(
                        "Analysis scope",
                        options=["Full trace", "Element region", "Manual window"],
                        index=1 if file_type == 'spe' else 0,
                        key=f'xps_scope_mode_{trace_key}',
                        help=xps_t("scope_help")
                    )

                    scope_id = "full_trace"
                    scope_label = "Full trace"
                    scope_min = energy_min
                    scope_max = energy_max
                    selected_element_label = ""

                    if scope_mode == "Element region":
                        element_options = [item['label'] for item in selected_reference_peaks] or ref_options
                        selected_element_label = st.selectbox(
                            "Element / core level",
                            options=element_options,
                            key=f'xps_scope_element_{trace_key}',
                            help="Choose an element peak and create an element-centered fitting/export window."
                        )
                        ref_lookup = {
                            str(item['label']): float(item['energy_eV'])
                            for item in reference_library
                            if 'label' in item and 'energy_eV' in item
                        }
                        center_e = float(ref_lookup.get(selected_element_label, 0.5 * (energy_min + energy_max)))
                        half_width = st.slider(
                            "Element window half-width (eV)",
                            min_value=0.5,
                            max_value=40.0,
                            value=8.0,
                            step=0.1,
                            key=f"xps_scope_halfwidth_{trace_key}_{stable_token(selected_element_label, 32)}",
                            help="Window half-width around the selected element peak."
                        )
                        scope_min = max(energy_min, center_e - half_width)
                        scope_max = min(energy_max, center_e + half_width)
                        scope_id = f"element_{stable_token(selected_element_label, 32)}"
                        scope_label = f"{selected_element_label} ({scope_min:.2f} - {scope_max:.2f} eV)"

                    elif scope_mode == "Manual window":
                        manual_window = st.slider(
                            "Manual analysis window (eV)",
                            min_value=energy_min,
                            max_value=energy_max,
                            value=(energy_min, energy_max),
                            step=0.05,
                            key=f'xps_scope_manual_window_{trace_key}',
                            help="Custom region for local fitting/export."
                        )
                        scope_min = float(min(manual_window))
                        scope_max = float(max(manual_window))
                        scope_id = "manual_window"
                        scope_label = f"Manual window ({scope_min:.2f} - {scope_max:.2f} eV)"

                    scope_mask = (energy >= scope_min) & (energy <= scope_max)
                    scope_points = int(np.sum(scope_mask))
                    scope_ready_detection = scope_points >= 10
                    scope_ready_fit = scope_points >= 20

                    analysis_energy = energy[scope_mask] if scope_points >= 2 else energy.copy()
                    analysis_intensity = intensity[scope_mask] if scope_points >= 2 else intensity.copy()
                    if analysis_energy.size == 0 or analysis_intensity.size == 0:
                        analysis_energy = energy.copy()
                        analysis_intensity = intensity.copy()
                        scope_ready_detection = False
                        scope_ready_fit = False

                    scope_min_used = float(np.min(analysis_energy))
                    scope_max_used = float(np.max(analysis_energy))
                    st.caption(
                        f"Active scope: {scope_label} | Points: {analysis_energy.size}"
                    )
                    if not scope_ready_detection:
                        st.warning(xps_t("scope_points_warning"))

            analysis_control_token = f"{filename}::{trace_name}::{scope_id}"
            analysis_control_key = short_hash(analysis_control_token, 12)
            analysis_fit_token = f"{analysis_control_token}::{scope_min_used:.3f}:{scope_max_used:.3f}"
            analysis_fit_key = short_hash(analysis_fit_token, 18)
            scope_slug = stable_token(selected_element_label or scope_id, max_len=32)

            with left_panel:
                with st.container(border=True):
                    st.subheader("Step 3: Background & Fitting ⚙️")
                    sl1, sl2 = st.columns(2)
                    with sl1:
                        bg_edge_pct = st.slider(
                            "Background edge (%)",
                            5, 30, 12, 1,
                            key=f'xps_bg_edge_pct_{analysis_control_key}',
                            help="Edge fraction for linear baseline."
                        )
                        smoothing_window = st.slider(
                            "Smoothing window",
                            5, 51, 5, 2,
                            key=f'xps_smooth_window_{analysis_control_key}',
                            help="Savitzky-Golay filtering window."
                        )
                        min_peak_distance = st.slider(
                            "Min peak distance",
                            1, 80, 10, 1,
                            key=f'xps_peak_distance_{analysis_control_key}',
                            help="Minimum points between detected peaks."
                        )
                    with sl2:
                        prominence_pct = st.slider(
                            "Peak prominence (%)",
                            1, 40, 4, 1,
                            key=f'xps_prominence_pct_{analysis_control_key}',
                            help="Peak prominence threshold."
                        )
                        max_peaks = st.slider(
                            "Max peaks",
                            1, 8, 4, 1,
                            key=f'xps_max_peaks_{analysis_control_key}',
                            help="Maximum detected peaks."
                        )
                        match_tolerance = st.slider(
                            "Element match (eV)",
                            0.2, 3.0, 1.2, 0.1,
                            key=f'xps_match_tolerance_{analysis_control_key}',
                            help="Tolerance to match references."
                        )

                    fit_profile = st.selectbox(
                        "Peak profile",
                        ["Gaussian", "Lorentzian"],
                        key=f'xps_peak_profile_{analysis_control_key}',
                        help="Line-shape model."
                    )

                    if scope_ready_detection:
                        fit_range_key = f'xps_fit_range_{analysis_control_key}'
                        fit_range_default = (scope_min_used, scope_max_used)
                        previous_range = st.session_state.get(fit_range_key)
                        if isinstance(previous_range, (tuple, list)) and len(previous_range) == 2:
                            prev_low = float(np.clip(min(previous_range), scope_min_used, scope_max_used))
                            prev_high = float(np.clip(max(previous_range), scope_min_used, scope_max_used))
                            if prev_high > prev_low:
                                fit_range_default = (prev_low, prev_high)

                        fit_range = st.slider(
                            "Fitting range (eV)",
                            min_value=scope_min_used,
                            max_value=scope_max_used,
                            value=fit_range_default,
                            step=0.05,
                            key=fit_range_key,
                            help="Peaks inside this range are fitted."
                        )
                    else:
                        fit_range = (scope_min_used, scope_max_used)

                    if scope_ready_detection:
                        detection = detect_spectrum_peaks(
                            energy=analysis_energy,
                            intensity=analysis_intensity,
                            edge_fraction=bg_edge_pct / 100.0,
                            smoothing_window=smoothing_window,
                            prominence_ratio=prominence_pct / 100.0,
                            min_distance_points=min_peak_distance,
                            max_peaks=max_peaks
                        )
                        detected_peaks = detection['peaks']
                    else:
                        detected_peaks = pd.DataFrame(
                            columns=["index", "center_eV", "height", "prominence", "raw_intensity"]
                        )

                    fit_min = float(min(fit_range))
                    fit_max = float(max(fit_range))
                    fit_range_peaks = detected_peaks[
                        (detected_peaks['center_eV'] >= fit_min) & (detected_peaks['center_eV'] <= fit_max)
                    ].copy()

                    if detected_peaks.empty:
                        st.warning("No peaks detected. Lower prominence or increase smoothing.")
                    else:
                        st.info(f"Detected {len(detected_peaks)} peaks. In fit range: {len(fit_range_peaks)}")

                    run_fit = st.button(
                        "🎯 Run Peak Fitting",
                        type="primary",
                        use_container_width=True,
                        key=f'xps_run_fit_{analysis_control_key}',
                        disabled=(fit_range_peaks.empty or not scope_ready_fit)
                    )

                    if run_fit:
                        try:
                            fit_output = fit_spectrum_peaks(
                                energy=analysis_energy,
                                intensity=analysis_intensity,
                                initial_centers=fit_range_peaks['center_eV'].to_numpy(dtype=float),
                                profile=fit_profile,
                                fit_range=fit_range,
                                edge_fraction=bg_edge_pct / 100.0
                            )
                            fit_output['peak_table'] = match_fitted_peaks(
                                peak_table=fit_output['peak_table'],
                                selected_reference_peaks=selected_reference_peaks,
                                tolerance_eV=match_tolerance
                            )
                            st.session_state.xps_fit_store[analysis_fit_key] = {
                                "fit_result": fit_output,
                                "file_name": filename,
                                "trace_name": trace_name,
                                "scope_label": scope_label,
                                "scope_id": scope_id,
                                "scope_min": scope_min_used,
                                "scope_max": scope_max_used,
                                "fit_min": fit_min,
                                "fit_max": fit_max,
                                "updated_at": datetime.utcnow().isoformat(timespec="seconds"),
                            }
                            st.success("Peak fitting completed!")
                        except Exception as exc:
                            st.session_state.xps_fit_store.pop(analysis_fit_key, None)
                            st.error(f"Peak fitting failed: {exc}")

            current_xps_entry = st.session_state.xps_fit_store.get(analysis_fit_key)
            current_xps_fit = current_xps_entry["fit_result"] if current_xps_entry else None

            with right_panel:
                with st.container():
                    st.subheader("Spectrum Visualization", divider="gray")
                    with st.expander("👁️ Display Options", expanded=False):
                        dc1, dc2, dc3, dc4 = st.columns(4)
                        with dc1:
                            show_detected_markers = st.checkbox("Detected peaks", value=True, key=f'xps_show_detected_{analysis_control_key}')
                            show_ref_lines = st.checkbox("Reference lines", value=True, key=f'xps_show_ref_lines_{analysis_control_key}')
                        with dc2:
                            show_fit_curve = st.checkbox("Fitting curve", value=current_xps_fit is not None, key=f'xps_show_fit_curve_{analysis_control_key}')
                            show_fit_components = st.checkbox("Fit components", value=current_xps_fit is not None, key=f'xps_show_fit_components_{analysis_control_key}')
                        with dc3:
                            show_baseline_line = st.checkbox("Baseline", value=current_xps_fit is not None, key=f'xps_show_baseline_{analysis_control_key}')
                            show_ref_labels = st.checkbox("Reference labels", value=True, key=f'xps_show_ref_labels_{analysis_control_key}')
                        with dc4:
                            show_peak_name_labels = st.checkbox("Peak names", value=False, key=f'xps_show_peak_names_{analysis_control_key}')
                            focus_to_fit_range = st.checkbox("Zoom to fit range", value=(scope_mode != "Full trace"), key=f'xps_zoom_fit_range_{analysis_control_key}')

                        style_col1, style_col2 = st.columns(2)
                        with style_col1:
                            font_size_level = st.select_slider("Chart text size", options=["Small", "Medium", "Large", "XL"], value="Large", key=f'xps_font_size_{analysis_control_key}')
                        with style_col2:
                            line_weight_level = st.select_slider("Spectrum line weight", options=["Fine", "Normal", "Bold"], value="Fine", key=f'xps_line_weight_{analysis_control_key}')

                    font_size_map = {"Small": 12, "Medium": 14, "Large": 16, "XL": 18}
                    line_weight_map = {"Fine": 0.8, "Normal": 1.0, "Bold": 1.25}
                    base_font_size = font_size_map[font_size_level]
                    line_width_scale = line_weight_map[line_weight_level]

                    figure_title = f"XPS Spectrum: {filename} [{trace_name}]"
                    if scope_mode != "Full trace":
                        figure_title = f"{figure_title} - {scope_label}"

                    fig_xps = create_xps_spectrum_figure(
                        energy=analysis_energy,
                        intensity=analysis_intensity,
                        detected_peaks=detected_peaks,
                        fit_result=current_xps_fit if show_fit_curve else None,
                        reference_peaks=selected_reference_peaks,
                        title=figure_title,
                        show_detected_peaks=show_detected_markers,
                        show_fit_curve=show_fit_curve,
                        show_baseline=show_baseline_line,
                        show_components=show_fit_components,
                        show_reference_lines=show_ref_lines,
                        show_reference_labels=show_ref_labels,
                        show_peak_name_labels=show_peak_name_labels,
                        include_legend=True,
                        base_font_size=base_font_size,
                        title_font_size=base_font_size + 8,
                        figure_height=650,
                        focus_range=fit_range if focus_to_fit_range else None,
                        line_width_scale=line_width_scale,
                    )
                    st.plotly_chart(fig_xps, use_container_width=True)

                if current_xps_fit is not None:
                    st.subheader("Fitted Peak Parameters")
                    colm1, colm2, colm3 = st.columns(3)
                    with colm1:
                        st.metric("Profile", current_xps_fit['profile'])
                    with colm2:
                        st.metric("R²", f"{current_xps_fit['r_squared']:.4f}")
                    with colm3:
                        st.metric("RMSE", f"{current_xps_fit['rmse']:.3f}")

                    peak_table = current_xps_fit['peak_table'].copy()
                    st.dataframe(peak_table, use_container_width=True)

            with left_panel:
                with st.container(border=True):
                    st.subheader("Step 4: Export & Downloads 📥")
                    if current_xps_fit is not None:
                        peak_csv = peak_table.to_csv(index=False).encode('utf-8')
                        download_name = f"{filename.rsplit('.', 1)[0]}_{stable_token(trace_name, 40)}_{scope_slug}_peak_fit.csv"
                        st.download_button(
                            label="⬇️ Download Peak Fit CSV",
                            data=peak_csv, file_name=download_name, mime="text/csv",
                            use_container_width=True, key=f'xps_download_peak_fit_{analysis_fit_key}'
                        )

                    export_format = st.selectbox("Format", ["PNG", "SVG", "PDF"], key=f'xps_export_format_{analysis_control_key}')
                    export_size_preset = st.selectbox("Size preset", ["Presentation (1900x1100)", "Journal Wide (2400x1400)", "Square (1600x1600)"], key=f'xps_export_size_{analysis_control_key}')
                    export_font_level = st.select_slider("Text size", ["Medium", "Large", "XL"], value="Large", key=f'xps_export_font_{analysis_control_key}')
                    export_dpi = st.selectbox("PNG DPI", [150, 300, 600], index=1, key=f'xps_export_dpi_{analysis_control_key}')

                    size_map = {"Presentation (1900x1100)": (1900, 1100), "Journal Wide (2400x1400)": (2400, 1400), "Square (1600x1600)": (1600, 1600)}
                    export_width, export_height = size_map[export_size_preset]
                    export_font_map = {"Medium": 14, "Large": 16, "XL": 18}
                    export_base_font = export_font_map[export_font_level]

                    export_line_weight = st.select_slider("Export line weight", ["Fine", "Normal", "Bold"], value=line_weight_level, key=f'xps_export_line_weight_{analysis_control_key}')
                    export_zoom_fit = st.checkbox("Export zoom to fit range", value=(scope_mode != "Full trace"), key=f'xps_export_zoom_fit_{analysis_control_key}')
                    export_line_width_scale = line_weight_map[export_line_weight]

                    with st.expander("Advanced Export Options", expanded=False):
                        eo1, eo2 = st.columns(2)
                        with eo1:
                            export_include_fit = st.checkbox("Include fitting", value=current_xps_fit is not None, disabled=current_xps_fit is None, key=f'xps_exp_fit_{analysis_control_key}')
                            export_include_detected = st.checkbox("Include detected peaks", value=True, key=f'xps_exp_det_{analysis_control_key}')
                            export_include_ref_lines = st.checkbox("Include ref lines", value=True, key=f'xps_exp_refl_{analysis_control_key}')
                            export_include_legend = st.checkbox("Include legend", value=True, key=f'xps_exp_leg_{analysis_control_key}')
                        with eo2:
                            export_include_peak_names = st.checkbox("Include peak names", value=False, disabled=(current_xps_fit is None or not export_include_fit), key=f'xps_exp_pkn_{analysis_control_key}')
                            export_include_components = st.checkbox("Include fit components", value=True, disabled=(current_xps_fit is None or not export_include_fit), key=f'xps_exp_comp_{analysis_control_key}')
                            export_include_ref_labels = st.checkbox("Include ref labels", value=True, key=f'xps_exp_reflb_{analysis_control_key}')
                            export_include_baseline = st.checkbox("Include baseline", value=True, disabled=(current_xps_fit is None or not export_include_fit), key=f'xps_exp_base_{analysis_control_key}')

                        export_title = st.text_input("Export title", value=f"{filename.rsplit('.', 1)[0]} - {trace_name} - {scope_label}", key=f'xps_export_title_{analysis_control_key}')

                    export_bytes_key = f'xps_export_bytes_{analysis_fit_key}'
                    export_name_key = f'xps_export_name_{analysis_fit_key}'
                    if export_bytes_key not in st.session_state:
                        st.session_state[export_bytes_key] = None
                    if export_name_key not in st.session_state:
                        st.session_state[export_name_key] = ""

                    if st.button("🎨 Prepare Export Figure", type="primary", use_container_width=True, key=f'xps_prepare_export_{analysis_control_key}'):
                        try:
                            export_fig = create_xps_spectrum_figure(
                                energy=analysis_energy, intensity=analysis_intensity,
                                detected_peaks=detected_peaks if export_include_detected else None,
                                fit_result=current_xps_fit if export_include_fit else None,
                                reference_peaks=selected_reference_peaks,
                                title=export_title.strip(),
                                show_detected_peaks=export_include_detected, show_fit_curve=export_include_fit,
                                show_baseline=export_include_fit and export_include_baseline,
                                show_components=export_include_fit and export_include_components,
                                show_reference_lines=export_include_ref_lines, show_reference_labels=export_include_ref_labels,
                                show_peak_name_labels=export_include_fit and export_include_peak_names,
                                include_legend=export_include_legend, base_font_size=export_base_font,
                                title_font_size=export_base_font + 10, figure_height=max(700, int(export_height * 0.62)),
                                focus_range=fit_range if export_zoom_fit else None, line_width_scale=export_line_width_scale,
                            )
                            img_bytes = export_xps_figure_bytes(
                                fig=export_fig, export_format=export_format.lower(),
                                width_px=export_width, height_px=export_height, dpi=int(export_dpi),
                            )
                            export_filename = f"{filename.rsplit('.', 1)[0]}_{stable_token(trace_name, 40)}_{scope_slug}_figure.{export_format.lower()}"
                            st.session_state[export_bytes_key] = img_bytes
                            st.session_state[export_name_key] = export_filename
                            st.success("Export figure prepared.")
                        except Exception as exc:
                            st.error(f"Figure export failed: {exc}")

                    if st.session_state.get(export_bytes_key):
                        mime_map = {"png": "image/png", "svg": "image/svg+xml", "pdf": "application/pdf"}
                        file_ext = st.session_state[export_name_key].split(".")[-1].lower()
                        st.download_button(
                            label=f"⬇️ Download {st.session_state[export_name_key]}",
                            data=st.session_state[export_bytes_key], file_name=st.session_state[export_name_key],
                            mime=mime_map.get(file_ext, "application/octet-stream"), use_container_width=True,
                            key=f'xps_download_figure_{analysis_fit_key}'
                        )

            with left_panel:
                if file_type == 'spe' and len(all_spectra) > 1:
                    st.markdown("---")
                    if scope_mode != "Full trace":
                        st.info("Batch fitting is enabled in Full trace mode.")
                    else:
                        if st.button("⚡ Run Batch Fit For All Traces", type="primary", use_container_width=True, key='xps_run_batch_fit'):
                            batch_rows = []
                            issues = []
                            for spec in all_spectra:
                                spec_trace_name = str(spec.get('trace_name', 'unknown_trace'))
                                try:
                                    spec_detection = detect_spectrum_peaks(
                                        energy=spec['energy'], intensity=spec['intensity'],
                                        edge_fraction=bg_edge_pct / 100.0, smoothing_window=smoothing_window,
                                        prominence_ratio=prominence_pct / 100.0, min_distance_points=min_peak_distance, max_peaks=max_peaks
                                    )
                                    spec_peaks = spec_detection['peaks']
                                    if spec_peaks.empty:
                                        issues.append(f"{spec_trace_name}: no peaks")
                                        continue

                                    spec_fit_min = max(float(np.min(spec['energy'])), float(min(fit_range)))
                                    spec_fit_max = min(float(np.max(spec['energy'])), float(max(fit_range)))
                                    if spec_fit_min >= spec_fit_max: continue

                                    spec_fit_peaks = spec_peaks[(spec_peaks['center_eV'] >= spec_fit_min) & (spec_peaks['center_eV'] <= spec_fit_max)]
                                    if spec_fit_peaks.empty: continue

                                    spec_fit = fit_spectrum_peaks(
                                        energy=spec['energy'], intensity=spec['intensity'],
                                        initial_centers=spec_fit_peaks['center_eV'].to_numpy(dtype=float),
                                        profile=fit_profile, fit_range=(spec_fit_min, spec_fit_max), edge_fraction=bg_edge_pct / 100.0
                                    )
                                    spec_table = match_fitted_peaks(peak_table=spec_fit['peak_table'], selected_reference_peaks=selected_reference_peaks, tolerance_eV=match_tolerance)
                                    spec_table.insert(0, 'trace', spec_trace_name)
                                    spec_table.insert(1, 'r_squared', spec_fit['r_squared'])
                                    spec_table.insert(2, 'rmse', spec_fit['rmse'])
                                    spec_table.insert(3, 'profile', fit_profile)
                                    batch_rows.append(spec_table)
                                except Exception as exc:
                                    issues.append(f"{spec_trace_name}: {exc}")

                            if batch_rows:
                                batch_df = pd.concat(batch_rows, ignore_index=True)
                                st.session_state.xps_batch_fit_result = {'table': batch_df, 'issues': issues}
                                st.session_state.xps_batch_fit_file = filename
                                st.success(f"Batch fit completed: {batch_df['trace'].nunique()} trace(s).")
                            else:
                                st.session_state.xps_batch_fit_result = None
                                st.session_state.xps_batch_fit_file = ""
                                st.error("Batch fit failed.")

                        if st.session_state.xps_batch_fit_file == filename and st.session_state.xps_batch_fit_result is not None:
                            batch_result = st.session_state.xps_batch_fit_result
                            st.dataframe(batch_result['table'], use_container_width=True)
                            if batch_result.get('issues'):
                                with st.expander("Batch Fit Issues"):
                                    for issue in batch_result['issues']:
                                        st.warning(issue)
                            batch_csv = batch_result['table'].to_csv(index=False).encode('utf-8')
                            st.download_button("⬇️ Download Batch Fit Summary", data=batch_csv, file_name=f"{filename.rsplit('.', 1)[0]}_batch_peak_fit.csv", mime="text/csv", use_container_width=True, key='xps_download_batch_fit')

# ============================================================================
# TAB 7: UPS SPECTRUM
# ============================================================================

with tab7:

    UPS_UI_LANGUAGE = "en"

    UPS_TEXT = {
        "intro": {
            "en": (
                "Upload a UPS spectrum (CSV) to determine the **work function** via the "
                "secondary electron cutoff (SEC) and valence band edge analysis."
            ),
            "zh": (
                "上传UPS光谱（CSV）文件，通过二次电子截止（SEC）和价带边缘分析确定**功函数**。"
            ),
        },
        "upload_required": {
            "en": "⬆️ Upload a UPS spectrum CSV to begin analysis.",
            "zh": "⬆️ 请上传UPS光谱CSV文件以开始分析。",
        },
        "upload_help": {
            "en": "CSV with energy (binding or kinetic) and intensity columns.",
            "zh": "包含能量（结合能或动能）和强度列的CSV文件。",
        },
        "workflow_title": {
            "en": "📋 UPS Analysis Workflow",
            "zh": "📋 UPS分析工作流程",
        },
        "workflow_body": {
            "en": (
                "1. **Upload** a UPS spectrum CSV file\n"
                "2. **Select** photon source and energy scale\n"
                "3. **Adjust** SEC cutoff and VBM edge regions\n"
                "4. **Run** work function calculation\n"
                "5. **Export** results and figures"
            ),
            "zh": (
                "1. **上传** UPS光谱CSV文件\n"
                "2. **选择** 光子源和能量标度\n"
                "3. **调整** SEC截止和VBM边缘区域\n"
                "4. **运行** 功函数计算\n"
                "5. **导出** 结果和图形"
            ),
        },
    }

    def ups_t(key: str) -> str:
        entry = UPS_TEXT.get(key, {})
        return str(entry.get(UPS_UI_LANGUAGE) or entry.get("en") or key)

    st.header("🔮 UPS Spectrum Analyzer")
    st.markdown(ups_t("intro"))

    col_left, col_right = st.columns([1.2, 2.0], gap="large")
    ups_left = col_left.container(height=860, border=False, key="ups_left_scroll_panel")
    ups_right = col_right.container(height=860, border=False, key="ups_right_scroll_panel")

    # ── Step 1: Data Import ──────────────────────────────────────────────
    with ups_left:
        with st.expander(ups_t("workflow_title"), expanded=False):
            st.markdown(ups_t("workflow_body"))

        with st.container(border=True):
            st.subheader("Step 1: Data Import 📁")

            with st.expander("🛠️ Developer / AI Agent Auto-Loader", expanded=False):
                ups_dev_path = st.text_input(
                    "Absolute path to UPS .csv file:",
                    key="ups_dev_path_input",
                )
                uc1, uc2 = st.columns(2)
                with uc1:
                    if st.button("Auto-Load File", key="ups_btn_autoload", use_container_width=True):
                        if os.path.exists(ups_dev_path):
                            st.session_state['ups_dev_autoload_path'] = ups_dev_path
                        else:
                            st.error("File not found.")
                with uc2:
                    if st.button("Clear loaded file", key="ups_btn_clear_autoload", use_container_width=True):
                        st.session_state.pop('ups_dev_autoload_path', None)

            ups_dev_file = st.session_state.get('ups_dev_autoload_path', '')

            uploaded_ups = st.file_uploader(
                "Upload UPS spectrum file",
                type=['csv'],
                key='ups_upload_csv',
                help=ups_t("upload_help"),
            )

    ups_has_file = (uploaded_ups is not None) or (ups_dev_file and os.path.exists(ups_dev_file))

    if not ups_has_file:
        with ups_left:
            st.info(ups_t("upload_required"))
    else:
        try:
            if uploaded_ups is not None:
                ups_file_content = uploaded_ups.getvalue()
                ups_filename = uploaded_ups.name
            else:
                with open(ups_dev_file, 'rb') as f:
                    ups_file_content = f.read()
                ups_filename = os.path.basename(ups_dev_file)

            ups_payload = load_uploaded_ups_file(
                file_content=ups_file_content,
                filename=ups_filename,
            )
        except Exception as exc:
            with ups_left:
                st.error(f"Failed to load UPS file: {exc}")
        else:
            ups_data = ups_payload['spectra'][0]
            ups_energy = ups_data['energy']
            ups_intensity = ups_data['intensity']
            ups_trace_name = str(ups_data.get('trace_name', ups_filename))

            with ups_left:
                st.caption(
                    f"File: {ups_filename} | Points: {ups_data['n_points']} | "
                    f"Range: {ups_data['energy_min']:.2f} – {ups_data['energy_max']:.2f} eV"
                )

            # ── Step 2: Source & Energy Scale ────────────────────────────
            with ups_left:
                with st.container(border=True):
                    st.subheader("Step 2: Photon Source & Scale ☀️")

                    sc1, sc2 = st.columns(2)
                    with sc1:
                        source_label = st.selectbox(
                            "Photon source",
                            options=list(PHOTON_SOURCES.keys()),
                            index=0,
                            key="ups_photon_source",
                            help="UV photon source energy.",
                        )
                    with sc2:
                        energy_mode = st.selectbox(
                            "Energy scale",
                            options=["Binding Energy", "Kinetic Energy"],
                            index=0,
                            key="ups_energy_mode",
                            help="Energy axis convention of your data.",
                        )

                    if source_label == "Custom":
                        photon_energy = st.number_input(
                            "Custom photon energy (eV)",
                            min_value=1.0, max_value=200.0,
                            value=21.22, step=0.01,
                            key="ups_custom_hv",
                        )
                    else:
                        photon_energy = PHOTON_SOURCES[source_label]
                        st.caption(f"hν = {photon_energy:.2f} eV")

                    energy_mode_key = "binding" if "Binding" in energy_mode else "kinetic"

            # ── Step 3: Edge Regions ─────────────────────────────────────
            e_min = float(np.min(ups_energy))
            e_max = float(np.max(ups_energy))
            e_span = e_max - e_min

            with ups_left:
                with st.container(border=True):
                    st.subheader("Step 3: Edge Regions 🔍")

                    st.markdown("**Secondary Electron Cutoff (SEC)**")
                    if energy_mode_key == "binding":
                        sec_default_min = e_min
                        sec_default_max = min(e_min + 0.30 * e_span, e_max)
                    else:
                        sec_default_min = e_min
                        sec_default_max = min(e_min + 0.30 * e_span, e_max)

                    cutoff_range = st.slider(
                        "SEC region (eV)",
                        min_value=e_min, max_value=e_max,
                        value=(sec_default_min, sec_default_max),
                        step=0.05,
                        key="ups_cutoff_range",
                        help="Select the energy region containing the secondary electron cutoff edge.",
                    )

                    st.markdown("**Valence Band / Fermi Edge**")
                    if energy_mode_key == "binding":
                        vb_default_max = e_max
                        vb_default_min = max(e_max - 0.30 * e_span, e_min)
                    else:
                        vb_default_max = e_max
                        vb_default_min = max(e_max - 0.30 * e_span, e_min)

                    vbm_range = st.slider(
                        "VBM / Fermi region (eV)",
                        min_value=e_min, max_value=e_max,
                        value=(vb_default_min, vb_default_max),
                        step=0.05,
                        key="ups_vbm_range",
                        help="Select the energy region containing the valence band edge or Fermi level.",
                    )

                    smoothing = st.slider(
                        "Smoothing window",
                        5, 51, 11, 2,
                        key="ups_smooth_window",
                        help="Savitzky-Golay smoothing for edge detection.",
                    )

            # ── Step 4: Run Analysis ─────────────────────────────────────
            with ups_left:
                with st.container(border=True):
                    st.subheader("Step 4: Work Function Analysis ⚙️")

                    run_ups = st.button(
                        "🎯 Calculate Work Function",
                        type="primary",
                        use_container_width=True,
                        key="ups_run_analysis",
                    )

                    if run_ups:
                        try:
                            cutoff_res = detect_secondary_cutoff(
                                energy=ups_energy,
                                intensity=ups_intensity,
                                cutoff_region=tuple(cutoff_range),
                                smoothing_window=smoothing,
                            )
                            vbm_res = detect_valence_band_edge(
                                energy=ups_energy,
                                intensity=ups_intensity,
                                vb_region=tuple(vbm_range),
                                smoothing_window=smoothing,
                            )
                            wf_res = calculate_work_function(
                                photon_energy=photon_energy,
                                cutoff_eV=cutoff_res['cutoff_eV'],
                                vbm_eV=vbm_res['vbm_eV'],
                                energy_mode=energy_mode_key,
                            )
                            st.session_state.ups_cutoff_result = cutoff_res
                            st.session_state.ups_vbm_result = vbm_res
                            st.session_state.ups_wf_result = wf_res
                            st.success(f"Work function: **{wf_res['work_function_eV']:.3f} eV**")
                        except Exception as exc:
                            st.session_state.ups_cutoff_result = None
                            st.session_state.ups_vbm_result = None
                            st.session_state.ups_wf_result = None
                            st.error(f"Analysis failed: {exc}")

                    # Show stored results summary
                    if st.session_state.ups_wf_result is not None:
                        wfr = st.session_state.ups_wf_result
                        m1, m2, m3 = st.columns(3)
                        with m1:
                            st.metric("Φ (eV)", f"{wfr['work_function_eV']:.3f}")
                        with m2:
                            st.metric("Width (eV)", f"{wfr['spectrum_width_eV']:.3f}")
                        with m3:
                            st.metric("IP (eV)", f"{wfr['ionization_potential_eV']:.3f}")

                        st.caption(
                            f"SEC cutoff = {wfr['cutoff_eV']:.3f} eV  |  "
                            f"VBM = {wfr['vbm_eV']:.3f} eV  |  "
                            f"hν = {wfr['photon_energy_eV']:.2f} eV"
                        )

            # ── Visualization (right panel) ──────────────────────────────
            with ups_right:
                with st.container():
                    st.subheader("Spectrum Visualization", divider="gray")

                    with st.expander("👁️ Display Options", expanded=False):
                        dc1, dc2, dc3 = st.columns(3)
                        with dc1:
                            ups_show_cutoff = st.checkbox("SEC cutoff fit", value=True, key="ups_show_cutoff_fit")
                        with dc2:
                            ups_show_vbm = st.checkbox("VBM/Fermi fit", value=True, key="ups_show_vbm_fit")
                        with dc3:
                            ups_show_annot = st.checkbox("Annotations", value=True, key="ups_show_annotations")

                        s1, s2 = st.columns(2)
                        with s1:
                            ups_font_level = st.select_slider(
                                "Chart text size",
                                options=["Small", "Medium", "Large", "XL"],
                                value="Large", key="ups_font_size",
                            )
                        with s2:
                            ups_lw_level = st.select_slider(
                                "Line weight",
                                options=["Fine", "Normal", "Bold"],
                                value="Fine", key="ups_line_weight",
                            )

                    ups_font_map = {"Small": 12, "Medium": 14, "Large": 16, "XL": 18}
                    ups_lw_map = {"Fine": 0.8, "Normal": 1.0, "Bold": 1.25}
                    ups_base_font = ups_font_map[ups_font_level]
                    ups_line_ws = ups_lw_map[ups_lw_level]

                    fig_ups = create_ups_spectrum_figure(
                        energy=ups_energy,
                        intensity=ups_intensity,
                        cutoff_result=st.session_state.ups_cutoff_result if ups_show_cutoff else None,
                        vbm_result=st.session_state.ups_vbm_result if ups_show_vbm else None,
                        wf_result=st.session_state.ups_wf_result if ups_show_annot else None,
                        title=f"UPS Spectrum: {ups_filename}",
                        show_cutoff_fit=ups_show_cutoff,
                        show_vbm_fit=ups_show_vbm,
                        show_annotations=ups_show_annot,
                        energy_mode=energy_mode_key,
                        base_font_size=ups_base_font,
                        title_font_size=ups_base_font + 8,
                        figure_height=650,
                        line_width_scale=ups_line_ws,
                    )
                    st.plotly_chart(fig_ups, use_container_width=True)

                # Results table
                if st.session_state.ups_wf_result is not None:
                    wfr = st.session_state.ups_wf_result
                    st.subheader("Analysis Results")
                    results_df = pd.DataFrame([{
                        "Parameter": "Work Function (Φ)",
                        "Value (eV)": f"{wfr['work_function_eV']:.4f}",
                    }, {
                        "Parameter": "SEC Cutoff",
                        "Value (eV)": f"{wfr['cutoff_eV']:.4f}",
                    }, {
                        "Parameter": "VBM / Fermi Edge",
                        "Value (eV)": f"{wfr['vbm_eV']:.4f}",
                    }, {
                        "Parameter": "Spectrum Width",
                        "Value (eV)": f"{wfr['spectrum_width_eV']:.4f}",
                    }, {
                        "Parameter": "Ionization Potential",
                        "Value (eV)": f"{wfr['ionization_potential_eV']:.4f}",
                    }, {
                        "Parameter": "Photon Energy (hν)",
                        "Value (eV)": f"{wfr['photon_energy_eV']:.2f}",
                    }])
                    st.dataframe(results_df, use_container_width=True, hide_index=True)

            # ── Step 5: Export ────────────────────────────────────────────
            with ups_left:
                with st.container(border=True):
                    st.subheader("Step 5: Export & Downloads 📥")

                    if st.session_state.ups_wf_result is not None:
                        wfr = st.session_state.ups_wf_result
                        csv_data = (
                            "parameter,value_eV\n"
                            f"work_function,{wfr['work_function_eV']:.6f}\n"
                            f"sec_cutoff,{wfr['cutoff_eV']:.6f}\n"
                            f"vbm_fermi,{wfr['vbm_eV']:.6f}\n"
                            f"spectrum_width,{wfr['spectrum_width_eV']:.6f}\n"
                            f"ionization_potential,{wfr['ionization_potential_eV']:.6f}\n"
                            f"photon_energy,{wfr['photon_energy_eV']:.6f}\n"
                        ).encode('utf-8')
                        st.download_button(
                            label="⬇️ Download Results CSV",
                            data=csv_data,
                            file_name=f"{ups_filename.rsplit('.', 1)[0]}_ups_results.csv",
                            mime="text/csv",
                            use_container_width=True,
                            key="ups_download_results",
                        )

                    ups_export_fmt = st.selectbox("Format", ["PNG", "SVG", "PDF"], key="ups_export_format")
                    ups_export_preset = st.selectbox(
                        "Size preset",
                        ["Presentation (1900x1100)", "Journal Wide (2400x1400)", "Square (1600x1600)"],
                        key="ups_export_size",
                    )
                    ups_export_dpi = st.selectbox("PNG DPI", [150, 300, 600], index=1, key="ups_export_dpi")

                    ups_size_map = {
                        "Presentation (1900x1100)": (1900, 1100),
                        "Journal Wide (2400x1400)": (2400, 1400),
                        "Square (1600x1600)": (1600, 1600),
                    }
                    ups_ew, ups_eh = ups_size_map[ups_export_preset]

                    if 'ups_export_bytes' not in st.session_state:
                        st.session_state.ups_export_bytes = None
                    if 'ups_export_name' not in st.session_state:
                        st.session_state.ups_export_name = ""

                    if st.button("🎨 Prepare Export Figure", type="primary", use_container_width=True, key="ups_prepare_export"):
                        try:
                            export_fig = create_ups_spectrum_figure(
                                energy=ups_energy, intensity=ups_intensity,
                                cutoff_result=st.session_state.ups_cutoff_result,
                                vbm_result=st.session_state.ups_vbm_result,
                                wf_result=st.session_state.ups_wf_result,
                                title=f"UPS: {ups_filename.rsplit('.', 1)[0]}",
                                show_cutoff_fit=True, show_vbm_fit=True,
                                show_annotations=True,
                                energy_mode=energy_mode_key,
                                base_font_size=16, title_font_size=26,
                                figure_height=max(700, int(ups_eh * 0.62)),
                                line_width_scale=ups_line_ws,
                            )
                            img_bytes = export_ups_figure_bytes(
                                fig=export_fig,
                                export_format=ups_export_fmt.lower(),
                                width_px=ups_ew, height_px=ups_eh,
                                dpi=int(ups_export_dpi),
                            )
                            export_name = f"{ups_filename.rsplit('.', 1)[0]}_ups_figure.{ups_export_fmt.lower()}"
                            st.session_state.ups_export_bytes = img_bytes
                            st.session_state.ups_export_name = export_name
                            st.success("Export figure prepared.")
                        except Exception as exc:
                            st.error(f"Export failed: {exc}")

                    if st.session_state.get('ups_export_bytes'):
                        mime_map = {"png": "image/png", "svg": "image/svg+xml", "pdf": "application/pdf"}
                        ext = st.session_state.ups_export_name.split(".")[-1].lower()
                        st.download_button(
                            label=f"⬇️ Download {st.session_state.ups_export_name}",
                            data=st.session_state.ups_export_bytes,
                            file_name=st.session_state.ups_export_name,
                            mime=mime_map.get(ext, "application/octet-stream"),
                            use_container_width=True,
                            key="ups_download_figure",
                        )

# ============================================================================

with tab8:
    st.markdown("""
    ## About This Tool

    This interactive visualization tool explores the physics of two-dimensional electron gases (2DEG)
    at oxide surfaces, connecting surface band bending to observable quantities in photoemission
    spectroscopy (XPS/UPS).

    ### Three Physical Models

    - **M1 (Triangular)**: Constant electric field approximation
    - **M2 (Fang-Howard)**: Self-consistent variational approach with adaptive width
    - **M3 (Parabolic)**: Linearly decaying field (harmonic potential)

    ### Key Relationships

    1. **Gauss's Law**: nₛ = (ε·Eₛ)/q
    2. **Band Bending**: Φₛ relates to surface field Eₛ
    3. **Work Function**: ΔWF = -Φₛ + ΔΦ_dip
    4. **XPS Shift**: ΔE_CL = weighted average of band bending potential

    ### Model Comparison

    The slope of ΔWF vs nₛ differs between models:
    - **M1**: slope = -q·W/ε
    - **M3**: slope = -q·W/(2ε) (exactly half of M1)
    - **M2**: intermediate, field-dependent

    ### New Features (v2.0)

    #### 🔬 Experiment Comparison
    - Import experimental CSV data (UPS/XPS measurements)
    - Automatic parameter fitting (W, η optimization)
    - Residual analysis and goodness-of-fit diagnostics
    - Theory-experiment overlay plots

    #### 📤 Publication Export
    - High-quality figure export (SVG/PNG/PDF)
    - Journal-specific styles (Nature, Science, ACS, Grayscale)
    - Customizable sizes and fonts
    - Ready for journal submission

    #### 🧪 Beta Features
    - Upcoming features in development
    - Self-consistent diagnostics
    - Experimental guidance mode
    - Uncertainty analysis

    ### Bug Fixes (v2.1)

    #### 🐛 Fixed Critical Issues
    - **Unit Conversion Bug**: Corrected density conversion in `ns_to_display()` function
    - **ns Range Fix**: Updated Figure 2 density range to realistic 0.2-2.0 × 10¹³ cm⁻²
    - **Adsorbate Validation**: Added input validation and warnings for unusual dipole shifts
    - **Plot Sanity Checks**: Added guardrail checks to catch unit conversion errors before plotting

    ### Bug Fixes (v2.1.1)

    #### 🐛 Fixed UI and Visualization Issues
    - **Comparison Plot**: Now displays theoretical curve alongside experimental data with η comparison
    - **m* Uncertainty Band**: Fixed visualization - now properly shows shaded region for m* = 0.30-0.35 m₀
    - **M2 Model UI**: Cleaned up Additional Plots tab - no longer shows empty "Quantum Subband Energy Levels" title for Fang-Howard model

    ### Bug Fixes (v2.1.2)

    #### 🐛 Fixed Experiment Comparison Critical Issues
    - **Sample Data Physics**: Fixed sample data generation to produce physically reasonable In₂O₃ vacuum annealing data
      - Both ΔWF and ΔE_CL now correctly show negative values (desorption scenario)
      - Positive correlation with η ≈ 0.85 (previously showed unphysical negative correlation)
      - Includes realistic desorption kinetics with temperature
    - **Data Validation**: Added comprehensive physics validation for experimental data
      - Detects negative correlations (sign convention errors)
      - Warns about unusual η values (<0.5 or >0.95)
      - Checks for insufficient data range
      - Validates sign consistency between ΔWF and ΔE_CL
    - **Publication Export**: Fixed Figure 3 (comparison plot) to include theory curve
      - Now properly displays blue theory line alongside experimental data
      - Shows red experimental fit line with η and R² values
      - Includes annotation box with fit statistics and theory comparison
      - Added zero-crossing reference lines

    ### Bug Fixes (v2.1.4)

    #### 🐛 Fixed Parameter Fitting Critical Issue
    - **Parameter Fitting Curve Direction**: Fixed inverted fitted curves in "Experiment Comparison" tab
      - Fitted curves now correctly match the trend of experimental data
      - Corrected sign convention: `Delta_CL = eta * Delta_WF` (direct proportional, not inverse)
      - Test verification: η = 0.860 (expected ~0.85), R² = 0.9909, trends match ✓

    ### Usage Tips

    1. Use the sidebar to adjust parameters
    2. Compare different models using "Add to Compare"
    3. **NEW**: Import experimental data for fitting with automatic validation
    4. **NEW**: Export publication-quality figures with theory curves
    5. Adsorbates shift ΔWF without changing slope
    6. **NEW**: Physics validation warns about data quality issues

    ### Use Case: In₂O₃ Annealing Experiment

    This tool is optimized for analyzing in-situ annealing experiments:
    - **Material**: In₂O₃ thin films (100 nm typical)
    - **Experiment**: UHV annealing (25-400°C)
    - **Measurements**: UPS (work function) + XPS (core level shifts)
    - **Goal**: Extract η factor and depletion width W
    - **Expected**: ΔWF < 0, ΔE_CL < 0, η ≈ 0.7-0.9

    ### References

    - Fang & Howard, Phys. Rev. B **13**, 1546 (1966)
    - Copie et al., Adv. Mater. **29**, 1604112 (2017)
    - Salvinelli et al., ACS Appl. Mater. Interfaces **10**, 25941 (2018)

    ---
    ### New Features (v2.2)

    #### 🔮 UPS Spectrum Analyzer
    - Upload UPS spectrum CSV files for work function analysis
    - Automatic secondary electron cutoff (SEC) detection via linear extrapolation
    - Valence band maximum (VBM) / Fermi edge detection
    - Work function calculation: Φ = hν − spectrum width
    - Support for He I, He II, Ne I, Ne II, and custom photon sources
    - Both binding energy and kinetic energy scale support
    - Interactive visualization with edge fit overlays
    - Export results (CSV) and publication-quality figures (PNG/SVG/PDF)

    **Version**: 2.2.0
    **Updated**: March 2026
    **Framework**: Python + Streamlit + Plotly + Matplotlib
    **GitHub**: [2DEG-S-P-toy](https://github.com/aaronderek/2DEG-S-P-toy)
    """)

# ============================================================================
# COMPARISON FUNCTIONALITY
# ============================================================================

if add_comparison:
    # Generate comparison curve data
    colors = ['red', 'green', 'orange', 'purple']
    color_idx = len(st.session_state.comparison_curves) % len(colors)

    comparison_data = {
        'name': f"{model_type} (W={W_nm:.1f}nm)",
        'model_type': model_type,
        'color': colors[color_idx],
        'Phi_s_range': Phi_s_range,
        'ns_array': ns_array,
        'ns_range': ns_range,
        'Delta_WF_array': Delta_WF_array,
        'with_adsorbates': show_adsorbates,
        'Delta_WF_with_ads': Delta_WF_array + Delta_Phi_dip if show_adsorbates else None,
        'params': {
            'm_star': m_star_ratio,
            'epsilon_r': epsilon_r,
            'W_nm': W_nm
        }
    }

    st.session_state.comparison_curves.append(comparison_data)
    st.sidebar.success(f"Added: {comparison_data['name']}")

if clear_comparison:
    st.session_state.comparison_curves = []
    st.sidebar.info("Comparison curves cleared")

# Display comparison list
if st.session_state.comparison_curves:
    st.sidebar.markdown("### Comparison Curves")
    for i, curve in enumerate(st.session_state.comparison_curves):
        st.sidebar.markdown(f"- {curve['name']}")

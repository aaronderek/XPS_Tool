# Changelog

## Version 2.1.5 - 2025-11-12

### Critical Bug Fixes 🐛

#### Potential Energy Profile Direction Corrected (CRITICAL)
- **Issue**: Theoretical η calculation gave incorrect values (0.29-0.57 instead of expected 0.35-0.55)
- **Root Cause**:
  - M1 (Triangular) and M2 (Fang-Howard) models used wrong potential direction
  - Code: `V(z) = Q × Es × z` (increases from 0 to Φs)
  - Physical: `V(z) = Q × Φs × (1 - z/W)` (decreases from Φs to 0)
  - M3 (Parabolic) was correct: `V(z) = Q × Φs × (1 - z/W)²`
- **Physical Context**:
  - Potential energy V(z) should be maximum at surface (z=0) and zero at bulk (z=W)
  - Electron experiences a potential barrier that decreases with depth
  - XPS samples this potential with exponential weighting w(z) = exp(-z/λ)
- **Fix**:
  - M1: Changed to `V(z) = Q × Φs × (1 - z/W)` (triangular profile, correct direction)
  - M2: Changed to `V(z) = Q × Φs × (1 - z/W_eff)` (self-consistent triangular)
  - Added detailed comments explaining electric potential vs. potential energy
- **Verification**:
  - All models now give V(z=0) = Φs ✓
  - All models now give V(z=W) ≈ 0 ✓
  - M1: η ≈ 0.50 (triangular averages to 50% of uniform potential) ✓
  - M2: η ≈ 0.43 (similar to triangular) ✓
  - M3: η ≈ 0.37 (parabolic averages to 33% of uniform potential) ✓
- **Impact**:
  - Theoretical η now physically reasonable and model-dependent
  - Users can distinguish between potential models based on experimental data
  - η values match theoretical expectations for each potential shape
- **Files Modified**: `models/triangular.py`, `models/fang_howard.py`

### User Interface Improvements 🎨

#### Publication Export Figure Enhancements
- **Improvements**:
  - **Annotation Box**: Moved from right to left side to avoid overlapping theory curve
  - **Transparency**: Increased annotation box transparency (alpha=0.98) for better curve visibility
  - **Legend Simplification**: Removed redundant experimental fit line (red dashed)
    - Before: 4 items (Theory, Experiment, Exp. fit, Best fit)
    - After: 3 items (Theory (initial params), Experimental data, Best fit)
  - **Label Clarity**:
    - "Theory" → "Theory (initial params)" (clarifies this uses input W, λ parameters)
    - "Experiment" → "Experimental data"
    - "Fitted" → "Best fit" (clearer terminology)
  - **X-axis Range**: Dynamic adjustment to minimize empty space
    - Auto-calculates range as [min(ΔWF) - 0.02, max(ΔWF) + 0.02]
  - **Visual Enhancements**:
    - Larger scatter markers: 100 → 120 (better visibility)
    - Thicker lines: 2.5 → 2.8 (theory), 3.0 (best fit)
    - Best fit now uses solid line (more important than initial theory)
    - Larger axis labels: 14 → 15 pt
    - Enhanced legend with shadow and rounded corners
- **Result**: Cleaner, more professional figures suitable for publication
- **Files Modified**: `utils/publication_export.py`

### Documentation & Testing 📝

#### Comprehensive Test Suite Added
- Added `test_eta_theory_debug.py`: Basic η calculation verification
- Added `test_theory_curve_calculation.py`: Full curve generation simulation
- Added `test_all_models_eta.py`: Comparison of all three models
- Added `test_m2_potential_bug.py`: Detailed M2 potential profile analysis
- Added `test_fix_validation.py`: Comprehensive validation of all fixes
- All tests confirm physical correctness of η calculations

## Version 2.1.4 - 2025-11-12

### Critical Bug Fix 🐛

#### Parameter Fitting Curve Direction Corrected (CRITICAL)
- **Issue**: Parameter fitting output curve showed opposite trend compared to sample data and theory
- **Root Cause**:
  - Incorrect sign in `calculate_theory_Delta_CL()`: used `-eta * Delta_WF` instead of `eta * Delta_WF`
  - Also incorrect sign extraction in `linear_fit_eta()`: used `eta_exp = -slope` instead of `eta_exp = slope`
  - This caused fitted curves to be inverted (positive when should be negative, and vice versa)
- **Physical Context**:
  - Correct relationship: `Delta_CL = eta * Delta_WF` (direct proportional)
  - Both quantities follow same sign convention (both negative in vacuum annealing, both positive in oxidation)
  - Previous code incorrectly assumed inverse relationship
- **Fix**:
  - Changed `Delta_CL_theory = -eta * Delta_WF_exp` → `Delta_CL_theory = eta * Delta_WF_exp`
  - Changed `eta_exp = -slope` → `eta_exp = slope`
  - Updated function docstrings to clarify correct relationship
- **Verification**:
  - Test with sample data: η = 0.860 (expected ~0.85) ✓
  - Correlation with experiment: R² = 0.9909 ✓
  - Curve trends now match correctly ✓
- **Impact**: Parameter fitting now produces physically correct curves that match experimental data
- **Files Modified**: `utils/fitting.py`

## Version 2.1.2 - 2025-11-12

### Critical Bug Fixes 🐛

#### Sample Data Physics Corrected (CRITICAL)
- **Issue**: Sample data showed unphysical negative correlation between ΔWF and ΔE_CL
- **Root Cause**:
  - Generated positive ΔWF values (work function increase)
  - Applied `In3d = In3d_base - Delta_CL`, causing negative ΔE_CL
  - Result: ΔWF positive, ΔE_CL negative → negative correlation, η = -0.850
- **Physical Context**:
  - In₂O₃ vacuum annealing: water desorption → surface potential increases
  - Should give: ΔWF < 0 (work function decreases), ΔE_CL < 0 (binding energy decreases)
  - Expected: Positive correlation with η ≈ 0.7-0.9
- **Fix**:
  - Rewrote `create_sample_data()` with correct physics
  - Simulates realistic desorption kinetics (sigmoid activation at 100-200°C)
  - Includes surface potential evolution (0.2 → 0.5 eV) and dipole layer effects
  - Both ΔWF and ΔE_CL now correctly negative with η = 0.85
  - Added noise simulation (±0.01 eV) for realism
- **Impact**: Sample data now demonstrates correct physics, suitable for teaching and testing
- **Files Modified**: `utils/experiment_data.py`

#### Data Validation System Added (NEW FEATURE)
- **Purpose**: Detect common experimental data issues and warn users
- **Implementation**: New function `validate_experimental_data_physics()`
- **Checks Performed**:
  1. **Correlation Check**: Detects negative correlation (sign convention errors)
  2. **η Range Check**: Warns if η < 0.5 or η > 0.95 (outside typical TCO range)
  3. **Data Range Check**: Warns if ΔWF range < 0.1 eV (insufficient for fitting)
  4. **Sign Consistency**: Validates ΔWF and ΔE_CL have same sign trend
- **User Experience**:
  - Warnings displayed in "Data Quality Checks" section after data import
  - Clear explanations of physical meaning and possible causes
  - Actionable suggestions for improvement
- **Files Modified**: `utils/experiment_data.py`, `app.py`

#### Publication Export Missing Theory Curve (CRITICAL)
- **Issue**: Exported Figure 3 (comparison plot) only showed experimental scatter points
- **Root Cause**: `create_publication_comparison_figure()` didn't accept model parameters
- **Impact**: Exported figures incomplete, not suitable for publication
- **Fix**:
  - Added `model`, `lambda_nm`, `theta_deg` parameters to function
  - Function now calculates theory curve from model (same as interactive plot)
  - Generates theoretical ΔE_CL vs ΔWF curve via XPS integration
  - Added experimental linear fit line (red dashed)
  - Added annotation box with η_exp, η_theory, R², and comparison
  - Added zero-crossing reference lines (x=0, y=0)
- **Result**: Exported figures now match interactive display quality
- **Files Modified**: `utils/publication_export.py`, `app.py`

### Verification
- ✅ Sample data shows η = 0.85 ± 0.01 with R² > 0.99
- ✅ Both ΔWF and ΔE_CL negative (vacuum annealing scenario)
- ✅ Data validation correctly identifies issues in problematic datasets
- ✅ Publication export includes blue theory line, red scatter points, and red fit line
- ✅ Annotation box shows η comparison and statistical metrics

### Files Modified
- `utils/experiment_data.py`: Rewrote sample data generation, added validation
- `utils/publication_export.py`: Added theory curve calculation to export function
- `app.py`: Integrated validation display, updated export function call, updated version to 2.1.2
- `README.md`: Added v2.1.2 to version history

### Documentation Updates
- Updated in-app "About" section with v2.1.2 fixes
- Updated `README.md` version history
- Updated `CHANGELOG.md` (this file)

---

## Version 2.1.1 - 2025-11-12

### Bug Fixes 🐛

#### Comparison Plot Missing Theory Curve (CRITICAL)
- **Issue**: Theory-experiment comparison plot only showed experimental data points
- **Root Cause**: `create_comparison_CL_vs_WF_plot()` was called with `theory_data=None`, and the function didn't calculate theoretical curve from model
- **Impact**: Users couldn't visually compare theory predictions with experimental data
- **Fix**:
  - Modified `create_comparison_CL_vs_WF_plot()` to accept model parameters
  - Function now calculates theoretical ΔE_CL vs ΔWF curve from model
  - Added experimental linear fit line (red dashed)
  - Added annotation showing η_exp vs η_theory comparison
- **Files Modified**: `ui/plots.py`, `app.py`

#### m* Uncertainty Band Not Displayed (HIGH)
- **Issue**: "Show m* uncertainty band" checkbox had no effect
- **Root Cause**: Checkbox value `show_uncertainty` was not passed to plotting functions
- **Impact**: Users couldn't visualize the impact of effective mass uncertainty (m* = 0.30-0.35 m₀)
- **Fix**:
  - Added `show_uncertainty` and `model_obj` parameters to `create_ns_vs_Phi_s_plot()`
  - Added `show_uncertainty` and `model_obj` parameters to `create_Delta_WF_vs_ns_plot()`
  - Functions now calculate upper/lower bounds with different m* values
  - Display shaded region using `fill='toself'` in Plotly
- **Files Modified**: `ui/plots.py`, `app.py`

#### M2 Model Shows Empty Subband Levels Title (MEDIUM)
- **Issue**: M2 (Fang-Howard) model displayed "Quantum Subband Energy Levels" title with no content below
- **Root Cause**: UI logic showed title before checking if model has subband methods
- **Background**: M2 is a variational approximation that doesn't solve for discrete subband levels (E₀, E₁, E₂...)
- **Fix**:
  - Moved title inside conditional check for `hasattr(model, 'get_subband_energies')`
  - Added informative message for M2 explaining why subbands aren't shown
  - Enhanced M2's n(z) display with key parameters (b, W_eff, ⟨z⟩)
- **Files Modified**: `app.py`

### Verification
- ✅ Comparison plot now shows blue theory curve + red experimental points + red dashed fit
- ✅ η_theory and η_exp values displayed in annotation box
- ✅ m* uncertainty band appears as shaded region when checkbox enabled
- ✅ M2 model no longer shows empty subband title

### Files Modified
- `ui/plots.py`: Updated plot functions with uncertainty bands and theory curve calculation
- `app.py`: Updated plot function calls and M2 subband display logic

## Version 2.1 - 2025-11-12

### Critical Bug Fixes 🐛

#### Unit Conversion Bug Fixed
- **Issue**: `ns_to_display()` function in `physics/units.py` had incorrect density conversion
- **Root Cause**: Used `CM2_TO_M2` (area conversion) instead of inverse for density conversion
- **Impact**: Sheet density displayed values were off by factor of 10⁸
- **Fix**: Corrected to use proper density conversion (1 m⁻² = 10⁻⁴ cm⁻², inverse of area)
- **Files Modified**: `physics/units.py`

#### ns Range Correction in Figure 2
- **Issue**: Figure 2 (ΔWF vs nₛ) used incorrect range `np.linspace(1e12, 2e13, 100)` m⁻²
- **Problem**: This corresponds to 0.00001-0.0002 × 10¹³ cm⁻² (far too small!)
- **Expected**: Typical 2DEG density range should be 0.2-2.0 × 10¹³ cm⁻²
- **Fix**: Updated to `np.linspace(2e16, 2e17, 100)` m⁻² (= 0.2-2.0 × 10¹³ cm⁻²)
- **Impact**: Curves now display proper slopes instead of appearing nearly horizontal
- **Files Modified**: `app.py` line 267-269

#### Adsorbate Input Validation
- **Added**: Input validation in `calculate_dipole_from_coverage()`
  - Clamps coverage θ to [0, 1] range
  - Warns if N_site outside typical range (10¹³-10¹⁶ cm⁻²)
  - Warns if μ⊥ outside typical range (0-10 Debye)
  - Warns if |ΔΦ_dip| > 2 eV (unusually large)
- **Purpose**: Catch unit conversion errors before they corrupt plots
- **Files Modified**: `physics/xps.py`

#### Plot Sanity Checks
- **Added**: Guardrail checks in `create_Delta_WF_vs_ns_plot()`
  - Warns if |ΔWF| > 5 eV (likely unit error)
  - Warns if nₛ display values outside 0.001-1000 × 10¹³ cm⁻² range
- **Purpose**: Prevent nonsensical plots from reaching user
- **Files Modified**: `ui/plots.py`

### Verification
- **Test Suite**: `test_adsorbate_bug.py` now passes all checks
  - ✅ nₛ range: 0.2 to 2.0 × 10¹³ cm⁻²
  - ✅ ΔΦ_dip: -0.283 eV (expected value)
  - ✅ Slope: -0.603 eV/(10¹³ cm⁻²) for M1 with εᵣ=9, W=3nm
  - ✅ Parallel curves: "with ads" and "no ads" have identical slopes
  - ✅ Vertical shift: Exactly equals ΔΦ_dip

### Files Modified
- `physics/units.py`: Fixed density conversion logic
- `physics/xps.py`: Added input validation and warnings
- `ui/plots.py`: Added sanity checks for plot data
- `app.py`: Updated ns_range and version to 2.1

### Documentation Updates
- Updated in-app "About" section with v2.1 bug fix summary
- Updated `CHANGELOG.md` (this file)
- Updated `README.md` with v2.1 information

## Version 2.0 - 2025-11-11

### Major New Features

#### 🔬 Experiment Comparison (High Priority)
- **CSV Data Import**: Support for two formats:
  - Format A: Raw measurements (T_degC, WF_eV, CL_eV)
  - Format B: Processed data (Delta_WF_eV, Delta_CL_eV)
- **Data Validation**: Automatic format detection and data quality checks
- **Sample Data**: Built-in sample In₂O₃ annealing data for testing

#### 🎯 Theory-Experiment Fitting (High Priority)
- **Automatic Parameter Fitting**: Optimize W (depletion width) and η (XPS sampling factor)
- **Scipy-based Optimization**: L-BFGS-B algorithm for robust convergence
- **Goodness-of-Fit Metrics**: R², RMSE, residual analysis
- **Fitting Diagnostics**:
  - Residuals vs ΔWF plot
  - Residuals histogram
  - Predicted vs actual scatter
  - Statistical summary table

#### 📊 New Visualization Plots
- **ΔE_CL vs ΔWF Comparison Plot**: Overlay experimental and fitted data
- **Annealing Trajectory Plot**: Show evolution with temperature (if T data available)
- **Comprehensive Residual Analysis**: 4-panel diagnostic plot

#### 📤 Publication-Quality Export (High Priority)
- **High-Quality Figure Export**:
  - SVG (vector graphics, infinite resolution)
  - PNG (raster, 300/600/1200 DPI)
  - PDF (document format)
- **Journal-Specific Styles**:
  - Nature (warm colors, grid)
  - Science (cool colors, no grid)
  - ACS (standard colors)
  - Grayscale (B&W with different markers)
- **Customizable Settings**:
  - Figure size presets (single/double column, custom)
  - Font size control
  - Line width adjustment
- **Matplotlib Backend**: High-quality rendering for publication

#### 🧪 Beta Features Tab
- Placeholder for upcoming features:
  - Self-consistent S-P diagnostic
  - Experimental guidance mode
  - Uncertainty propagation analysis
  - Material database
  - Annealing trajectory animation

### Technical Improvements

#### New Modules
- `utils/experiment_data.py`: Data import and validation
- `utils/fitting.py`: Parameter optimization and linear regression
- `utils/publication_export.py`: High-quality figure generation with journal styles

#### Enhanced UI
- 6 tabs (was 4): Core Figures, Additional Plots, Experiment Comparison, Publication Export, Beta Features, About
- Better organization and workflow
- Improved user feedback and error messages

#### Dependencies
- Added `matplotlib>=3.7.0` for publication-quality exports
- Existing `scipy>=1.11.0` now used for fitting

### Files Added
- `utils/experiment_data.py`
- `utils/fitting.py`
- `utils/publication_export.py`
- `sample_data_In2O3.csv`
- `CHANGELOG.md`

### Files Modified
- `app.py`: Major refactor with new tabs and features
- `ui/plots.py`: Added new plotting functions
- `requirements.txt`: Added matplotlib dependency

## Version 1.0 - 2025-11-XX

Initial release with:
- Three physical models (M1-Triangular, M2-Fang-Howard, M3-Parabolic)
- Core figures (ns vs Φs, ΔWF vs ns)
- Additional plots (V(z), n(z), w(z))
- Parameter controls
- Model comparison
- Basic data export (CSV, JSON)

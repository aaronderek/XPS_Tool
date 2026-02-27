# 2DEG Surface Electron Gas Visualization Tool

Interactive tool for exploring two-dimensional electron gas (2DEG) physics at oxide surfaces, connecting surface band bending to observable quantities in X-ray photoelectron spectroscopy (XPS/UPS).

![2DEG Visualization](https://img.shields.io/badge/physics-2DEG-blue) ![Python](https://img.shields.io/badge/python-3.8+-green) ![Streamlit](https://img.shields.io/badge/streamlit-1.28+-red)

**Version**: v2.1.4 (2025-11-14)

## Features

### 🎯 Core Functionality
- **Two Essential Figures**:
  - **Figure 1**: Sheet Density (nₛ) vs Surface Potential (Φₛ)
  - **Figure 2**: Work Function Change (ΔWF) vs Sheet Density (nₛ)

- **Three Physical Models**:
  - **M1 (Triangular)**: Constant electric field approximation
  - **M2 (Fang-Howard)**: Self-consistent variational approach
  - **M3 (Parabolic)**: Linearly decaying field (harmonic potential)

- **XPS Modeling**:
  - Core level shift calculations
  - Escape depth effects (λ, θ dependence)
  - Adsorbate dipole layer effects

### 🔬 Experiment Comparison (NEW in v2.0)
- **CSV Data Import**: Load experimental UPS/XPS measurements
  - Format A: Raw data (T, WF, core level)
  - Format B: Processed data (ΔWF, ΔE_CL)
- **Automatic Parameter Fitting**: Optimize W and η to match experimental data
- **Residual Analysis**: Comprehensive goodness-of-fit diagnostics (R², RMSE, residual plots)
- **Theory-Experiment Overlay**: Direct visual comparison
- **Annealing Trajectory**: Visualize temperature-dependent evolution

### 📤 Publication Export (NEW in v2.0)
- **High-Quality Figures**:
  - SVG (vector graphics for editing)
  - PNG (300/600/1200 DPI)
  - PDF (document-ready)
- **Journal Styles**:
  - Nature (warm colors, grid)
  - Science (cool colors, clean)
  - ACS (standard colors)
  - Grayscale (B&W with markers)
- **Customizable**: Size presets, fonts, line widths

### 🧲 XPS Spectrum Analyzer (NEW)
- **Direct `.spe` Upload**: Auto-convert PHI `.spe` traces to analysis-ready spectra using the same conversion workflow as `convert_spe_to_csv.py`
- **Peak Detection & Fitting**: Automatic peak finding + multi-peak fitting (Gaussian/Lorentzian)
- **Element Annotation**: Built-in reference lines (C 1s, O 1s, In 3d, etc.) and fitted-peak matching
- **Batch Mode for SPE**: One-click fitting for all traces in a single `.spe` file with CSV summary export

### 🔧 Interactive Controls
- Real-time parameter adjustment with sliders
- Model comparison (overlay up to 4 curves)
- Adsorbate effects (coverage, dipole moment)
- Additional visualizations: V(z), n(z), w(z) profiles

### 💾 Export Options
- **Figures**: PNG, SVG, PDF (publication-quality)
- **Data**: CSV with all parameters
- **Configuration**: JSON for reproducibility

## Installation

### Requirements
- Python 3.8 or higher (3.11 recommended)
- pip package manager
- git (for cloning)
- (Recommended) virtual environment via `venv`

### Setup (Windows, macOS, Linux)
1. **Clone the repository**
   ```bash
   git clone https://github.com/yourusername/2DEG-S-P-toy.git
   cd 2DEG-S-P-toy
   ```
2. **Create and activate a virtual environment (recommended)**
   - Windows:
     ```bash
     py -3.11 -m venv .venv
     .venv\Scripts\activate
     ```
   - macOS/Linux:
     ```bash
     python3 -m venv .venv
     source .venv/bin/activate
     ```
3. **Install dependencies**
   ```bash
   python -m pip install --upgrade pip
   python -m pip install -r requirements.txt
   ```
4. **Run the app**
   - Windows: `run_app.bat` or `python -m streamlit run app.py`
   - macOS/Linux: `bash run.sh` or `streamlit run app.py`
   - The UI opens at `http://localhost:8501`

### Dependency summary
```
streamlit>=1.28.0
plotly>=5.17.0
matplotlib>=3.7.0
numpy>=1.24.0
scipy>=1.11.0
pandas>=2.0.0
kaleido>=0.2.1
```

## Usage

### Running the Application

Start the Streamlit app:
```bash
streamlit run app.py
```

The application will open in your default web browser at `http://localhost:8501`.

### Basic Workflow

1. **Select Model**: Choose M1, M2, or M3 from the sidebar
2. **Adjust Parameters**:
   - Material: m*/m₀, εᵣ, T
   - Surface: Φₛ, W
   - XPS: λ, θ
3. **View Results**: Observe real-time updates in both figures
4. **Compare Models**: Click "Add to Compare" to overlay multiple curves
5. **Export Data**: Use the Export tab to save figures and data

### Example Use Cases

#### 1. Analyze In₂O₃ Annealing Experiment (NEW)
```
1. Go to "Experiment Comparison" tab
2. Click "Load Sample Data" (or upload your CSV)
3. Review data preview and estimated η
4. Enable "Fit η" and click "Run Fitting"
5. View R² and residual analysis
6. Export publication figure in "Publication Export" tab
```

#### 2. Compare Models for Different Depletion Widths
```
1. Set W = 2 nm, select M1, click "Add to Compare"
2. Set W = 3 nm, select M1, click "Add to Compare"
3. Set W = 5 nm, select M1, click "Add to Compare"
4. Observe linear scaling of nₛ with 1/W
```

#### 3. Explore Adsorbate Effects
```
1. Set Φₛ = 0.4 eV
2. Check "Show with adsorbates"
3. Adjust coverage θ and dipole μ⊥
4. Observe vertical shift in Figure 2 (slope unchanged)
```

#### 4. Generate Publication Figures (NEW)
```
1. Configure your desired plot in "Core Figures"
2. Go to "Publication Export" tab
3. Select journal style (e.g., "Nature")
4. Choose format (SVG for vector graphics)
5. Select size preset (e.g., "Double column")
6. Click "Generate Publication Figures"
7. Download and import into Illustrator/Inkscape
```

#### 5. Analyze Raw PHI XPS `.spe` Files (NEW)
```
1. Go to "XPS Spectrum" tab
2. Drag in a .spe file (or .csv)
3. For .spe files, pick a converted trace (e.g., C1s/O1s/In3d)
4. Tune detection parameters and click "Run Peak Fitting"
5. Inspect fitted peaks + element matches, then export peak-fit CSV
6. (Optional) Run "Batch Fit For All Traces" for one-click summary
```

## Project Structure

```
2DEG-S-P-toy/
├── app.py                      # Main Streamlit application
├── models/
│   ├── __init__.py
│   ├── triangular.py          # M1: Triangular model
│   ├── fang_howard.py         # M2: Fang-Howard model
│   └── parabolic.py           # M3: Parabolic model
├── physics/
│   ├── __init__.py
│   ├── constants.py           # Physical constants
│   ├── units.py               # Unit conversions
│   └── xps.py                 # XPS modeling
├── ui/
│   ├── __init__.py
│   └── plots.py               # Plotly plotting functions
├── utils/
│   ├── __init__.py
│   ├── export.py              # Export utilities
│   ├── experiment_data.py     # NEW: Data import & validation
│   ├── fitting.py             # NEW: Parameter optimization
│   ├── xps_spectrum.py        # NEW: SPE/CSV conversion, peak fitting, XPS annotation
│   └── publication_export.py  # NEW: High-quality figure export
├── tests/
│   └── test_models.py         # Unit tests
├── sample_data_In2O3.csv      # NEW: Example experimental data
├── CHANGELOG.md               # NEW: Version history
├── README.md
└── requirements.txt
```

## Physics Background

### Key Relationships

1. **Gauss's Law** (universal):
   ```
   nₛ = (ε·Eₛ) / q
   ```

2. **Model-Specific**:
   - **M1**: Φₛ = Eₛ·W  →  nₛ = (ε/q)·(Φₛ/W)
   - **M2**: Self-consistent iteration with b = (12m*qEₛ/ℏ²)^(1/3)
   - **M3**: Eₛ = 2Φₛ/W  →  nₛ = (2ε/q)·(Φₛ/W)

3. **Work Function**:
   ```
   ΔWF = -Φₛ + ΔΦ_dip
   ```

4. **XPS Core Level Shift**:
   ```
   ΔE_CL = -∫ w(z)·V(z) dz
   w(z) = exp(-z/(λ·cosθ)) / (λ·cosθ)
   ```

### Model Comparison

| Property | M1 | M2 | M3 |
|----------|----|----|----|
| Field profile | Constant | Exponential decay | Linear decay |
| nₛ at Φₛ=0.3eV, W=3nm | 5.3×10¹³ cm⁻² | ~7×10¹³ cm⁻² | 10.6×10¹³ cm⁻² |
| ΔWF slope | -q·W/ε | Intermediate | -q·W/(2ε) |
| Physical regime | Low accumulation | General | High accumulation |

**Key Insight**: M3 gives exactly **2× the nₛ** and **½ the slope** compared to M1!

## Testing

Run validation tests:
```bash
python tests/test_models.py
```

### Test Cases

The test suite validates:
1. ✅ M1 model numerical accuracy
2. ✅ M3 model numerical accuracy
3. ✅ M3/M1 ratio exactly equals 2.0
4. ✅ Adsorbate shift (constant offset)
5. ✅ M2 self-consistent convergence
6. ✅ Zero potential boundary condition

Expected output:
```
=== Test 1: M1 Model Basic ===
Result: ns = 5.285 × 10¹³ cm⁻²
Expected: ns ≈ 5.3 × 10¹³ cm⁻²
✓ Test 1 passed!
...
Results: 6 passed, 0 failed
```

## Default Parameters

```python
m*/m₀ = 0.32          # Effective mass (typical for SrTiO₃)
εᵣ = 9                # Relative permittivity
T = 300 K             # Temperature
W = 3.0 nm            # Depletion width
Φₛ = 0.4 eV           # Surface potential
λ = 1.8 nm            # XPS mean free path
θ = 0°                # XPS detection angle (normal emission)
```

## Scientific References

1. **Fang-Howard Model**:
   - Fang & Howard, *Phys. Rev. B* **13**, 1546 (1966)

2. **2DEG at Oxide Surfaces**:
   - Ohtomo & Hwang, *Nature* **427**, 423 (2004)
   - Copie et al., *Adv. Mater.* **29**, 1604112 (2017)

3. **XPS and Work Function**:
   - Salvinelli et al., *ACS Appl. Mater. Interfaces* **10**, 25941 (2018)
   - Dudy et al., *Adv. Mater.* **28**, 7443 (2016)

## Troubleshooting

### Common Issues

**Issue**: `ModuleNotFoundError: No module named 'streamlit'`
- **Solution**: Run `pip install -r requirements.txt`

**Issue**: Figures not exporting to PNG
- **Solution**: Install kaleido: `pip install kaleido`

**Issue**: M2 model not converging
- **Solution**: Try reducing Φₛ or increasing W (avoid extreme field strengths)

**Issue**: Slow performance
- **Solution**: Reduce number of comparison curves (max 4 recommended)

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## License

MIT License - see LICENSE file for details.

## Authors

- **Aaron** - Initial development
- **Claude (Anthropic)** - Code generation assistance

## Acknowledgments

- Physics models based on established semiconductor theory
- UI framework powered by Streamlit
- Visualizations created with Plotly

## Version History

- **v2.1.4** (2025-11-14): Compatibility and documentation update
  - Fixed XPS integration for NumPy 2.4+ (np.trapezoid) to resolve M1/M3 crashes
  - Clarified installation, virtual environment, and run steps across Windows/macOS/Linux
  - See [CHANGELOG.md](CHANGELOG.md) for details

- **v2.1.2** (2025-11-12): Experiment comparison critical fixes
  - 🐛 **Fixed**: Sample data now shows physically reasonable In₂O₃ vacuum annealing behavior
    - Both ΔWF and ΔE_CL correctly negative (desorption scenario)
    - Positive correlation with η ≈ 0.85 (previously unphysical negative correlation)
  - ✨ **New**: Comprehensive physics validation for experimental data
    - Detects negative correlations (sign convention errors)
    - Warns about unusual η values or insufficient data range
  - 🐛 **Fixed**: Publication export Figure 3 now includes theory curve
    - Blue theory line with proper XPS calculation
    - Red experimental fit line with statistics
    - Annotation box with η comparison
  - See [CHANGELOG.md](CHANGELOG.md) for details

- **v2.1.1** (2025-11-12): UI and visualization bug fixes
  - 🐛 **Fixed**: Comparison plot now shows theoretical curve with η comparison
  - 🐛 **Fixed**: m* uncertainty band visualization now working correctly
  - 🐛 **Fixed**: M2 model no longer shows empty subband levels title
  - ✨ Enhanced M2 Additional Plots with b, W_eff, ⟨z⟩ parameters
  - See [CHANGELOG.md](CHANGELOG.md) for details

- **v2.1.0** (2025-11-12): Critical bug fix release
  - 🐛 Fixed unit conversion error in density display
  - 🐛 Corrected Figure 2 ns range (was 10⁸× too small!)
  - ✅ Added input validation for adsorbate parameters
  - ✅ Added plot sanity checks to catch unit errors
  - ✅ All validation tests now pass
  - See [CHANGELOG.md](CHANGELOG.md) for details

- **v2.0.0** (2025-11-11): Major feature update
  - ✨ Experiment comparison with CSV data import
  - 🎯 Automatic parameter fitting (W, η optimization)
  - 📊 Comprehensive residual analysis
  - 📤 Publication-quality figure export (SVG/PNG/PDF)
  - 🎨 Journal-specific styles (Nature, Science, ACS, Grayscale)
  - 🧪 Beta features tab for upcoming features
  - See [CHANGELOG.md](CHANGELOG.md) for details

- **v1.0.0** (2025-11-11): Initial release
  - Three physical models (M1, M2, M3)
  - Interactive Streamlit interface
  - XPS modeling and adsorbate effects
  - Export functionality
  - Comprehensive test suite

---

**Contact**: For questions or feedback, please open an issue on GitHub.

**Citation**: If you use this tool in your research, please cite:
```bibtex
@software{2deg_visualization_2025,
  title = {2DEG Surface Electron Gas Visualization Tool},
  author = {Aaron},
  year = {2025},
  url = {https://github.com/yourusername/2DEG-S-P-toy}
}
```

# Changelog

All notable changes to this repository are documented in this file.

## 2026-02-27 (Workspace / XPS Tab6 Update)

### Changed

- Refactored XPS Spectrum `tab6` in `app.py` from single-column into a two-column workflow layout.
  - Left: step-based controls (`Data Import`, `Region & Elements`, `Background & Fitting`, `Export & Downloads`)
  - Right: spectrum visualization and fit table
- Fixed manual upload parsing bug:
  - `filename = filename` -> `filename = uploaded_xps.name`
  - Resolved runtime error: `Failed to parse XPS file: name 'filename' is not defined`
- Replaced unstable sticky CSS/JS behavior with Streamlit official scroll containers:
  - `left_panel = col_left.container(height=860, ...)`
  - `right_panel = col_right.container(height=860, ...)`
- Removed failed global auto-height CSS override that caused severe layout compression.
- Updated default peak detection parameters for better XPS detail:
  - smoothing window `11 -> 5`
  - peak prominence `% 6 -> 4`

### Added

- Developer / AI Agent Auto-Loader in Step 1 for absolute-path local `.spe` loading.

### Known Tradeoff

- Independent pane scrolling currently uses fixed height `860`, which may need tuning for unusual viewport sizes.

## 2025-11-12 (Imported from Main Program History)

### v2.1.5

- Corrected potential profile direction for M1/M2 theoretical `η` calculation.
- Improved publication export visuals (annotation placement, legend cleanup, dynamic x-range).
- Added comprehensive fix-validation test scripts.

### v2.1.4

- Fixed fitting curve direction sign issue in `utils/fitting.py`.
- Corrected relationship to `ΔE_CL = η * ΔWF` for comparison/fitting workflow.

### v2.1.2

- Reworked sample data physics to align with expected annealing behavior.
- Added experimental data physics validation checks.
- Fixed publication comparison export to include theory curve.

### v2.1.1 / v2.1

- Added uncertainty/theory-comparison improvements in plotting.
- Fixed density unit conversion and range issues.
- Added XPS/adsorbate input sanity checks.

---

For full historical detail, also see:

- `XPS Tool/Main Program/CHANGELOG.md`

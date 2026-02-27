# XPS Tool Carbon Redesign Commit Notes

Date: 2026-02-27  
Branch: `tiger/updateddesign`

## Scope

This commit contains a visual refresh of the Streamlit app using Carbon-style design patterns, plus run/startup reliability improvements for Codex environments.

## Main Changes

1. Carbon theme foundation added:
- `.streamlit/config.toml`
- `carbon_streamlit.css`
- `carbon_theme.py`

2. App UI upgraded:
- Injected Carbon CSS on startup.
- Added sidebar brand block.
- Added hero summary section and KPI metrics row.
- Unified tab labels and comparison curve colors.
- Kept original physics/business logic intact.

3. Plot system restyled:
- Added shared Carbon palette and Plotly layout helpers.
- Unified chart colorway, grid, typography, spacing, and legend style.
- Applied to core physics plots, comparison plots, residual plots, and trajectory plots.

4. XPS spectrum visuals aligned:
- Introduced Carbon color tokens for spectrum rendering.
- Updated markers/lines/annotations/background/grid styling for consistency.

5. Run script hardened:
- `run.sh` now auto-detects Python (`python3`/`python`) and repairs `.venv` when needed.
- Installs dependencies via venv interpreter.
- Starts Streamlit with explicit host/port.
- Auto-opens browser by default (`AUTO_OPEN_BROWSER=0` to disable).

## Files Included

- `XPS Tool/Main Program/app.py`
- `XPS Tool/Main Program/run.sh`
- `XPS Tool/Main Program/ui/plots.py`
- `XPS Tool/Main Program/utils/xps_spectrum.py`
- `XPS Tool/Main Program/.streamlit/config.toml`
- `XPS Tool/Main Program/carbon_streamlit.css`
- `XPS Tool/Main Program/carbon_theme.py`
- `XPS Tool/Main Program/CARBON_REDESIGN_COMMIT_NOTES.md`

## Validation Done

- `python3 -m py_compile app.py ui/plots.py utils/xps_spectrum.py carbon_theme.py`
- Carbon audit script rerun: high-priority missing-theme issues resolved.

## Notes

- No merge to `main` is performed in this commit.
- Untracked local `.codex/` workspace files are intentionally excluded.

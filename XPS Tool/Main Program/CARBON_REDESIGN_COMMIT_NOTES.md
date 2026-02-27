# XPS Tool Carbon Redesign Commit Notes

Date: 2026-02-27  
Branch: `tiger/updateddesign`

## Scope

This branch applies a Carbon-style visual refresh to the Streamlit app and then performs a follow-up UI cleanup pass based on design feedback.

## Main Changes

1. Carbon theme foundation
- Added `.streamlit/config.toml` theme tokens.
- Added `carbon_streamlit.css` for component-level styling.
- Added `carbon_theme.py` and injects CSS from `app.py`.

2. App-level visual refactor
- Sidebar brand block introduced and simplified.
- Hero block retained with compact metadata chips.
- Removed redundant marketing copy from sidebar and hero.
- Removed the global KPI metric strip that appeared under hero across all tabs.

3. Tabs and navigation behavior
- Kept Streamlit native tab highlight as the single active indicator.
- Removed extra static tab underline styling to avoid dual-indicator visuals.

4. Expander and panel hierarchy cleanup
- XPS tab step sections switched from `st.container(border=True)` to plain `st.container()`.
- Expander wrappers flattened to reduce nested "box inside rounded box" effect.
- Border radius normalized to square corners for cards, expanders, and upload dropzones.

5. Color consistency tuning
- Unified slider accents and thumb-value color with Carbon primary blue.
- Unified checkbox/radio accent color with Carbon primary blue.
- Updated adsorbate curve color mapping to use the same app primary token.

6. Plot styling alignment
- Shared Carbon palette and Plotly layout helpers applied in plotting modules.
- Grid, typography, legend style, and chart colorway aligned across main plot surfaces.

## Files Included (Redesign + Cleanup)

- `XPS Tool/Main Program/app.py`
- `XPS Tool/Main Program/carbon_streamlit.css`
- `XPS Tool/Main Program/carbon_theme.py`
- `XPS Tool/Main Program/.streamlit/config.toml`
- `XPS Tool/Main Program/ui/plots.py`
- `XPS Tool/Main Program/utils/xps_spectrum.py`
- `XPS Tool/Main Program/run.sh`
- `XPS Tool/Main Program/CARBON_REDESIGN_COMMIT_NOTES.md`

## Validation Done

- `python3 -m py_compile app.py`
- Streamlit smoke run on localhost with clean startup and shutdown.
- Carbon audit rerun confirms no high-priority theme integration gaps.

## Notes

- Changes are committed on branch `tiger/updateddesign`.
- No merge to `main` is performed.

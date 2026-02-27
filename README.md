# XPS Tool

An integrated Streamlit toolkit for:

- 2DEG surface electron gas modeling (M1/M2/M3)
- XPS/UPS observable analysis (`ΔWF`, `ΔE_CL`, `η`)
- PHI `.spe` / `.csv` spectrum parsing, peak detection, fitting, and export

This repository focuses on the **software**. Large datasets and archive materials are kept locally and excluded from Git tracking.

## What It Can Do

- Interactive 2DEG physics visualization and model comparison
- Theory vs experiment fitting from imported CSV data
- Publication-ready figure export (PNG/SVG/PDF styles)
- XPS Spectrum Analyzer:
  - direct `.spe` upload and conversion fallback chain
  - element-region/manual-window fitting
  - per-scope export and batch fitting summaries

## Repository Layout

```text
XPS_Tool/
├── README.md
├── CHANGELOG.md
├── Launch_XPS_Tool.command
├── XPS Tool/
│   ├── Main Program/
│   │   ├── app.py
│   │   ├── models/
│   │   ├── physics/
│   │   ├── ui/
│   │   ├── utils/
│   │   ├── tests/
│   │   └── requirements.txt
│   └── XPS Plugin/
│       ├── convert_spe_to_csv.py
│       ├── run_conversion.py
│       ├── generate_xps_graphs.py
│       └── merge_xps_to_excel.py
└── LICENSE
```

Local-only (ignored by `.gitignore`):

- `Data/`
- `Archive/`
- `.runtime/`

## Quick Start

### macOS one-click launch

From repo root:

```bash
./Launch_XPS_Tool.command
```

App URL:

```text
http://127.0.0.1:8501
```

### Manual run

```bash
cd "XPS Tool/Main Program"
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m streamlit run app.py --server.address 127.0.0.1 --server.port 8501
```

### Stop background app

```bash
kill "$(cat '.runtime/streamlit.pid')"
```

## Notes

- Paths include spaces (`XPS Tool/...`), so quote paths in shell commands.
- For current updates, see [CHANGELOG.md](CHANGELOG.md).
- Legacy in-program docs remain available:
  - `XPS Tool/Main Program/README.md`
  - `XPS Tool/Main Program/CHANGELOG.md`

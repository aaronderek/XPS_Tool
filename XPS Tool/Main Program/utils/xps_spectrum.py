"""
XPS spectrum utilities for loading, peak detection, fitting, and plotting.
"""

from __future__ import annotations

import io
import os
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path
from typing import Any, Dict, List, Sequence

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, savgol_filter


ENERGY_COLUMN_CANDIDATES = (
    "e",
    "energy",
    "binding_energy",
    "binding energy",
    "be",
    "x",
)

INTENSITY_COLUMN_CANDIDATES = (
    "y",
    "intensity",
    "counts",
    "count",
    "cps",
    "signal",
)


ELEMENT_REFERENCE_PEAKS = [
    {"label": "C 1s", "energy_eV": 284.8},
    {"label": "O 1s", "energy_eV": 530.1},
    {"label": "N 1s", "energy_eV": 399.8},
    {"label": "Si 2p", "energy_eV": 99.4},
    {"label": "Si 2s", "energy_eV": 150.5},
    {"label": "In 3d5/2", "energy_eV": 444.7},
    {"label": "In 3d3/2", "energy_eV": 452.3},
    {"label": "In 4d", "energy_eV": 17.8},
    {"label": "Ti 2p3/2", "energy_eV": 458.6},
    {"label": "Ti 2p1/2", "energy_eV": 464.3},
    {"label": "Sn 3d5/2", "energy_eV": 486.6},
    {"label": "Sn 3d3/2", "energy_eV": 495.0},
    {"label": "Zn 2p3/2", "energy_eV": 1021.8},
    {"label": "Zn 2p1/2", "energy_eV": 1044.9},
    {"label": "S 2p", "energy_eV": 164.0},
    {"label": "S 2s", "energy_eV": 226.0},
    {"label": "Al 2p", "energy_eV": 74.4},
    {"label": "Ga 3d", "energy_eV": 20.1},
]


def _find_best_column(columns: Sequence[str], candidates: Sequence[str]) -> str | None:
    normalized_to_original = {c.strip().lower(): c for c in columns}

    for cand in candidates:
        if cand in normalized_to_original:
            return normalized_to_original[cand]

    for normalized, original in normalized_to_original.items():
        for cand in candidates:
            if cand in normalized:
                return original

    return None


def _flatten_dataframe_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if isinstance(out.columns, pd.MultiIndex):
        new_cols = []
        for col in out.columns.to_flat_index():
            if isinstance(col, tuple):
                parts = [str(x) for x in col if str(x) not in ("", "None")]
                new_cols.append("_".join(parts) if parts else "col")
            else:
                new_cols.append(str(col))
        out.columns = new_cols
    return out


def _spectrum_from_dataframe(df: pd.DataFrame, filename: str = "") -> Dict[str, Any]:
    out = _flatten_dataframe_columns(df)
    if isinstance(out.index, pd.MultiIndex) or out.index.name is not None:
        out = out.reset_index()

    energy_col = _find_best_column(out.columns, ENERGY_COLUMN_CANDIDATES)
    intensity_col = _find_best_column(out.columns, INTENSITY_COLUMN_CANDIDATES)

    if energy_col is None or intensity_col is None:
        numeric_cols = [c for c in out.columns if pd.api.types.is_numeric_dtype(out[c])]
        if len(numeric_cols) >= 2:
            energy_col = numeric_cols[0]
            intensity_col = numeric_cols[1]

    if energy_col is None or intensity_col is None:
        raise ValueError(
            "Unable to detect energy/intensity columns. Expected columns similar to E,y."
        )

    energy = pd.to_numeric(out[energy_col], errors="coerce").to_numpy(dtype=float)
    intensity = pd.to_numeric(out[intensity_col], errors="coerce").to_numpy(dtype=float)

    valid = np.isfinite(energy) & np.isfinite(intensity)
    energy = energy[valid]
    intensity = intensity[valid]

    if energy.size < 20:
        raise ValueError("Not enough valid points for XPS analysis (need at least 20).")

    sort_idx = np.argsort(energy)
    energy = energy[sort_idx]
    intensity = intensity[sort_idx]

    return {
        "filename": filename,
        "energy_col": energy_col,
        "intensity_col": intensity_col,
        "energy": energy,
        "intensity": intensity,
        "n_points": int(energy.size),
        "energy_min": float(np.min(energy)),
        "energy_max": float(np.max(energy)),
    }


def _find_converter_context(start: Path) -> tuple[Path, Path] | None:
    for candidate in [start] + list(start.parents):
        new_script = candidate / "XPS Tool" / "XPS Plugin" / "convert_spe_to_csv.py"
        if new_script.exists():
            return candidate, (candidate / "XPS Tool" / "XPS Plugin")

        legacy_script = candidate / "XPS data" / "convert_spe_to_csv.py"
        if legacy_script.exists():
            return candidate, (candidate / "XPS data")
    return None


def _local_camel_to_snake(s: str) -> str:
    import re

    s = re.sub(r"(.)([A-Z][a-z]+)", r"\1_\2", s)
    return re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", s).lower()


def _local_parse_header_lines(lines: Sequence[bytes]) -> tuple[Dict[str, Any], int]:
    sofh_idx = None
    eofh_idx = None
    for i, line in enumerate(lines):
        token = line.strip()
        if token == b"SOFH":
            sofh_idx = i
        elif token == b"EOFH":
            eofh_idx = i
            break

    if sofh_idx is None or eofh_idx is None or eofh_idx <= sofh_idx:
        raise ValueError("Unable to locate SOFH/EOFH header markers in .spe file.")

    header: Dict[str, Any] = {}
    for line in lines[sofh_idx + 1 : eofh_idx]:
        if b":" not in line:
            continue
        key_raw, value_raw = line.split(b":", 1)
        key = _local_camel_to_snake(key_raw.decode(errors="ignore").strip())
        value = value_raw.decode(errors="ignore").strip()
        if key in header:
            if isinstance(header[key], list):
                header[key].append(value)
            else:
                header[key] = [header[key], value]
        else:
            header[key] = value

    return header, eofh_idx


def _local_parse_trace_defs(header: Dict[str, Any]) -> List[Dict[str, Any]]:
    defs_raw = header.get("spectral_reg_def")
    if defs_raw is None:
        raise ValueError("Header missing 'spectral_reg_def'.")

    if isinstance(defs_raw, str):
        defs = [defs_raw]
    else:
        defs = list(defs_raw)

    out = []
    for reg_def in defs:
        parts = str(reg_def).split()
        if len(parts) < 13:
            continue
        out.append(
            {
                "trace_number": int(parts[0]),
                "name": parts[2],
                "num_datapoints": int(parts[4]),
                "start": float(parts[6]),
                "stop": float(parts[7]),
            }
        )

    if not out:
        raise ValueError("No valid spectral_reg_def entries parsed.")
    return out


def _extract_traces_with_local_phi_parser(input_path: Path) -> List[Dict[str, Any]]:
    """
    Minimal local parser for PHI .spe files.
    Handles common non-standard files where:
    - header uses CRLF
    - trace_header.end_of_data is 0 (sequential data layout)
    """
    spe_bytes = input_path.read_bytes()
    lines = spe_bytes.splitlines(keepends=True)
    header, eofh_idx = _local_parse_header_lines(lines)
    trace_defs = _local_parse_trace_defs(header)

    binary_data = b"".join(lines[eofh_idx + 1 :])

    data_header_dtype = np.dtype(
        [
            ("group", "<u4"),
            ("num_traces", "<u4"),
            ("trace_header_size", "<u4"),
            ("data_header_size", "<u4"),
        ]
    )
    trace_header_dtype = np.dtype(
        [
            ("trace_number", "<u4"),
            ("bool_01", "<u4"),
            ("bool_02", "<u4"),
            ("trace_number_again", "<u4"),
            ("bool_03", "<u4"),
            ("num_datapoints", "<u4"),
            ("bool_04", "<u4"),
            ("bool_05", "<u4"),
            ("string_01", "|S4"),
            ("string_02", "|S4"),
            ("string_03", "|S4"),
            ("int_02", "<u4"),
            ("string_04", "|S4"),
            ("string_05", "|S4"),
            ("y_unit", "|S4"),
            ("int_05", "<u4"),
            ("int_06", "<u4"),
            ("int_07", "<u4"),
            ("data_dtype", "|S4"),
            ("num_data_bytes", "<u4"),
            ("num_datapoints_tot", "<u4"),
            ("int_10", "<u4"),
            ("int_11", "<u4"),
            ("end_of_data", "<u4"),
        ]
    )

    if len(binary_data) < data_header_dtype.itemsize:
        raise ValueError("Binary section too small for PHI data header.")

    data_header = np.frombuffer(binary_data, dtype=data_header_dtype, count=1)[0]
    num_traces = int(data_header["num_traces"])
    data_header_size = int(data_header["data_header_size"])
    trace_header_size = int(data_header["trace_header_size"])

    if num_traces <= 0:
        raise ValueError("PHI parser: num_traces <= 0.")

    trace_headers = np.frombuffer(
        binary_data,
        offset=data_header_size,
        dtype=trace_header_dtype,
        count=num_traces,
    )

    defs_by_number = {int(d["trace_number"]): d for d in trace_defs}
    sequential_offset = data_header_size + trace_header_size

    spectra: List[Dict[str, Any]] = []
    for th in trace_headers:
        trace_number = int(th["trace_number"])
        trace_def = defs_by_number.get(trace_number)
        if trace_def is None:
            continue

        num_datapoints = int(th["num_datapoints"])
        num_data_bytes = int(th["num_data_bytes"])
        end_of_data = int(th["end_of_data"])

        raw_dtype = th["data_dtype"].decode(errors="ignore").strip().strip("\x00")
        if raw_dtype not in ("f4", "f8"):
            raw_dtype = "f4"
        data_dtype = np.dtype(raw_dtype)

        if (
            end_of_data > 0
            and end_of_data >= num_data_bytes
            and end_of_data <= len(binary_data)
        ):
            data_offset = end_of_data - num_data_bytes
        else:
            data_offset = sequential_offset

        if data_offset < 0 or data_offset + num_data_bytes > len(binary_data):
            raise ValueError(
                f"Trace {trace_number}: invalid data offset/size "
                f"(offset={data_offset}, bytes={num_data_bytes}, total={len(binary_data)})."
            )

        datapoints = np.frombuffer(
            binary_data,
            offset=data_offset,
            dtype=data_dtype,
            count=num_datapoints,
        )
        sequential_offset = data_offset + num_data_bytes

        energies = np.linspace(
            trace_def["start"],
            trace_def["stop"],
            num_datapoints,
            endpoint=True,
        )
        df = pd.DataFrame({"E": energies, "y": datapoints})
        spectrum = _spectrum_from_dataframe(
            df, filename=f"{input_path.stem}_{trace_def['name']}.csv"
        )
        spectrum["trace_name"] = str(trace_def["name"])
        spectra.append(spectrum)

    if not spectra:
        raise ValueError("Local PHI parser produced no traces.")
    return spectra


def _extract_traces_with_yadg(input_path: Path) -> List[Dict[str, Any]]:
    from yadg.extractors import extract

    ds = extract(filetype="phi.spe", path=str(input_path))
    traces: List[Dict[str, Any]] = []

    for trace_name in ds.keys():
        trace_node = ds[trace_name]
        if hasattr(trace_node, "ds"):
            trace_df = trace_node.ds.to_dataframe()
        elif hasattr(trace_node, "to_dataframe"):
            trace_df = trace_node.to_dataframe()
        else:
            continue

        spectrum = _spectrum_from_dataframe(
            trace_df, filename=f"{input_path.stem}_{trace_name}.csv"
        )
        spectrum["trace_name"] = str(trace_name)
        traces.append(spectrum)

    return traces


def _extract_traces_by_script(input_path: Path) -> tuple[List[Dict[str, Any]], str]:
    here = Path(__file__).resolve()
    converter_context = _find_converter_context(here)
    if converter_context is None:
        raise RuntimeError(
            "Could not locate converter script in either "
            "'XPS Tool/XPS Plugin' or legacy 'XPS data'."
        )
    _, xps_plugin_dir = converter_context

    converter_script = xps_plugin_dir / "convert_spe_to_csv.py"
    if not converter_script.exists():
        raise RuntimeError(f"Converter script not found: {converter_script}")

    py_candidates: List[str] = []

    for pattern in [
        ".venv*/Scripts/python.exe",
        ".venv*/Scripts/python*.exe",
        ".venv*/bin/python",
        ".venv*/bin/python3",
    ]:
        for candidate in sorted(xps_plugin_dir.glob(pattern)):
            py_candidates.append(str(candidate))

    py_candidates.append(str(xps_plugin_dir / ".venv311" / "Scripts" / "python.exe"))
    py_candidates.append(sys.executable)
    py_candidates.extend(["python3", "python"])

    # Keep order while deduplicating
    py_candidates = list(dict.fromkeys(py_candidates))

    run_log_lines = []
    with tempfile.TemporaryDirectory(prefix="xps_spe_convert_") as tmp_dir:
        tmp_path = Path(tmp_dir)
        local_spe = tmp_path / input_path.name
        local_spe.write_bytes(input_path.read_bytes())

        success = False
        for py_exe in py_candidates:
            if not py_exe:
                continue
            if os.path.isabs(py_exe) and not os.path.exists(py_exe):
                continue
            try:
                result = subprocess.run(
                    [py_exe, str(converter_script)],
                    cwd=str(tmp_path),
                    capture_output=True,
                    text=True,
                    timeout=180,
                )
            except Exception as exc:
                run_log_lines.append(f"[{py_exe}] failed to execute: {exc}")
                continue

            run_log_lines.append(f"[{py_exe}] returncode={result.returncode}")
            if result.stdout:
                run_log_lines.append(result.stdout.strip())
            if result.stderr:
                run_log_lines.append(result.stderr.strip())

            csv_files = sorted(tmp_path.glob(f"{local_spe.stem}_*.csv"))
            if result.returncode == 0 and csv_files:
                success = True
                break

        if not success:
            raise RuntimeError(
                "Fallback conversion script failed. "
                "Ensure patched yadg env is available in "
                f"'{xps_plugin_dir / '.venv311'}'."
            )

        traces: List[Dict[str, Any]] = []
        prefix = f"{local_spe.stem}_"
        for csv_file in sorted(tmp_path.glob(f"{local_spe.stem}_*.csv")):
            spec = load_xps_spectrum_data(csv_file.read_bytes(), filename=csv_file.name)
            trace_name = csv_file.stem[len(prefix):] if csv_file.stem.startswith(prefix) else csv_file.stem
            spec["trace_name"] = trace_name
            traces.append(spec)

    return traces, "\n".join([x for x in run_log_lines if x])


def extract_spe_traces(file_content: bytes, filename: str) -> Dict[str, Any]:
    """
    Convert a .spe upload into one or more trace spectra.

    Strategy:
    1) Use in-process yadg.extractors.extract (same method as convert_spe_to_csv.py)
    2) Fallback to running convert_spe_to_csv.py from the local plugin directory
    """
    logs: List[str] = []

    if not filename.lower().endswith(".spe"):
        raise ValueError("Input file is not .spe format.")

    with tempfile.TemporaryDirectory(prefix="xps_spe_upload_") as tmp_dir:
        tmp_path = Path(tmp_dir)
        spe_path = tmp_path / filename
        spe_path.write_bytes(file_content)

        try:
            traces = _extract_traces_with_yadg(spe_path)
            if traces:
                logs.append("Converted with in-process yadg.extractors.extract.")
                return {"spectra": traces, "logs": "\n".join(logs), "method": "yadg.extract"}
            logs.append("In-process yadg returned no traces.")
        except Exception as exc:
            logs.append(f"In-process yadg conversion failed: {exc}")

        try:
            traces = _extract_traces_with_local_phi_parser(spe_path)
            if traces:
                logs.append("Converted with built-in local PHI parser fallback.")
                return {"spectra": traces, "logs": "\n".join(logs), "method": "local_phi_parser"}
            logs.append("Local PHI parser returned no traces.")
        except Exception as exc:
            logs.append(f"Local PHI parser conversion failed: {exc}")

        try:
            traces, script_log = _extract_traces_by_script(spe_path)
            if script_log:
                logs.append(script_log)
            return {"spectra": traces, "logs": "\n".join(logs), "method": "convert_spe_to_csv.py"}
        except Exception as exc:
            logs.append(f"Fallback conversion script failed: {exc}")
            raise RuntimeError("\n".join(logs))


def load_uploaded_xps_file(file_content: bytes, filename: str) -> Dict[str, Any]:
    """
    Load CSV directly or auto-convert .spe into spectra traces.
    """
    lower = (filename or "").lower()

    if lower.endswith(".csv"):
        spec = load_xps_spectrum_data(file_content, filename=filename)
        spec["trace_name"] = Path(filename).stem
        return {
            "file_type": "csv",
            "spectra": [spec],
            "logs": "Loaded CSV spectrum directly.",
            "method": "csv",
        }

    if lower.endswith(".spe"):
        spe_result = extract_spe_traces(file_content=file_content, filename=filename)
        return {
            "file_type": "spe",
            "spectra": spe_result["spectra"],
            "logs": spe_result.get("logs", ""),
            "method": spe_result.get("method", "spe"),
        }

    raise ValueError("Unsupported file format. Please upload .csv or .spe.")


def build_spectra_zip_bytes(spectra: Sequence[Dict[str, Any]], base_name: str = "xps") -> bytes:
    """
    Package converted spectra into a ZIP containing one CSV per trace.
    """
    bio = io.BytesIO()
    with zipfile.ZipFile(bio, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
        for idx, spec in enumerate(spectra):
            trace = str(spec.get("trace_name") or f"trace_{idx+1}")
            safe_trace = "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in trace)
            csv_name = f"{base_name}_{safe_trace}.csv"
            df = pd.DataFrame({"E": spec["energy"], "y": spec["intensity"]})
            zf.writestr(csv_name, df.to_csv(index=False))
    return bio.getvalue()


def load_xps_spectrum_data(file_content: bytes, filename: str = "") -> Dict[str, Any]:
    """
    Load and sanitize XPS spectrum data from CSV bytes.
    """
    df = pd.read_csv(io.BytesIO(file_content))
    if df.empty:
        raise ValueError("CSV is empty.")
    return _spectrum_from_dataframe(df, filename=filename)


def _safe_savgol(y: np.ndarray, window_length: int) -> np.ndarray:
    if y.size < 7:
        return y.copy()

    max_window = y.size - 1 if y.size % 2 == 0 else y.size
    win = int(min(window_length, max_window))
    if win % 2 == 0:
        win -= 1
    if win < 5:
        return y.copy()

    polyorder = min(3, win - 2)
    return savgol_filter(y, window_length=win, polyorder=polyorder)


def estimate_linear_background(
    energy: np.ndarray,
    intensity: np.ndarray,
    edge_fraction: float = 0.12,
) -> Dict[str, Any]:
    """
    Estimate a linear background using both spectrum edges.
    """
    if energy.size < 10:
        raise ValueError("Need at least 10 points to estimate background.")

    edge_fraction = float(np.clip(edge_fraction, 0.05, 0.45))
    n_edge = max(5, int(edge_fraction * energy.size))
    n_edge = min(n_edge, energy.size // 2)

    x_edges = np.concatenate([energy[:n_edge], energy[-n_edge:]])
    y_edges = np.concatenate([intensity[:n_edge], intensity[-n_edge:]])

    slope, intercept = np.polyfit(x_edges, y_edges, deg=1)
    baseline = slope * energy + intercept

    return {
        "baseline": baseline,
        "slope": float(slope),
        "intercept": float(intercept),
        "n_edge": int(n_edge),
    }


def detect_spectrum_peaks(
    energy: np.ndarray,
    intensity: np.ndarray,
    edge_fraction: float = 0.12,
    smoothing_window: int = 11,
    prominence_ratio: float = 0.06,
    min_distance_points: int = 10,
    max_peaks: int = 4,
) -> Dict[str, Any]:
    """
    Detect prominent XPS peaks after linear background subtraction.
    """
    background = estimate_linear_background(energy, intensity, edge_fraction=edge_fraction)
    y_corr = intensity - background["baseline"]
    y_corr_smoothed = _safe_savgol(y_corr, window_length=smoothing_window)

    signal = np.clip(y_corr_smoothed, a_min=0.0, a_max=None)
    span = float(np.ptp(signal))
    if span <= 0:
        return {
            "background": background,
            "corrected": y_corr,
            "smoothed": y_corr_smoothed,
            "peaks": pd.DataFrame(
                columns=["index", "center_eV", "height", "prominence", "raw_intensity"]
            ),
        }

    prominence = max(1e-12, span * float(np.clip(prominence_ratio, 0.005, 0.8)))
    distance = max(1, int(min_distance_points))
    peak_indices, props = find_peaks(signal, prominence=prominence, distance=distance)

    if peak_indices.size == 0:
        peak_df = pd.DataFrame(
            columns=["index", "center_eV", "height", "prominence", "raw_intensity"]
        )
    else:
        peak_df = pd.DataFrame(
            {
                "index": peak_indices.astype(int),
                "center_eV": energy[peak_indices],
                "height": signal[peak_indices],
                "prominence": props["prominences"],
                "raw_intensity": intensity[peak_indices],
            }
        )
        peak_df = peak_df.sort_values("prominence", ascending=False).head(int(max_peaks))
        peak_df = peak_df.sort_values("center_eV", ascending=False).reset_index(drop=True)

    return {
        "background": background,
        "corrected": y_corr,
        "smoothed": y_corr_smoothed,
        "peaks": peak_df,
    }


def _gaussian(x: np.ndarray, amp: float, center: float, width: float) -> np.ndarray:
    sigma = max(width, 1e-6)
    return amp * np.exp(-0.5 * ((x - center) / sigma) ** 2)


def _lorentzian(x: np.ndarray, amp: float, center: float, width: float) -> np.ndarray:
    gamma = max(width, 1e-6)
    return amp * (gamma**2) / ((x - center) ** 2 + gamma**2)


def _get_profile(name: str):
    return _lorentzian if name.lower().startswith("lorentz") else _gaussian


def fit_spectrum_peaks(
    energy: np.ndarray,
    intensity: np.ndarray,
    initial_centers: Sequence[float],
    profile: str = "Gaussian",
    fit_range: Sequence[float] | None = None,
    edge_fraction: float = 0.12,
) -> Dict[str, Any]:
    """
    Fit multi-peak XPS spectrum using Gaussian or Lorentzian components.
    """
    if len(initial_centers) == 0:
        raise ValueError("No initial peak centers provided for fitting.")

    x = np.asarray(energy, dtype=float)
    y = np.asarray(intensity, dtype=float)

    if fit_range is not None:
        fit_min = min(float(fit_range[0]), float(fit_range[1]))
        fit_max = max(float(fit_range[0]), float(fit_range[1]))
        mask = (x >= fit_min) & (x <= fit_max)
        x = x[mask]
        y = y[mask]
        if x.size < 20:
            raise ValueError("Selected fit range is too narrow (need at least 20 points).")

    bg = estimate_linear_background(x, y, edge_fraction=edge_fraction)
    y_corr = np.clip(y - bg["baseline"], a_min=0.0, a_max=None)
    profile_fn = _get_profile(profile)

    centers = np.asarray(initial_centers, dtype=float)
    centers = centers[np.isfinite(centers)]
    x_min = float(np.min(x))
    x_max = float(np.max(x))
    centers = centers[(centers >= x_min) & (centers <= x_max)]

    if centers.size == 0:
        # Fallback to the strongest point in the selected fit range.
        centers = np.asarray([float(x[np.argmax(y_corr)])], dtype=float)

    centers = np.sort(np.unique(np.round(centers, decimals=6)))
    n_peaks = centers.size

    p0 = [bg["slope"], bg["intercept"]]
    lower = [-np.inf, -np.inf]
    upper = [np.inf, np.inf]

    y_span = max(float(np.ptp(y_corr)), 1.0)
    x_span = max(x_max - x_min, 1e-6)
    width_guess = max(0.05, min(0.8, 0.01 * x_span))
    width_upper = max(0.08, min(3.0, 0.5 * x_span))
    center_window = max(0.2, min(1.5, 0.3 * x_span))

    for center in centers:
        idx = int(np.argmin(np.abs(x - center)))
        amp_guess = max(y_corr[idx], 0.05 * y_span)
        c_low = max(x_min, center - center_window)
        c_high = min(x_max, center + center_window)
        if c_high <= c_low:
            c_low = max(x_min, center - 0.05)
            c_high = min(x_max, center + 0.05)
        if c_high <= c_low:
            c_low = x_min
            c_high = x_max

        amp_guess = float(np.clip(amp_guess, 0.0, 10.0 * y_span))
        center_guess = float(np.clip(center, c_low, c_high))
        width_guess_i = float(np.clip(width_guess, 0.05, width_upper))

        p0.extend([amp_guess, center_guess, width_guess_i])
        lower.extend([0.0, c_low, 0.05])
        upper.extend([10.0 * y_span, c_high, width_upper])

    def model(x_input: np.ndarray, *params: float) -> np.ndarray:
        slope, intercept = params[0], params[1]
        out = slope * x_input + intercept
        for i in range(n_peaks):
            amp = params[2 + 3 * i]
            center = params[3 + 3 * i]
            width = params[4 + 3 * i]
            out = out + profile_fn(x_input, amp, center, width)
        return out

    # Final guard: ensure initial values are within bounds.
    p0 = np.asarray(p0, dtype=float)
    lower_arr = np.asarray(lower, dtype=float)
    upper_arr = np.asarray(upper, dtype=float)
    eps = 1e-9
    p0 = np.minimum(np.maximum(p0, lower_arr + eps), upper_arr - eps)

    popt, _ = curve_fit(
        model,
        x,
        y,
        p0=p0,
        bounds=(lower_arr, upper_arr),
        maxfev=30000,
    )

    y_fit = model(x, *popt)
    y_baseline = popt[0] * x + popt[1]
    residuals = y - y_fit

    ss_res = float(np.sum(residuals**2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    rmse = float(np.sqrt(np.mean(residuals**2)))

    component_rows = []
    components = []
    for i in range(n_peaks):
        amp = float(popt[2 + 3 * i])
        center = float(popt[3 + 3 * i])
        width = float(popt[4 + 3 * i])
        comp_y = profile_fn(x, amp, center, width)
        components.append(
            {
                "label": f"Peak {i + 1}",
                "x": x,
                "y": comp_y,
                "amplitude": amp,
                "center_eV": center,
                "width_eV": width,
            }
        )

        if profile_fn is _gaussian:
            area = amp * width * np.sqrt(2.0 * np.pi)
        else:
            area = amp * np.pi * width

        component_rows.append(
            {
                "peak": f"Peak {i + 1}",
                "center_eV": center,
                "amplitude": amp,
                "width_eV": width,
                "area_au_eV": float(area),
            }
        )

    peak_table = pd.DataFrame(component_rows).sort_values("center_eV", ascending=False).reset_index(drop=True)

    return {
        "profile": profile,
        "x": x,
        "y": y,
        "fit_y": y_fit,
        "baseline_y": y_baseline,
        "residuals": residuals,
        "r_squared": r_squared,
        "rmse": rmse,
        "components": components,
        "peak_table": peak_table,
        "params": popt,
    }


def get_element_reference_library() -> List[Dict[str, float | str]]:
    """
    Return built-in element reference peak positions.
    """
    return ELEMENT_REFERENCE_PEAKS.copy()


def suggest_reference_labels(filename: str) -> List[str]:
    """
    Suggest element lines from filename hints.
    """
    name = (filename or "").lower()
    suggestions: List[str] = []
    token_map = {
        "c1s": "C 1s",
        "o1s": "O 1s",
        "in3d5": "In 3d5/2",
        "in3d3": "In 3d3/2",
        "in3d": "In 3d5/2",
        "si2p": "Si 2p",
        "si2s": "Si 2s",
        "n1s": "N 1s",
        "s2p": "S 2p",
        "ti2p": "Ti 2p3/2",
        "sn3d": "Sn 3d5/2",
        "zn2p": "Zn 2p3/2",
    }

    for token, label in token_map.items():
        if token in name and label not in suggestions:
            suggestions.append(label)

    if not suggestions:
        suggestions = ["C 1s", "O 1s", "In 3d5/2", "In 3d3/2"]
    return suggestions


def match_fitted_peaks(
    peak_table: pd.DataFrame,
    selected_reference_peaks: Sequence[Dict[str, Any]],
    tolerance_eV: float = 1.2,
) -> pd.DataFrame:
    """
    Match fitted peak positions to selected element reference lines.
    """
    if peak_table.empty:
        out = peak_table.copy()
        out["matched_peak"] = ""
        out["offset_eV"] = np.nan
        return out

    out = peak_table.copy()
    out["matched_peak"] = ""
    out["offset_eV"] = np.nan

    if not selected_reference_peaks:
        return out

    refs = [
        (str(item["label"]), float(item["energy_eV"]))
        for item in selected_reference_peaks
        if "label" in item and "energy_eV" in item
    ]
    if not refs:
        return out

    tolerance = float(max(tolerance_eV, 0.01))
    for idx, row in out.iterrows():
        center = float(row["center_eV"])
        best_label = ""
        best_offset = np.nan
        best_abs = float("inf")
        for label, ref_energy in refs:
            delta = center - ref_energy
            if abs(delta) < best_abs and abs(delta) <= tolerance:
                best_abs = abs(delta)
                best_label = label
                best_offset = delta
        out.at[idx, "matched_peak"] = best_label
        out.at[idx, "offset_eV"] = best_offset

    return out


def create_xps_spectrum_figure(
    energy: np.ndarray,
    intensity: np.ndarray,
    detected_peaks: pd.DataFrame | None = None,
    fit_result: Dict[str, Any] | None = None,
    reference_peaks: Sequence[Dict[str, Any]] | None = None,
    title: str = "XPS Spectrum",
    show_detected_peaks: bool = True,
    show_fit_curve: bool = True,
    show_baseline: bool = True,
    show_components: bool = True,
    show_reference_lines: bool = True,
    show_reference_labels: bool = True,
    show_peak_name_labels: bool = False,
    include_legend: bool = True,
    base_font_size: int = 16,
    title_font_size: int = 24,
    figure_height: int = 760,
    focus_range: Sequence[float] | None = None,
    line_width_scale: float = 0.85,
) -> go.Figure:
    """
    Create an annotated XPS spectrum figure.
    """
    fig = go.Figure()
    width_scale = float(np.clip(line_width_scale, 0.5, 2.0))

    fig.add_trace(
        go.Scatter(
            x=energy,
            y=intensity,
            mode="lines",
            line=dict(color="#1f2a44", width=1.6 * width_scale),
            name="Raw spectrum",
            hovertemplate="E: %{x:.3f} eV<br>I: %{y:.3f}<extra></extra>",
        )
    )

    y_min = float(np.min(intensity))
    y_max = float(np.max(intensity))
    y_span = max(y_max - y_min, 1e-9)

    if show_detected_peaks and detected_peaks is not None and not detected_peaks.empty:
        peak_x = detected_peaks["center_eV"].to_numpy(dtype=float)
        peak_y = np.interp(peak_x, energy, intensity)
        fig.add_trace(
            go.Scatter(
                x=peak_x,
                y=peak_y,
                mode="markers",
                marker=dict(color="#d94801", size=max(6.0, 9.0 * width_scale), symbol="x"),
                name="Detected peaks",
                hovertemplate="Detected peak<br>E: %{x:.3f} eV<extra></extra>",
            )
        )

    if fit_result is not None and show_fit_curve:
        x_fit = fit_result["x"]
        if show_baseline:
            fig.add_trace(
                go.Scatter(
                    x=x_fit,
                    y=fit_result["baseline_y"],
                    mode="lines",
                    line=dict(color="#7f7f7f", width=1.2 * width_scale, dash="dot"),
                    name="Baseline",
                    hovertemplate="Baseline<br>E: %{x:.3f} eV<br>I: %{y:.3f}<extra></extra>",
                )
            )
        fig.add_trace(
            go.Scatter(
                x=x_fit,
                y=fit_result["fit_y"],
                mode="lines",
                line=dict(color="#006d2c", width=2.1 * width_scale),
                name=f"Fit ({fit_result['profile']})",
                hovertemplate="Fit<br>E: %{x:.3f} eV<br>I: %{y:.3f}<extra></extra>",
            )
        )

        if show_components:
            for component in fit_result["components"]:
                fig.add_trace(
                    go.Scatter(
                        x=component["x"],
                        y=component["y"] + fit_result["baseline_y"],
                        mode="lines",
                        line=dict(width=1.0 * width_scale, dash="dash"),
                        name=component["label"],
                        opacity=0.75,
                        hovertemplate=(
                            f"{component['label']}<br>"
                            "E: %{x:.3f} eV<br>I: %{y:.3f}<extra></extra>"
                        ),
                    )
                )

        if show_peak_name_labels and "peak_table" in fit_result and not fit_result["peak_table"].empty:
            peak_table = fit_result["peak_table"]
            fit_y = fit_result["fit_y"]
            for _, row in peak_table.head(16).iterrows():
                center = float(row["center_eV"])
                y_here = float(np.interp(center, x_fit, fit_y))
                matched = str(row.get("matched_peak", "")).strip()
                label = matched if matched else str(row.get("peak", "Peak"))
                fig.add_annotation(
                    x=center,
                    y=y_here + 0.03 * y_span,
                    text=label,
                    showarrow=True,
                    arrowhead=2,
                    arrowsize=1,
                    arrowwidth=1,
                    arrowcolor="rgba(60,60,60,0.7)",
                    font=dict(size=max(11, int(base_font_size * 0.85)), color="#202020"),
                    bgcolor="rgba(255,255,255,0.72)",
                    bordercolor="rgba(120,120,120,0.45)",
                    borderwidth=1,
                )

    if show_reference_lines and reference_peaks:
        ref_count = 0
        for item in reference_peaks:
            if "energy_eV" not in item:
                continue
            x_ref = float(item["energy_eV"])
            if x_ref < float(np.min(energy)) or x_ref > float(np.max(energy)):
                continue
            fig.add_vline(
                x=x_ref,
                line_width=max(0.8, 1.0 * width_scale),
                line_dash="dot",
                line_color="rgba(217, 95, 2, 0.55)",
            )
            if show_reference_labels and ref_count < 12:
                fig.add_annotation(
                    x=x_ref,
                    y=y_min + 0.96 * y_span,
                    text=str(item.get("label", "")),
                    showarrow=False,
                    textangle=-90,
                    font=dict(size=max(10, int(base_font_size * 0.9)), color="#7f2704"),
                    yanchor="top",
                )
            ref_count += 1

    fig.update_layout(
        title=dict(text=title, font=dict(size=title_font_size)),
        template="plotly_white",
        hovermode="closest",
        xaxis_title="Binding Energy (eV)",
        yaxis_title="Intensity (a.u.)",
        font=dict(size=base_font_size),
        hoverlabel=dict(font=dict(size=max(11, base_font_size - 1))),
        height=int(figure_height),
        margin=dict(l=90, r=30, t=100, b=90),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0.0,
            font=dict(size=max(10, base_font_size - 1)),
            itemwidth=60,
        ),
        showlegend=include_legend,
    )
    fig.update_xaxes(
        autorange="reversed",
        showgrid=True,
        gridcolor="LightGray",
        title_font=dict(size=max(12, base_font_size + 4)),
        tickfont=dict(size=max(10, base_font_size)),
        automargin=True,
    )
    fig.update_yaxes(
        showgrid=True,
        gridcolor="LightGray",
        title_font=dict(size=max(12, base_font_size + 4)),
        tickfont=dict(size=max(10, base_font_size)),
        automargin=True,
    )

    if focus_range is not None and len(focus_range) == 2:
        x0 = float(min(focus_range))
        x1 = float(max(focus_range))
        if x1 > x0:
            fig.update_xaxes(range=[x1, x0])

    return fig


def export_xps_figure_bytes(
    fig: go.Figure,
    export_format: str = "png",
    width_px: int = 2000,
    height_px: int = 1200,
    dpi: int = 300,
) -> bytes:
    """
    Export a Plotly XPS figure to image bytes.
    """
    fmt = str(export_format).strip().lower().replace(".", "")
    if fmt not in {"png", "svg", "pdf"}:
        raise ValueError("Unsupported export format. Use png/svg/pdf.")

    width_px = int(max(800, width_px))
    height_px = int(max(500, height_px))

    render_kwargs: Dict[str, Any] = {
        "format": fmt,
        "width": width_px,
        "height": height_px,
    }

    if fmt == "png":
        scale = max(1.0, float(dpi) / 120.0)
        render_kwargs["scale"] = min(scale, 5.0)

    try:
        return fig.to_image(**render_kwargs)
    except Exception as exc:
        # Fallback path: render with matplotlib to avoid Chrome/kaleido runtime issues.
        try:
            return _export_plotly_figure_with_matplotlib(
                fig=fig,
                export_format=fmt,
                width_px=width_px,
                height_px=height_px,
                dpi=dpi,
            )
        except Exception as fallback_exc:
            raise RuntimeError(
                "Figure export failed with both Plotly and Matplotlib fallback. "
                f"Plotly error: {exc} | Fallback error: {fallback_exc}"
            ) from fallback_exc


def _plotly_color_to_mpl(color: Any) -> Any:
    if color is None:
        return None
    if not isinstance(color, str):
        return color

    c = color.strip()
    if c.startswith("rgba(") and c.endswith(")"):
        vals = [v.strip() for v in c[5:-1].split(",")]
        if len(vals) == 4:
            r, g, b = [float(vals[i]) / 255.0 for i in range(3)]
            a = float(vals[3])
            return (r, g, b, a)
    if c.startswith("rgb(") and c.endswith(")"):
        vals = [v.strip() for v in c[4:-1].split(",")]
        if len(vals) == 3:
            r, g, b = [float(vals[i]) / 255.0 for i in range(3)]
            return (r, g, b)
    return c


def _dash_to_linestyle(dash: Any) -> str:
    if not dash:
        return "-"
    d = str(dash).lower()
    if d == "dot":
        return ":"
    if d == "dash":
        return "--"
    if d == "dashdot":
        return "-."
    return "-"


def _export_plotly_figure_with_matplotlib(
    fig: go.Figure,
    export_format: str,
    width_px: int,
    height_px: int,
    dpi: int,
) -> bytes:
    import matplotlib.pyplot as plt

    fmt = export_format.lower()
    if fmt not in {"png", "svg", "pdf"}:
        raise ValueError("Matplotlib fallback: unsupported format.")

    figsize = (max(400, width_px) / max(72, dpi), max(300, height_px) / max(72, dpi))
    mfig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # Draw traces
    for trace in fig.data:
        x = np.asarray(trace.x, dtype=float) if trace.x is not None else None
        y = np.asarray(trace.y, dtype=float) if trace.y is not None else None
        if x is None or y is None or x.size == 0 or y.size == 0:
            continue

        mode = str(trace.mode or "lines")
        name = str(trace.name or "")
        line = getattr(trace, "line", None)
        marker = getattr(trace, "marker", None)

        line_color = _plotly_color_to_mpl(getattr(line, "color", None))
        line_width = float(getattr(line, "width", 1.8) or 1.8)
        line_style = _dash_to_linestyle(getattr(line, "dash", None))

        marker_symbol = "o"
        marker_size = 6.0
        marker_color = _plotly_color_to_mpl(getattr(marker, "color", line_color))
        if marker is not None:
            m_symbol = str(getattr(marker, "symbol", "") or "")
            if "x" in m_symbol:
                marker_symbol = "x"
            marker_size = float(getattr(marker, "size", 6) or 6)

        label = name if name else None
        if "lines" in mode:
            ax.plot(
                x,
                y,
                linestyle=line_style,
                linewidth=line_width,
                color=line_color,
                label=label,
                alpha=float(getattr(trace, "opacity", 1.0) or 1.0),
            )
            label = None
        if "markers" in mode:
            ax.plot(
                x,
                y,
                linestyle="None",
                marker=marker_symbol,
                markersize=max(3.0, marker_size * 0.45),
                color=marker_color,
                label=label,
                alpha=float(getattr(trace, "opacity", 1.0) or 1.0),
            )

    # Draw vertical reference lines from layout shapes
    for shape in list(getattr(fig.layout, "shapes", []) or []):
        if getattr(shape, "type", "") != "line":
            continue
        x0 = getattr(shape, "x0", None)
        x1 = getattr(shape, "x1", None)
        if x0 is None or x1 is None or float(x0) != float(x1):
            continue
        line = getattr(shape, "line", None)
        ax.axvline(
            float(x0),
            color=_plotly_color_to_mpl(getattr(line, "color", "#999999")),
            linewidth=float(getattr(line, "width", 1.0) or 1.0),
            linestyle=_dash_to_linestyle(getattr(line, "dash", "dot")),
            alpha=0.8,
        )

    # Titles and labels
    title_text = ""
    if getattr(fig.layout, "title", None) is not None:
        title_text = str(getattr(fig.layout.title, "text", "") or "")
    if title_text:
        ax.set_title(title_text, fontsize=18, fontweight="bold", pad=12)

    xaxis_obj = getattr(fig.layout, "xaxis", None)
    yaxis_obj = getattr(fig.layout, "yaxis", None)
    x_title_obj = getattr(xaxis_obj, "title", None)
    y_title_obj = getattr(yaxis_obj, "title", None)
    x_title = str(getattr(x_title_obj, "text", "") or "Binding Energy (eV)")
    y_title = str(getattr(y_title_obj, "text", "") or "Intensity (a.u.)")
    ax.set_xlabel(x_title, fontsize=16, fontweight="bold")
    ax.set_ylabel(y_title, fontsize=16, fontweight="bold")

    ax.tick_params(axis="both", labelsize=13, direction="in", length=6, width=1.4, top=True, right=True)
    ax.minorticks_on()
    ax.grid(True, alpha=0.25, linewidth=0.8)

    # Axis direction and ranges
    xaxis = getattr(fig.layout, "xaxis", None)
    if xaxis is not None:
        x_range = getattr(xaxis, "range", None)
        if x_range and len(x_range) == 2:
            x0, x1 = float(x_range[0]), float(x_range[1])
            ax.set_xlim(min(x0, x1), max(x0, x1))
        if str(getattr(xaxis, "autorange", "")).lower() == "reversed":
            ax.invert_xaxis()

    yaxis = getattr(fig.layout, "yaxis", None)
    if yaxis is not None:
        y_range = getattr(yaxis, "range", None)
        if y_range and len(y_range) == 2:
            y0, y1 = float(y_range[0]), float(y_range[1])
            ax.set_ylim(min(y0, y1), max(y0, y1))

    # Text annotations
    for ann in list(getattr(fig.layout, "annotations", []) or []):
        x = getattr(ann, "x", None)
        y = getattr(ann, "y", None)
        text = str(getattr(ann, "text", "") or "")
        if x is None or y is None or not text:
            continue
        rotation = float(getattr(ann, "textangle", 0) or 0)
        ax.text(
            float(x),
            float(y),
            text,
            fontsize=11,
            rotation=rotation,
            ha="center",
            va="center",
            color="#3a3a3a",
            alpha=0.92,
        )

    if bool(getattr(fig.layout, "showlegend", True)):
        handles, labels = ax.get_legend_handles_labels()
        if labels:
            ax.legend(
                loc="upper center",
                bbox_to_anchor=(0.5, 1.02),
                ncol=min(6, max(1, len(labels))),
                fontsize=11,
                frameon=False,
            )

    out = io.BytesIO()
    mfig.savefig(out, format=fmt, dpi=dpi, bbox_inches="tight")
    plt.close(mfig)
    return out.getvalue()

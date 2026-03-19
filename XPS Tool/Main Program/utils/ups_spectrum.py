"""
UPS spectrum utilities for loading, work function analysis, and plotting.

UPS (Ultraviolet Photoelectron Spectroscopy) uses UV photon sources
(He I = 21.22 eV, He II = 40.81 eV) to probe the valence band and
measure the work function via the secondary electron cutoff.
"""

from __future__ import annotations

import io
from typing import Any, Dict, List, Sequence

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.signal import savgol_filter


# ── Photon source presets ─────────────────────────────────────────────────
PHOTON_SOURCES: Dict[str, float] = {
    "He I (21.22 eV)": 21.22,
    "He II (40.81 eV)": 40.81,
    "Ne I (16.67 eV)": 16.67,
    "Ne II (26.91 eV)": 26.91,
    "Custom": 0.0,
}

# ── Column detection candidates ──────────────────────────────────────────
ENERGY_COLUMN_CANDIDATES = (
    "e", "energy", "binding_energy", "binding energy", "be",
    "kinetic_energy", "kinetic energy", "ke", "x",
)
INTENSITY_COLUMN_CANDIDATES = (
    "y", "intensity", "counts", "count", "cps", "signal",
)


# ── Helpers ──────────────────────────────────────────────────────────────

def _find_best_column(columns: Sequence[str], candidates: Sequence[str]) -> str | None:
    normalized = {c.strip().lower(): c for c in columns}
    for cand in candidates:
        if cand in normalized:
            return normalized[cand]
    for norm, orig in normalized.items():
        for cand in candidates:
            if cand in norm:
                return orig
    return None


def _flatten_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if isinstance(out.columns, pd.MultiIndex):
        out.columns = [
            "_".join(str(x) for x in col if str(x) not in ("", "None"))
            if isinstance(col, tuple) else str(col)
            for col in out.columns.to_flat_index()
        ]
    return out


def _safe_savgol(y: np.ndarray, window_length: int) -> np.ndarray:
    if y.size < 7:
        return y.copy()
    max_win = y.size - 1 if y.size % 2 == 0 else y.size
    win = int(min(window_length, max_win))
    if win % 2 == 0:
        win -= 1
    if win < 5:
        return y.copy()
    polyorder = min(3, win - 2)
    return savgol_filter(y, window_length=win, polyorder=polyorder)


# ── Data loading ─────────────────────────────────────────────────────────

def load_ups_spectrum_data(file_content: bytes, filename: str = "") -> Dict[str, Any]:
    """Load and validate a UPS spectrum from CSV bytes."""
    df = pd.read_csv(io.BytesIO(file_content))
    if df.empty:
        raise ValueError("CSV file is empty.")
    return _spectrum_from_dataframe(df, filename=filename)


def _spectrum_from_dataframe(df: pd.DataFrame, filename: str = "") -> Dict[str, Any]:
    out = _flatten_columns(df)
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
        raise ValueError("Cannot detect energy/intensity columns in CSV.")

    energy = pd.to_numeric(out[energy_col], errors="coerce").to_numpy(dtype=float)
    intensity = pd.to_numeric(out[intensity_col], errors="coerce").to_numpy(dtype=float)

    valid = np.isfinite(energy) & np.isfinite(intensity)
    energy = energy[valid]
    intensity = intensity[valid]

    if energy.size < 10:
        raise ValueError("Not enough valid points for UPS analysis (need >= 10).")

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


def load_uploaded_ups_file(file_content: bytes, filename: str) -> Dict[str, Any]:
    """Load a UPS CSV file into the standard payload format."""
    lower = (filename or "").lower()
    if not lower.endswith(".csv"):
        raise ValueError("UPS panel currently supports .csv files only.")

    spec = load_ups_spectrum_data(file_content, filename=filename)
    spec["trace_name"] = filename.rsplit(".", 1)[0] if "." in filename else filename
    return {
        "file_type": "csv",
        "spectra": [spec],
        "logs": "Loaded UPS CSV spectrum.",
        "method": "csv",
    }


# ── Edge detection via linear extrapolation ──────────────────────────────

def _linear_edge_fit(
    energy: np.ndarray,
    intensity: np.ndarray,
    region_min: float,
    region_max: float,
) -> Dict[str, Any]:
    """
    Fit the steepest part of an edge region with a line and extrapolate
    the x-intercept (intensity = 0).

    Returns slope, intercept, x-intercept, and the fit segment.
    """
    mask = (energy >= region_min) & (energy <= region_max)
    e_seg = energy[mask]
    i_seg = intensity[mask]

    if e_seg.size < 4:
        raise ValueError(
            f"Too few points in edge region [{region_min:.2f}, {region_max:.2f}]."
        )

    # Smooth, then find the steepest gradient window
    i_smooth = _safe_savgol(i_seg, min(11, max(5, e_seg.size // 3)))
    grad = np.gradient(i_smooth, e_seg)

    # For cutoff: steepest rising edge (large positive gradient)
    # We pick the window with maximum |gradient|
    abs_grad = np.abs(grad)
    # Use a rolling-window of ~30 % of the segment to find the linear region
    win = max(4, int(0.3 * e_seg.size))
    best_start = 0
    best_sum = 0.0
    for start in range(0, max(1, e_seg.size - win + 1)):
        s = float(np.sum(abs_grad[start:start + win]))
        if s > best_sum:
            best_sum = s
            best_start = start

    fit_slice = slice(best_start, best_start + win)
    e_fit = e_seg[fit_slice]
    i_fit = i_seg[fit_slice]

    if e_fit.size < 3:
        e_fit = e_seg
        i_fit = i_seg

    slope, intercept = np.polyfit(e_fit, i_fit, 1)

    # x-intercept where extrapolated line meets y=0
    if abs(slope) < 1e-15:
        x_intercept = float(e_fit[0])
    else:
        x_intercept = float(-intercept / slope)

    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "x_intercept": x_intercept,
        "fit_energy": e_fit,
        "fit_intensity": i_fit,
    }


def detect_secondary_cutoff(
    energy: np.ndarray,
    intensity: np.ndarray,
    cutoff_region: tuple[float, float] | None = None,
    smoothing_window: int = 11,
) -> Dict[str, Any]:
    """
    Detect the secondary electron cutoff (SEC) edge.

    The SEC is typically at the low-kinetic-energy (high-binding-energy) side
    of the UPS spectrum. The work function onset is found by linear
    extrapolation of the leading edge to zero intensity.
    """
    if cutoff_region is not None:
        rmin, rmax = float(min(cutoff_region)), float(max(cutoff_region))
    else:
        # Auto: use the lowest 30 % of the energy range
        e_span = float(np.ptp(energy))
        rmin = float(np.min(energy))
        rmax = rmin + 0.30 * e_span

    fit = _linear_edge_fit(energy, intensity, rmin, rmax)
    return {
        "cutoff_eV": fit["x_intercept"],
        "fit": fit,
        "region": (rmin, rmax),
    }


def detect_valence_band_edge(
    energy: np.ndarray,
    intensity: np.ndarray,
    vb_region: tuple[float, float] | None = None,
    smoothing_window: int = 11,
) -> Dict[str, Any]:
    """
    Detect the valence band maximum (VBM) / Fermi edge by linear
    extrapolation of the high-energy trailing edge to zero intensity.
    """
    if vb_region is not None:
        rmin, rmax = float(min(vb_region)), float(max(vb_region))
    else:
        # Auto: use the highest 30 % of the energy range
        e_span = float(np.ptp(energy))
        rmax = float(np.max(energy))
        rmin = rmax - 0.30 * e_span

    fit = _linear_edge_fit(energy, intensity, rmin, rmax)
    return {
        "vbm_eV": fit["x_intercept"],
        "fit": fit,
        "region": (rmin, rmax),
    }


def calculate_work_function(
    photon_energy: float,
    cutoff_eV: float,
    vbm_eV: float,
    energy_mode: str = "binding",
) -> Dict[str, Any]:
    """
    Calculate the work function from UPS spectrum width.

    For binding energy scale:
        Spectrum width = E_cutoff - E_Fermi (where Fermi ~ 0 eV)
        Work function  = hν - spectrum_width

    For kinetic energy scale:
        Work function = hν - (KE_max - KE_min)
        where KE_min is SEC and KE_max is Fermi edge.

    Also computes the ionization potential (IP):
        IP = hν - (E_cutoff - E_VBM)   [binding energy]
    """
    if energy_mode == "kinetic":
        # KE: cutoff is low-KE edge, vbm is high-KE Fermi edge
        spectrum_width = abs(vbm_eV - cutoff_eV)
        work_function = photon_energy - spectrum_width
        ionization_potential = photon_energy - (vbm_eV - cutoff_eV)
    else:
        # Binding energy: cutoff at high BE, Fermi at 0
        spectrum_width = abs(cutoff_eV - vbm_eV)
        work_function = photon_energy - spectrum_width
        ionization_potential = photon_energy - spectrum_width + abs(vbm_eV)

    return {
        "work_function_eV": float(work_function),
        "spectrum_width_eV": float(spectrum_width),
        "ionization_potential_eV": float(ionization_potential),
        "photon_energy_eV": float(photon_energy),
        "cutoff_eV": float(cutoff_eV),
        "vbm_eV": float(vbm_eV),
        "energy_mode": energy_mode,
    }


# ── Visualization ────────────────────────────────────────────────────────

def create_ups_spectrum_figure(
    energy: np.ndarray,
    intensity: np.ndarray,
    cutoff_result: Dict[str, Any] | None = None,
    vbm_result: Dict[str, Any] | None = None,
    wf_result: Dict[str, Any] | None = None,
    title: str = "UPS Spectrum",
    show_cutoff_fit: bool = True,
    show_vbm_fit: bool = True,
    show_annotations: bool = True,
    energy_mode: str = "binding",
    base_font_size: int = 16,
    title_font_size: int = 24,
    figure_height: int = 700,
    focus_range: Sequence[float] | None = None,
    line_width_scale: float = 0.85,
) -> go.Figure:
    """Create an annotated UPS spectrum figure with edge fits."""
    fig = go.Figure()
    ws = float(np.clip(line_width_scale, 0.5, 2.0))

    # Raw spectrum
    fig.add_trace(go.Scatter(
        x=energy, y=intensity,
        mode="lines",
        line=dict(color="#1f2a44", width=1.6 * ws),
        name="UPS spectrum",
        hovertemplate="E: %{x:.3f} eV<br>I: %{y:.3f}<extra></extra>",
    ))

    y_min = float(np.min(intensity))
    y_max = float(np.max(intensity))
    y_span = max(y_max - y_min, 1e-9)

    # SEC cutoff linear fit
    if show_cutoff_fit and cutoff_result is not None:
        fit = cutoff_result["fit"]
        e_fit = fit["fit_energy"]
        slope, intercept = fit["slope"], fit["intercept"]
        # Extend the line to y=0
        x_int = fit["x_intercept"]
        e_line = np.array([x_int, float(e_fit[-1]) if slope > 0 else float(e_fit[0])])
        i_line = slope * e_line + intercept

        fig.add_trace(go.Scatter(
            x=e_line, y=i_line,
            mode="lines",
            line=dict(color="#d62728", width=2.0 * ws, dash="dash"),
            name="SEC cutoff fit",
        ))
        # Mark cutoff point
        fig.add_trace(go.Scatter(
            x=[x_int], y=[0],
            mode="markers",
            marker=dict(color="#d62728", size=10 * ws, symbol="diamond"),
            name=f"Cutoff = {x_int:.2f} eV",
        ))
        if show_annotations:
            fig.add_vline(
                x=x_int,
                line_width=max(0.8, 1.0 * ws),
                line_dash="dot",
                line_color="rgba(214, 39, 40, 0.5)",
            )
            fig.add_annotation(
                x=x_int, y=y_min + 0.10 * y_span,
                text=f"SEC = {x_int:.2f} eV",
                showarrow=True, arrowhead=2,
                font=dict(size=max(11, int(base_font_size * 0.85)), color="#d62728"),
                bgcolor="rgba(255,255,255,0.8)",
            )

    # VBM / Fermi edge linear fit
    if show_vbm_fit and vbm_result is not None:
        fit = vbm_result["fit"]
        e_fit = fit["fit_energy"]
        slope, intercept = fit["slope"], fit["intercept"]
        x_int = fit["x_intercept"]
        e_line = np.array([float(e_fit[0]) if slope < 0 else float(e_fit[-1]), x_int])
        i_line = slope * e_line + intercept

        fig.add_trace(go.Scatter(
            x=e_line, y=i_line,
            mode="lines",
            line=dict(color="#2ca02c", width=2.0 * ws, dash="dash"),
            name="VBM/Fermi fit",
        ))
        fig.add_trace(go.Scatter(
            x=[x_int], y=[0],
            mode="markers",
            marker=dict(color="#2ca02c", size=10 * ws, symbol="diamond"),
            name=f"VBM = {x_int:.2f} eV",
        ))
        if show_annotations:
            fig.add_vline(
                x=x_int,
                line_width=max(0.8, 1.0 * ws),
                line_dash="dot",
                line_color="rgba(44, 160, 44, 0.5)",
            )
            fig.add_annotation(
                x=x_int, y=y_min + 0.10 * y_span,
                text=f"VBM = {x_int:.2f} eV",
                showarrow=True, arrowhead=2,
                font=dict(size=max(11, int(base_font_size * 0.85)), color="#2ca02c"),
                bgcolor="rgba(255,255,255,0.8)",
            )

    # Work function annotation box
    if show_annotations and wf_result is not None:
        wf = wf_result["work_function_eV"]
        sw = wf_result["spectrum_width_eV"]
        ip = wf_result["ionization_potential_eV"]
        hv = wf_result["photon_energy_eV"]
        annotation_text = (
            f"<b>Work Function</b><br>"
            f"Φ = {wf:.3f} eV<br>"
            f"hν = {hv:.2f} eV<br>"
            f"Width = {sw:.3f} eV<br>"
            f"IP = {ip:.3f} eV"
        )
        fig.add_annotation(
            x=0.98, y=0.95, xref="paper", yref="paper",
            text=annotation_text,
            showarrow=False,
            font=dict(size=max(12, base_font_size - 1), family="monospace"),
            bgcolor="rgba(255,255,255,0.92)",
            bordercolor="rgba(100,100,100,0.6)",
            borderwidth=1,
            borderpad=8,
            align="left",
            xanchor="right", yanchor="top",
        )

    x_label = "Kinetic Energy (eV)" if energy_mode == "kinetic" else "Binding Energy (eV)"
    fig.update_layout(
        title=dict(text=title, font=dict(size=title_font_size)),
        template="plotly_white",
        hovermode="closest",
        xaxis_title=x_label,
        yaxis_title="Intensity (a.u.)",
        font=dict(size=base_font_size),
        height=int(figure_height),
        margin=dict(l=90, r=30, t=100, b=90),
        legend=dict(
            orientation="h", yanchor="bottom", y=1.02,
            xanchor="left", x=0.0,
            font=dict(size=max(10, base_font_size - 1)),
        ),
        showlegend=True,
    )

    reversed_x = (energy_mode == "binding")
    fig.update_xaxes(
        autorange="reversed" if reversed_x else True,
        showgrid=True, gridcolor="LightGray",
        title_font=dict(size=max(12, base_font_size + 4)),
        tickfont=dict(size=max(10, base_font_size)),
    )
    fig.update_yaxes(
        showgrid=True, gridcolor="LightGray",
        title_font=dict(size=max(12, base_font_size + 4)),
        tickfont=dict(size=max(10, base_font_size)),
    )

    if focus_range is not None and len(focus_range) == 2:
        x0, x1 = float(min(focus_range)), float(max(focus_range))
        if x1 > x0:
            if reversed_x:
                fig.update_xaxes(range=[x1, x0])
            else:
                fig.update_xaxes(range=[x0, x1])

    return fig


def export_ups_figure_bytes(
    fig: go.Figure,
    export_format: str = "png",
    width_px: int = 2000,
    height_px: int = 1200,
    dpi: int = 300,
) -> bytes:
    """Export a UPS Plotly figure to image bytes (png/svg/pdf)."""
    fmt = str(export_format).strip().lower().replace(".", "")
    if fmt not in {"png", "svg", "pdf"}:
        raise ValueError("Unsupported format. Use png/svg/pdf.")

    width_px = int(max(800, width_px))
    height_px = int(max(500, height_px))
    kwargs: Dict[str, Any] = {"format": fmt, "width": width_px, "height": height_px}

    if fmt == "png":
        scale = max(1.0, float(dpi) / 120.0)
        kwargs["scale"] = min(scale, 5.0)

    return fig.to_image(**kwargs)

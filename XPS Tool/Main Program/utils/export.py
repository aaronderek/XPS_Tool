"""
Export utilities for figures and data.
"""

import json
import pandas as pd
import numpy as np
from datetime import datetime
import io
import base64


def export_to_csv(data_dict, filename):
    """
    Export curve data to CSV file.

    Parameters
    ----------
    data_dict : dict
        Dictionary with curve data
        Expected format: {
            'curves': [
                {
                    'name': 'curve_name',
                    'x': [x1, x2, ...],
                    'y': [y1, y2, ...],
                    'params': {...}
                },
                ...
            ]
        }
    filename : str
        Output filename

    Returns
    -------
    str
        Path to saved file
    """
    rows = []

    for curve in data_dict.get('curves', []):
        name = curve['name']
        x_data = curve['x']
        y_data = curve['y']
        params = curve.get('params', {})

        for x_val, y_val in zip(x_data, y_data):
            row = {
                'curve_name': name,
                'x': x_val,
                'y': y_val,
            }
            # Add parameters
            row.update(params)
            rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(filename, index=False)

    return filename


def export_parameters_to_json(params_dict, curves_data, filename):
    """
    Export parameters and curve data to JSON for reproducibility.

    Parameters
    ----------
    params_dict : dict
        Dictionary of parameters
    curves_data : list
        List of curve data
    filename : str
        Output filename

    Returns
    -------
    str
        Path to saved file
    """
    export_data = {
        'timestamp': datetime.now().isoformat(),
        'parameters': params_dict,
        'curves': []
    }

    # Convert numpy arrays to lists for JSON serialization
    for curve in curves_data:
        curve_data = {
            'name': curve['name'],
            'type': curve.get('type', 'unknown'),
            'data': []
        }

        x = curve['x']
        y = curve['y']

        # Convert to lists if numpy arrays
        if isinstance(x, np.ndarray):
            x = x.tolist()
        if isinstance(y, np.ndarray):
            y = y.tolist()

        for x_val, y_val in zip(x, y):
            curve_data['data'].append([x_val, y_val])

        export_data['curves'].append(curve_data)

    with open(filename, 'w') as f:
        json.dump(export_data, f, indent=2)

    return filename


def load_parameters_from_json(filename):
    """
    Load parameters from JSON file.

    Parameters
    ----------
    filename : str
        Input filename

    Returns
    -------
    dict
        Dictionary with parameters and curve data
    """
    with open(filename, 'r') as f:
        data = json.load(f)

    return data


def fig_to_base64(fig, format='png'):
    """
    Convert matplotlib or plotly figure to base64 string.

    Parameters
    ----------
    fig : matplotlib.figure.Figure or plotly.graph_objects.Figure
        Figure to convert
    format : str
        Output format ('png', 'svg', 'jpeg')

    Returns
    -------
    str
        Base64 encoded string
    """
    buf = io.BytesIO()

    # Check if it's a plotly figure
    if hasattr(fig, 'write_image'):
        fig.write_image(buf, format=format)
    else:
        # Matplotlib figure
        fig.savefig(buf, format=format, dpi=300, bbox_inches='tight')

    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode()

    return img_base64


def create_download_link(data, filename, file_format='csv'):
    """
    Create a download link for Streamlit.

    Parameters
    ----------
    data : bytes or str
        Data to download
    filename : str
        Suggested filename
    file_format : str
        File format ('csv', 'json', 'png', 'svg')

    Returns
    -------
    str
        HTML download link
    """
    if isinstance(data, str):
        data = data.encode()

    b64 = base64.b64encode(data).decode()

    mime_types = {
        'csv': 'text/csv',
        'json': 'application/json',
        'png': 'image/png',
        'svg': 'image/svg+xml'
    }

    mime_type = mime_types.get(file_format, 'application/octet-stream')

    href = f'<a href="data:{mime_type};base64,{b64}" download="{filename}">Download {filename}</a>'

    return href

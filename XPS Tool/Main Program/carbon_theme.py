"""Helpers for injecting Carbon-inspired CSS in Streamlit apps."""

from pathlib import Path
import streamlit as st


def inject_carbon_css(css_path: str = "carbon_streamlit.css") -> None:
    css_file = Path(css_path)
    if not css_file.exists():
        raise FileNotFoundError(f"CSS file not found: {css_file}")
    st.markdown(f"<style>{css_file.read_text(encoding='utf-8')}</style>", unsafe_allow_html=True)

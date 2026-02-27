#!/usr/bin/env bash
# Codex-friendly run script: install deps (if needed) and start Streamlit.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

resolve_system_python() {
  if command -v python3 >/dev/null 2>&1; then
    command -v python3
    return
  fi
  if command -v python >/dev/null 2>&1; then
    command -v python
    return
  fi
  echo ""
}

echo "========================================="
echo "2DEG Visualization Tool"
echo "========================================="
echo "Project dir: ${SCRIPT_DIR}"
echo ""

if [ ! -d ".venv" ] || { [ ! -x ".venv/bin/python" ] && [ ! -x ".venv/bin/python3" ]; }; then
  PY_BOOTSTRAP="$(resolve_system_python)"
  if [ -z "${PY_BOOTSTRAP}" ]; then
    echo "[error] Python is not available in PATH. Please install Python 3."
    exit 1
  fi
  echo "[setup] Creating or repairing virtual environment (.venv)"
  "${PY_BOOTSTRAP}" -m venv .venv
fi

# shellcheck source=/dev/null
source .venv/bin/activate

if [ -x ".venv/bin/python" ]; then
  PY_BIN=".venv/bin/python"
elif [ -x ".venv/bin/python3" ]; then
  PY_BIN=".venv/bin/python3"
else
  echo "[error] .venv exists but has no Python binary (.venv/bin/python or .venv/bin/python3)."
  exit 1
fi

echo "[setup] Installing dependencies from requirements.txt"
"${PY_BIN}" -m ensurepip --upgrade >/dev/null 2>&1 || true
"${PY_BIN}" -m pip install --upgrade pip
"${PY_BIN}" -m pip install -r requirements.txt

echo ""
echo "[run] Starting Streamlit application..."
APP_URL="http://127.0.0.1:8501"
echo "Open your browser at: ${APP_URL}"
echo ""

open_browser() {
  local url="$1"
  if command -v open >/dev/null 2>&1; then
    open "$url" >/dev/null 2>&1 && return 0
  fi
  if command -v xdg-open >/dev/null 2>&1; then
    xdg-open "$url" >/dev/null 2>&1 && return 0
  fi
  if command -v cmd.exe >/dev/null 2>&1; then
    cmd.exe /c start "" "$url" >/dev/null 2>&1 && return 0
  fi
  return 1
}

if [ "${AUTO_OPEN_BROWSER:-1}" = "1" ]; then
  (
    sleep 2
    if open_browser "${APP_URL}"; then
      echo "[run] Browser opened: ${APP_URL}"
    else
      echo "[run] Could not auto-open browser. Open manually: ${APP_URL}"
    fi
  ) &
fi

exec "${PY_BIN}" -m streamlit run app.py --server.address 127.0.0.1 --server.port 8501

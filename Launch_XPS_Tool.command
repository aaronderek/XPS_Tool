#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
APP_DIR="$ROOT_DIR/XPS Tool/Main Program"
STATE_DIR="$ROOT_DIR/.runtime"
PID_FILE="$STATE_DIR/streamlit.pid"
LOG_FILE="$STATE_DIR/streamlit.log"
APP_URL="http://127.0.0.1:8501"

mkdir -p "$STATE_DIR"

if [ ! -d "$APP_DIR" ]; then
  echo "Error: App directory not found: $APP_DIR"
  read -r -p "Press Enter to close..."
  exit 1
fi

cd "$APP_DIR"

echo "========================================="
echo "XPS Tool"
echo "========================================="

if [ -f "$PID_FILE" ]; then
  old_pid="$(cat "$PID_FILE" 2>/dev/null || true)"
  if [ -n "${old_pid:-}" ] && kill -0 "$old_pid" 2>/dev/null; then
    echo "Streamlit is already running (PID: $old_pid)"
    echo "Open: $APP_URL"
    open "$APP_URL" >/dev/null 2>&1 || true
    exit 0
  fi
fi

PY_CMD=""
if [ -x ".venv/bin/python" ]; then
  PY_CMD=".venv/bin/python"
elif command -v python3 >/dev/null 2>&1; then
  PY_CMD="python3"
fi

if [ -z "$PY_CMD" ]; then
  echo "Error: No suitable Python runtime found."
  read -r -p "Press Enter to close..."
  exit 1
fi

echo "Starting Streamlit in background..."
nohup "$PY_CMD" -m streamlit run app.py --server.address 127.0.0.1 --server.port 8501 >"$LOG_FILE" 2>&1 &
new_pid=$!
echo "$new_pid" >"$PID_FILE"

sleep 2
if ! kill -0 "$new_pid" 2>/dev/null; then
  echo "Failed to start Streamlit. Last logs:"
  tail -n 60 "$LOG_FILE" || true
  read -r -p "Press Enter to close..."
  exit 1
fi

echo "Started (PID: $new_pid)"
echo "Open: $APP_URL"
echo "Log:  $LOG_FILE"
open "$APP_URL" >/dev/null 2>&1 || true

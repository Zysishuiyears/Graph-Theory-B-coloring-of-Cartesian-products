#!/usr/bin/env bash
set -euo pipefail

if [[ ! -f "README.md" || ! -d "src/active" ]]; then
    echo "Run this script from the repository root." >&2
    exit 1
fi

PYTHON_BIN="${PYTHON_BIN:-}"
if [[ -z "$PYTHON_BIN" ]]; then
    if command -v python >/dev/null 2>&1; then
        PYTHON_BIN="python"
    elif command -v python3 >/dev/null 2>&1; then
        PYTHON_BIN="python3"
    else
        echo "Neither python nor python3 was found in PATH." >&2
        exit 1
    fi
fi

CMSAT_BIN="${CMSAT_PATH:-cryptominisat5}"

echo "repo_root=$(pwd)"
echo "branch=$(git branch --show-current)"
echo "python=$PYTHON_BIN"
"$PYTHON_BIN" --version
echo "git=$(command -v git || true)"
git --version

echo "tmux=$(command -v tmux || true)"
if ! command -v tmux >/dev/null 2>&1; then
    echo "warning: tmux not found; long-running jobs should not be started yet." >&2
fi

if command -v "$CMSAT_BIN" >/dev/null 2>&1; then
    echo "cmsat=$CMSAT_BIN"
    version_out="$($CMSAT_BIN --version 2>&1 || true)"
    printf '%s\n' "$version_out" | sed -n '1,2p'
else
    echo "cmsat=$CMSAT_BIN"
    echo "CryptoMiniSat not found in PATH or CMSAT_PATH." >&2
    exit 1
fi

"$PYTHON_BIN" - <<'PY'
from src.active.general_product_decide import decide_existence_B
from src.active.general_product_search import search_best
print("python_imports=ok")
print("decide_symbol=", decide_existence_B.__name__)
print("search_symbol=", search_best.__name__)
PY

echo "environment_check=ok"

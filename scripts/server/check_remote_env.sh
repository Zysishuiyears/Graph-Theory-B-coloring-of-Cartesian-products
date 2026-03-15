#!/usr/bin/env bash
set -euo pipefail

if [[ ! -f "README.md" || ! -d "src/active" ]]; then
    echo "Run this script from the repository root." >&2
    exit 1
fi

CMSAT_BIN="${CMSAT_PATH:-cryptominisat5}"

echo "repo_root=$(pwd)"
echo "python=$(command -v python || true)"
python --version
echo "git=$(command -v git || true)"
git --version
echo "tmux=$(command -v tmux || true)"

if command -v "$CMSAT_BIN" >/dev/null 2>&1; then
    echo "cmsat=$CMSAT_BIN"
    "$CMSAT_BIN" --version | head -n 2
else
    echo "cmsat=$CMSAT_BIN"
    echo "CryptoMiniSat not found in PATH or CMSAT_PATH." >&2
    exit 1
fi

python - <<'PY'
from src.active.general_product_decide import decide_existence_B
from src.active.general_product_search import search_best
print("python_imports=ok")
print("decide_symbol=", decide_existence_B.__name__)
print("search_symbol=", search_best.__name__)
PY

echo "environment_check=ok"

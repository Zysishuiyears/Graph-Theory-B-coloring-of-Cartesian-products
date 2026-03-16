#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage: bash scripts/server/run_decide.sh '[4,6,6,10]' '[]' 8 [dump_one_solution]" >&2
    exit 1
fi

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

CYCLES_EXPR="$1"
PATHS_EXPR="$2"
K_VALUE="$3"
DUMP_ONE_SOLUTION="${4:-1}"
CPU_COUNT="$(nproc)"
DEFAULT_WORKERS="$(( CPU_COUNT / 3 ))"
if (( DEFAULT_WORKERS < 1 )); then
    DEFAULT_WORKERS=1
elif (( DEFAULT_WORKERS > 4 )); then
    DEFAULT_WORKERS=4
fi

export CMSAT_PATH="${CMSAT_PATH:-cryptominisat5}"
export BCOLOR_RESULTS_ROOT="${BCOLOR_RESULTS_ROOT:-results/runs}"
export CMSAT_SOLVER_MODE="${CMSAT_SOLVER_MODE:-portfolio}"
export CMSAT_THREADS="${CMSAT_THREADS:-$CPU_COUNT}"
export CMSAT_PORTFOLIO_WORKERS="${CMSAT_PORTFOLIO_WORKERS:-$DEFAULT_WORKERS}"
export BCOLOR_SYMMETRY_BREAKING="${BCOLOR_SYMMETRY_BREAKING:-1}"

if ! command -v "$CMSAT_PATH" >/dev/null 2>&1; then
    echo "Solver not found: $CMSAT_PATH" >&2
    exit 1
fi

echo "PYTHON_BIN=$PYTHON_BIN"
echo "CMSAT_PATH=$CMSAT_PATH"
echo "BCOLOR_RESULTS_ROOT=$BCOLOR_RESULTS_ROOT"
echo "CMSAT_SOLVER_MODE=$CMSAT_SOLVER_MODE"
echo "CMSAT_THREADS=$CMSAT_THREADS"
echo "CMSAT_PORTFOLIO_WORKERS=$CMSAT_PORTFOLIO_WORKERS"
echo "BCOLOR_SYMMETRY_BREAKING=$BCOLOR_SYMMETRY_BREAKING"
echo "cycles=$CYCLES_EXPR paths=$PATHS_EXPR k=$K_VALUE dump=$DUMP_ONE_SOLUTION"

"$PYTHON_BIN" - "$CYCLES_EXPR" "$PATHS_EXPR" "$K_VALUE" "$DUMP_ONE_SOLUTION" <<'PY'
import ast
import sys

from src.active.general_product_decide import decide_existence_B

cycles = ast.literal_eval(sys.argv[1])
paths = ast.literal_eval(sys.argv[2])
k = int(sys.argv[3])
dump_one_solution = bool(int(sys.argv[4]))

result = decide_existence_B(
    cycles=cycles,
    paths=paths,
    k=k,
    dump_one_solution=dump_one_solution,
)
print("decide_result=", result)
PY

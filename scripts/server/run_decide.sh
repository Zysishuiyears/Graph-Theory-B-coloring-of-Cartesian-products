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

CYCLES_EXPR="$1"
PATHS_EXPR="$2"
K_VALUE="$3"
DUMP_ONE_SOLUTION="${4:-1}"

export CMSAT_PATH="${CMSAT_PATH:-cryptominisat5}"
export BCOLOR_RESULTS_ROOT="${BCOLOR_RESULTS_ROOT:-results/runs}"
export CMSAT_SOLVER_MODE="${CMSAT_SOLVER_MODE:-portfolio}"
export CMSAT_THREADS="${CMSAT_THREADS:-$(nproc)}"
export CMSAT_PORTFOLIO_WORKERS="${CMSAT_PORTFOLIO_WORKERS:-4}"
export BCOLOR_SYMMETRY_BREAKING="${BCOLOR_SYMMETRY_BREAKING:-1}"

echo "CMSAT_PATH=$CMSAT_PATH"
echo "BCOLOR_RESULTS_ROOT=$BCOLOR_RESULTS_ROOT"
echo "CMSAT_SOLVER_MODE=$CMSAT_SOLVER_MODE"
echo "CMSAT_THREADS=$CMSAT_THREADS"
echo "CMSAT_PORTFOLIO_WORKERS=$CMSAT_PORTFOLIO_WORKERS"
echo "BCOLOR_SYMMETRY_BREAKING=$BCOLOR_SYMMETRY_BREAKING"
echo "cycles=$CYCLES_EXPR paths=$PATHS_EXPR k=$K_VALUE dump=$DUMP_ONE_SOLUTION"

python - "$CYCLES_EXPR" "$PATHS_EXPR" "$K_VALUE" "$DUMP_ONE_SOLUTION" <<'PY'
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

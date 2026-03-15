#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 5 ]]; then
    echo "Usage: bash scripts/server/run_decide_tmux.sh '[4,6,6,10]' '[]' 8 [session_name] [dump_one_solution]" >&2
    exit 1
fi

if [[ ! -f "README.md" || ! -d "src/active" ]]; then
    echo "Run this script from the repository root." >&2
    exit 1
fi

if ! command -v tmux >/dev/null 2>&1; then
    echo "tmux is required but not installed." >&2
    exit 1
fi

CYCLES_EXPR="$1"
PATHS_EXPR="$2"
K_VALUE="$3"
SESSION_NAME="${4:-decide-$(date +%Y%m%d-%H%M%S)}"
DUMP_ONE_SOLUTION="${5:-1}"

mkdir -p results/runs/_remote_logs
LOG_FILE="results/runs/_remote_logs/${SESSION_NAME}.log"
REPO_ROOT="$(pwd)"

CMD=$(printf "cd %q && bash scripts/server/run_decide.sh %q %q %q %q 2>&1 | tee -a %q" \
    "$REPO_ROOT" "$CYCLES_EXPR" "$PATHS_EXPR" "$K_VALUE" "$DUMP_ONE_SOLUTION" "$LOG_FILE")

tmux new-session -d -s "$SESSION_NAME" "$CMD"

echo "tmux_session=$SESSION_NAME"
echo "log_file=$LOG_FILE"
echo "attach_cmd=tmux attach -t $SESSION_NAME"

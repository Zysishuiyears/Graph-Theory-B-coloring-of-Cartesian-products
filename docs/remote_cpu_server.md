# Remote CPU Server Workflow

This project can be run efficiently on a remote Linux CPU server through VS Code Remote - SSH.

The active runtime targets are:

- `src/active/general_product_decide.py`
- `src/active/general_product_search.py`

For the current stage, the recommended remote workload is the decision core only.

## 1. Recommended Setup

- Server OS: Ubuntu 22.04
- Python: 3.10 in a dedicated conda environment
- Solver: `cryptominisat5`
- Long-running job manager: `tmux`
- Repo path: use an ASCII-only path such as `~/projects/bcoloring`

## 2. VS Code Remote - SSH

Install these VS Code extensions locally:

- `Remote - SSH`
- `Python`

Add the server to your local SSH config:

```sshconfig
Host bcolor-server
    HostName <server-ip-or-domain>
    User <your-user>
    Port 22
```

Then connect through:

- `Ctrl+Shift+P`
- `Remote-SSH: Connect to Host`
- choose `bcolor-server`

After the remote window opens, open the repository folder on the server rather than editing local files through a mapped path.

## 3. Clone and Checkout

Run these commands on the server terminal:

```bash
mkdir -p ~/projects
cd ~/projects
git clone https://github.com/Zysishuiyears/Graph-Theory-B-coloring-of-Cartesian-products.git bcoloring
cd bcoloring
git checkout codex/split-decide-search
```

## 4. Python Environment

If conda is available:

```bash
conda create -n bcolor python=3.10 -y
conda activate bcolor
```

The current active code only needs Python standard library plus the external `cryptominisat5` executable.

## 5. Solver Check

First check whether the solver is already available:

```bash
which cryptominisat5
cryptominisat5 --version
```

If it is not installed system-wide, set `CMSAT_PATH` to the full binary path after you install or place the executable on the server.

The upstream references are:

- VS Code Remote - SSH: https://code.visualstudio.com/docs/remote/ssh
- CryptoMiniSat repo: https://github.com/msoos/cryptominisat
- CryptoMiniSat releases: https://github.com/msoos/cryptominisat/releases

## 6. Runtime Environment Variables

Recommended defaults on the server:

```bash
export CMSAT_PATH="${CMSAT_PATH:-cryptominisat5}"
export BCOLOR_RESULTS_ROOT="${BCOLOR_RESULTS_ROOT:-results/runs}"
export CMSAT_SOLVER_MODE="${CMSAT_SOLVER_MODE:-portfolio}"
export CMSAT_THREADS="${CMSAT_THREADS:-$(nproc)}"
export CMSAT_PORTFOLIO_WORKERS="${CMSAT_PORTFOLIO_WORKERS:-4}"
export BCOLOR_SYMMETRY_BREAKING="${BCOLOR_SYMMETRY_BREAKING:-1}"
```

## 7. Helper Scripts

This repository provides server-side helper scripts:

- `scripts/server/check_remote_env.sh`
- `scripts/server/run_decide.sh`
- `scripts/server/run_decide_tmux.sh`

Run the environment check from the repo root:

```bash
bash scripts/server/check_remote_env.sh
```

Run a foreground decision task:

```bash
bash scripts/server/run_decide.sh '[4,6,6,10]' '[]' 8
```

Run the same task inside `tmux`:

```bash
bash scripts/server/run_decide_tmux.sh '[4,6,6,10]' '[]' 8 c46610-k8
```

## 8. Why tmux

Do not bind long SAT jobs to the VS Code UI session alone.

Use `tmux` so the task survives:

- VS Code reload
- local network interruption
- SSH disconnect

Useful commands:

```bash
tmux ls
tmux attach -t c46610-k8
tmux kill-session -t c46610-k8
```

## 9. Results

Default outputs go to `results/runs/` on the server.

Recommended workflow:

- keep all raw runs on the server
- manually inspect `one_solution.txt` and `solver_summary.txt`
- only copy important outputs back to local machine
- only curate selected results into `results/snapshots/`

You can download files through:

- VS Code remote file explorer
- `scp`
- `rsync`

Example:

```bash
scp <user>@<host>:~/projects/bcoloring/results/runs/<run-dir>/solver_summary.txt .
```

## 10. Minimal Smoke Test

After setup, run:

```bash
python -c "from src.active.general_product_decide import decide_existence_B; decide_existence_B(cycles=[4], paths=[], k=4, dump_one_solution=False)"
```

Expected outcome:

- the import succeeds
- CryptoMiniSat starts normally
- a new run directory is created under `results/runs/`

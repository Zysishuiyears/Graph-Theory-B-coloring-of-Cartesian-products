# Bcoloring

本仓库用于研究笛卡尔积图上的 B-染色（B-coloring）与 B-色数（B-chromatic number）的 SAT 建模求解。项目采用参数驱动流程：在核心脚本中调整笛卡尔积图的图因子长度参数和数量参数（例如奇圈因子$C_{2n+1}$、路因子$P_m$），以及颜色参数（使用$k$种颜色进行染色），调用 SAT 求解器后，将结果输出到 `results/runs/`，再人工分析解文件。

## 项目简介

- 研究主题：笛卡尔积图 B-染色的可行性判定与染色模式搜索。
- 核心方法：将约束编码为 CNF，调用 CryptoMiniSat5 求解。
- 当前主线代码：`src/active/general_product_search.py`。
- 兼容入口：`code/general_product_search.py`（转发到主线脚本，便于旧命令继续使用）。

## Project Status

- 主维护入口：`src/active/general_product_search.py`。
- `code/general_product_search.py` 仅作为兼容入口壳，不承载主维护逻辑。
- `src/legacy/` 与 `scripts/legacy/` 为历史归档，默认不作为当前维护目标。

## Quick Start

### 1) 前置条件

- Python 3.10+
- [CryptoMiniSat5](https://github.com/msoos/cryptominisat)
  - 方式 A：可执行文件已在系统 `PATH` 中
  - 方式 B：通过环境变量 `CMSAT_PATH` 显式指定路径

### 2) 运行默认示例

```bash
python src/active/general_product_search.py
```

默认按脚本 `__main__` 中的示例参数运行。常规输出通常写入 `results/runs/`，重要结果可手工整理到 `results/snapshots/`。

### 3) 运行自定义最小示例

```bash
python -c "from src.active.general_product_search import decide_existence_B; decide_existence_B(cycles=[11], paths=[5], k=4)"
```

可通过环境变量覆盖默认配置：

```bash
# Windows PowerShell
$env:CMSAT_PATH='E:\\cryptominisat5-win\\cryptominisat5.exe'
$env:BCOLOR_RESULTS_ROOT='results/runs'
```

```bash
# Linux/macOS
export CMSAT_PATH=/path/to/cryptominisat5
export BCOLOR_RESULTS_ROOT=results/runs
```

## 依赖说明

### 核心运行依赖（主线必需）

- Python 标准库（`os`、`subprocess`、`itertools`、`typing`、`datetime`、`random` 等）
- CryptoMiniSat5

### Legacy/可视化脚本可选依赖（仅历史脚本可能使用）

- `networkx`
- `matplotlib`
- `numpy`
- `pandas`
- `sympy`
- `pulp`
- `pycryptosat`

## 输出与结果约定

- `results/runs/`：原始运行输出与历史结果归档（默认忽略）。
- `results/snapshots/`：人工筛选的关键结果快照（允许跟踪）。
- `results/runs/legacy_2026/_early_outputs/`：早期结构化之前遗留的 CNF 产物归档。

## 仓库结构

```text
Bcoloring/
  src/
    active/                     # 当前主线研究代码
    legacy/                     # 历史版本代码归档
  code/
    general_product_search.py   # 兼容入口壳（转发到 src/active）
  results/
    runs/                       # 原始实验输出（ignored）
    snapshots/                  # 关键结果快照（tracked）
  scripts/
    legacy/                     # 历史辅助脚本
  docs/
    legacy/                     # 历史文档
  README.md
  .gitignore
```

## 运行配置

- `CMSAT_PATH`：可选，指定 SAT 求解器路径。默认 `cryptominisat5`（从 PATH 解析）。
- `BCOLOR_RESULTS_ROOT`：可选，指定结果输出根目录。默认 `results/runs`。

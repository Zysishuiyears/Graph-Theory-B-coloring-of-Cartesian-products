# Bcoloring

本仓库用于研究笛卡尔积图上的 B-染色问题，并用 SAT 建模与求解器完成可满足性判定和候选染色搜索。

当前工作流已经拆分为两套核心代码：
- 判定核心：`src/active/general_product_decide.py`
- 搜索核心：`src/active/general_product_search.py`
- 共享基础层：`src/active/general_product_core.py`

## Project Status

- 主维护判定入口：`src/active/general_product_decide.py`
- 主维护搜索入口：`src/active/general_product_search.py`
- 兼容入口：
  - `code/general_product_decide.py`
  - `code/general_product_search.py`
- `src/legacy/` 与 `scripts/legacy/` 为历史归档，不作为当前主维护目标。

## Quick Start

### 1) 前置条件

- Python 3.10+
- [CryptoMiniSat5](https://github.com/msoos/cryptominisat)
  - 已在系统 `PATH` 中，或
  - 通过环境变量 `CMSAT_PATH` 显式指定

### 2) 运行判定示例

```bash
python src/active/general_product_decide.py
```

或：

```bash
python -c "from src.active.general_product_decide import decide_existence_B; decide_existence_B(cycles=[4,6], paths=[], k=5)"
```

### 3) 运行搜索示例

```bash
python -c "from src.active.general_product_search import search_best; search_best(cycles=[4], paths=[], k=4)"
```

常规运行输出默认写入 `results/runs/`。重要结果可手工筛选后整理到 `results/snapshots/`。

### 4) 常用环境变量

```bash
# Windows PowerShell
$env:CMSAT_PATH='E:\\cryptominisat5-win\\cryptominisat5.exe'
$env:BCOLOR_RESULTS_ROOT='results/runs'
$env:CMSAT_SOLVER_MODE='portfolio'
$env:CMSAT_PORTFOLIO_WORKERS='4'
$env:BCOLOR_SYMMETRY_BREAKING='1'
```

## Remote CPU Server

如果你准备把判定任务迁到远程 Linux CPU 服务器运行，见：

- `docs/remote_cpu_server.md`

仓库内已经附带远程辅助脚本：

- `scripts/server/check_remote_env.sh`
- `scripts/server/run_decide.sh`
- `scripts/server/run_decide_tmux.sh`

## 仓库结构

```text
Bcoloring/
  src/
    active/
      general_product_core.py      # 图结构、CNF、DIMACS、单次 solver 接口
      general_product_decide.py    # 判定核心
      general_product_search.py    # 搜索核心
    legacy/                        # 历史版本归档
  code/
    general_product_decide.py      # 判定兼容入口壳
    general_product_search.py      # 搜索兼容入口壳
  results/
    runs/                          # 原始运行输出（ignored）
    snapshots/                     # 关键结果快照（tracked）
  scripts/
    legacy/                        # 历史辅助脚本
  docs/
    legacy/                        # 历史文档
```

## 模块职责

### 判定模块

`src/active/general_product_decide.py` 只负责：
- 构造判定 CNF
- 可选颜色对称性破缺
- 单实例 / portfolio 求解
- 输出 `instance.cnf`、`one_solution.txt`、`solver_summary.txt`

主入口函数：

```python
from src.active.general_product_decide import decide_existence_B
```

### 搜索模块

`src/active/general_product_search.py` 只负责：
- 枚举 SAT 解
- blocking clause 去重
- 自同构 canonical 去同构
- 按冻结评分准则选择当前“最佳”方案

当前搜索评分准则固定为：

```text
S = W_LAYER * L + W_TRANS * T - W_UNIF * U
```

其中：
- `L`：层规律度，按每个维度的线性拟合命中比例平均
- `T`：平移对称性比例，只统计周期维平移保持不变的比例
- `U`：颜色分布不均匀惩罚

主入口函数：

```python
from src.active.general_product_search import search_best
```

## 依赖说明

### 主线必需

- Python 标准库
- CryptoMiniSat5

### Legacy 或可视化脚本可选

- `networkx`
- `matplotlib`
- `numpy`
- `pandas`
- `sympy`
- `pulp`
- `pycryptosat`


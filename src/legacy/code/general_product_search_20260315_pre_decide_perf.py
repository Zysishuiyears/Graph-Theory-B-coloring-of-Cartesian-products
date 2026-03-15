import os
import subprocess
from datetime import datetime
from typing import List, Tuple, Optional, Set, Dict
from itertools import product
import random


# ============================================================
# =============== 用户配置区：你只需要改这里 ==================
# ============================================================

USE_SOLVER = "cmsat"
# Environment override is supported: CMSAT_PATH=<path-to-cryptominisat5>
CMSAT_PATH = os.getenv("CMSAT_PATH", "cryptominisat5")
CMSAT_THREADS = min(12, os.cpu_count() or 1)
CMSAT_VERB = 0

# 约束开关：是否对每个“基本 square(4-cycle)”要求同色至多一次
USE_C4_DISTINCT = True

# 搜索停止条件
MAX_SOLVER_CALLS = 10000
MAX_NO_NEW_CLASS = 500
MAX_NO_BEST_IMPROVE = 2000

# 评分权重（默认给一个“通用版本”：平移对称 + 均匀性 + 维度线性规律）
W_LAYER = 1.0
W_TRANS = 1.0
W_UNIF  = 0.1

# 高维自同构太多时，canonical 只采样这么多个 shift（flip 仍全枚举）
MAX_SHIFT_SAMPLES = 2000
RANDOM_SEED = 0

# Environment override is supported: BCOLOR_RESULTS_ROOT=<output-root>
RESULTS_ROOT = os.getenv("BCOLOR_RESULTS_ROOT", os.path.join("results", "runs"))


# ============================================================
# ===================== 图结构：任意(圈/路)笛卡尔积 =====================
# ============================================================

def var_color(e: int, c: int, k: int) -> int:
    # x_{e,c} in [1..E*k]
    return e * k + c + 1

def var_rep(e: int, c: int, E: int, k: int) -> int:
    # rep_{e,c} in [E*k+1..2*E*k]
    return E * k + e * k + c + 1


class ProductGraph:
    """
    G = □_{d=0..D-1} (C_{L[d]} if periodic[d] else P_{L[d]})
    顶点：坐标 tuple (x0,...,x_{D-1})
    边：在某一维 +1 的无向边（周期维 wrap）
    """
    def __init__(self, sizes: List[int], periodic: List[bool]):
        assert len(sizes) == len(periodic)
        self.sizes = sizes
        self.periodic = periodic
        self.D = len(sizes)

        # 生成顶点列表（用于遍历）
        self.vertices = list(product(*[range(L) for L in sizes]))
        self.VN = len(self.vertices)

        # 边列表与查找表
        self.edges: List[Tuple[Tuple[int, ...], Tuple[int, ...]]] = []
        self.edge_of_uv: Dict[Tuple[Tuple[int, ...], Tuple[int, ...]], int] = {}

        self._build_edges()

        self.E = len(self.edges)

        # 每个顶点 incident 的边
        self.incident: Dict[Tuple[int, ...], List[int]] = {v: [] for v in self.vertices}
        for eid, (u, v) in enumerate(self.edges):
            self.incident[u].append(eid)
            self.incident[v].append(eid)

        # line-graph 邻接（共享端点）
        self.neigh = [set() for _ in range(self.E)]
        for v in self.vertices:
            inc = self.incident[v]
            for i in range(len(inc)):
                for j in range(i + 1, len(inc)):
                    e1, e2 = inc[i], inc[j]
                    self.neigh[e1].add(e2)
                    self.neigh[e2].add(e1)
        self.neigh = [list(s) for s in self.neigh]

        # square(4-cycle) 列表（可选使用）
        self.C4s = self._enumerate_squares()

    def _canon_uv(self, a: Tuple[int, ...], b: Tuple[int, ...]):
        return (a, b) if a <= b else (b, a)

    def _step(self, v: Tuple[int, ...], dim: int, delta: int) -> Optional[Tuple[int, ...]]:
        L = self.sizes[dim]
        x = v[dim]
        y = x + delta
        if self.periodic[dim]:
            y %= L
        else:
            if not (0 <= y < L):
                return None
        vv = list(v)
        vv[dim] = y
        return tuple(vv)

    def _build_edges(self):
        # 只枚举“+1方向”的边，避免重复
        for v in self.vertices:
            for d in range(self.D):
                L = self.sizes[d]
                if self.periodic[d]:
                    # 周期维：每个顶点都有一条 +1 边
                    w = self._step(v, d, +1)
                    assert w is not None
                else:
                    # 路维：只有 x < L-1 才有 +1
                    if v[d] >= L - 1:
                        continue
                    w = self._step(v, d, +1)
                    assert w is not None

                key = self._canon_uv(v, w)
                if key in self.edge_of_uv:
                    # 理论上不会发生（因为我们只用+1）
                    continue
                eid = len(self.edges)
                self.edges.append((key[0], key[1]))
                self.edge_of_uv[key] = eid

    def _enumerate_squares(self) -> List[List[int]]:
        """
        枚举所有“基本 square”：
          选两维 a<b，选 base 顶点 v00，使得在 a、b 维都能 +1（周期维总能，路维需 <L-1）
          v10=v00+ea, v01=v00+eb, v11=v00+ea+eb
          边：v00-v10, v00-v01, v01-v11, v10-v11
        """
        C4s = []
        for a in range(self.D):
            for b in range(a + 1, self.D):
                for v00 in self.vertices:
                    # 检查 a,b 维是否能 +1
                    if (not self.periodic[a]) and (v00[a] >= self.sizes[a] - 1):
                        continue
                    if (not self.periodic[b]) and (v00[b] >= self.sizes[b] - 1):
                        continue

                    v10 = self._step(v00, a, +1)
                    v01 = self._step(v00, b, +1)
                    v11 = self._step(v10, b, +1)  # type: ignore
                    if v10 is None or v01 is None or v11 is None:
                        continue

                    e_a0 = self.edge_of_uv[self._canon_uv(v00, v10)]
                    e_b0 = self.edge_of_uv[self._canon_uv(v00, v01)]
                    e_a1 = self.edge_of_uv[self._canon_uv(v01, v11)]
                    e_b1 = self.edge_of_uv[self._canon_uv(v10, v11)]
                    C4s.append([e_a0, e_b0, e_a1, e_b1])
        # add 1D C4 cycles for periodic dimensions with length 4
        for d in range(self.D):
            if not (self.periodic[d] and self.sizes[d] == 4):
                continue
            other_ranges = [range(self.sizes[i]) for i in range(self.D) if i != d]
            for fixed in product(*other_ranges):
                base = [0] * self.D
                idx = 0
                for i in range(self.D):
                    if i == d:
                        continue
                    base[i] = fixed[idx]
                    idx += 1
                verts = []
                for t in range(4):
                    vv = list(base)
                    vv[d] = t
                    verts.append(tuple(vv))
                edges = []
                for t in range(4):
                    u = verts[t]
                    v = verts[(t + 1) % 4]
                    edges.append(self.edge_of_uv[self._canon_uv(u, v)])
                C4s.append(edges)
        return C4s


# ============================================================
# =================== CNF 生成（B-边染色） ====================
# ============================================================
def build_cnf_proper_plus_rainbowC4(graph: ProductGraph, k: int) -> Tuple[int, List[List[int]]]:
    """
    变量：x_{e,c} (e=0..E-1, c=0..k-1)
    约束：
      1) 每条边恰好一种颜色
      2) 任意 incident edges 不同色（proper edge-coloring）
      3) 任意 4-cycle 的四条边两两不同色（rainbow C4）
    """
    E = graph.E
    clauses: List[List[int]] = []

    # 1) 每边恰一色
    for e in range(E):
        clauses.append([var_color(e, c, k) for c in range(k)])  # at least one
        for c1 in range(k):
            for c2 in range(c1 + 1, k):
                clauses.append([-var_color(e, c1, k), -var_color(e, c2, k)])  # at most one

    # 2) incident 边异色
    for v in graph.vertices:
        inc = graph.incident[v]
        for c in range(k):
            for i in range(len(inc)):
                for j in range(i + 1, len(inc)):
                    clauses.append([-var_color(inc[i], c, k), -var_color(inc[j], c, k)])

    # 3) 4-cycle rainbow：四条边两两不同色
    # 你原来写的“对每个颜色 at-most-one”也能推出四条边全不同（因为每种颜色最多出现一次）
    # 这里用“成对禁止同色”更直接：
    for cyc in graph.C4s:
        for i in range(4):
            for j in range(i + 1, 4):
                for c in range(k):
                    clauses.append([-var_color(cyc[i], c, k), -var_color(cyc[j], c, k)])

    nvars = E * k
    return nvars, clauses

def write_dimacs(nvars: int, clauses: List[List[int]], path: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"p cnf {nvars} {len(clauses)}\n")
        for cls in clauses:
            f.write(" ".join(map(str, cls)) + " 0\n")


def run_cmsat(cnf_path: str) -> Optional[List[int]]:
    cmd = [
        CMSAT_PATH,
        "-t", str(CMSAT_THREADS),
        "--verb", str(CMSAT_VERB),
        cnf_path
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    out = proc.stdout + "\n" + proc.stderr
    if "UNSAT" in out:
        return None

    model = []
    for line in out.splitlines():
        if line.startswith("v ") or line.startswith("V "):
            for t in line.split()[1:]:
                if t == "0":
                    continue
                try:
                    model.append(int(t))
                except:
                    pass
    return model


def decode_model_to_edge_colors(model: List[int], E: int, k: int) -> List[int]:
    max_color_var = E * k
    pos = set(v for v in model if 0 < v <= max_color_var)

    edge_color = [-1] * E
    for e in range(E):
        for c in range(k):
            if var_color(e, c, k) in pos:
                edge_color[e] = c
                break
        if edge_color[e] < 0:
            raise RuntimeError("decode failed: some edge has no true color.")
    return edge_color


# ============================================================
# ============== 去同构：颜色置换 +（周期维平移）+ 反射 ==============
# ============================================================

def canonical_relabel_key(colors: List[int]) -> Tuple[int, ...]:
    mp = {}
    nxt = 0
    out = []
    for x in colors:
        if x not in mp:
            mp[x] = nxt
            nxt += 1
        out.append(mp[x])
    return tuple(out)


def transform_vertex(v: Tuple[int, ...], sizes: List[int], periodic: List[bool],
                     shift: List[int], flip: List[bool]) -> Tuple[int, ...]:
    """
    flip:
      - 周期维：x -> (-x) mod L
      - 路维：x -> (L-1-x)
    shift:
      - 只对周期维生效：x -> (x + shift[d]) mod L
      - 路维 shift 应为 0（否则不是自同构）
    """
    vv = list(v)
    for d in range(len(vv)):
        L = sizes[d]
        x = vv[d]
        if flip[d]:
            if periodic[d]:
                x = (-x) % L
            else:
                x = (L - 1 - x)
        if periodic[d]:
            x = (x + shift[d]) % L
        vv[d] = x
    return tuple(vv)


def apply_aut_to_coloring(graph: ProductGraph, edge_colors: List[int],
                          shift: List[int], flip: List[bool]) -> List[int]:
    """
    把 edge_colors 经过自同构变换，得到新 coloring（按 graph.edges 的固定顺序输出）
    """
    E = graph.E
    new_colors = [-1] * E
    for eid, (u, v) in enumerate(graph.edges):
        u2 = transform_vertex(u, graph.sizes, graph.periodic, shift, flip)
        v2 = transform_vertex(v, graph.sizes, graph.periodic, shift, flip)
        key = (u2, v2) if u2 <= v2 else (v2, u2)
        eid2 = graph.edge_of_uv[key]
        new_colors[eid2] = edge_colors[eid]
    assert all(c >= 0 for c in new_colors)
    return new_colors


def canonical_key_under_aut_color(graph: ProductGraph, edge_colors: List[int]) -> Tuple[int, ...]:
    """
    自同构枚举：周期维平移 × 全维反射（不做维度交换）
    若周期维平移总数太大，则随机采样 MAX_SHIFT_SAMPLES 个。
    """
    random.seed(RANDOM_SEED)

    cycle_dims = [d for d, per in enumerate(graph.periodic) if per]
    D = graph.D

    # 生成 shift 列表
    total_shifts = 1
    for d in cycle_dims:
        total_shifts *= graph.sizes[d]

    def iter_shifts():
        if total_shifts <= MAX_SHIFT_SAMPLES:
            # 全枚举
            for vals in product(*[range(graph.sizes[d]) for d in cycle_dims]):
                shift = [0] * D
                for idx, d in enumerate(cycle_dims):
                    shift[d] = vals[idx]
                yield shift
        else:
            # 采样
            for _ in range(MAX_SHIFT_SAMPLES):
                shift = [0] * D
                for d in cycle_dims:
                    shift[d] = random.randrange(graph.sizes[d])
                yield shift

    best = None
    for shift in iter_shifts():
        for flip_mask in product([False, True], repeat=D):
            flip = list(flip_mask)
            colors2 = apply_aut_to_coloring(graph, edge_colors, shift, flip)
            key = canonical_relabel_key(colors2)
            if best is None or key < best:
                best = key
    assert best is not None
    return best

def decide_existence_B(cycles: List[int], paths: List[int], k: int = 5, dump_one_solution: bool = True):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    name = "x".join([f"C{c}" for c in cycles] + [f"P{p}" for p in paths]) + f"_k{k}"
    out_dir = os.path.join(RESULTS_ROOT, f"{stamp}_{name}_DECIDE_B")
    os.makedirs(out_dir, exist_ok=True)

    sizes = cycles + paths
    periodic = [True] * len(cycles) + [False] * len(paths)
    graph = ProductGraph(sizes, periodic)

    max_deg = max(len(graph.incident[v]) for v in graph.vertices) if graph.vertices else 0
    print(f"Graph: D={graph.D}, |V|={graph.VN}, |E|={graph.E}, maxdeg={max_deg}, #C4={len(graph.C4s)}")
    print(f"Decide existence for: proper edge-coloring + rainbow 4-cycles, k={k}")

    nvars, clauses = build_cnf_proper_plus_rainbowC4(graph, k)
    cnf_path = os.path.join(out_dir, "instance.cnf")
    write_dimacs(nvars, clauses, cnf_path)

    model = run_cmsat(cnf_path)
    if model is None:
        print("RESULT: UNSAT (不存在满足你两条定义的 k-色方案)")
        return False, out_dir, None

    print("RESULT: SAT (存在满足你两条定义的 k-色方案)")

    if not dump_one_solution:
        return True, out_dir, None

    edge_colors = decode_model_to_edge_colors(model, graph.E, k)
    sol_path = os.path.join(out_dir, "one_solution.txt")
    with open(sol_path, "w", encoding="utf-8") as f:
        f.write(f"# SAT solution: sizes={graph.sizes}, periodic={graph.periodic}, k={k}\n")
        f.write("# eid : u -- v : color(1..k)\n")
        for eid, (u, v) in enumerate(graph.edges):
            f.write(f"{eid} : {u} -- {v} : {edge_colors[eid] + 1}\n")

    print("One solution dumped to:", sol_path)
    return True, out_dir, sol_path


# ============================================================
# ======================== 通用评分函数（可自行替换） ======================
# ============================================================

def uniform_penalty(edge_colors: List[int], k: int) -> float:
    cnt = [0] * k
    for c in edge_colors:
        cnt[c] += 1
    ideal = len(edge_colors) / k
    return sum((x - ideal) ** 2 for x in cnt)

def translation_symmetry_score(graph: ProductGraph, edge_colors: List[int]) -> float:
    """
    只统计“周期维平移”能保持 coloring 不变的比例。
    """
    cycle_dims = [d for d, per in enumerate(graph.periodic) if per]
    if not cycle_dims:
        return 0.0

    D = graph.D
    total = 0
    good = 0

    for vals in product(*[range(graph.sizes[d]) for d in cycle_dims]):
        if all(v == 0 for v in vals):
            continue
        shift = [0] * D
        for idx, d in enumerate(cycle_dims):
            shift[d] = vals[idx]
        colors2 = apply_aut_to_coloring(graph, edge_colors, shift, [False] * D)
        total += 1
        if colors2 == edge_colors:
            good += 1

    return good / total if total > 0 else 0.0

def layer_regular_score_generic(graph: ProductGraph, edge_colors: List[int], k: int) -> float:
    """
    通用“线性层规律”评分（给你一个可用默认）：
    对每个维度 d，把沿 d 的边按“固定其它坐标”的每条线分组，
    在该线中尝试拟合：color(t) = a*t + b (mod k)，取匹配最多的比例。
    最后对所有维度取平均。
    """
    # 先为每条边记录：它属于哪个维度 d，以及它在线上的位置 t（取较小端点的坐标当 t）
    # 通过比较 u,v 哪一维不同得到维度 d
    per_dim_lines: Dict[int, Dict[Tuple[int, ...], List[Tuple[int, int]]]] = {d: {} for d in range(graph.D)}
    # line_key: 去掉维度 d 的坐标，剩余坐标作为 key
    for eid, (u, v) in enumerate(graph.edges):
        diff = [i for i in range(graph.D) if u[i] != v[i]]
        if len(diff) != 1:
            continue
        d = diff[0]
        # t 用 min(u[d],v[d]) 对路维OK；周期维只是一个约定（足够做“局部线性”指标）
        t = min(u[d], v[d])
        line_key = tuple(u[i] for i in range(graph.D) if i != d)
        per_dim_lines[d].setdefault(line_key, []).append((t, edge_colors[eid]))

    scores = []
    for d in range(graph.D):
        lines = per_dim_lines[d].values()
        if not lines:
            continue
        hit_sum = 0
        len_sum = 0
        for line in lines:
            # line: list of (t, color)
            line_sorted = sorted(line, key=lambda x: x[0])
            ts = [x[0] for x in line_sorted]
            cs = [x[1] for x in line_sorted]
            Lline = len(cs)
            if Lline == 0:
                continue

            best = 0
            for a in range(k):
                for b in range(k):
                    hits = 0
                    for idx in range(Lline):
                        if cs[idx] == (a * ts[idx] + b) % k:
                            hits += 1
                    if hits > best:
                        best = hits
            hit_sum += best
            len_sum += Lline
        scores.append(hit_sum / len_sum if len_sum else 0.0)

    return sum(scores) / len(scores) if scores else 0.0

def score_coloring(graph: ProductGraph, edge_colors: List[int], k: int) -> Tuple[float, float, float, float]:
    L = layer_regular_score_generic(graph, edge_colors, k)
    T = translation_symmetry_score(graph, edge_colors)
    U = uniform_penalty(edge_colors, k)
    S = W_LAYER * L + W_TRANS * T - W_UNIF * U
    return S, L, T, U


# ============================================================
# =================== blocking clause（避免同一解重复） =================
# ============================================================

def make_blocking_clause(edge_colors: List[int], k: int) -> List[int]:
    # OR_e ¬x_{e, color[e]}
    return [-var_color(e, edge_colors[e], k) for e in range(len(edge_colors))]


# ============================================================
# =============================== 主函数：搜索型停止找最优 =============================
# ============================================================

def search_best(cycles: List[int], paths: List[int], k: int = 5):
    """
    cycles: [c1,c2,...] 表示 C_c1 □ C_c2 □ ...
    paths : [p1,p2,...] 表示 P_p1 □ P_p2 □ ...
    """
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    name = "x".join([f"C{c}" for c in cycles] + [f"P{p}" for p in paths]) + f"_k{k}"
    out_dir = os.path.join(RESULTS_ROOT, f"{stamp}_{name}")
    os.makedirs(out_dir, exist_ok=True)

    sizes = cycles + paths
    periodic = [True] * len(cycles) + [False] * len(paths)

    graph = ProductGraph(sizes, periodic)
    max_deg = max(len(graph.incident[v]) for v in graph.vertices) if graph.vertices else 0
    print(f"Graph: D={graph.D}, |V|={graph.VN}, |E|={graph.E}, maxdeg={max_deg}, squares={len(graph.C4s)}")

    nvars, base_clauses = build_base_cnf(graph, k, use_c4_distinct=USE_C4_DISTINCT)

    blocking_clauses: List[List[int]] = []
    seen_keys: Set[Tuple[int, ...]] = set()

    best = None  # (score, L, T, U, edge_colors)
    no_new = 0
    no_improve = 0

    cnf_path = os.path.join(out_dir, "working.cnf")

    for call_id in range(1, MAX_SOLVER_CALLS + 1):
        all_clauses = base_clauses + blocking_clauses
        write_dimacs(nvars, all_clauses, cnf_path)

        model = run_cmsat(cnf_path)
        if model is None:
            print(f"[{call_id}] UNSAT：解空间已被 blocking 穷尽。")
            break

        edge_colors = decode_model_to_edge_colors(model, graph.E, k)

        # 先 block 这个具体解
        blocking_clauses.append(make_blocking_clause(edge_colors, k))

        # canonical 去同构
        key = canonical_key_under_aut_color(graph, edge_colors)
        if key in seen_keys:
            no_new += 1
            no_improve += 1
        else:
            seen_keys.add(key)
            no_new = 0

            score, L, T, U = score_coloring(graph, edge_colors, k)
            if best is None or score > best[0]:
                best = (score, L, T, U, edge_colors)
                no_improve = 0
                print(f"[{call_id}] NEW BEST  S={score:.4f}  L={L:.3f}  T={T:.3f}  U={U:.1f}  unique={len(seen_keys)}")
            else:
                no_improve += 1

        if no_new >= MAX_NO_NEW_CLASS:
            print(f"[{call_id}] Stop: 连续 {MAX_NO_NEW_CLASS} 次没有新等价类。unique={len(seen_keys)}")
            break
        if no_improve >= MAX_NO_BEST_IMPROVE:
            print(f"[{call_id}] Stop: 连续 {MAX_NO_BEST_IMPROVE} 次 best 不提升。unique={len(seen_keys)}")
            break

        if call_id % 200 == 0:
            print(f"[{call_id}] progress: unique={len(seen_keys)} blocking={len(blocking_clauses)} no_new={no_new} no_improve={no_improve}")

    if best is None:
        print("没有找到任何可行解（SAT 直接 UNSAT）。")
        return None

    score, L, T, U, best_colors = best

    # 输出 best_solution.txt（输出每条边的端点与颜色）
    sol_path = os.path.join(out_dir, "best_solution.txt")
    with open(sol_path, "w", encoding="utf-8") as f:
        f.write(f"# Best S={score:.6f}, LayerReg={L:.6f}, TransSym={T:.6f}, Penalty={U:.3f}\n")
        f.write(f"# Graph sizes={graph.sizes}, periodic={graph.periodic}, E={graph.E}\n")
        f.write("# Format: eid : u -> v : color(1..k)\n")
        for eid, (u, v) in enumerate(graph.edges):
            f.write(f"{eid} : {u} -- {v} : {best_colors[eid] + 1}\n")

    print("\n==== DONE ====")
    print("Results dir:", out_dir)
    print("Best solution:", sol_path)
    return out_dir, best


if __name__ == "__main__":
    # 示例：C3 □ P7 （就是你原来“圈乘路”的典型）
    # search_best(cycles=[7], paths=[3], k=4)

    # 再举个高维示例：C3 □ C5 □ P4
    # search_best(cycles=[3, 5], paths=[4], k=6)
    # C3 □ P7
    decide_existence_B(cycles=[4,6,6,6], paths=[], k=8)

    # C3 □ C5 □ P4
    # decide_existence_B(cycles=[3,5], paths=[4], k=6)


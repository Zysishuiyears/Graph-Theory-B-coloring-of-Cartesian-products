import os
import subprocess
from datetime import datetime
from typing import List, Tuple, Optional, Set

import matplotlib.pyplot as plt
import networkx as nx


# ============================================================
# =============== 用户配置区：你只需要改这里 ==================
# ============================================================

USE_SOLVER = "cmsat"
CMSAT_PATH = r"E:\cryptominisat5-win\cryptominisat5.exe"

# 约束开关：是否要求每个 4-cycle “同色最多一次”（你的旧代码里的 C4 约束）
USE_C4_DISTINCT = True

# 搜索停止条件（搜索型停止）
MAX_SOLVER_CALLS = 200000          # 最多调用 solver 的次数（上限保护）
MAX_NO_NEW_CLASS = 5000            # 连续这么多次没有出现“新等价类”就停
MAX_NO_BEST_IMPROVE = 20000        # 连续这么多次 best 分数不提升就停

# 评分权重（按你之前的“层规律 + 平移对称 + 均匀性”）
W_LAYER = 1.0
W_TRANS = 0.5
W_UNIF  = 0.01

# 结果目录根
RESULTS_ROOT = "results"


# ============================================================
# ===================== 图结构：C_m □ C_n =====================
# ============================================================

def num_edges(m, n): return 2 * m * n

def edge_id_H(i, j, m, n):
    """水平边 H(i,j): (i,j)-(i,j+1)"""
    return i * n + j

def edge_id_V(i, j, m, n):
    """垂直边 V(i,j): (i,j)-(i+1,j)"""
    return m * n + i * n + j

def incident_edges_at_vertex(i, j, m, n):
    """顶点 (i,j) incident 的 4 条边"""
    return [
        edge_id_H(i, j, m, n),
        edge_id_H(i, (j - 1) % n, m, n),
        edge_id_V(i, j, m, n),
        edge_id_V((i - 1) % m, j, m, n),
    ]

def build_edge_neighbors(m, n):
    """line-graph 邻接：共享端点则相邻"""
    E = num_edges(m, n)
    neigh = [set() for _ in range(E)]
    for i in range(m):
        for j in range(n):
            inc = incident_edges_at_vertex(i, j, m, n)
            for a in range(len(inc)):
                for b in range(a + 1, len(inc)):
                    e1, e2 = inc[a], inc[b]
                    neigh[e1].add(e2)
                    neigh[e2].add(e1)
    return [list(s) for s in neigh]

def enumerate_C4_edge_sets(m, n):
    """
    枚举 torus 上的 mn 个基本 4-cycle（一个格子一个）
    cycle: (i,j)-(i,j+1)-(i+1,j+1)-(i+1,j)-(i,j)
    """
    C4s = []
    for i in range(m):
        for j in range(n):
            e_left  = edge_id_V(i, j, m, n)                 # (i,j)-(i+1,j)
            e_bot   = edge_id_H(i, j, m, n)                 # (i,j)-(i,j+1)
            e_right = edge_id_V(i, (j+1) % n, m, n)         # (i,j+1)-(i+1,j+1)
            e_top   = edge_id_H((i+1) % m, j, m, n)         # (i+1,j)-(i+1,j+1)
            C4s.append([e_left, e_bot, e_right, e_top])
    return C4s


# ============================================================
# ===================== SAT 变量编号 ==========================
# ============================================================

def var_color(e, c, k):
    """x_{e,c} : 边 e 取颜色 c（c=0..k-1），范围 1..E*k"""
    return e * k + c + 1

def var_rep(e, c, E, k):
    """rep_{e,c} : 边 e 作为颜色 c 的 B-边代表，范围 E*k+1..2*E*k"""
    return E * k + e * k + c + 1


# ============================================================
# =================== CNF 生成（B-边染色） ====================
# ============================================================

def build_base_cnf(m, n, k, use_c4_distinct=True) -> Tuple[int, List[List[int]]]:
    """
    base CNF（不要加 symmetry breaking；我们做去同构是在 Python 层完成）
    约束：
      - 每边恰一色
      - incident 异色（proper edge-coloring）
      - (可选) 每个 4-cycle 对每个颜色 at-most-one（你旧代码那条）
      - B-coloring 支配：每色至少一条 B-边 rep(e,c)，且 rep(e,c) -> 邻域覆盖所有其他色
    """
    E = num_edges(m, n)
    neigh = build_edge_neighbors(m, n)
    clauses: List[List[int]] = []

    # 1) 每边至少一色 + 至多一色
    for e in range(E):
        clauses.append([var_color(e, c, k) for c in range(k)])  # atleast1
        for c1 in range(k):
            for c2 in range(c1 + 1, k):
                clauses.append([-var_color(e, c1, k), -var_color(e, c2, k)])

    # 2) incident 异色
    for i in range(m):
        for j in range(n):
            inc = incident_edges_at_vertex(i, j, m, n)
            for c in range(k):
                for a in range(len(inc)):
                    for b in range(a + 1, len(inc)):
                        clauses.append([-var_color(inc[a], c, k), -var_color(inc[b], c, k)])

    # 3) (可选) 每个 4-cycle 四色（对每个颜色至多一次）
    if use_c4_distinct:
        C4s = enumerate_C4_edge_sets(m, n)
        for cyc in C4s:
            for c in range(k):
                for a in range(4):
                    for b in range(a + 1, 4):
                        clauses.append([-var_color(cyc[a], c, k), -var_color(cyc[b], c, k)])

    # 4) B-coloring 支配 via rep(e,c)
    # (a) 每种颜色至少一个代表
    for c in range(k):
        clauses.append([var_rep(e, c, E, k) for e in range(E)])

    # (b) rep(e,c) -> color(e,c)
    for e in range(E):
        for c in range(k):
            clauses.append([-var_rep(e, c, E, k), var_color(e, c, k)])

    # (c) rep(e,c) -> 对每个 d!=c，邻边中至少有一个颜色 d
    for e in range(E):
        N = neigh[e]
        for c in range(k):
            for d in range(k):
                if d == c:
                    continue
                clause = [-var_rep(e, c, E, k)]
                clause += [var_color(f, d, k) for f in N]
                clauses.append(clause)

    nvars = 2 * E * k
    return nvars, clauses


def write_dimacs(nvars: int, clauses: List[List[int]], path: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"p cnf {nvars} {len(clauses)}\n")
        for cls in clauses:
            f.write(" ".join(map(str, cls)) + " 0\n")


def run_cmsat(cnf_path: str) -> Optional[List[int]]:
    proc = subprocess.run([CMSAT_PATH, cnf_path], capture_output=True, text=True)
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


def decode_model_to_HV(model: List[int], m: int, n: int, k: int) -> Tuple[List[List[int]], List[List[int]]]:
    """只解码颜色变量（不看 rep），返回 H,V 颜色矩阵（0..k-1）"""
    E = num_edges(m, n)
    max_color_var = E * k
    pos = set(v for v in model if 0 < v <= max_color_var)

    edge_color = [None] * E
    for e in range(E):
        for c in range(k):
            v = var_color(e, c, k)
            if v in pos:
                edge_color[e] = c
                break
        if edge_color[e] is None:
            raise RuntimeError("Model decode failed: some edge has no color true.")

    H = [[0]*n for _ in range(m)]
    V = [[0]*n for _ in range(m)]
    for i in range(m):
        for j in range(n):
            H[i][j] = edge_color[edge_id_H(i, j, m, n)]
            V[i][j] = edge_color[edge_id_V(i, j, m, n)]
    return H, V


# ============================================================
# ================= 去同构：颜色置换 + 平移 + 翻转 =================
# ============================================================

def canonical_relabel_key_from_HV(H, V, m, n):
    """
    颜色置换去除：按固定扫描顺序，把“首次出现的颜色”依次重编号 0,1,2,...
    返回扁平 tuple key（长度 2mn）
    """
    mp = {}
    nxt = 0

    def relabel(x):
        nonlocal nxt
        if x not in mp:
            mp[x] = nxt
            nxt += 1
        return mp[x]

    seq = []
    for i in range(m):
        for j in range(n):
            seq.append(relabel(H[i][j]))
    for i in range(m):
        for j in range(n):
            seq.append(relabel(V[i][j]))
    return tuple(seq)


def shift_flip_HV(H, V, m, n, a, b, flip_i, flip_j):
    """
    对 H,V 做 D_m×D_n 的自同构：先反射再平移（不交换因子）
    这是在矩阵上直接做，速度快。
    """
    H2 = [[0]*n for _ in range(m)]
    V2 = [[0]*n for _ in range(m)]
    for i in range(m):
        for j in range(n):
            ii = (-i) % m if flip_i else i
            jj = (-j) % n if flip_j else j
            i2 = (ii + a) % m
            j2 = (jj + b) % n
            H2[i2][j2] = H[i][j]
            V2[i2][j2] = V[i][j]
    return H2, V2


def canonical_key_under_aut_color(H, V, m, n):
    """
    在所有 4mn 个自同构下取 canonical（字典序最小）：
      - 每个自同构变换后先做颜色 canonical relabel
    """
    best = None
    for a in range(m):
        for b in range(n):
            for fi in (False, True):
                for fj in (False, True):
                    H2, V2 = shift_flip_HV(H, V, m, n, a, b, fi, fj)
                    key = canonical_relabel_key_from_HV(H2, V2, m, n)
                    if best is None or key < best:
                        best = key
    return best


# ============================================================
# ======================== 评分函数（规律指标） ======================
# ============================================================

def horizontal_layer_score(H, k):
    m = len(H)
    n = len(H[0])
    best_over_a = 0.0
    for a in range(1, k):
        total_hits = 0
        for i in range(m):
            row = H[i]
            best_row = 0
            for b in range(k):
                hits = 0
                for j in range(n):
                    if row[j] == (a*j + b) % k:
                        hits += 1
                if hits > best_row:
                    best_row = hits
            total_hits += best_row
        best_over_a = max(best_over_a, total_hits / (m*n))
    return best_over_a

def vertical_layer_score(V, k):
    m = len(V)
    n = len(V[0])
    best_over_a = 0.0
    for a in range(1, k):
        total_hits = 0
        for j in range(n):
            best_col = 0
            for b in range(k):
                hits = 0
                for i in range(m):
                    if V[i][j] == (a*i + b) % k:
                        hits += 1
                if hits > best_col:
                    best_col = hits
            total_hits += best_col
        best_over_a = max(best_over_a, total_hits / (m*n))
    return best_over_a

def layer_regular_score(H, V, k):
    return 0.5 * (horizontal_layer_score(H, k) + vertical_layer_score(V, k))

def count_edge_periods(H, V):
    m, n = len(H), len(H[0])
    cnt = 0
    for a in range(m):
        for b in range(n):
            if a == 0 and b == 0:
                continue
            ok = True
            for i in range(m):
                if not ok:
                    break
                for j in range(n):
                    if H[i][j] != H[(i+a) % m][(j+b) % n]:
                        ok = False
                        break
                    if V[i][j] != V[(i+a) % m][(j+b) % n]:
                        ok = False
                        break
            if ok:
                cnt += 1
    return cnt

def translation_symmetry_score(H, V):
    m, n = len(H), len(H[0])
    maxP = m*n - 1
    if maxP <= 0:
        return 0.0
    return count_edge_periods(H, V) / maxP

def uniform_penalty_edges(H, V, k):
    m, n = len(H), len(H[0])
    cnt = [0]*k
    for i in range(m):
        for j in range(n):
            cnt[H[i][j]] += 1
            cnt[V[i][j]] += 1
    ideal = (2*m*n)/k
    return sum((x - ideal)**2 for x in cnt)

def score_edge_coloring(H, V, k):
    L = layer_regular_score(H, V, k)
    T = translation_symmetry_score(H, V)
    U = uniform_penalty_edges(H, V, k)
    score = W_LAYER*L + W_TRANS*T - W_UNIF*U
    return score, L, T, U


# ============================================================
# =================== blocking clause（避免同一解重复） =================
# ============================================================

def make_blocking_clause_for_HV(H, V, m, n, k):
    """
    阻止“完全同一边染色”再次出现：
      OR_e  ¬x_{e, color_in_solution(e)}
    """
    clause = []
    for i in range(m):
        for j in range(n):
            eH = edge_id_H(i, j, m, n)
            clause.append(-var_color(eH, H[i][j], k))
    for i in range(m):
        for j in range(n):
            eV = edge_id_V(i, j, m, n)
            clause.append(-var_color(eV, V[i][j], k))
    return clause


# ============================================================
# ========================= 画图：边染色 + 边标签 =======================
# ============================================================

def draw_graph_edge_coloring(H, V, m, n, k, path):
    G = nx.Graph()
    for i in range(m):
        for j in range(n):
            G.add_node((i, j))

    for i in range(m):
        for j in range(n):
            G.add_edge((i, j), (i, (j+1) % n), etype="H", i=i, j=j)
            G.add_edge((i, j), ((i+1) % m, j), etype="V", i=i, j=j)

    # 网格布局（你要圆形布局可改：pos = nx.circular_layout(G)）
    pos = {(i, j): (j, -i) for i in range(m) for j in range(n)}
    cmap = plt.get_cmap("tab10")

    edges = []
    edge_colors = []
    edge_labels = {}

    for u, v, data in G.edges(data=True):
        i = data["i"]; j = data["j"]
        if data["etype"] == "H":
            c = H[i][j]
        else:
            c = V[i][j]
        edges.append((u, v))
        edge_colors.append(cmap(c % 10))
        edge_labels[(u, v)] = str(c + 1)

    plt.figure(figsize=(8, 7))
    nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=edge_colors, width=2.2)

    nx.draw_networkx_nodes(
        G, pos,
        node_color="white",
        edgecolors="black",
        linewidths=0.8,
        node_size=520
    )

    # 顶点标签
    labels = {(i, j): f"({i},{j})" for i in range(m) for j in range(n)}
    nx.draw_networkx_labels(G, pos, labels, font_size=8)

    # 边标签（颜色编号）
    nx.draw_networkx_edge_labels(
        G, pos,
        edge_labels=edge_labels,
        font_size=8,
        font_color="black",
        label_pos=0.5
    )

    from matplotlib.patches import Patch
    legend_handles = [Patch(color=cmap(c % 10), label=f"color {c+1}") for c in range(k)]
    plt.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.)

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()


# ============================================================
# =============================== 主函数：搜索型停止找最优 =============================
# ============================================================

def search_best(m: int, n: int, k: int = 5):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = os.path.join(RESULTS_ROOT, f"{stamp}_C{m}xC{n}_k{k}")
    os.makedirs(out_dir, exist_ok=True)

    nvars, base_clauses = build_base_cnf(m, n, k, use_c4_distinct=USE_C4_DISTINCT)

    blocking_clauses: List[List[int]] = []
    seen_keys: Set[Tuple[int, ...]] = set()

    best = None  # (score, L, T, U, H, V)
    no_new = 0
    no_improve = 0

    cnf_path = os.path.join(out_dir, "working.cnf")

    for call_id in range(1, MAX_SOLVER_CALLS + 1):
        # 写 CNF（base + blocking）
        all_clauses = base_clauses + blocking_clauses
        write_dimacs(nvars, all_clauses, cnf_path)

        model = run_cmsat(cnf_path)
        if model is None:
            print(f"[{call_id}] UNSAT，解空间已被穷尽（在当前 blocking 下）。")
            break

        H, V = decode_model_to_HV(model, m, n, k)

        # 不管接不接受，都先 block 这个具体解（避免 solver 重复给同一个）
        blocking_clauses.append(make_blocking_clause_for_HV(H, V, m, n, k))

        # canonical 去同构（颜色置换 + 平移 + 翻转）
        key = canonical_key_under_aut_color(H, V, m, n)
        if key in seen_keys:
            no_new += 1
            no_improve += 1
        else:
            seen_keys.add(key)
            no_new = 0

            score, L, T, U = score_edge_coloring(H, V, k)
            if best is None or score > best[0]:
                best = (score, L, T, U, H, V)
                no_improve = 0
                print(f"[{call_id}] NEW BEST  score={score:.4f}  L={L:.3f}  T={T:.3f}  U={U:.1f}  unique={len(seen_keys)}")
            else:
                no_improve += 1

        # 停止条件（搜索型）
        if no_new >= MAX_NO_NEW_CLASS:
            print(f"[{call_id}] Stop: 连续 {MAX_NO_NEW_CLASS} 次没有新等价类。unique={len(seen_keys)}")
            break
        if no_improve >= MAX_NO_BEST_IMPROVE:
            print(f"[{call_id}] Stop: 连续 {MAX_NO_BEST_IMPROVE} 次 best 不提升。unique={len(seen_keys)}")
            break

        if call_id % 200 == 0:
            print(f"[{call_id}] progress: unique={len(seen_keys)}  blocking={len(blocking_clauses)}  no_new={no_new}  no_improve={no_improve}")

    # 输出最优方案
    if best is None:
        print("没有找到任何可行解（SAT 直接 UNSAT）。")
        return None

    score, L, T, U, Hbest, Vbest = best

    # 写 best_solution.txt
    sol_path = os.path.join(out_dir, "best_solution.txt")
    with open(sol_path, "w", encoding="utf-8") as f:
        f.write(f"# Best score={score:.6f}, LayerReg={L:.6f}, TransSym={T:.6f}, Penalty={U:.3f}\n")
        f.write("# H(i,j): (i,j)-(i,j+1)\n")
        for i in range(m):
            f.write(" ".join(str(Hbest[i][j] + 1) for j in range(n)) + "\n")
        f.write("\n# V(i,j): (i,j)-(i+1,j)\n")
        for i in range(m):
            f.write(" ".join(str(Vbest[i][j] + 1) for j in range(n)) + "\n")

    # 画 best 图
    img_path = os.path.join(out_dir, "best_B_edge_coloring_network.png")
    draw_graph_edge_coloring(Hbest, Vbest, m, n, k, img_path)

    print("\n==== DONE ====")
    print("Results dir:", out_dir)
    print("Best solution:", sol_path)
    print("Best figure:", img_path)

    return out_dir, best


if __name__ == "__main__":
    # 你也可以直接在这里手动跑
    # 示例：C3 x C7, k=5
    search_best(3, 7, 5)

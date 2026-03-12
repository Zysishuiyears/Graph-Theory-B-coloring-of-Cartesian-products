import os
import subprocess
from datetime import datetime
import matplotlib.pyplot as plt
import networkx as nx

########################################
#          Solver 配置
########################################
USE_SOLVER = "cmsat"
CMSAT_PATH = r"E:\cryptominisat5-win\cryptominisat5.exe"

BASE_CNF = os.path.join("out.cnf")
WORKING_CNF = os.path.join("working.cnf")
SYM_BREAK_COLOR_PERM = True

########################################
#     边编号 & 变量编号（颜色 + B-边）
########################################
# 边：
#  H(i,j): (i,j) -- (i,j+1)
#  V(i,j): (i,j) -- (i+1,j)
#  i in [0,m), j in [0,n)
#
# 边总数 E = 2*m*n
#
# SAT 变量：
#  color(e,c): 边 e 取颜色 c
#  rep(e,c):   e 被选为颜色 c 的 B-边
#
# 编号：
#  id(color(e,c)) = e*k + c + 1           ∈ [1, E*k]
#  id(rep(e,c))   = E*k + e*k + c + 1     ∈ [E*k+1, 2*E*k]


def edge_id_H(i, j, m, n):
    return i * n + j


def edge_id_V(i, j, m, n):
    return m * n + i * n + j


def num_edges(m, n):
    return 2 * m * n


def var_color(e, c, k):
    return e * k + c + 1


def var_rep(e, c, E, k):
    return E * k + e * k + c + 1


########################################
#   顶点 incident 边 & 边邻接关系 (line graph)
########################################
def incident_edges_at_vertex(i, j, m, n):
    """返回顶点 (i,j) incident 的边编号列表"""
    ids = []
    # 水平 (i,j)-(i,j+1) 和 (i,j-1)-(i,j)
    ids.append(edge_id_H(i, j, m, n))
    ids.append(edge_id_H(i, (j - 1) % n, m, n))
    # 垂直 (i,j)-(i+1,j) 和 (i-1,j)-(i,j)
    ids.append(edge_id_V(i, j, m, n))
    ids.append(edge_id_V((i - 1) % m, j, m, n))
    return ids


def build_edge_neighbors(m, n):
    """
    line graph 中的邻接：两条边共享一个顶点即为邻接。
    neighbors[e] = list of 边 f
    """
    E = num_edges(m, n)
    neighbors = [set() for _ in range(E)]

    for i in range(m):
        for j in range(n):
            inc = incident_edges_at_vertex(i, j, m, n)
            for a in range(len(inc)):
                for b in range(a + 1, len(inc)):
                    e1, e2 = inc[a], inc[b]
                    neighbors[e1].add(e2)
                    neighbors[e2].add(e1)

    return [list(s) for s in neighbors]

def enumerate_C4_edge_sets(m, n):
    """
    返回所有 4-圈对应的边编号列表列表。
    每个 4-圈给出一个长度为 4 的列表 [e1,e2,e3,e4]。
    这里按 C_m □ C_n 的“方格”来枚举：
      顶点: (i,j) -- (i,j+1) -- (i+1,j+1) -- (i+1,j) -- (i,j)
      边:   V(i,j), H(i+1,j), V(i,(j+1)), H(i,j)
    全部下标都 mod m / mod n 处理，所以是 torus。
    """
    C4s = []
    for i in range(m):
        for j in range(n):
            # 四条边分别是：
            e1 = edge_id_V(i, j, m, n)            # (i,j)-(i+1,j)
            e2 = edge_id_H((i+1) % m, j, m, n)    # (i+1,j)-(i+1,j+1)
            e3 = edge_id_V(i, (j+1) % n, m, n)    # (i,j+1)-(i+1,j+1)
            e4 = edge_id_H(i, j, m, n)            # (i,j)-(i,j+1)
            C4s.append([e1, e2, e3, e4])
    return C4s


########################################
#        写入带 B-coloring 约束的 CNF
########################################
def write_base_cnf(m, n, k):
    E = num_edges(m, n)
    neighbors = build_edge_neighbors(m, n)
    clauses = []

    # 1. 每条边至少一种颜色
    for e in range(E):
        clauses.append([var_color(e, c, k) for c in range(k)])

    # 2. 每条边至多一种颜色
    for e in range(E):
        for c1 in range(k):
            for c2 in range(c1 + 1, k):
                clauses.append([
                    -var_color(e, c1, k),
                    -var_color(e, c2, k)
                ])
    # ------- 3. proper edge-coloring：同一顶点 incident 的边不能同色 -------
    for i in range(m):
        for j in range(n):
            inc = incident_edges_at_vertex(i, j, m, n)
            for a in range(len(inc)):
                for b in range(a + 1, len(inc)):
                    e1, e2 = inc[a], inc[b]
                    for c in range(k):
                        clauses.append([
                            -var_color(e1, c, k),
                            -var_color(e2, c, k)
                        ])

    # ------- 3.5  每个 4-圈同一颜色最多出现一次（对应你旧代码里的 C4s 约束） -------
    C4s = enumerate_C4_edge_sets(m, n)
    for cyc in C4s:            # cyc 是这 4 条边的编号
        for c in range(k):     # 颜色 0..k-1
            # 相当于 atmost1_pairwise([var_color(e,c,k) for e in cyc])
            for i1 in range(4):
                for i2 in range(i1 + 1, 4):
                    e1, e2 = cyc[i1], cyc[i2]
                    clauses.append([
                        -var_color(e1, c, k),
                        -var_color(e2, c, k)
                    ])

        # ------- 3.8 颜色置换对称破坏（推荐开启） -------
    if SYM_BREAK_COLOR_PERM:
        # 基准顶点 v0 = (0,0) 的四条 incident 边（固定顺序）
        # 右：H(0,0)           (0,0)-(0,1)
        # 左：H(0,n-1)         (0,n-1)-(0,0)
        # 下：V(0,0)           (0,0)-(1,0)
        # 上：V(m-1,0)         (m-1,0)-(0,0)
        e_right = edge_id_H(0, 0, m, n)
        e_left  = edge_id_H(0, (n - 1) % n, m, n)
        e_down  = edge_id_V(0, 0, m, n)
        e_up    = edge_id_V((m - 1) % m, 0, m, n)

        ordered_edges = [e_right, e_left, e_down, e_up]

        # 给这 4 条边分别固定为颜色 1,2,3,4（代码里是 0,1,2,3）
        # 若 k < 4，就只固定前 k 条
        for c0, e in enumerate(ordered_edges[:min(k, 4)]):
            clauses.append([var_color(e, c0, k)])


    # 4. B-coloring 支配约束
    # rep(e,c): e 是颜色 c 的 B-边
    # (a) 每种颜色至少有一条 B-边:  OR_e rep(e,c)
    for c in range(k):
        clauses.append([var_rep(e, c, E, k) for e in range(E)])

    # (b) rep(e,c) -> color(e,c)
    for e in range(E):
        for c in range(k):
            clauses.append([
                -var_rep(e, c, E, k),
                var_color(e, c, k)
            ])

    # (c) rep(e,c) -> 对每个 d!=c，都存在邻边 f 颜色为 d
    #     rep(e,c) -> OR_{f∈N(e)} color(f,d)
    for e in range(E):
        N_e = neighbors[e]
        if not N_e:
            continue
        for c in range(k):
            for d in range(k):
                if d == c:
                    continue
                clause = [-var_rep(e, c, E, k)]
                for f in N_e:
                    clause.append(var_color(f, d, k))
                clauses.append(clause)

    num_vars = 2 * E * k
    with open(BASE_CNF, "w") as f:
        f.write(f"p cnf {num_vars} {len(clauses)}\n")
        for cls in clauses:
            f.write(" ".join(map(str, cls)) + " 0\n")

    # 初始化 working.cnf
    os.system(f"copy {BASE_CNF} {WORKING_CNF} >nul")


########################################
#             调用 solver
########################################
def run_cmsat(path):
    cmd = [CMSAT_PATH, path]
    res = subprocess.run(cmd, capture_output=True, text=True)

    if "UNSAT" in res.stdout:
        return None

    model = []
    for line in res.stdout.splitlines():
        if line.startswith("v "):
            for x in line.split()[1:]:
                if x != "0":
                    model.append(int(x))
    return model


########################################
#          解码模型 -> 边颜色 H,V
########################################
def decode_model_to_edge_colors(model, m, n, k):
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

    H = [[None] * n for _ in range(m)]
    V = [[None] * n for _ in range(m)]
    for i in range(m):
        for j in range(n):
            eH = edge_id_H(i, j, m, n)
            eV = edge_id_V(i, j, m, n)
            H[i][j] = edge_color[eH]
            V[i][j] = edge_color[eV]
    return H, V


########################################
#      画 “图二风格” 的 B-边染色网络图
########################################
def draw_graph_edge_coloring(H, V, m, n, k, path):
    """
    顶点：白色圆点 + 黑边框，标签 (i,j)
    边：依据 H,V 颜色上色，并在边的中点标注颜色数字（1..k）
    右边：颜色 1..k 图注
    """
    G = nx.Graph()
    for i in range(m):
        for j in range(n):
            G.add_node((i, j))

    # 添加水平 / 垂直边
    for i in range(m):
        for j in range(n):
            G.add_edge((i, j), (i, (j+1) % n), etype="H", i=i, j=j)
            G.add_edge((i, j), ((i+1) % m, j), etype="V", i=i, j=j)

    # 顶点位置（网格）
    # pos = {(i, j): (j, -i) for i in range(m) for j in range(n)}
    pos = nx.circular_layout(G)
    cmap = plt.get_cmap("tab10")

    # 收集边颜色
    edges = []
    edge_colors = []
    edge_labels = {}  # <---- 我们要添加的字典，存储边上显示的文本

    for u, v, data in G.edges(data=True):
        etype = data["etype"]
        i = data["i"]
        j = data["j"]

        if etype == "H":
            c = H[i][j]
        else:
            c = V[i][j]

        edges.append((u, v))
        edge_colors.append(cmap(c % 10))

        # ⭐⭐ 在边上加入颜色数字（1..k）
        edge_labels[(u, v)] = str(c + 1)

    # 绘图
    plt.figure(figsize=(7, 7))

    # 画彩色边
    nx.draw_networkx_edges(
        G, pos,
        edgelist=edges,
        edge_color=edge_colors,
        width=2.2
    )

    # 画顶点
    nx.draw_networkx_nodes(
        G, pos,
        node_color="white",
        edgecolors="black",
        linewidths=0.8,
        node_size=500
    )

    # 顶点标签
    nx.draw_networkx_labels(
        G, pos,
        labels={(i, j): f"({i},{j})" for i in range(m) for j in range(n)},
        font_size=8
    )

    # ⭐⭐⭐ 最重要：画边的文字标签（颜色编号）
    nx.draw_networkx_edge_labels(
        G, pos,
        edge_labels=edge_labels,
        font_size=8,
        font_color="black",
        label_pos=0.5   # 边的中点
    )

    # 图例
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(color=cmap(c % 10), label=f"color {c+1}")
        for c in range(k)
    ]
    plt.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(1.05, 1.0)
    )

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()

########################################
#       评分：周期性 + 颜色均匀性（边）
########################################
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
                    if H[i][j] != H[(i + a) % m][(j + b) % n]:
                        ok = False
                        break
                    if V[i][j] != V[(i + a) % m][(j + b) % n]:
                        ok = False
                        break
            if ok:
                cnt += 1
    return cnt


def uniform_penalty_edges(H, V, k):
    m, n = len(H), len(H[0])
    counts = [0] * k
    for i in range(m):
        for j in range(n):
            counts[H[i][j]] += 1
            counts[V[i][j]] += 1
    ideal = (2 * m * n) / k
    return sum((c - ideal) ** 2 for c in counts)

def horizontal_layer_score(H, k):
    """
    水平层规律度：在所有步长 a 中选一个，使
    （在该 a 下，各层分别选最优平移 b 的总命中比例）最大。
    返回值在 [0,1] 之间。
    """
    m = len(H)
    n = len(H[0])
    if m == 0 or n == 0:
        return 0.0

    best_over_a = 0.0
    for a in range(1, k):  # a=0 对应常数序列，一般没意义，跳过
        total_hits = 0
        for i in range(m):
            row = H[i]
            best_row_hits = 0
            for b in range(k):
                hits = 0
                for j in range(n):
                    if row[j] == (a * j + b) % k:
                        hits += 1
                if hits > best_row_hits:
                    best_row_hits = hits
            total_hits += best_row_hits
        score_a = total_hits / (m * n)
        if score_a > best_over_a:
            best_over_a = score_a
    return best_over_a


def vertical_layer_score(V, k):
    """
    竖直层规律度：类似于 horizontal_layer_score，只是换成纵向的层。
    """
    m = len(V)
    n = len(V[0])
    if m == 0 or n == 0:
        return 0.0

    best_over_a = 0.0
    for a in range(1, k):
        total_hits = 0
        for j in range(n):
            # 这一列是一条“竖直层”
            best_col_hits = 0
            for b in range(k):
                hits = 0
                for i in range(m):
                    if V[i][j] == (a * i + b) % k:
                        hits += 1
                if hits > best_col_hits:
                    best_col_hits = hits
            total_hits += best_col_hits
        score_a = total_hits / (m * n)
        if score_a > best_over_a:
            best_over_a = score_a
    return best_over_a


def layer_regular_score(H, V, k):
    """
    综合水平 + 竖直的层规律度，简单取平均。
    范围大致在 [0,1]。
    """
    h = horizontal_layer_score(H, k)
    v = vertical_layer_score(V, k)
    return (h + v) / 2.0



# def score_edge_coloring(H, V, k, lam=0.01):
#     P = count_edge_periods(H, V)
#     U = uniform_penalty_edges(H, V, k)
#     return P - lam * U, P, U
def score_edge_coloring(H, V, k, lam=0.01):
    """
    新的评分：
        score = layer_regular_score(H,V,k) - lam * uniform_penalty_edges(H,V,k)

    返回：
        score: 综合评分（越大越好）
        L    : 层规律度（0~1，越大越规则）
        U    : 颜色均匀性惩罚（越小越好）
    """
    L = layer_regular_score(H, V, k)
    U = uniform_penalty_edges(H, V, k)
    score = L - lam * U
    return score, L, U


########################################
#     blocking clause（只针对颜色变量）
########################################
def shift_HV(H, V, a, b):
    """
    将边染色 (H,V) 做平移 (a,b)，得到新的 (Hs,Vs)。
    平移含义：把原来坐标 (i-a, j-b) 的颜色搬到 (i,j)。
    """
    m, n = len(H), len(H[0])
    Hs = [[None]*n for _ in range(m)]
    Vs = [[None]*n for _ in range(m)]
    for i in range(m):
        for j in range(n):
            Hs[i][j] = H[(i - a) % m][(j - b) % n]
            Vs[i][j] = V[(i - a) % m][(j - b) % n]
    return Hs, Vs


def add_translation_orbit_blocking(H, V, m, n, k):
    """
    把当前解的所有平移副本都 block 掉（包括 (0,0) 本身也行，反正也要 block）。
    为避免重复（解有周期时平移副本可能相同），用集合去重。
    """
    seen = set()

    for a in range(m):
        for b in range(n):
            Hs, Vs = shift_HV(H, V, a, b)

            # 用 tuple 哈希去重：把 H、V 拉平成一条序列
            key = (tuple(x for row in Hs for x in row),
                   tuple(x for row in Vs for x in row))
            if key in seen:
                continue
            seen.add(key)

            clause = []

            # block 水平边
            for i in range(m):
                for j in range(n):
                    eH = edge_id_H(i, j, m, n)
                    clause.append(-var_color(eH, Hs[i][j], k))

            # block 垂直边
            for i in range(m):
                for j in range(n):
                    eV = edge_id_V(i, j, m, n)
                    clause.append(-var_color(eV, Vs[i][j], k))

            with open(WORKING_CNF, "a") as f:
                f.write(" ".join(map(str, clause)) + " 0\n")

def add_blocking_clause(H, V, m, n, k):
    clause = []

    for i in range(m):
        for j in range(n):
            eH = edge_id_H(i, j, m, n)
            cH = H[i][j]
            clause.append(-var_color(eH, cH, k))

    for i in range(m):
        for j in range(n):
            eV = edge_id_V(i, j, m, n)
            cV = V[i][j]
            clause.append(-var_color(eV, cV, k))

    with open(WORKING_CNF, "a") as f:
        f.write(" ".join(map(str, clause)) + " 0\n")


########################################
#              主搜索函数
########################################
def search_best(m, n, k=5, max_solutions=60):

    write_base_cnf(m, n, k)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = f"results/{stamp}_C{m}xC{n}_k{k}"
    os.makedirs(out_dir, exist_ok=True)

    best_H = None
    best_V = None
    best_info = None

    for s in range(1, max_solutions + 1):
        model = run_cmsat(WORKING_CNF)
        if model is None:
            print("UNSAT，停止枚举。")
            break

        H, V = decode_model_to_edge_colors(model, m, n, k)
        # score, P, U = score_edge_coloring(H, V, k)

        # print(f"[解 {s}] Score={score:.3f}, Periods={P}, Penalty={U:.1f}")   

        score, L, U = score_edge_coloring(H, V, k)
        print(f"[解 {s}] Score={score:.3f}, LayerReg={L:.3f}, Penalty={U:.1f}")


        # 保存方案文本：先 H 后 V
        sol_path = os.path.join(out_dir, f"solution_{s}.txt")
        with open(sol_path, "w") as f:
            f.write("# 水平边 H(i,j): (i,j)-(i,j+1)\n")
            for i in range(m):
                f.write(" ".join(str(H[i][j] + 1) for j in range(n)) + "\n")
            f.write("\n# 垂直边 V(i,j): (i,j)-(i+1,j)\n")
            for i in range(m):
                f.write(" ".join(str(V[i][j] + 1) for j in range(n)) + "\n")

        if best_info is None or score > best_info[0]:
            best_H, best_V = H, V
            # best_info = (score, P, U, sol_path)
            best_info = (score, L, U, sol_path)

        # add_blocking_clause(H, V, m, n, k)
        add_translation_orbit_blocking(H, V, m, n, k)


    if best_H is not None:
        img_path = os.path.join(out_dir, "best_B_edge_coloring_network.png")
        draw_graph_edge_coloring(best_H, best_V, m, n, k, img_path)
        print("最优 B-边染色网络图已输出到：", img_path)

    return (best_H, best_V), best_info


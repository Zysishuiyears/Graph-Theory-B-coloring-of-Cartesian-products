import itertools
import networkx as nx
from pulp import LpProblem, LpVariable, lpSum, LpMaximize, LpBinary, PULP_CBC_CMD

def cartesian_cycle_graph(m, n):
    """生成 G = C_m □ C_n"""
    G = nx.Graph()
    for i in range(m):
        for j in range(n):
            G.add_node((i, j))
            G.add_edge((i, j), ((i + 1) % m, j))   # 水平环
            G.add_edge((i, j), (i, (j + 1) % n))   # 垂直环
    return G

import networkx as nx

def expand_line_graph(G: nx.Graph) -> nx.Graph:
    """
    构造扩张线图 L^+(G)，假设 G = C_m □ C_n，
    顶点是 G 的边，两类邻接：
      1) 在 G 中共端点的两条边
      2) 在同一个 4-圈中成对边的两条边
    """

    # --- 1. 先把 G 的边做一个规范化列表（端点排序） ---
    edges = []
    for u, v in G.edges():
        if u > v:
            u, v = v, u
        edges.append((u, v))

    mE = len(edges)
    edge_index = {e: idx for idx, e in enumerate(edges)}

    L = nx.Graph()
    L.add_nodes_from(range(mE))   # 0..|E|-1 对应每一条边

    # --- 2. 普通线图邻接：共端点的边在 L 里相邻 ---
    for i in range(mE):
        u1, v1 = edges[i]
        for j in range(i + 1, mE):
            u2, v2 = edges[j]
            # 四种共端点情形
            if u1 == u2 or u1 == v2 or v1 == u2 or v1 == v2:
                L.add_edge(i, j)

    # --- 3. 扩张部分：每个 4-圈里的对边在 L 里相邻 ---
    # 利用 C_m □ C_n 的坐标来枚举每个 4-圈
    rows = sorted({i for i, _ in G.nodes()})
    cols = sorted({j for _, j in G.nodes()})
    m = len(rows)
    n = len(cols)

    def norm_edge(a, b):
        """端点排序后作为 key"""
        return (a, b) if a < b else (b, a)

    for i in rows:
        for j in cols:
            # 当前小方块的 4 个顶点
            v00 = (i, j)
            v10 = ((i + 1) % m, j)
            v11 = ((i + 1) % m, (j + 1) % n)
            v01 = (i, (j + 1) % n)

            # 这个 4-圈的四条边
            e1 = norm_edge(v00, v10)  # 底边
            e2 = norm_edge(v10, v11)  # 右边
            e3 = norm_edge(v11, v01)  # 上边
            e4 = norm_edge(v01, v00)  # 左边

            id1 = edge_index[e1]
            id2 = edge_index[e2]
            id3 = edge_index[e3]
            id4 = edge_index[e4]

            # 对边连边： (e1, e3) 和 (e2, e4)
            L.add_edge(id1, id3)
            L.add_edge(id2, id4)

    return L


def max_independent_set_exact(L):
    """ILP 精确求解 α(L)"""
    prob = LpProblem("MaxIndependentSet", LpMaximize)
    x = {v: LpVariable(f"x_{v}", cat=LpBinary) for v in L.nodes()}
    # 目标
    prob += lpSum(x[v] for v in L.nodes())
    # 相邻约束
    for u, v in L.edges():
        prob += x[u] + x[v] <= 1
    # 求解
    solver = PULP_CBC_CMD(msg=False)
    prob.solve(solver)
    value = sum(var.value() for var in x.values())
    return int(value)

if __name__ == "__main__":
    m, n = 3, 7   # 你可以改成 3x5, 5x5 等
    G = cartesian_cycle_graph(m, n)
    Lplus = expand_line_graph(G)
    print(f"|V(L+)| = {Lplus.number_of_nodes()}, |E(L+)| = {Lplus.number_of_edges()}")
    alpha = max_independent_set_exact(Lplus)
    print(f"α(L^+(C_{m}□C_{n})) = {alpha}")



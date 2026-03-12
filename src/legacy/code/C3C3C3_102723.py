#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
C3xC3xC3 (C_3^3) 7-color B-coloring via SAT (CryptoMiniSat5)
"""

import argparse
import itertools
import os
import subprocess
import sys
from collections import defaultdict
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  # 触发 3D 投影
import pandas as pd
import networkx as nx

# =========================
#  Utilities
# =========================

def mod3(x): return x % 3

def vert_id(i, j, k):
    return (i % 3) * 9 + (j % 3) * 3 + (k % 3)

def vert_tuple(v):
    i = v // 9
    r = v % 9
    j = r // 3
    k = r % 3
    return (i, j, k)

# =========================
#  Build C3^3 graph  (FIXED)
# =========================

def build_graph_c3_c3_c3():
    """
    Return:
      V: [0..26]
      edges: list of (u,v) with u<v, |E| should be 81
      inc_edges: dict v -> [edge ids]
      axes_edge_ids: {'X': [...], 'Y': [...], 'Z': [...]}
    """
    V = list(range(27))
    edges = []
    inc_edges = defaultdict(list)
    axes_edge_ids = {'X': [], 'Y': [], 'Z': []}

    seen = set()  # store (min(u,v), max(u,v)) to avoid duplicates

    def add_edge(u, v, axis_key):
        a, b = (u, v) if u < v else (v, u)
        key = (a, b)
        if a == b or key in seen:
            return
        seen.add(key)
        eid = len(edges)
        edges.append(key)
        inc_edges[a].append(eid)
        inc_edges[b].append(eid)
        axes_edge_ids[axis_key].append(eid)

    # enumerate all +X, +Y, +Z from every vertex; rely on 'seen' to dedup
    for i in range(3):
        for j in range(3):
            for k in range(3):
                u = vert_id(i, j, k)
                add_edge(u, vert_id(mod3(i+1), j, k), 'X')
                add_edge(u, vert_id(i, mod3(j+1), k), 'Y')
                add_edge(u, vert_id(i, j, mod3(k+1)), 'Z')

    assert len(V) == 27
    assert len(edges) == 81, f"Expect 81 edges, got {len(edges)}"
    return V, edges, inc_edges, axes_edge_ids

# =========================
#  Enumerate all 4-cycles
# =========================

def enumerate_all_4cycles(edges, inc_edges):
    edge_index = {}
    for eid, (u, v) in enumerate(edges):
        edge_index[(u, v)] = eid
        edge_index[(v, u)] = eid

    def edge_id(u, v): return edge_index[(u, v)]

    C4s = []

    # XY planes (fix k)
    for k in range(3):
        for i in range(3):
            for j in range(3):
                a = vert_id(i, j, k)
                b = vert_id(mod3(i+1), j, k)
                c = vert_id(mod3(i+1), mod3(j+1), k)
                d = vert_id(i, mod3(j+1), k)
                C4s.append([edge_id(a,b), edge_id(b,c), edge_id(c,d), edge_id(d,a)])

    # YZ planes (fix i)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                a = vert_id(i, j, k)
                b = vert_id(i, mod3(j+1), k)
                c = vert_id(i, mod3(j+1), mod3(k+1))
                d = vert_id(i, j, mod3(k+1))
                C4s.append([edge_id(a,b), edge_id(b,c), edge_id(c,d), edge_id(d,a)])

    # ZX planes (fix j)
    for j in range(3):
        for i in range(3):
            for k in range(3):
                a = vert_id(i, j, k)
                b = vert_id(mod3(i+1), j, k)
                c = vert_id(mod3(i+1), j, mod3(k+1))
                d = vert_id(i, j, mod3(k+1))
                C4s.append([edge_id(a,b), edge_id(b,c), edge_id(c,d), edge_id(d,a)])

    # unique (defensive; should already be unique)
    seen = set()
    uniq = []
    for cyc in C4s:
        t = tuple(cyc)
        if t not in seen:
            seen.add(t)
            uniq.append(cyc)
    assert len(uniq) == 81, f"Expect 81 C4 cycles, got {len(uniq)}"
    return uniq

# =========================
#  CNF builder
# =========================

class CNFBuilder:
    def __init__(self, nvars):
        self.nvars = nvars
        self.clauses = []

    def add_clause(self, lits):
        assert all(1 <= abs(x) <= self.nvars for x in lits)
        self.clauses.append(list(lits))

    def add_atleast_one(self, lits):
        self.add_clause(lits)

    def add_atmost_one_pairwise(self, lits):
        for i in range(len(lits)):
            for j in range(i+1, len(lits)):
                self.add_clause([-lits[i], -lits[j]])

    def add_exactly_one_pairwise(self, lits):
        self.add_atleast_one(lits)
        self.add_atmost_one_pairwise(lits)

    def to_dimacs(self):
        out = [f"p cnf {self.nvars} {len(self.clauses)}"]
        out += [" ".join(map(str, c)) + " 0" for c in self.clauses]
        return "\n".join(out)

def var_id(eid, color, K): return eid * K + color

# =========================
#  Build model
# =========================

def build_bcolor_cnf(K=7, symbreak=True):
    V, edges, inc_edges, axes_edge_ids = build_graph_c3_c3_c3()
    C4s = enumerate_all_4cycles(edges, inc_edges)

    n_edges = len(edges)
    nvars = n_edges * K
    cnf = CNFBuilder(nvars=nvars)

    # (A) each edge exactly one color
    for eid in range(n_edges):
        lits = [var_id(eid, c, K) for c in range(1, K+1)]
        cnf.add_exactly_one_pairwise(lits)

    # (B) proper edge-coloring at each vertex
    for v in V:
        Es = inc_edges[v]  # 6 incident edges
        for c in range(1, K+1):
            cnf.add_atmost_one_pairwise([var_id(e, c, K) for e in Es])

    # (C) every 4-cycle rainbow
    for cyc in C4s:
        for c in range(1, K+1):
            cnf.add_atmost_one_pairwise([var_id(e, c, K) for e in cyc])

    # (SB) symmetry breaking: fix 6 edges at v0 to colors 1..6
    if symbreak:
        v0 = vert_id(0, 0, 0)
        eids = inc_edges[v0]

        def edge_axis(eid):
            u, w = edges[eid]
            i1,j1,k1 = vert_tuple(u)
            i2,j2,k2 = vert_tuple(w)
            di,dj,dk = ((i2-i1)%3, (j2-j1)%3, (k2-k1)%3)
            if di in (1,2) and dj==0 and dk==0: return 'X'
            if dj in (1,2) and di==0 and dk==0: return 'Y'
            if dk in (1,2) and di==0 and dj==0: return 'Z'
            raise RuntimeError("Bad edge axis")

        groups = {'X':[], 'Y':[], 'Z':[]}
        for e in eids:
            groups[edge_axis(e)].append(e)

        ordered = []
        for ax in ['X','Y','Z']:
            tmp = []
            for e in groups[ax]:
                u, w = edges[e]
                nb = w if u == v0 else u
                tmp.append((nb, e))
            tmp.sort()
            ordered += [e for _,e in tmp]

        for fixed_color, eid in enumerate(ordered, start=1):
            if fixed_color > min(6, K): break
            cnf.add_clause([var_id(eid, fixed_color, K)])

    return cnf, (V, edges, inc_edges, C4s), K

# =========================
#  Solvers
# =========================

def solve_with_pycryptosat(cnf_builder):
    try:
        from pycryptosat import Solver
    except ImportError:
        print("[!] pip install pycryptosat", file=sys.stderr)
        return None, None
    s = Solver()
    for c in cnf_builder.clauses:
        s.add_clause(c)
    sat, model = s.solve()
    if not sat: return False, None
    assign = {}
    if isinstance(model, dict):
        assign = model
    else:
        for v in range(1, cnf_builder.nvars+1):
            assign[v] = bool(model[v])
    return True, assign

def write_dimacs(cnf_builder, path):
    if os.path.dirname(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(cnf_builder.to_dimacs())

def parse_cmsat_output(text, nvars):
    sat = None
    vals = []
    for line in text.splitlines():
        s = line.strip()
        if not s: continue
        if s.startswith('s '):
            if 'UNSAT' in s: sat = False
            elif 'SAT' in s: sat = True
        elif s.startswith('v ') or s.startswith('V '):
            for p in s.split()[1:]:
                if p == '0': continue
                try: vals.append(int(p))
                except: pass
    if sat is None:
        if "UNSAT" in text: sat = False
        elif "SAT" in text: sat = True
        else: return None, None
    if not sat: return False, None
    assign = {v: False for v in range(1, nvars+1)}
    for lit in vals:
        if lit == 0: continue
        v = abs(lit)
        if v <= nvars: assign[v] = (lit > 0)
    return True, assign

def solve_with_cmsat_dimacs(cnf_builder, cmsat_path="cryptominisat5", cnf_path="out.cnf"):
    write_dimacs(cnf_builder, cnf_path)
    try:
        proc = subprocess.run([cmsat_path, cnf_path],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except FileNotFoundError:
        print(f"[!] Solver not found: {cmsat_path}", file=sys.stderr)
        return None, None
    out = proc.stdout + "\n" + proc.stderr
    return parse_cmsat_output(out, cnf_builder.nvars)

# =========================
#  Decode & verify
# =========================

def var_to_color_list(assignment, n_edges, K):
    colors = [-1]*n_edges
    for eid in range(n_edges):
        picks = []
        for c in range(1, K+1):
            if assignment.get(eid*K + c, False):
                picks.append(c)
        if len(picks) != 1:
            raise RuntimeError(f"Edge {eid} has {len(picks)} colors: {picks}")
        colors[eid] = picks[0]
    return colors

def verify_solution(colors, V, edges, inc_edges, C4s):
    # vertex proper
    for v in V:
        seen = set()
        for e in inc_edges[v]:
            c = colors[e]
            if c in seen:
                raise AssertionError(f"Vertex {v} has duplicate color {c}")
            seen.add(c)
    # rainbow C4
    for cyc in C4s:
        cols = [colors[e] for e in cyc]
        if len(set(cols)) != 4:
            raise AssertionError(f"C4 {cyc} not rainbow: {cols}")
        

def plot_flat_networkx_from_csv(edges_csv_path, layout="spring",
                                node_size=250, edge_width=2, K=7):
    """
    从 edges.csv 读取 (u,v,color)，用 NetworkX 平面图画出来。
    关键点：显式传 edgelist，并让 edge_color 与 edgelist 一一对应。
    """
    # 1) 读 CSV：要求有列 Source, Target, Color（或 color）
    import pandas as pd
    df = pd.read_csv(edges_csv_path)
    # 兼容不同列名
    ucol = "Source" if "Source" in df.columns else "source"
    vcol = "Target" if "Target" in df.columns else "target"
    ccol = "Color"  if "Color"  in df.columns else ("color" if "color" in df.columns else "col")

    # 2) 按 CSV 原顺序构造 edgelist 和 colorlist（这是我们的“真顺序”）
    edgelist = list(zip(df[ucol].tolist(), df[vcol].tolist()))
    color_ids = df[ccol].astype(int).tolist()

    # 3) 建图（顺序不重要），但画的时候要用 edgelist=edgelist
    G = nx.Graph()
    G.add_edges_from(edgelist)

    # 4) 选择布局
    if layout == "spring":
        pos = nx.spring_layout(G, seed=42)
    elif layout == "kamada_kawai":
        pos = nx.kamada_kawai_layout(G)
    elif layout == "circular":
        pos = nx.circular_layout(G)
    elif layout == "shell":
        pos = nx.shell_layout(G)
    else:
        pos = nx.spring_layout(G, seed=42)

    # 5) 颜色表（1..7）
    palette = ["#E41A1C", "#377EB8", "#4DAF4A",
               "#984EA3", "#FF7F00", "#FFFF33", "#A65628"]
    edge_colors = [palette[(c-1) % len(palette)] for c in color_ids]

    # 6) 画图 —— 一定要把 edgelist 和 edge_color 同时传入
    plt.figure(figsize=(8, 7))
    nx.draw_networkx_nodes(G, pos, node_color="white", edgecolors="black",
                           node_size=node_size, linewidths=0.8)
    nx.draw_networkx_edges(G, pos,
                           edgelist=edgelist,          # ✅ 固定顺序
                           edge_color=edge_colors,     # ✅ 与之对齐
                           width=edge_width)
    nx.draw_networkx_labels(G, pos, font_size=9, font_color="black")

    # 图例
    import matplotlib.lines as mlines
    handles = [mlines.Line2D([], [], color=palette[i], lw=3, label=f"color {i+1}")
               for i in range(K)]
    plt.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0))

    plt.axis("off")
    plt.tight_layout()
    plt.title("C3 x C3 x C3  —  7-color B-coloring (flat view)")
    plt.show()


def plot_flat_networkx(edges, colors, K=7, layout="spring", node_size=250, edge_width=2):
    """
    用 NetworkX 平面方式画图（无 xyz 坐标）。
    edges : [(u,v), ...]
    colors : list[int]  每条边的颜色编号 (1..K)
    layout : "spring" (默认) / "circular" / "kamada_kawai" / "shell"
    """
    G = nx.Graph()
    for eid, (u, v) in enumerate(edges):
        G.add_edge(u, v, color=colors[eid])

    # 不同布局选项
    if layout == "spring":
        pos = nx.spring_layout(G, seed=42)   # 弹簧布局（自动拉伸）
    elif layout == "circular":
        pos = nx.circular_layout(G)
    elif layout == "kamada_kawai":
        pos = nx.kamada_kawai_layout(G)
    elif layout == "shell":
        pos = nx.shell_layout(G)
    else:
        pos = nx.spring_layout(G)

    # 颜色映射表（7 色）
    palette = ["#E41A1C", "#377EB8", "#4DAF4A",
               "#984EA3", "#FF7F00", "#FFFF33", "#A65628"]

    edge_colors = [palette[(c-1) % len(palette)] for c in colors]

    plt.figure(figsize=(8, 7))
    nx.draw_networkx_nodes(G, pos, node_color="white", edgecolors="black",
                           node_size=node_size, linewidths=0.8)
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_width)
    nx.draw_networkx_labels(G, pos, font_size=9, font_color="black")

    # 生成简单图例
    import matplotlib.lines as mlines
    handles = [mlines.Line2D([], [], color=palette[i], lw=3, label=f"color {i+1}")
               for i in range(K)]
    plt.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0))

    plt.axis("off")
    plt.tight_layout()
    plt.title(f"C3 x C3 x C3  —  {K}-color B-coloring (flat view)")
    plt.show()


def print_solution_summary(colors, edges, sample=24):
    print("\n=== Solution summary ===")
    print(f"Total edges: {len(edges)}")
    for eid in range(min(sample, len(edges))):
        print(f"  edge {eid}: {edges[eid]} -> color {colors[eid]}")
    hist = defaultdict(int)
    for c in colors: hist[c] += 1
    print("Color histogram:")
    for c in sorted(hist):
        print(f"  color {c}: {hist[c]} edges")


def export_for_gephi(V, edges, colors, out_prefix="C3x3x3_K7"):
    """导出节点表 nodes.csv 与边表 edges.csv，可直接导入 Gephi。"""
    # === 节点表 ===
    node_rows = []
    for v in V:
        i,j,k = vert_tuple(v)
        node_rows.append({"Id": v, "Label": f"({i},{j},{k})", "x": i, "y": j, "z": k})
    df_nodes = pd.DataFrame(node_rows)
    nodes_path = f"{out_prefix}_nodes.csv"
    df_nodes.to_csv(nodes_path, index=False)

    # === 边表 ===
    # 预设 7 种颜色（对应 tab10 调色板或你想用的七色方案）
    palette_hex = ["#E41A1C", "#377EB8", "#4DAF4A",
                   "#984EA3", "#FF7F00", "#FFFF33", "#A65628"]
    edge_rows = []
    for eid, (u, v) in enumerate(edges):
        c = colors[eid]
        edge_rows.append({
            "Id": eid,
            "Source": u,
            "Target": v,
            "Type": "Undirected",
            "Color": c,
            "ColorRGB": palette_hex[c-1]
        })
    df_edges = pd.DataFrame(edge_rows)
    edges_path = f"{out_prefix}_edges.csv"
    df_edges.to_csv(edges_path, index=False)

    print(f"[Gephi export] Wrote {nodes_path} and {edges_path}")


# =========================
#  Main
# =========================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=7)
    parser.add_argument("--solver", choices=["pycryptosat","cmsat"], default="pycryptosat")
    parser.add_argument("--cnf", type=str, default="out.cnf")
    parser.add_argument("--cmsat_path", type=str, default="cryptominisat5")
    parser.add_argument("--symbreak", action="store_true")
    parser.add_argument("--no-symbreak", dest="symbreak", action="store_false")
    parser.add_argument("--plot", action="store_true",
                    help="求解后绘制 3D 可视化")
    parser.add_argument("--annotate", action="store_true",
                    help="在每个顶点旁标注编号")
    parser.add_argument("--save_png", type=str, default=None,
                    help="若提供路径（如 C3x3x3_K7.png），保存可视化图片")

    parser.set_defaults(symbreak=True)
    args = parser.parse_args()

    cnf, meta, K = build_bcolor_cnf(K=args.k, symbreak=args.symbreak)
    V, edges, inc_edges, C4s = meta
    print(f"[Info] Vars={cnf.nvars}, Clauses={len(cnf.clauses)}, |V|=27, |E|={len(edges)}, #C4={len(C4s)}, symbreak={args.symbreak}")

    if args.solver == "pycryptosat":
        sat, assignment = solve_with_pycryptosat(cnf)
    else:
        sat, assignment = solve_with_cmsat_dimacs(cnf, cmsat_path=args.cmsat_path, cnf_path=args.cnf)

    if sat is None:
        print("[!] Solver failed to run.")
        sys.exit(2)
    if not sat:
        print(f"[Result] UNSAT for K={K}")
        sys.exit(1)

    print("[Result] SAT")
    colors = var_to_color_list(assignment, len(edges), K)
    verify_solution(colors, V, edges, inc_edges, C4s)
    print("[Check] Verified OK.")
    print_solution_summary(colors, edges)

    sol_path = f"solution_C3x3x3_K{K}.txt"
    with open(sol_path, "w") as f:
        for eid,(u,v) in enumerate(edges):
            f.write(f"{eid}\t{u}\t{v}\t{colors[eid]}\n")
    print(f"[Saved] {sol_path}")

    # === 导出给 Gephi ===
    prefix = f"C3x3x3_K{K}"
    export_for_gephi(V, edges, colors, out_prefix=prefix)

    # === 平面 NetworkX 可视化 ===
    # plot_flat_networkx(edges, colors, K=K, layout="spring")
    plot_flat_networkx_from_csv("C3x3x3_K7_edges.csv", layout="kamada_kawai")


if __name__ == "__main__":
    main()
    

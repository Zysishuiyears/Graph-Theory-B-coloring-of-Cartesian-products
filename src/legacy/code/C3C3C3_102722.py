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
        
def _vertex_positions_3d():
    """返回顶点坐标字典 pos[v] = (x,y,z)，使用 vert_tuple(v) 的自然坐标 (0,1,2)."""
    return {v: vert_tuple(v) for v in range(27)}

def _set_axes_equal(ax):
    """让 3D 轴比例相等（避免拉伸变形）。"""
    xs = [0, 2]; ys = [0, 2]; zs = [0, 2]
    ax.set_box_aspect([max(xs)-min(xs), max(ys)-min(ys), max(zs)-min(zs)])
    # 也可以手动设置范围：
    ax.set_xlim( -0.3, 2.3 )
    ax.set_ylim( -0.3, 2.3 )
    ax.set_zlim( -0.3, 2.3 )

def plot_c3x3x3_solution(V, edges, colors, annotate=False, save_path=None,
                         linewidth=2.0, node_size=40, show=True):
    """
    将求解得到的边染色可视化为 3D 图。
    - V, edges: 构图结果
    - colors:  len(edges) 的颜色编号列表（1..K）
    - annotate: 是否在点旁标 v 的编号
    - save_path: 若不为 None，则保存到该路径（如 'C3x3x3_K7.png'）
    """
    pos = _vertex_positions_3d()

    # 给 1..7 号颜色准备可区分的颜色表（tab10 的前 7 种）
    base = plt.cm.tab10.colors
    palette = [base[i % len(base)] for i in range(7)]  # idx 0..6 对应 color 1..7

    fig = plt.figure(figsize=(7.5, 7.0))
    ax = fig.add_subplot(111, projection='3d')

    # 画顶点
    xs = [pos[v][0] for v in V]
    ys = [pos[v][1] for v in V]
    zs = [pos[v][2] for v in V]
    ax.scatter(xs, ys, zs, s=node_size, depthshade=True, alpha=0.9)

    if annotate:
        for v in V:
            x, y, z = pos[v]
            ax.text(x+0.05, y+0.05, z+0.05, str(v), fontsize=8, alpha=0.8)

    # 画边（按颜色）
    for eid, (u, w) in enumerate(edges):
        c_idx = colors[eid] - 1  # 颜色编号 1..7 -> palette 下标 0..6
        col = palette[c_idx]
        x = [pos[u][0], pos[w][0]]
        y = [pos[u][1], pos[w][1]]
        z = [pos[u][2], pos[w][2]]
        ax.plot(x, y, z, linewidth=linewidth, color=col, alpha=0.95)

    # 设置坐标与观感
    ax.set_xlabel('i')
    ax.set_ylabel('j')
    ax.set_zlabel('k')
    _set_axes_equal(ax)
    ax.view_init(elev=20, azim=35)
    ax.set_title('C3 x C3 x C3  —  7-color B-coloring')

    # 手动做个简单的图例（颜色1..7）
    import matplotlib.lines as mlines
    handles = []
    for k in range(1, 8):
        handles.append(mlines.Line2D([], [], color=palette[k-1], lw=3, label=f'color {k}'))
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.02, 1.0))

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=220, bbox_inches='tight')
        print(f"[Saved figure] {save_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)


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


        # 可视化（按需）
    if args.plot:
        png = args.save_png if args.save_png else f"C3x3x3_K{K}.png"
        plot_c3x3x3_solution(V, edges, colors,
                             annotate=args.annotate,
                             save_path=png,
                             linewidth=2.0, node_size=40, show=True)


if __name__ == "__main__":
    main()
    


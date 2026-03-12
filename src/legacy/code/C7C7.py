# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 22:48:06 2025

@author: JZX
"""

# encode_c3xc7_sat.py
# 生成 C3 x C7 的 5-color B-coloring 的 SAT CNF（DIMACS），并包含解析器
# 用法:
#   python encode_c3xc7_sat.py c7 c7 5
# 结果:
#   生成文件 c7xc7_k5.cnf，可用 minisat/cminisat/Glucose 求解
#   若求解器给出可满足解，使用脚本解析出边->color 映射

import sys
from itertools import combinations
import numpy as np
import csv

def mk_vars(m, n, k):
    # assign each undirected edge an index and then variable x_{e,c}
    V = [(i,j) for i in range(m) for j in range(n)]
    edges = []
    edge_index = {}
    for i in range(m):
        for j in range(n):
            u=(i,j); v=((i+1)%m,j)
            e = tuple(sorted((u,v)))
            if e not in edge_index:
                edge_index[e] = len(edges)
                edges.append(e)
    for i in range(m):
        for j in range(n):
            u=(i,j); v=(i,(j+1)%n)
            e = tuple(sorted((u,v)))
            if e not in edge_index:
                edge_index[e] = len(edges)
                edges.append(e)
    # variable numbering 1..N
    var_of = {}
    cur = 1
    for ei,e in enumerate(edges):
        for c in range(1,k+1):
            var_of[(ei,c)] = cur
            cur += 1
    return V, edges, edge_index, var_of, cur-1




def exactly_one_clause(vars_list):
    # returns CNF clauses for exactly-one constraint on vars_list (list of ints)
    # encoding: at-least-one (OR) + pairwise at-most-one (pairwise)
    clauses = []
    clauses.append(" ".join(map(str, vars_list)) + " 0")
    for a,b in combinations(vars_list,2):
        clauses.append(f"-{a} -{b} 0")
    return clauses

def build_cnf(m,n,k, symmetry_break=True):
    V, edges, edge_index, var_of, nvars = mk_vars(m,n,k)
    clauses = []
    # each edge exactly one color
    for ei in range(len(edges)):
        vars_list = [var_of[(ei,c)] for c in range(1,k+1)]
        clauses += exactly_one_clause(vars_list)
    # adjacent-edge constraints: edges that share a vertex must have different colors
    # build incidence: for each vertex, collect incident edge indices
    incident = {}
    for v in V:
        incident[v]=[]
    for ei,e in enumerate(edges):
        u,vv = e
        incident[u].append(ei); incident[vv].append(ei)
    # for each vertex, any two incident edges cannot have same color
    for v in V:
        inc = incident[v]
        for a,b in combinations(inc,2):
            for c in range(1,k+1):
                clauses.append(f"-{var_of[(a,c)]} -{var_of[(b,c)]} 0")
    # 4-cycle all-different: for each 4-cycle, pairwise inequality among its four edges
    # enumerate 4-cycles inherent to Cartesian product
    for i in range(m):
        for j in range(n):
            a=(i,j); b=((i+1)%m,j); cnode=((i+1)%m,(j+1)%n); d=(i,(j+1)%n)
            e1 = edge_index[tuple(sorted((a,b)))]
            e2 = edge_index[tuple(sorted((b,cnode)))]
            e3 = edge_index[tuple(sorted((cnode,d)))]
            e4 = edge_index[tuple(sorted((d,a)))]
            quartet = [e1,e2,e3,e4]
            # require all four colors pairwise distinct -> for each pair (p,q) and color c: not both have color c
            for p,q in combinations(quartet,2):
                for c in range(1,k+1):
                    clauses.append(f"-{var_of[(p,c)]} -{var_of[(q,c)]} 0")
    # symmetry breaking: fix some edges' colors to reduce search
    if symmetry_break:
        # find the three edges of the C3 triangle at column 0 (horizontal edges):
        # edges correspond to (0,0)-(1,0), (1,0)-(2,0), (0,0)-(2,0) — note ordering
        # also fix vertical (0,0)-(0,1) to color k-1 (if k>=4) to reduce symmetry
        def find_edge_index(edges, a,b):
            e=tuple(sorted((a,b)))
            return edges.index(e)
        try:
            e01 = find_edge_index(edges,(0,0),(1,0))
            e12 = find_edge_index(edges,(1,0),(2,0))
            e02 = find_edge_index(edges,(0,0),(2,0))
            clauses.append(f"{var_of[(e01,1)]} 0")
            clauses.append(f"{var_of[(e12,2)]} 0")
            clauses.append(f"{var_of[(e02,3)]} 0")
            if k>=4:
                ev = find_edge_index(edges,(0,0),(0,1))
                clauses.append(f"{var_of[(ev,k-1)]} 0")
        except ValueError:
            pass
    return nvars, len(clauses), clauses, V, edges, var_of

def write_dimacs(filename, nvars, nclauses, clauses):
    with open(filename,"w") as f:
        f.write(f"p cnf {nvars} {nclauses}\n")
        for cl in clauses:
            f.write(cl + "\n")
            
def write_edges_csv(m, n, edge_index, filename):
    """
    导出每条边的有向规范：V(i,j): (i,j)->(i+1,j), H(i,j): (i,j)->(i,j+1)
    按照 mk_vars 的顺序（先竖直、后水平），用 edge_index 查到 cnf 的边索引 ei。
    CSV 列：ei,kind,i,j,i2,j2
    """
    with open(filename, "w", encoding="utf-8") as f:
        f.write("ei,kind,i,j,i2,j2\n")
        # 竖直边 (V): (i,j) -> (i+1 mod m, j)
        for i in range(m):
            for j in range(n):
                a = (i, j)
                b = ((i+1) % m, j)
                ei = edge_index[tuple(sorted((a, b)))]
                f.write(f"{ei},V,{i},{j},{(i+1)%m},{j}\n")
        # 水平边 (H): (i,j) -> (i, j+1 mod n)
        for i in range(m):
            for j in range(n):
                a = (i, j)
                b = (i, (j+1) % n)
                ei = edge_index[tuple(sorted((a, b)))]
                f.write(f"{ei},H,{i},{j},{i},{(j+1)%n}\n")


def decode_model_to_HV(m, n, k, out_lines, edges_csv, var_of):
    """
    读取 out.txt 的模型 + edges.csv 的映射，复原 H(i,j), V(i,j)。
    返回 H, V 两个 m×n 矩阵（int）。
    """
    # 收集正文字（true 变量）
    model_true = set()
    for line in out_lines:
        line = line.strip()
        if not line or line[0] in "cs":
            continue
        if line in ("SAT", "UNSAT"):
            continue
        for tok in line.split():
            v = int(tok)
            if v > 0:
                model_true.add(v)

    # 反查字典 var -> (ei,c)
    inv = {v: (ei, c) for (ei, c), v in var_of.items()}

    # 读 edges.csv
    rows = {}
    with open(edges_csv, newline="", encoding="utf-8") as f:
        rdr = csv.DictReader(f)
        for r in rdr:
            ei = int(r["ei"])
            rows[ei] = r  # kind, i, j, i2, j2

    # 构建 H, V
    H = np.zeros((m, n), dtype=int)
    V = np.zeros((m, n), dtype=int)

    # 遍历所有 true 变量，给对应 (i,j) 填色
    for v in model_true:
        if v not in inv:  # 不是颜色变量（极少见），跳过
            continue
        ei, c = inv[v]
        r = rows[ei]
        kind = r["kind"]
        i = int(r["i"])
        j = int(r["j"])
        if kind == "H":
            H[i, j] = c
        else:
            V[i, j] = c

    return H, V


def parse_minisat_output(output_lines, var_of, edges, k):
    # output_lines from minisat: SAT/UNSAT then model list
    # model list is line of integers; positive means var=true
    model = []
    for line in output_lines:
        line=line.strip()
        if line=="" or line.startswith("c"): continue
        if line=="SAT":
            continue
        if line=="UNSAT":
            return None
        # else model numbers
        parts=line.split()
        for p in parts:
            v=int(p)
            if v>0:
                model.append(v)
    # invert var_of
    inv = {v: (ei,c) for (ei,c),v in var_of.items()}
    edge_color = {}
    for v in model:
        if v in inv:
            ei,c = inv[v]
            edge_color[edges[ei]] = c
    return edge_color

if __name__ == "__main__":
    # call: python encode_c3xc7_sat.py
    m = 7
    n = 7
    k = 5
    nvars, nclauses, clauses, V, edges, var_of = build_cnf(m, n, k, symmetry_break=True)

    fn = f"c{m}xc{n}_k{k}.cnf"
    write_dimacs(fn, nvars, nclauses, clauses)
    print("Wrote DIMACS to", fn)
    print("Variables:", nvars, "Clauses:", nclauses)

    # === 新增：导出 edges.csv，描述 ei -> (kind,i,j,i2,j2) ===
    edges_csv = f"c{m}xc{n}_edges.csv"
    # 注意：build_cnf 返回了 edge_index；我们从中重建它
    # 这里简单复用 mk_vars 一次，或把 build_cnf 的 edge_index 也返回
    _, _, edge_index, _, _ = mk_vars(m, n, k)
    write_edges_csv(m, n, edge_index, edges_csv)
    print("Wrote edges mapping to", edges_csv)

    print("Now run a SAT solver, e.g.:")
    print(f"  minisat {fn} out.txt")
    print("Or with cryptominisat:")
    print(f"  cryptominisat5 {fn} > out.txt")
    print("After solver finishes, run this script again with 'parse' to decode the model.")

    # 可选：提供一个简单的解析入口
    if len(sys.argv) >= 2 and sys.argv[1] == "parse":
        out_file = "out.txt" if len(sys.argv) < 3 else sys.argv[2]
        with open(out_file, "r", encoding="utf-8") as f:
            out_lines = f.readlines()
        H, V = decode_model_to_HV(m, n, k, out_lines, edges_csv, var_of)

        # 打印/检验
        print("H matrix:")
        for i in range(m):
            print(" ".join(map(str, H[i])))

        print("V matrix:")
        for i in range(m):
            print(" ".join(map(str, V[i])))

        # 简单校验：顶点四边互异 + C4 四边互异
        def check(H, V):
            # 顶点
            for i in range(m):
                for j in range(n):
                    colors = [H[i,j], H[i,(j-1)%n], V[i,j], V[(i-1)%m, j]]
                    if len(set(colors)) != 4:
                        return False, ("vertex", i, j, colors)
            # C4
            for i in range(m):
                for j in range(n):
                    e1 = H[i,j]
                    e2 = V[i,(j+1)%n]
                    e3 = H[(i+1)%m, j]
                    e4 = V[i,j]
                    if len(set([e1,e2,e3,e4])) != 4:
                        return False, ("face", i, j, [e1,e2,e3,e4])
            return True, None

        ok, info = check(H, V)
        print("VALID B-coloring?" , ok, info if not ok else "")


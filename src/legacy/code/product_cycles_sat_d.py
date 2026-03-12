#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cartesian product of cycles ☐_i C_{n_i}  edge B-coloring via SAT
- Proper edge-coloring (incident edges pairwise different)
- Every 4-cycle rainbow
- Parameters edited IN-CODE (see CONFIG below)
- Solver: cryptominisat5 (exe) or pycryptosat
"""

import os, subprocess
from collections import defaultdict

# =========================
# ======== CONFIG =========
# =========================
DIMS = (3, 3, 3, 3, 3)          # 因子圈长 (如 (3,3,3) / (7,7))
K = 11                # 颜色数
USE_SOLVER = "cmsat"   # "cmsat" 调用 cryptominisat5 可执行程序, 或 "pycryptosat"
CMSAT_PATH = "E:\cryptominisat5-win\cryptominisat5.exe" # 若不在 PATH，请填完整路径
CNF_PATH = "out.cnf"   # 输出 CNF 文件路径
SYM_BREAK = True        # 启用对称破坏
EXPORT_GEPHI = True     # 导出 Gephi 节点/边 CSV
PLOT_FLAT = True        # 是否绘制 NetworkX 平面示意图

# =========================
# ===== 基础构造函数 ======
# =========================

def prod(it):
    from math import prod as _p
    return _p(it)

def tup_id(t, dims):
    v = 0; mul = 1
    for i in range(len(dims)-1, -1, -1):
        v += t[i]*mul
        mul *= dims[i]
    return v

def id_tup(v, dims):
    r = []
    for n in reversed(dims):
        r.append(v % n); v//=n
    return tuple(reversed(r))

def build_graph(dims):
    d = len(dims)
    Vn = prod(dims)
    V = list(range(Vn))
    edges = []
    inc = defaultdict(list)
    seen = set()
    for v in V:
        t = id_tup(v, dims)
        for ax in range(d):
            s = list(t); s[ax] = (s[ax]+1) % dims[ax]; s = tuple(s)
            w = tup_id(s, dims)
            a,b = (v,w) if v < w else (w,v)
            if (a,b) in seen: continue
            seen.add((a,b))
            eid = len(edges)
            edges.append((a,b))
            inc[a].append(eid)
            inc[b].append(eid)
    assert len(edges) == Vn*d, f"|E| expect {Vn*d}, got {len(edges)}"
    for v in V:
        assert len(inc[v]) == 2*d
    return V, edges, inc

def enumerate_all_C4(dims, edges):
    d = len(dims); Vn = prod(dims)
    edge_index = {}
    for eid,(u,v) in enumerate(edges):
        edge_index[(u,v)] = eid; edge_index[(v,u)] = eid
    def e(u,v): return edge_index[(u,v)]
    C4s = []
    for a in range(d):
        for b in range(a+1, d):
            for v in range(Vn):
                t  = id_tup(v, dims)
                t1 = list(t); t1[a]=(t1[a]+1)%dims[a]; t1=tuple(t1)
                t3 = list(t); t3[b]=(t3[b]+1)%dims[b]; t3=tuple(t3)
                t2 = list(t3); t2[a]=(t2[a]+1)%dims[a]; t2=tuple(t2)
                v0=v
                v1=tup_id(t1,dims)
                v3=tup_id(t3,dims)
                v2=tup_id(t2,dims)
                C4s.append([ e(v0,v1), e(v1,v2), e(v2,v3), e(v3,v0) ])
    return C4s

# =========================
# ====== CNF Builder ======
# =========================

class CNF:
    def __init__(self, nvars):
        self.nvars = nvars
        self.clauses = []
    def add(self, lits): self.clauses.append(list(lits))
    def atleast1(self, lits): self.add(lits)
    def atmost1_pairwise(self, lits):
        for i in range(len(lits)):
            for j in range(i+1,len(lits)):
                self.add([-lits[i], -lits[j]])
    def exactly1_pairwise(self, lits):
        self.atleast1(lits); self.atmost1_pairwise(lits)

def var_id(eid,c,K): return eid*K + c

def build_cnf(dims, K, sym_break=True):
    V, edges, inc = build_graph(dims)
    C4s = enumerate_all_C4(dims, edges)
    nvars = len(edges)*K
    cnf = CNF(nvars)
    # 每边恰一色
    for eid in range(len(edges)):
        lits = [var_id(eid,c,K) for c in range(1,K+1)]
        cnf.exactly1_pairwise(lits)
    # incident 异色
    for v in V:
        Es = inc[v]
        for c in range(1,K+1):
            cnf.atmost1_pairwise([var_id(e,c,K) for e in Es])
    # 每个 4-圈四色
    for cyc in C4s:
        for c in range(1,K+1):
            cnf.atmost1_pairwise([var_id(e,c,K) for e in cyc])
    # 对称破坏
    if sym_break:
        v0 = 0
        Es = inc[v0]
        tmp = []
        for e in Es:
            u,w = edges[e]
            nb = w if u==v0 else u
            tmp.append((nb,e))
        tmp.sort()
        ordered = [e for (_,e) in tmp]
        for kcol,eid in enumerate(ordered, start=1):
            if kcol>min(K,len(ordered)): break
            cnf.add([var_id(eid,kcol,K)])
    info = {"V":V,"edges":edges,"inc":inc,"C4s":C4s}
    return cnf, info

# =========================
# ====== Solver 调用 ======
# =========================

def write_dimacs(cnf, path):
    with open(path,"w") as f:
        f.write(f"p cnf {cnf.nvars} {len(cnf.clauses)}\n")
        for cls in cnf.clauses:
            f.write(" ".join(str(x) for x in cls) + " 0\n")

def parse_cmsat_output(text, nvars):
    sat=None; vals=[]
    for line in text.splitlines():
        s=line.strip()
        if s.startswith("s "):
            if "UNSAT" in s: sat=False
            elif "SAT" in s: sat=True
        elif s.startswith("v ") or s.startswith("V "):
            for p in s.split()[1:]:
                if p=="0": continue
                try: vals.append(int(p))
                except: pass
    if sat is None:
        if "UNSAT" in text: sat=False
        elif "SAT" in text: sat=True
        else: return None,None
    if not sat: return False,None
    assign={v:False for v in range(1,nvars+1)}
    for lit in vals:
        if lit==0: continue
        v=abs(lit)
        if v<=nvars: assign[v]=(lit>0)
    return True,assign

def solve_with_cmsat(cnf, path=CMSAT_PATH, cnf_path=CNF_PATH):
    write_dimacs(cnf, cnf_path)
    try:
        proc = subprocess.run([path, cnf_path],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except FileNotFoundError:
        print(f"[!] solver not found: {path}")
        return None,None
    out = proc.stdout + "\n" + proc.stderr
    return parse_cmsat_output(out, cnf.nvars)

def solve_with_pycryptosat(cnf):
    from pycryptosat import Solver
    s = Solver()
    for cls in cnf.clauses: s.add_clause(cls)
    sat,model = s.solve()
    if not sat: return False,None
    assign={}
    if isinstance(model,dict): assign=model
    else:
        for v in range(1,cnf.nvars+1):
            assign[v]=bool(model[v])
    return True,assign

# =========================
# ====== 后处理 ==========
# =========================

# —— 生成任意 K 个可区分颜色（返回 hex 字符串列表）——
def make_palette(K, name="tab20"):
    """
    生成长度为 K 的颜色列表（十六进制字符串）。
    默认使用 matplotlib 的离散 colormap（tab20/tab20b/tab20c），
    若 K > 20 则退化到 hsv 均匀取样。
    """
    import matplotlib
    import matplotlib.cm as cm

    def rgba_to_hex(rgba):
        r,g,b,a = rgba
        return "#{:02X}{:02X}{:02X}".format(int(r*255), int(g*255), int(b*255))

    if K <= 10:
        cmap = cm.get_cmap("tab10", K)
        return [rgba_to_hex(cmap(i)) for i in range(K)]
    elif K <= 20:
        cmap = cm.get_cmap(name, K)  # "tab20"
        return [rgba_to_hex(cmap(i)) for i in range(K)]
    else:
        # K 很大时用 HSV 等角采样
        cmap = cm.get_cmap("hsv", K)
        return [rgba_to_hex(cmap(i)) for i in range(K)]


def assignment_to_colors(assign, edges, K):
    colors=[-1]*len(edges)
    for eid in range(len(edges)):
        picks=[c for c in range(1,K+1) if assign.get(var_id(eid,c,K),False)]
        if len(picks)!=1:
            raise RuntimeError(f"edge {eid} wrong colors {picks}")
        colors[eid]=picks[0]
    return colors

def verify_solution(colors, info):
    V,edges,inc,C4s = info["V"],info["edges"],info["inc"],info["C4s"]
    for v in V:
        seen=set()
        for e in inc[v]:
            c=colors[e]
            if c in seen:
                raise AssertionError(f"vertex {v} duplicate color {c}")
            seen.add(c)
    for cyc in C4s:
        cols=[colors[e] for e in cyc]
        if len(set(cols))!=4:
            raise AssertionError(f"C4 not rainbow {cols}")

def export_gephi(V, edges, colors, dims, prefix):
    import pandas as pd
    nodes=[{"Id":v,"Label":str(id_tup(v,dims))} for v in V]
    pd.DataFrame(nodes).to_csv(f"{prefix}_nodes.csv",index=False)
    # palette=["#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628"]
    palette = make_palette(max(colors))  # 或 make_palette(K)

    rows=[]
    for eid,(u,v) in enumerate(edges):
        c=colors[eid]
        rows.append({"Id":eid,"Source":u,"Target":v,
                     "Type":"Undirected","Color":c,
                     "ColorRGB":palette[(c-1)%len(palette)]})
    pd.DataFrame(rows).to_csv(f"{prefix}_edges.csv",index=False)
    print(f"[Gephi] wrote {prefix}_nodes.csv / {prefix}_edges.csv")

def plot_flat(edges, colors, K):
    import networkx as nx, matplotlib.pyplot as plt, matplotlib.lines as mlines
    G=nx.Graph(); G.add_edges_from(edges)
    pos=nx.kamada_kawai_layout(G)

    palette = make_palette(K)
    edge_colors = [palette[(c-1) % len(palette)] for c in colors]

    # 图例：K 很大时图例会很挤，可按需裁剪或分栏
    max_legend = min(K, 20)  # 只展示前 20 个颜色的图例，避免遮挡
    handles = [mlines.Line2D([], [], color=palette[i], lw=3, label=f"color {i+1}")
           for i in range(max_legend)]

    # palette=["#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628"]
    # edge_colors=[palette[(c-1)%len(palette)] for c in colors]
    plt.figure(figsize=(8,7))
    nx.draw_networkx_nodes(G,pos,node_color="white",edgecolors="black",node_size=260)
    nx.draw_networkx_edges(G,pos,edgelist=edges,edge_color=edge_colors,width=2.2)
    nx.draw_networkx_labels(G,pos,font_size=9)
    # handles=[mlines.Line2D([],[],color=palette[i],lw=3,label=f"color {i+1}") for i in range(K)]
    plt.legend(handles=handles,loc="upper left",bbox_to_anchor=(1.02,1.0))
    plt.axis("off"); plt.tight_layout(); plt.show()

# =========================
# ========= 主程序 ========
# =========================

def main():
    dims=DIMS; Kc=K
    print(f"[Info] dims={dims}, K={Kc}")
    cnf,info = build_cnf(dims,Kc,sym_break=SYM_BREAK)
    print(f"[Build] Vars={cnf.nvars}, Clauses={len(cnf.clauses)}, |E|={len(info['edges'])}, |C4|={len(info['C4s'])}")
    if USE_SOLVER=="cmsat":
        sat,assign = solve_with_cmsat(cnf, CMSAT_PATH, CNF_PATH)
    else:
        sat,assign = solve_with_pycryptosat(cnf)
    if sat is None:
        print("[!] Solver failed."); return
    if not sat:
        print("[Result] UNSAT"); return
    print("[Result] SAT")
    colors=assignment_to_colors(assign,info["edges"],Kc)
    verify_solution(colors,info)
    print("[Check] Verified OK (proper + rainbow C4).")
    from collections import Counter
    cnt=Counter(colors)
    for c in sorted(cnt): print(f"  color {c}: {cnt[c]} edges")
    prefix=f"Cprod_{'x'.join(map(str,dims))}_K{Kc}"
    if EXPORT_GEPHI:
        export_gephi(info["V"],info["edges"],colors,dims,prefix)
    if PLOT_FLAT:
        plot_flat(info["edges"],colors,Kc)

if __name__=="__main__":
    main()

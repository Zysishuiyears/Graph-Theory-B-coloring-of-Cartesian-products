#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cartesian product of cycles ☐_i C_{n_i}  edge B-coloring via SAT  (Parallel Edition)
- Proper edge-coloring (incident edges pairwise different)
- Every 4-cycle rainbow
- Parameters edited IN-CODE (see CONFIG below)
- Solver: cryptominisat5 (exe) or pycryptosat
- Parallel CNF building via multiprocessing; solver threading configurable
"""

import os, subprocess
from collections import defaultdict
import csv
import math
from multiprocessing import Pool, cpu_count

# =========================
# ======== CONFIG =========
# =========================
DIMS = (3, 3, 3)          # 因子圈长 (如 (3,3,3) / (7,7))
K = 7                    # 颜色数
USE_SOLVER = "cmsat"      # "cmsat" 调用 cryptominisat5 可执行程序, 或 "pycryptosat"
CMSAT_PATH = r"E:\cryptominisat5-win\cryptominisat5.exe"  # 建议用原始字符串 r"..."
CNF_PATH = "out.cnf"      # 输出 CNF 文件路径
SYM_BREAK = True          # 启用对称破坏
EXPORT_GEPHI = True       # 导出 Gephi 节点/边 CSV
PLOT_FLAT = True          # 绘制总体 NetworkX 平面示意图

# 新增：输出“低一维所有副本”的图（例如 C3^3 -> 三张 C3^2）
PLOT_SLICES_LOWER = True
SLICE_AXIS = -1           # 默认沿“最后一维”切片；也可填 0/1/2... 指定维度

# ====== 并行相关（新增）======
N_WORKERS = max(1, min(cpu_count(), 8))  # Python 侧 CNF 构造的并行进程数
CMSAT_THREADS = max(1, min(cpu_count(), 8))  # 传给 CMS 的线程数（cmsat 与 pycryptosat）

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
    def add_many(self, list_of_lits_lists):
        self.clauses.extend(list_of_lits_lists)
    def atleast1(self, lits): self.add(lits)
    def atmost1_pairwise(self, lits):
        for i in range(len(lits)):
            li = lits[i]
            for j in range(i+1,len(lits)):
                self.add([-li, -lits[j]])
    def exactly1_pairwise(self, lits):
        self.atleast1(lits); self.atmost1_pairwise(lits)

def var_id(eid,c,K): return eid*K + c

# ---- 并行子任务：incident 约束 ----
def _build_incident_chunk(args):
    Es_chunk, K = args
    out = []
    # 对每个颜色：同一顶点所有入边两两不同色 => pairwise at-most-1
    for Es in Es_chunk:
        for c in range(1, K+1):
            if len(Es) <= 1: 
                continue
            for i in range(len(Es)):
                vi = var_id(Es[i], c, K)
                for j in range(i+1, len(Es)):
                    out.append([-vi, -var_id(Es[j], c, K)])
    return out

# ---- 并行子任务：C4 每色 at-most-1 ----
def _build_c4_chunk(args):
    C4_chunk, K = args
    out = []
    for cyc in C4_chunk:
        for c in range(1, K+1):
            # 4 条边，同色至多 1 条 => pairwise
            out.extend([[-var_id(cyc[i], c, K), -var_id(cyc[j], c, K)]
                        for i in range(4) for j in range(i+1, 4)])
    return out

def build_cnf_parallel(dims, K, sym_break=True, n_workers=1):
    """
    将 CNF 构造中的重活拆到多进程：
      - incident 异色约束
      - 每个 C4 的四色性（每色 at-most-1）
    其余部分（每边恰一色、对称破坏）保持串行添加。
    """
    V, edges, inc = build_graph(dims)
    C4s = enumerate_all_C4(dims, edges)
    nvars = len(edges)*K
    cnf = CNF(nvars)

    # 每边恰一色（exactly-1）——串行
    for eid in range(len(edges)):
        lits = [var_id(eid,c,K) for c in range(1,K+1)]
        cnf.exactly1_pairwise(lits)

    # incident 异色 —— 并行
    Es_all = [inc[v] for v in V]
    if n_workers <= 1 or len(Es_all) < 32:
        # 小实例不用并行
        for Es in Es_all:
            for c in range(1, K+1):
                cnf.atmost1_pairwise([var_id(e, c, K) for e in Es])
    else:
        chunk = math.ceil(len(Es_all) / n_workers)
        tasks = [(Es_all[i:i+chunk], K) for i in range(0, len(Es_all), chunk)]
        with Pool(processes=n_workers) as pool:
            for sub in pool.imap_unordered(_build_incident_chunk, tasks):
                cnf.add_many(sub)

    # 每个 4-圈四色 —— 并行
    if n_workers <= 1 or len(C4s) < 64:
        for cyc in C4s:
            for c in range(1, K+1):
                cnf.atmost1_pairwise([var_id(e, c, K) for e in cyc])
    else:
        chunk = math.ceil(len(C4s) / n_workers)
        tasks = [(C4s[i:i+chunk], K) for i in range(0, len(C4s), chunk)]
        with Pool(processes=n_workers) as pool:
            for sub in pool.imap_unordered(_build_c4_chunk, tasks):
                cnf.add_many(sub)

    # 对称破坏 —— 串行（保持确定性）
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

def solve_with_cmsat(cnf, path=CMSAT_PATH, cnf_path=CNF_PATH, threads=1):
    write_dimacs(cnf, cnf_path)
    args = [path, f"--threads={threads}", cnf_path]
    try:
        proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except FileNotFoundError:
        print(f"[!] solver not found: {path}")
        return None,None
    out = proc.stdout + "\n" + proc.stderr
    return parse_cmsat_output(out, cnf.nvars)

def solve_with_pycryptosat(cnf, threads=1):
    from pycryptosat import Solver
    s = Solver(threads=threads)
    for cls in cnf.clauses: s.add_clause(cls)
    sat,model = s.solve()
    if not sat: return False,None
    assign={}
    if isinstance(model,dict): assign=model
    else:
        for v in range(1,cnf.nvars+1):
            assign[v]=bool(model[v])
    return True,assign

## ======== 优化SAT求解结果 ======== 
def calculate_color_patterns(colors, dims, edges, K):
    """
    计算每一行或每一列的颜色模式并检测周期性。
    返回：周期模式得分（越接近1表示越规则）
    """
    color_patterns = defaultdict(list)
    for eid, color in enumerate(colors):
        u, v = edges[eid]
        
        # 根据 dims 维度，解包为 3 维坐标
        coords_u = id_tup(u, dims)
        coords_v = id_tup(v, dims)
        
        # 计算每行的颜色模式
        row_u, col_u, _ = coords_u  # 假设 dims 为 (3,3,3)，所以返回值是一个三维元组
        row_v, col_v, _ = coords_v  # 需要解包三个值
        
        color_patterns[row_u].append(color)
        color_patterns[row_v].append(color)
        
    # 检测周期性（例如：计算每列的颜色变化频率）
    row_pattern_similarity = []
    for row in color_patterns:
        pattern = color_patterns[row]
        pattern_similarity = sum([1 for i in range(len(pattern)-1) if pattern[i] == pattern[i+1]]) / len(pattern)
        row_pattern_similarity.append(pattern_similarity)
        
    # 返回一个“规则性”得分，越高越有规律
    return sum(row_pattern_similarity) / len(row_pattern_similarity) if row_pattern_similarity else 0
        
def optimize_for_periodicity(assign, dims, edges, K):
    """
    优化染色方案以增加其周期性：基于现有解调整，使其更具规则性。
    """
    colors = assignment_to_colors(assign, edges, K)
    initial_pattern_score = calculate_color_patterns(colors, dims, edges, K)
    
    # 基于初始模式得分，逐步尝试优化
    optimized_assign = assign.copy()
    for _ in range(10):  # 循环优化10次
        optimized_assign = adjust_for_pattern(optimized_assign, dims, edges, K)
        new_colors = assignment_to_colors(optimized_assign, edges, K)
        new_pattern_score = calculate_color_patterns(new_colors, dims, edges, K)
        
        # 如果优化后规则性得分更高，则更新
        if new_pattern_score > initial_pattern_score:
            colors = new_colors
            initial_pattern_score = new_pattern_score
    return optimized_assign

## == 调整解的规则性：实现一个调整函数，按照模式和周期性规则调整当前解 == 

def adjust_for_pattern(assign, dims, edges, K):
    """
    尝试调整解中的颜色模式，使其尽可能更具规律性。
    """
    colors = assignment_to_colors(assign, edges, K)
    new_assign = assign.copy()
    
    for i, color in enumerate(colors):
        # 如果某些颜色频繁重复，则调整其他部分的颜色
        if color == 1:  # 假设我们认为1号颜色很频繁，则调整1号色附近的颜色
            new_assign[i] = adjust_color(new_assign[i], dims, edges, K)
    
    return new_assign

def adjust_color(assign_value, dims, edges, K):
    """
    根据某些规则调整颜色方案（例如通过颜色模式平滑）
    """
    # 假设通过简单的递增或减少颜色
    return (assign_value + 1) % K + 1



# =========================
# ====== 后处理 ========== 
# =========================

def make_palette(K, name="tab20"):
    import matplotlib
    import matplotlib.cm as cm
    def rgba_to_hex(rgba):
        r,g,b,a = rgba
        return "#{:02X}{:02X}{:02X}".format(int(r*255), int(g*255), int(b*255))
    if K <= 10:
        cmap = cm.get_cmap("tab10", K)
        return [rgba_to_hex(cmap(i)) for i in range(K)]
    elif K <= 20:
        cmap = cm.get_cmap(name, K)
        return [rgba_to_hex(cmap(i)) for i in range(K)]
    else:
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
    palette = make_palette(max(colors))
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
    max_legend = min(K, 20)
    handles = [mlines.Line2D([], [], color=palette[i], lw=3, label=f"color {i+1}")
           for i in range(max_legend)]
    plt.figure(figsize=(8,7))
    nx.draw_networkx_nodes(G,pos,node_color="white",edgecolors="black",node_size=260)
    nx.draw_networkx_edges(G,pos,edgelist=edges,edge_color=edge_colors,width=2.2)
    nx.draw_networkx_labels(G,pos,font_size=9)
    plt.legend(handles=handles,loc="upper left",bbox_to_anchor=(1.02,1.0))
    plt.axis("off"); plt.tight_layout(); plt.show()

def plot_lower_dim_slices(dims, edges, colors, K, slice_axis=-1):
    import matplotlib.pyplot as plt, matplotlib.lines as mlines
    import numpy as np
    d = len(dims)
    if d < 3:
        print("[Slices] skipped: dims has no lower 1-D copies."); return
    if slice_axis < 0: slice_axis += d
    assert 0 <= slice_axis < d
    Vn = prod(dims)
    e2id = {}
    for eid,(u,v) in enumerate(edges):
        key = (u,v) if u<v else (v,u)
        e2id[key] = eid
    def col_of(u,v):
        key = (u,v) if u<v else (v,u)
        return colors[e2id[key]]
    remain_axes = [ax for ax in range(d) if ax != slice_axis]
    assert len(remain_axes) >= 2, "切片后至少需两维用于平面"
    a, b = remain_axes[:2]
    na, nb = dims[a], dims[b]
    palette = make_palette(K)
    def draw_periodic_curve(ax, p, q, color, box=(na, nb), bend=0.28, lw=2.4, alpha=0.98):
        na_, nb_ = box
        px, py = p; qx, qy = q
        choices_x = [qx- na_, qx, qx+ na_]
        choices_y = [qy- nb_, qy, qy+ nb_]
        best = None; best_dist = 1e9
        for xx in choices_x:
            for yy in choices_y:
                dist = (xx-px)**2 + (yy-py)**2
                if dist < best_dist:
                    best_dist = dist; best = (xx, yy)
        qx_u, qy_u = best
        vx, vy = (qx_u - px, qy_u - py)
        L = (vx**2 + vy**2)**0.5 + 1e-9
        nx_, ny_ = -vy / L, vx / L
        cx, cy = (px + qx_u) / 2 + bend * nx_, (py + qy_u) / 2 + bend * ny_
        t = np.linspace(0, 1, 60)
        Xu = (1-t)**2 * px + 2*(1-t)*t * cx + t**2 * qx_u
        Yu = (1-t)**2 * py + 2*(1-t)*t * cy + t**2 * qy_u
        X = np.mod(Xu, na_); Y = np.mod(Yu, nb_)
        seg_x, seg_y = [X[0]], [Y[0]]
        thr = max(na_, nb_) * 0.35
        for i in range(1, len(X)):
            dx = abs(X[i] - X[i-1]); dy = abs(Y[i] - Y[i-1])
            dx = min(dx, na_ - dx); dy = min(dy, nb_ - dy)
            if (dx*dx + dy*dy) ** 0.5 > thr:
                ax.plot(seg_x, seg_y, color=color, lw=lw, alpha=alpha, solid_capstyle="round")
                seg_x, seg_y = [X[i]], [Y[i]]
            else:
                seg_x.append(X[i]); seg_y.append(Y[i])
        if len(seg_x) > 1:
            ax.plot(seg_x, seg_y, color=color, lw=lw, alpha=alpha, solid_capstyle="round")
    for s in range(dims[slice_axis]):
        Vn = prod(dims)
        S = {v for v in range(Vn) if id_tup(v, dims)[slice_axis] == s}
        e_slice = [(u, v) for (u, v) in edges if (u in S and v in S)]
        base_pos = {}
        for v in S:
            t = id_tup(v, dims); base_pos[v] = (t[a], t[b])
        fig, ax = plt.subplots(figsize=(6.4, 6.2))
        for (u, v) in e_slice:
            cu = palette[(col_of(u, v) - 1) % K]
            p = base_pos[u]; q = base_pos[v]
            draw_periodic_curve(ax, p, q, cu, box=(na, nb), bend=0.26, lw=2.4, alpha=0.98)
        xs = [base_pos[v][0] for v in S]; ys = [base_pos[v][1] for v in S]
        ax.scatter(xs, ys, s=260, facecolors="white", edgecolors="black", zorder=3)
        for v in S:
            x, y = base_pos[v]
            ax.text(x, y, str(id_tup(v, dims)), fontsize=9, ha="center", va="center", zorder=4)
        ax.plot([0, na, na, 0, 0], [0, 0, nb, nb, 0], color="#999999", lw=0.8, ls="dotted", alpha=0.8)
        ax.set_xlim(-0.5, na + 0.5); ax.set_ylim(-0.5, nb + 0.5)
        ax.set_aspect("equal"); ax.set_axis_off()
        import matplotlib.lines as mlines
        handles = [mlines.Line2D([], [], color=palette[i], lw=3, label=f"color {i+1}")
                   for i in range(min(K, 12))]
        if handles:
            ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0))
        import matplotlib.pyplot as plt
        plt.tight_layout()
        plt.title(f"Slice at axis {slice_axis} = {s}  (periodic curves on torus)")
        plt.show()

def report_intercopy_colors(dims, edges, colors, slice_axis=2, export_csv=True, prefix="intercopy"):
    d = len(dims)
    assert 0 <= slice_axis < d
    Vn = prod(dims)
    e2id = {}
    for eid, (u, v) in enumerate(edges):
        key = (u, v) if u < v else (v, u)
        e2id[key] = eid
    def edge_color(u, v):
        key = (u, v) if u < v else (v, u)
        eid = e2id[key]
        return colors[eid], eid
    rows = []
    remain_axes = [ax for ax in range(d) if ax != slice_axis]
    ns = dims[slice_axis]
    for v in range(Vn):
        t = list(id_tup(v, dims))
        t2 = t[:]
        t2[slice_axis] = (t2[slice_axis] + 1) % ns
        u = tup_id(tuple(t2), dims)
        c, eid = edge_color(v, u)
        rows.append({
            "eid": eid,
            "color": c,
            "from": tuple(t),
            "to": tuple(t2),
            "slice_axis": slice_axis,
            "s": t[slice_axis],
        })
    if export_csv:
        fn = f"{prefix}_axis{slice_axis}_{'x'.join(map(str,dims))}.csv"
        with open(fn, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["eid", "color", "from", "to", "slice_axis", "s"])
            for r in rows:
                w.writerow([r["eid"], r["color"], r["from"], r["to"], r["slice_axis"], r["s"]])
        print(f"[InterCopy] wrote {fn} (total {len(rows)} edges)")
    if d == 3:
        a, b = remain_axes
        na, nb = dims[a], dims[b]
        for s in range(ns):
            grid = [[None for _ in range(nb)] for __ in range(na)]
            for r in rows:
                if r["s"] != s: continue
                t = r["from"]
                ia, ib = t[a], t[b]
                grid[ia][ib] = r["color"]
            print(f"\n[InterCopy axis={slice_axis} layer {s} → {(s+1)%ns}]  (shape {na}x{nb})")
            for ia in range(na):
                print("  " + " ".join(str(grid[ia][ib]) for ib in range(nb)))
    else:
        from collections import Counter
        cnt = Counter(r["color"] for r in rows)
        print(f"[InterCopy axis={slice_axis}] color histogram:", dict(sorted(cnt.items())))

# =========================
# ========= 主程序 ========
# =========================

# def main():
#     dims=DIMS; Kc=K
#     print(f"[Info] dims={dims}, K={Kc}, N_WORKERS={N_WORKERS}, CMSAT_THREADS={CMSAT_THREADS}")
#     cnf,info = build_cnf_parallel(dims,Kc,sym_break=SYM_BREAK,n_workers=N_WORKERS)
#     print(f"[Build] Vars={cnf.nvars}, Clauses={len(cnf.clauses)}, |E|={len(info['edges'])}, |C4|={len(info['C4s'])}")
#     if USE_SOLVER=="cmsat":
#         sat,assign = solve_with_cmsat(cnf, CMSAT_PATH, CNF_PATH, threads=CMSAT_THREADS)
#     else:
#         sat,assign = solve_with_pycryptosat(cnf, threads=CMSAT_THREADS)
#     if sat is None:
#         print("[!] Solver failed."); return
#     if not sat:
#         print("[Result] UNSAT"); return
#     print("[Result] SAT")
#     colors=assignment_to_colors(assign,info["edges"],Kc)
#     verify_solution(colors,info)
#     print("[Check] Verified OK (proper + rainbow C4).")
#     from collections import Counter
#     cnt=Counter(colors)
#     for c in sorted(cnt): print(f"  color {c}: {cnt[c]} edges")
#     prefix=f"Cprod_{'x'.join(map(str,dims))}_K{Kc}"
#     if EXPORT_GEPHI:
#         export_gephi(info["V"],info["edges"],colors,dims,prefix)
#     if PLOT_FLAT:
#         plot_flat(info["edges"],colors,Kc)
#     if PLOT_SLICES_LOWER:
#         slice_ax = (len(dims)-1) if SLICE_AXIS==-1 else SLICE_AXIS
#         plot_lower_dim_slices(dims, info["edges"], colors, Kc, slice_axis=slice_ax)
def main():
    dims=DIMS; Kc=K
    print(f"[Info] dims={dims}, K={Kc}, N_WORKERS={N_WORKERS}, CMSAT_THREADS={CMSAT_THREADS}")
    cnf, info = build_cnf_parallel(dims, Kc, sym_break=SYM_BREAK, n_workers=N_WORKERS)
    print(f"[Build] Vars={cnf.nvars}, Clauses={len(cnf.clauses)}, |E|={len(info['edges'])}, |C4|={len(info['C4s'])}")
    
    if USE_SOLVER == "cmsat":
        sat, assign = solve_with_cmsat(cnf, CMSAT_PATH, CNF_PATH, threads=CMSAT_THREADS)
    else:
        sat, assign = solve_with_pycryptosat(cnf, threads=CMSAT_THREADS)
        
    if sat is None:
        print("[!] Solver failed.")
        return
    if not sat:
        print("[Result] UNSAT")
        return
    print("[Result] SAT")
    
    # 优化解，使其具有更强的周期性
    optimized_assign = optimize_for_periodicity(assign, dims, info["edges"], Kc)
    
    colors = assignment_to_colors(optimized_assign, info["edges"], Kc)
    verify_solution(colors, info)
    print("[Check] Verified OK (proper + rainbow C4).")
    
    from collections import Counter
    cnt = Counter(colors)
    for c in sorted(cnt): 
        print(f"  color {c}: {cnt[c]} edges")
        
    prefix = f"Cprod_{'x'.join(map(str, dims))}_K{Kc}"
    if EXPORT_GEPHI:
        export_gephi(info["V"], info["edges"], colors, dims, prefix)
    if PLOT_FLAT:
        plot_flat(info["edges"], colors, Kc)
    if PLOT_SLICES_LOWER:
        slice_ax = (len(dims) - 1) if SLICE_AXIS == -1 else SLICE_AXIS
        plot_lower_dim_slices(dims, info["edges"], colors, Kc, slice_axis=slice_ax)


if __name__=="__main__":
    main()


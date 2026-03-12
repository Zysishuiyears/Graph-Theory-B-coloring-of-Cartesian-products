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
DIMS = (3, 3, 3)          # 因子圈长 (如 (3,3,3) / (7,7))
K = 7                     # 颜色数
USE_SOLVER = "cmsat"      # "cmsat" 调用 cryptominisat5 可执行程序, 或 "pycryptosat"
CMSAT_PATH = r"E:\cryptominisat5-win\cryptominisat5.exe"  # 建议用原始字符串 r"..."
CNF_PATH = "out.cnf"      # 输出 CNF 文件路径
SYM_BREAK = True          # 启用对称破坏
EXPORT_GEPHI = True       # 导出 Gephi 节点/边 CSV
PLOT_FLAT = True          # 绘制总体 NetworkX 平面示意图

# 新增：输出“低一维所有副本”的图（例如 C3^3 -> 三张 C3^2）
PLOT_SLICES_LOWER = True
SLICE_AXIS = -1          # 默认沿“最后一维”切片；也可填 0/1/2... 指定维度

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
    max_legend = min(K, 20)
    handles = [mlines.Line2D([], [], color=palette[i], lw=3, label=f"color {i+1}")
           for i in range(max_legend)]
    plt.figure(figsize=(8,7))
    nx.draw_networkx_nodes(G,pos,node_color="white",edgecolors="black",node_size=260)
    nx.draw_networkx_edges(G,pos,edgelist=edges,edge_color=edge_colors,width=2.2)
    nx.draw_networkx_labels(G,pos,font_size=9)
    plt.legend(handles=handles,loc="upper left",bbox_to_anchor=(1.02,1.0))
    plt.axis("off"); plt.tight_layout(); plt.show()

# === 新增：输出低一维所有副本的染色图 ===

# def plot_lower_dim_slices(dims, edges, colors, K, slice_axis=-1):
#     import networkx as nx, matplotlib.pyplot as plt, matplotlib.lines as mlines
#     d = len(dims)
#     if d < 3:
#         print("[Slices] skipped: dims has no lower 1-D copies."); return
#     if slice_axis < 0: slice_axis += d
#     assert 0 <= slice_axis < d

#     Vn = prod(dims)
#     # (u,v)->eid
#     e2id = {}
#     for eid, (u, v) in enumerate(edges):
#         key = (u, v) if u < v else (v, u)
#         e2id[key] = eid

#     def col_of(u, v):
#         key = (u, v) if u < v else (v, u)
#         return colors[e2id[key]]

#     remain_axes = [ax for ax in range(d) if ax != slice_axis]
#     assert len(remain_axes) >= 2, "切到 (d-1) 维后至少应有2维用于平面"
#     a, b = remain_axes[:2]           # 用这两维做平面坐标
#     na, nb = dims[a], dims[b]
#     palette = make_palette(K)

#     for s in range(dims[slice_axis]):
#         # 该切片的顶点 & 边
#         S = {v for v in range(Vn) if id_tup(v, dims)[slice_axis] == s}
#         e_slice = [(u, v) for (u, v) in edges if (u in S and v in S)]

#         # 规则网格原始坐标
#         base_pos = {}
#         for v in S:
#             t = id_tup(v, dims)
#             base_pos[v] = (t[a], t[b])

#         # 画布
#         fig, ax = plt.subplots(figsize=(6.2, 6.0))
#         # 节点（一次性画所有真实点）
#         ax.scatter([base_pos[v][0] for v in S], [base_pos[v][1] for v in S],
#                    s=260, facecolors="white", edgecolors="black", zorder=3)

#         # —— 画边：非回环边一段；回环边拆成两段并用“幽灵坐标”闭合 ——
#         for (u, v) in e_slice:
#             tu, tv = id_tup(u, dims), id_tup(v, dims)
#             cu = palette[(col_of(u, v) - 1) % K]

#             # 哪一维在此切片内变化（只会有一维）
#             if (tv[a] - tu[a]) % na != 0:
#                 changed = a; n_changed = na
#             else:
#                 changed = b; n_changed = nb

#             x1, y1 = base_pos[u]
#             x2, y2 = base_pos[v]

#             def draw_seg(p, q, color, style="solid", alpha=0.98, lw=2.4):
#                 ax.plot([p[0], q[0]], [p[1], q[1]], color=color, lw=lw,
#                         linestyle=style, alpha=alpha, zorder=1)

#             # 非回环（例如 0->1 或 1->2）：直接一段
#             wrap = False
#             if changed == a:
#                 # 回环当且仅当 (na-1 -> 0) 或 (0 -> na-1)
#                 wrap = (tu[a] == na - 1 and tv[a] == 0) or (tv[a] == na - 1 and tu[a] == 0)
#                 if wrap:
#                     # 拆成两段： (na-1, j) -> (na, j)  与  (0, j) -> (na, j)
#                     ghost = (na, y1) if tu[a] == na - 1 else (na, y2)
#                     if tu[a] == na - 1:  # u 在右边界
#                         draw_seg((x1, y1), ghost, cu, style="dashed", alpha=0.8, lw=2.0)
#                         draw_seg((0, y2),  ghost, cu, style="dashed", alpha=0.8, lw=2.0)
#                     else:                # v 在右边界
#                         draw_seg((x2, y2), ghost, cu, style="dashed", alpha=0.8, lw=2.0)
#                         draw_seg((0, y1),  ghost, cu, style="dashed", alpha=0.8, lw=2.0)
#                 else:
#                     draw_seg((x1, y1), (x2, y2), cu)
#             else:
#                 # changed == b
#                 wrap = (tu[b] == nb - 1 and tv[b] == 0) or (tv[b] == nb - 1 and tu[b] == 0)
#                 if wrap:
#                     ghost = (x1, nb) if tu[b] == nb - 1 else (x2, nb)
#                     if tu[b] == nb - 1:  # u 在上边界
#                         draw_seg((x1, y1), ghost, cu, style="dashed", alpha=0.8, lw=2.0)
#                         draw_seg((x2, 0),  ghost, cu, style="dashed", alpha=0.8, lw=2.0)
#                     else:                # v 在上边界
#                         draw_seg((x2, y2), ghost, cu, style="dashed", alpha=0.8, lw=2.0)
#                         draw_seg((x1, 0),  ghost, cu, style="dashed", alpha=0.8, lw=2.0)
#                 else:
#                     draw_seg((x1, y1), (x2, y2), cu)

#         # 标签
#         for v in S:
#             x, y = base_pos[v]
#             ax.text(x, y, str(id_tup(v, dims)), fontsize=9,
#                     ha="center", va="center", zorder=4)

#         # 轴范围：把“幽灵坐标”也放进去，刚好看到闭合效果
#         ax.set_xlim(-0.5, na + 0.5)
#         ax.set_ylim(-0.5, nb + 0.5)

#         # 图例
#         handles = [mlines.Line2D([], [], color=palette[i], lw=3, label=f"color {i+1}")
#                    for i in range(min(K, 12))]
#         if handles:
#             ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0))

#         ax.set_aspect("equal")
#         ax.set_axis_off()
#         plt.tight_layout()
#         plt.title(f"Slice at axis {slice_axis} = {s}  (lower-dim copy; dashed = wrap)")
#         plt.show()

def plot_lower_dim_slices(dims, edges, colors, K, slice_axis=-1):
    """
    低一维切片：回环边用“连续曲线（周期化贝塞尔）”跨边界再进入，
    视觉上是 C_{n_a} x C_{n_b} 的闭合网格而非 P_n。
    """
    import matplotlib.pyplot as plt, matplotlib.lines as mlines
    import numpy as np

    d = len(dims)
    if d < 3:
        print("[Slices] skipped: dims has no lower 1-D copies."); return
    if slice_axis < 0: slice_axis += d
    assert 0 <= slice_axis < d

    Vn = prod(dims)

    # (u,v)->eid
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

    # —— 工具：画“周期化贝塞尔曲线” —— #
    def draw_periodic_curve(ax, p, q, color, box=(na, nb), bend=0.28, lw=2.4, alpha=0.98):
        """
        p,q 在基本方块坐标系内（0..na, 0..nb）。
        若需要回环，内部自动把 q 平移 ±box 以取最近路径，
        在展开平面生成二次贝塞尔 B(t)，再按模 box 投影回基本方块。
        期间如遇跨边界跳变，自动分段绘制，从视觉上是一条连续曲线。
        """
        na_, nb_ = box
        px, py = p
        qx, qy = q

        # 选择展开平面中的 q'：对每个维度，允许平移 ±box 以最小化距离
        choices_x = [qx- na_, qx, qx+ na_]
        choices_y = [qy- nb_, qy, qy+ nb_]
        best = None; best_dist = 1e9
        for xx in choices_x:
            for yy in choices_y:
                dist = (xx-px)**2 + (yy-py)**2
                if dist < best_dist:
                    best_dist = dist; best = (xx, yy)
        qx_u, qy_u = best  # 展开平面中的 q'

        # 选一个“弯曲方向”的控制点（与直线法向）：
        vx, vy = (qx_u - px, qy_u - py)
        L = (vx**2 + vy**2)**0.5 + 1e-9
        nx, ny = -vy / L, vx / L  # 法向单位向量
        cx, cy = (px + qx_u) / 2 + bend * nx, (py + qy_u) / 2 + bend * ny

        # 采样贝塞尔曲线（展开平面）
        t = np.linspace(0, 1, 60)
        Xu = (1-t)**2 * px + 2*(1-t)*t * cx + t**2 * qx_u
        Yu = (1-t)**2 * py + 2*(1-t)*t * cy + t**2 * qy_u

        # 投影回基本方块，并在“跳变”处自动分段
        X = np.mod(Xu, na_)
        Y = np.mod(Yu, nb_)
        seg_x, seg_y = [X[0]], [Y[0]]
        # 阈值：若相邻采样点投影后距离 > box/2，视为跨边界，开新段
        thr = max(na_, nb_) * 0.35
        for i in range(1, len(X)):
            dx = abs(X[i] - X[i-1])
            dy = abs(Y[i] - Y[i-1])
            # 考虑环上最近距离
            dx = min(dx, na_ - dx)
            dy = min(dy, nb_ - dy)
            if (dx*dx + dy*dy) ** 0.5 > thr:
                ax.plot(seg_x, seg_y, color=color, lw=lw, alpha=alpha, solid_capstyle="round")
                seg_x, seg_y = [X[i]], [Y[i]]
            else:
                seg_x.append(X[i]); seg_y.append(Y[i])
        if len(seg_x) > 1:
            ax.plot(seg_x, seg_y, color=color, lw=lw, alpha=alpha, solid_capstyle="round")

    for s in range(dims[slice_axis]):
        # 该切片的顶点与边
        S = {v for v in range(Vn) if id_tup(v, dims)[slice_axis] == s}
        e_slice = [(u, v) for (u, v) in edges if (u in S and v in S)]

        # 基本方块中的网格坐标
        base_pos = {}
        for v in S:
            t = id_tup(v, dims)
            base_pos[v] = (t[a], t[b])

        fig, ax = plt.subplots(figsize=(6.4, 6.2))

        # 先画边：所有边都用曲线（非回环也能轻微弯曲，便于风格统一）
        for (u, v) in e_slice:
            cu = palette[(col_of(u, v) - 1) % K]
            p = base_pos[u]; q = base_pos[v]
            draw_periodic_curve(ax, p, q, cu, box=(na, nb), bend=0.26, lw=2.4, alpha=0.98)

        # 再画节点与标签
        xs = [base_pos[v][0] for v in S]; ys = [base_pos[v][1] for v in S]
        ax.scatter(xs, ys, s=260, facecolors="white", edgecolors="black", zorder=3)
        for v in S:
            x, y = base_pos[v]
            ax.text(x, y, str(id_tup(v, dims)), fontsize=9, ha="center", va="center", zorder=4)

        # 边界虚线框，仅做参照（可去掉）
        ax.plot([0, na, na, 0, 0], [0, 0, nb, nb, 0], color="#999999", lw=0.8, ls="dotted", alpha=0.8)

        # 轴范围覆盖到边界（不额外加幽灵格）
        ax.set_xlim(-0.5, na + 0.5)
        ax.set_ylim(-0.5, nb + 0.5)
        ax.set_aspect("equal")
        ax.set_axis_off()

        # 图例
        handles = [mlines.Line2D([], [], color=palette[i], lw=3, label=f"color {i+1}")
                   for i in range(min(K, 12))]
        if handles:
            ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0))

        plt.tight_layout()
        plt.title(f"Slice at axis {slice_axis} = {s}  (periodic curves on torus)")
        plt.show()


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
    if PLOT_SLICES_LOWER:
        plot_lower_dim_slices(dims, info["edges"], colors, Kc,
                              slice_axis = (len(dims)-1) if SLICE_AXIS==-1 else SLICE_AXIS)

if __name__=="__main__":
    main()

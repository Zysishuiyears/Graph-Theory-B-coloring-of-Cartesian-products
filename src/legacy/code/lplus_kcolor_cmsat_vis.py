#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
鏋勯€犲湀涔樺湀绗涘崱灏旂Н鍥?G = C_{n1} x ... x C_{nd} 鐨勬墿寮犵嚎鍥?L^+(G)锛?
鐢?SAT (cryptominisat5) 姹備竴涓?k-鐐规煋鑹诧紝骞跺鍑猴細
- L^+(G) 鐨?Gephi CSV
- 鍘熷浘 G 鐨?Gephi CSV
- L^+(G) 椤剁偣搴︾粺璁?
- L^+(G) 鐨?PNG 鍥撅紙鑺傜偣鎸夐鑹蹭笂鑹诧級

鎺ㄨ崘鏀惧湪: code/lplus_kcolor_cmsat.py
"""

import os
import subprocess
import tempfile
from itertools import product
from collections import Counter
from datetime import datetime

import matplotlib.pyplot as plt
import networkx as nx
import matplotlib as mpl

# =============== 鍦ㄨ繖閲屾敼鍙傛暟 ===============

# 绗涘崱灏旂Н鍦堢殑闀垮害锛屾瘮濡?[3, 5] 灏辨槸 C3 x C5
LENGTHS = [3, 3, 3, 3]

# 鐐规煋鑹茬殑棰滆壊鏁帮紙浣犺娴嬭瘯鐨?k锛?
K = 9

# 瀹為獙鍚嶅瓧锛岀敤鏉ュ尯鍒嗕笉鍚屽疄楠?
EXPERIMENT_NAME = "C3xC3xC3xC3_k9"

# 鏄惁鎶?CNF 鏂囦欢淇濆瓨鍒?results/runs/legacy_2026/_early_outputs/ 鐩綍锛圱rue/False锛?
WRITE_CNF_FILE = True
CNF_FILENAME = "C3xC3xC3xC3_k9.cnf"   # 鍙湪 WRITE_CNF_FILE=True 鏃剁敤

# SAT 姹傝В鍣ㄨ缃細鐩墠灏辨槸鐢?cryptominisat5
USE_SOLVER = "cmsat"
CMSAT_PATH = r"E:\cryptominisat5-win\cryptominisat5.exe"  # 鎸変綘鑷繁鐨勮矾寰勬敼


# =============== 璺緞璁剧疆锛氬熀浜庨」鐩粨鏋?===============

# 褰撳墠鑴氭湰鎵€鍦ㄧ洰褰曪細.../code
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
# 椤圭洰鏍圭洰褰曪細鍋囪鏄?code 鐨勪笂涓€灞?
ROOT_DIR = os.path.abspath(os.path.join(_THIS_DIR, ".."))
DATA_DIR = os.path.join(ROOT_DIR, "results", "runs", "legacy_2026", "_early_outputs")
RESULTS_BASE_DIR = os.path.join(ROOT_DIR, "results")

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_BASE_DIR, exist_ok=True)

# 甯︽椂闂存埑鐨勫疄楠岀洰褰?
_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
RESULTS_DIR = os.path.join(RESULTS_BASE_DIR, f"{_timestamp}_{EXPERIMENT_NAME}")
os.makedirs(RESULTS_DIR, exist_ok=True)


# ----------------------------------------------------------------------
# 1. 鏋勯€?C_{n1} x ... x C_{nd} 鍙婂叾鎵╁紶绾垮浘 L^+(G)
# ----------------------------------------------------------------------

def build_cycle_product_graph(lengths):
    """
    鏋勯€犵瑳鍗″皵绉浘 G = C_{n1} x ... x C_{nd}.

    杩斿洖
    ----
    vertices : list[tuple[int]]
        椤剁偣鍧愭爣鍒楄〃锛屾瘡涓《鐐规槸 d 缁村厓缁?(i1,...,id).
        椤剁偣缂栧彿涓哄叾鍦ㄥ垪琛ㄤ腑鐨勭储寮曪細0..|V|-1.

    edges : list[tuple[int,int]]
        鏃犲悜杈瑰垪琛?(u, v)锛屼繚璇?u < v.

    edge_index : dict[(int,int) -> int]
        浠?(min(u,v), max(u,v)) 鏄犲皠鍒拌竟缂栧彿 e_id.
    """
    d = len(lengths)
    vertices = list(product(*[range(n) for n in lengths]))
    vid = {v: i for i, v in enumerate(vertices)}

    edges = []
    edge_index = {}
    for v in vertices:
        v_id = vid[v]
        for dim in range(d):
            w = list(v)
            w[dim] = (w[dim] + 1) % lengths[dim]
            w = tuple(w)
            w_id = vid[w]
            a, b = (v_id, w_id) if v_id < w_id else (w_id, v_id)
            key = (a, b)
            if key not in edge_index:
                edge_index[key] = len(edges)
                edges.append(key)

    return vertices, edges, edge_index


def build_expanded_line_graph(lengths):
    """
    鏋勯€犵瑳鍗″皵绉浘 G = C_{n1} x ... x C_{nd} 鐨勬墿寮犵嚎鍥?L^+(G).

    L^+(G) 鐨勯《鐐?= G 鐨勮竟锛?
    L^+(G) 鐨勮竟 = 绾垮浘 L(G) 鐨勮竟 + 鎵€鏈?4-鍦堜腑瀵硅竟涔嬮棿鐨勮竟銆?

    杩斿洖
    ----
    vertices : list[tuple[int]]
        G 鐨勯《鐐瑰潗鏍囧垪琛紙缂栧彿 0..|V|-1锛?

    edges_G : list[tuple[int,int]]
        G 鐨勮竟鍒楄〃锛堢紪鍙?0..|E|-1锛夛紝涔熸槸 L^+(G) 鐨勯《鐐归泦鍚堛€?

    adj_Lplus : list[set[int]]
        L^+(G) 鐨勯偦鎺ヨ〃锛岄暱搴︿负 |E|.
        adj_Lplus[e_id] 鏄笌 e_id 鐩搁偦鐨?e'_id 闆嗗悎銆?
    """
    vertices, edges_G, edge_index = build_cycle_product_graph(lengths)
    nV = len(vertices)
    nE = len(edges_G)
    d = len(lengths)

    adj_Lplus = [set() for _ in range(nE)]

    # 姣忎釜椤剁偣鍦?G 涓殑 incident edges
    incident = [[] for _ in range(nV)]
    for e_id, (u, v) in enumerate(edges_G):
        incident[u].append(e_id)
        incident[v].append(e_id)

    # 1) 绾垮浘 L(G) 閮ㄥ垎锛氬叡绔偣杈?-> 鐩搁偦
    for inc in incident:
        m = len(inc)
        for i in range(m):
            e1 = inc[i]
            for j in range(i + 1, m):
                e2 = inc[j]
                adj_Lplus[e1].add(e2)
                adj_Lplus[e2].add(e1)

    # 2) F(G) 閮ㄥ垎锛氭瘡涓?4-鍦堜腑鐨勫杈?-> 鐩搁偦
    vid = {v: i for i, v in enumerate(vertices)}

    def get_eid(a, b):
        key = (a, b) if a < b else (b, a)
        return edge_index[key]

    # 閬嶅巻姣忎釜椤剁偣 v0 鍜屾瘡涓€瀵圭淮搴?(i, j)锛屾瀯閫?4-鍦堬細
    for v0 in vertices:
        v0_id = vid[v0]
        for i in range(d):
            for j in range(i + 1, d):
                # v1 = v0 + e_i
                v1 = list(v0)
                v1[i] = (v1[i] + 1) % lengths[i]
                v1 = tuple(v1)
                # v3 = v0 + e_j
                v3 = list(v0)
                v3[j] = (v3[j] + 1) % lengths[j]
                v3 = tuple(v3)
                # v2 = v1 + e_j = v3 + e_i
                v2 = list(v1)
                v2[j] = (v2[j] + 1) % lengths[j]
                v2 = tuple(v2)

                v1_id = vid[v1]
                v2_id = vid[v2]
                v3_id = vid[v3]

                e_a = get_eid(v0_id, v1_id)
                e_b = get_eid(v1_id, v2_id)
                e_c = get_eid(v3_id, v2_id)
                e_d = get_eid(v0_id, v3_id)

                # 瀵硅竟瀵癸細(e_a, e_c), (e_b, e_d)
                adj_Lplus[e_a].add(e_c)
                adj_Lplus[e_c].add(e_a)
                adj_Lplus[e_b].add(e_d)
                adj_Lplus[e_d].add(e_b)

    return vertices, edges_G, adj_Lplus


# ----------------------------------------------------------------------
# 2. SAT 缂栫爜锛歀^+(G) 鐨?k-鐐规煋鑹?
# ----------------------------------------------------------------------

def build_k_coloring_cnf(adj, k):
    """
    缁欏畾鍥?H 鐨勯偦鎺ヨ〃 adj锛屾瀯閫?k-鐐规煋鑹茬殑 CNF 缂栫爜銆?

    鍙橀噺锛歺_{v,c} 琛ㄧず椤剁偣 v 鍙栭鑹?c (c = 0..k-1).
    缂栧彿锛歷ar(v,c) = v*k + c + 1  锛圖IMACS 浠?1 寮€濮嬶級銆?

    绾︽潫锛?
    1) 姣忎釜椤剁偣鑷冲皯涓€绉嶉鑹诧細
       (x_{v,0} 鈭?x_{v,1} 鈭?... 鈭?x_{v,k-1})
    2) 姣忎釜椤剁偣鑷冲涓€绉嶉鑹诧細
       瀵规墍鏈?0 <= c1 < c2 < k:
       (卢x_{v,c1} 鈭?卢x_{v,c2})
    3) 鐩搁偦椤剁偣涓嶈兘鍚岃壊锛?
       瀵规瘡鏉¤竟 (u,v)锛屽姣忎釜棰滆壊 c锛?
       (卢x_{u,c} 鈭?卢x_{v,c})

    杩斿洖
    ----
    clauses : list[list[int]]
    num_vars : int
    """
    n = len(adj)
    num_vars = n * k

    def var(v, c):
        return v * k + c + 1

    clauses = []

    # 1) 姣忎釜椤剁偣鑷冲皯涓€绉嶉鑹?
    for v in range(n):
        clauses.append([var(v, c) for c in range(k)])

    # 2) 姣忎釜椤剁偣鑷冲涓€绉嶉鑹?
    for v in range(n):
        for c1 in range(k):
            for c2 in range(c1 + 1, k):
                clauses.append([-var(v, c1), -var(v, c2)])

    # 3) 鐩搁偦椤剁偣涓嶈兘鍚岃壊
    for v in range(n):
        for u in adj[v]:
            if u > v:  # 鏃犲悜鍥撅紝鍙紪鐮佷竴娆?
                for c in range(k):
                    clauses.append([-var(v, c), -var(u, c)])

    return clauses, num_vars


def write_cnf(clauses, num_vars, path):
    """
    灏?CNF 瀛愬彞鍐欐垚 DIMACS 鏍煎紡鏂囦欢銆?
    """
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"p cnf {num_vars} {len(clauses)}\n")
        for cl in clauses:
            line = " ".join(str(lit) for lit in cl) + " 0\n"
            f.write(line)


# ----------------------------------------------------------------------
# 3. 璋冪敤 cryptominisat5 姹傝В & 瑙ｆ瀽妯″瀷
# ----------------------------------------------------------------------

def run_cmsat(cnf_path, cmsat_path=CMSAT_PATH):
    """
    璋冪敤 cryptominisat5 姹傝В CNF 鏂囦欢锛岃繑鍥炴ā鍨嬩腑鐨勬暣鍨嬫枃瀛楀垪琛ㄣ€?

    杩斿洖
    ----
    model_lits : list[int] 鎴?None
        鑻?SAT锛岃繑鍥炴ā鍨嬩腑鍑虹幇鐨勬枃瀛楋紙鍖呭惈姝ｈ礋锛夛紱鑻?UNSAT锛岃繑鍥?None銆?
    """
    if not os.path.isfile(cmsat_path):
        raise FileNotFoundError(f"鎵句笉鍒?cryptominisat5 鍙墽琛屾枃浠? {cmsat_path}")

    cmd = [cmsat_path, cnf_path]

    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    out = proc.stdout.splitlines()
    sat = None
    model_lits = []
    for line in out:
        line = line.strip()
        if not line:
            continue
        if line.startswith('s '):
            if 'UNSAT' in line:
                sat = False
            elif 'SAT' in line:
                sat = True
        elif line.startswith('v ') or line.startswith('V '):
            parts = line.split()[1:]
            for p in parts:
                if p == '0':
                    continue
                try:
                    model_lits.append(int(p))
                except ValueError:
                    pass

    if sat is False:
        return None
    if sat is None:
        raise RuntimeError("cryptominisat5 杈撳嚭鏈寘鍚?SAT/UNSAT 淇℃伅锛岃妫€鏌ョ増鏈拰鍙傛暟銆?)

    return model_lits


def decode_model_to_colors(model_lits, n_vertices, k):
    """
    灏?SAT 妯″瀷涓殑鏂囧瓧鍒楄〃瑙ｆ瀽涓烘瘡涓《鐐圭殑棰滆壊銆?

    缂栫爜绾﹀畾锛歷ar(v,c) = v*k + c + 1
    鍙湅姝ｆ枃瀛椼€?
    """
    colors = [-1] * n_vertices
    for lit in model_lits:
        if lit <= 0:
            continue
        idx = lit - 1
        v = idx // k
        c = idx % k
        if v < 0 or v >= n_vertices:
            continue
        if colors[v] != -1 and colors[v] != c:
            raise RuntimeError(f"妯″瀷瀵归《鐐?{v} 缁欎簡澶氫釜棰滆壊: {colors[v]}, {c}")
        colors[v] = c

    for v in range(n_vertices):
        if colors[v] == -1:
            raise RuntimeError(f"妯″瀷涓《鐐?{v} 娌℃湁鍒嗛厤棰滆壊銆?)

    return colors


# ----------------------------------------------------------------------
# 4. 妫€鏌ョ偣鏌撹壊姝ｇ‘鎬э紙鐩搁偦椤剁偣涓嶅悓鑹诧級
# ----------------------------------------------------------------------

def check_coloring_correct(adj, colors, k):
    """
    妫€鏌ョ偣鏌撹壊鏄惁姝ｇ‘锛堢浉閭婚《鐐逛笉鍚岃壊锛岄鑹插湪 0..k-1锛夈€?

    鑻ラ敊璇紝鎶涘嚭 RuntimeError銆?
    """
    n = len(adj)
    if len(colors) != n:
        raise RuntimeError(f"棰滆壊鏁扮粍闀垮害 {len(colors)} 涓庨《鐐规暟 {n} 涓嶇銆?)

    for v, c in enumerate(colors):
        if not (0 <= c < k):
            raise RuntimeError(f"椤剁偣 {v} 鎷垮埌浜嗛潪娉曢鑹?{c}銆?)

    for v in range(n):
        for u in adj[v]:
            if u > v and colors[u] == colors[v]:
                raise RuntimeError(
                    f"鏌撹壊閿欒锛氶《鐐?{v} 鍜?{u} 鐩搁偦浣嗗悓鑹?{colors[v]}銆?
                )


# ----------------------------------------------------------------------
# 5. 搴﹀垎甯冪粺璁?& 瀵煎嚭 Gephi CSV & 缁樺浘
# ----------------------------------------------------------------------

def summarize_degrees(adj, out_path):
    """
    缁熻 L^+(G) 鐨勯《鐐瑰害鎯呭喌锛屽苟鍐欏埌鏂囨湰鏂囦欢銆?

    鍐呭鍖呮嫭锛?
    - 椤剁偣鏁?
    - min / max / 骞冲潎搴︽暟
    - 鍚勪釜搴︽暟鍑虹幇娆℃暟
    """
    degrees = [len(nei) for nei in adj]
    n = len(degrees)
    min_deg = min(degrees)
    max_deg = max(degrees)
    avg_deg = sum(degrees) / n if n > 0 else 0.0
    counter = Counter(degrees)

    with open(out_path, "w", encoding="utf-8") as f:
        f.write(f"# L^+(G) 椤剁偣搴︾粺璁n")
        f.write(f"n_vertices = {n}\n")
        f.write(f"min_degree = {min_deg}\n")
        f.write(f"max_degree = {max_deg}\n")
        f.write(f"avg_degree = {avg_deg:.4f}\n")
        f.write("degree_counts:\n")
        for deg in sorted(counter.keys()):
            f.write(f"  degree {deg}: {counter[deg]}\n")


def export_lplus_gephi(vertices, edges_G, adj_Lplus, colors, out_prefix):
    """
    瀵煎嚭 L^+(G) 鐨勭偣鏌撹壊缁撴灉涓?Gephi 鍙鐨?CSV.

    杈撳嚭锛?
        <out_prefix>_Lplus_nodes.csv : id;label;color
        <out_prefix>_Lplus_edges.csv : source;target
    """
    nodes_path = out_prefix + "_Lplus_nodes.csv"
    edges_path = out_prefix + "_Lplus_edges.csv"

    with open(nodes_path, "w", encoding="utf-8") as f:
        f.write("id;label;color\n")
        for e_id, (u, v) in enumerate(edges_G):
            label = f"{vertices[u]}-{vertices[v]}"
            color = colors[e_id] if colors is not None else ""
            f.write(f"{e_id};\"{label}\";{color}\n")

    with open(edges_path, "w", encoding="utf-8") as f:
        f.write("source;target\n")
        n = len(adj_Lplus)
        for v in range(n):
            for u in adj_Lplus[v]:
                if u > v:
                    f.write(f"{v};{u}\n")


def export_original_gephi(vertices, edges_G, colors, out_prefix):
    """
    瀵煎嚭鍘熷浘 G 鐨勭粨鏋?+ B-杈规煋鑹诧紙鐢?L^+(G) 鐨勭偣鏌撹壊杞洖锛夈€?

    杈撳嚭锛?
        <out_prefix>_G_nodes.csv : id;label
        <out_prefix>_G_edges.csv : source;target;color
    """
    nodes_path = out_prefix + "_G_nodes.csv"
    edges_path = out_prefix + "_G_edges.csv"

    with open(nodes_path, "w", encoding="utf-8") as f:
        f.write("id;label\n")
        for v_id, coord in enumerate(vertices):
            f.write(f"{v_id};\"{coord}\"\n")

    with open(edges_path, "w", encoding="utf-8") as f:
        f.write("source;target;color\n")
        for e_id, (u, v) in enumerate(edges_G):
            color = colors[e_id] if colors is not None else ""
            f.write(f"{u};{v};{color}\n")


def plot_lplus_graph(adj_Lplus, colors, out_path, figsize=(8, 8)):
    """
    浣跨敤 networkx + matplotlib 灏?L^+(G) 鐢诲嚭鏉ュ苟淇濆瓨涓?PNG銆?

    - 鑺傜偣缂栧彿涓?0..|E|-1
    - 鑻ユ彁渚涗簡 colors锛屽垯鐢?colormap 鎸夐鑹叉秱鐐癸紝骞跺姞涓婇鑹插浘娉紙colorbar锛?
    """
    n = len(adj_Lplus)
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for v in range(n):
        for u in adj_Lplus[v]:
            if u > v:
                G.add_edge(v, u)

    # 寤虹珛 figure / axes锛屽悗闈?colorbar 浼氱敤鍒?ax
    fig, ax = plt.subplots(figsize=figsize)

    # spring layout锛屽浐瀹氶殢鏈虹瀛愭柟渚垮鐜?
    pos = nx.spring_layout(G, seed=42)

    if colors is not None and len(colors) == n:
        max_color = max(colors) if colors else 0
        # 鎺ㄨ崘鐨勮幏鍙?colormap 鐨勬柟寮?
        cmap = plt.get_cmap("tab10", max_color + 1)
        node_colors = [colors[v] for v in G.nodes()]

        # 鐢荤偣锛屽苟鎷垮埌 PathCollection 浣滀负 mappable
        nodes = nx.draw_networkx_nodes(
            G, pos,
            node_size=40,
            node_color=node_colors,
            cmap=cmap,
            ax=ax,
        )
        # 鐢昏竟
        nx.draw_networkx_edges(G, pos, alpha=0.3, width=0.5, ax=ax)

        # 鐩存帴鐢?nodes 浣滀负 mappable锛屽姞 colorbar锛屽苟鎸囧畾 ax
        cbar = fig.colorbar(nodes, ax=ax)
        cbar.set_ticks(range(0, max_color + 1))
        cbar.set_ticklabels([str(i) for i in range(0, max_color + 1)])
        cbar.set_label("vertex color (0..k-1)")
    else:
        # 鏃犻鑹蹭俊鎭椂鐢ㄧ粺涓€棰滆壊锛屼笉鐢?colorbar
        nx.draw_networkx_nodes(
            G, pos,
            node_size=40,
            node_color="lightgray",
            ax=ax,
        )
        nx.draw_networkx_edges(G, pos, alpha=0.3, width=0.5, ax=ax)

    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


# ----------------------------------------------------------------------
# 6. 涓绘祦绋嬶細鐩存帴鐢ㄤ笂闈㈢殑鍏ㄥ眬鍙傛暟
# ----------------------------------------------------------------------

def main():
    lengths = LENGTHS
    k = K

    print(f"[INFO] 鏋勯€?G = " + " x ".join(f"C{n}" for n in lengths))
    vertices, edges_G, adj_Lplus = build_expanded_line_graph(lengths)
    print(f"[INFO] |V(G)| = {len(vertices)}, |E(G)| = {len(edges_G)}")
    print(f"[INFO] L^+(G) 椤剁偣鏁?= {len(edges_G)} (姣忔潯杈瑰搴斾竴涓《鐐?")

    # 搴﹀垎甯冪粺璁?
    degree_summary_path = os.path.join(RESULTS_DIR, "Lplus_degree_summary.txt")
    summarize_degrees(adj_Lplus, degree_summary_path)
    print(f"[INFO] 椤剁偣搴︾粺璁″啓鍏? {degree_summary_path}")

    # 鏋勯€?CNF
    print(f"[INFO] 鏋勯€?L^+(G) 鐨?{k}-鐐规煋鑹?CNF 缂栫爜...")
    clauses, num_vars = build_k_coloring_cnf(adj_Lplus, k)
    print(f"[INFO] CNF: 鍙橀噺鏁?= {num_vars}, 瀛愬彞鏁?= {len(clauses)}")

    # CNF 鏂囦欢璺緞
    if WRITE_CNF_FILE and CNF_FILENAME:
        cnf_path = os.path.join(DATA_DIR, CNF_FILENAME)
    else:
        fd, cnf_path = tempfile.mkstemp(suffix=".cnf", prefix="lplus_")
        os.close(fd)

    write_cnf(clauses, num_vars, cnf_path)
    print(f"[INFO] CNF 鍐欏叆: {cnf_path}")

    colors = None

    if USE_SOLVER == "cmsat":
        print(f"[INFO] 璋冪敤 CryptoMiniSat5 姹傝В ({CMSAT_PATH})...")
        model_lits = run_cmsat(cnf_path, CMSAT_PATH)
        if model_lits is None:
            print(f"[INFO] 鏈壘鍒?{k}-鐐规煋鑹诧紙UNSAT锛夈€?)
        else:
            print("[INFO] SAT锛岃В鏋愭ā鍨嬩负椤剁偣棰滆壊...")
            colors = decode_model_to_colors(model_lits, n_vertices=len(adj_Lplus), k=k)
            print("[INFO] 妫€鏌ユ煋鑹叉纭€э紙鐩搁偦椤剁偣涓嶅悓鑹诧級...")
            check_coloring_correct(adj_Lplus, colors, k)
            print("[INFO] 鏌撹壊妫€鏌ラ€氳繃銆?)
    else:
        print("[INFO] 鏈惎鐢?SAT 瑙ｇ畻锛屽彧鐢熸垚 CNF銆?)

    # 瀵煎嚭 Gephi & 鐢诲浘
    if colors is not None:
        prefix = os.path.join(RESULTS_DIR, "graph")
        print(f"[INFO] 瀵煎嚭 Gephi CSV 鍒?{RESULTS_DIR} ...")
        export_lplus_gephi(vertices, edges_G, adj_Lplus, colors, prefix)
        export_original_gephi(vertices, edges_G, colors, prefix)
        print(f"[INFO] L^+(G) : graph_Lplus_nodes.csv / graph_Lplus_edges.csv")
        print(f"[INFO] 鍘熷浘 G : graph_G_nodes.csv / graph_G_edges.csv")

        # 鐢诲浘
        png_path = os.path.join(RESULTS_DIR, "Lplus_graph.png")
        print(f"[INFO] 缁樺埗 L^+(G) 鍥惧儚鍒?{png_path} ...")
        plot_lplus_graph(adj_Lplus, colors, png_path)
        plot_lplus_graph(adj_Lplus, colors, png_path)

    else:
        print("[INFO] 娌℃湁棰滆壊瑙ｏ紝璺宠繃 Gephi 瀵煎嚭鍜岀粯鍥俱€?)


if __name__ == "__main__":
    main()




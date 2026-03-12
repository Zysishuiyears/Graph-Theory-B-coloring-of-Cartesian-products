#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
鏋勯€?2-鍦堢瑳鍗″皵绉浘 G = C_m 鈻?C_n 鍙婂叾鎵╁紶绾垮浘 L^+(G)锛?
骞跺皢缁撴灉瀵煎嚭鍒?B-coloring 宸ョ▼鐨勬爣鍑嗙洰褰曠粨鏋勪腑銆?

宸ョ▼鐩綍缁撴瀯鍋囪涓猴細
  project_root/
    code/
    results/
    results/
    scripts/
    readme.txt

鏈剼鏈墍鍦ㄤ綅缃細project_root/scripts/extended_line_graph_cmn.py

杈撳嚭鐩綍锛?
  榛樿锛歱roject_root/results/YYYY-MM-DD_瀹為獙鍚?
    - Lplus_Cm_Cn.gexf       鎵╁紶绾垮浘 L^+(G)
    - base_Cm_Cn.gexf        鍘熷浘 G锛堜粎褰撴寚瀹?--export-base 鏃讹級

浣犱篃鍙互鐢?--out-dir 鑷繁鎸囧畾杈撳嚭鐩綍锛堢浉瀵硅矾寰勪細鎸夊伐绋嬫牴鐩綍瑙ｆ瀽锛夈€?
"""

import argparse
from datetime import date
from pathlib import Path

import networkx as nx


# ====================== 鍥炬瀯閫犻儴鍒?======================

def build_cartesian_cycles(m: int, n: int) -> nx.Graph:
    """
    鏋勯€?G = C_m 鈻?C_n锛岄《鐐圭敤 (i, j) 琛ㄧず锛? <= i < m, 0 <= j < n銆?
    杈规寜鐜粫鏂瑰紡杩炴帴锛堝嵆鍦?m脳n 鐨勭綉鏍间笂褰㈡垚鈥滅幆闈㈢綉鏍尖€濓級銆?
    """
    G = nx.Graph()

    # 娣诲姞椤剁偣骞堕檮涓婂潗鏍囧睘鎬э紙鏂逛究涔嬪悗鍦?Gephi 閲岃皟甯冨眬锛?
    for i in range(m):
        for j in range(n):
            G.add_node((i, j), x=float(i), y=float(j))

    # 娣诲姞姘村钩鍜屽瀭鐩磋竟锛堟敞鎰忕幆缁曪級
    for i in range(m):
        for j in range(n):
            u = (i, j)
            v_h = ((i + 1) % m, j)       # 姘村钩鐩搁偦
            v_v = (i, (j + 1) % n)       # 鍨傜洿鐩搁偦
            G.add_edge(u, v_h)
            G.add_edge(u, v_v)

    return G


def _norm_edge(u, v):
    """缁欏畾涓や釜绔偣锛岃繑鍥炰竴涓棤鍚戣竟鐨勮鑼冭〃绀?(u, v)锛堟寜瀛楀吀搴忔帓搴忥級銆?""
    return tuple(sorted((u, v)))


def extended_line_graph_Cmn(m: int, n: int) -> nx.Graph:
    """
    鏋勯€?C_m 鈻?C_n 鐨勬墿寮犵嚎鍥?L^+(G)銆?

    鑺傜偣锛氭瘡鏉¤竟 e 鈭?E(G)锛岀敤 (u, v) 鐨勮鑼冨寲鍏冪粍浣滀负鑺傜偣鍚嶃€?
    杈癸細
      - 绾垮浘閮ㄥ垎锛氬鏋?e, f 鍦?G 涓叡浜鐐癸紝鍒欏湪 L^+(G) 涓浉杩烇紱
      - 4-鍦堝杈归儴鍒嗭細鑻?e, f 鏄煇涓?4-cycle 涓殑瀵硅竟锛屽垯鍦?L^+(G) 涓浉杩炪€?
    """
    G = build_cartesian_cycles(m, n)
    Lp = nx.Graph()

    # 1. 涓烘瘡鏉¤竟寤轰竴涓《鐐?
    for u, v in G.edges():
        e = _norm_edge(u, v)
        (i1, j1), (i2, j2) = u, v

        # 鐢ㄤ腑鐐逛綔涓鸿妭鐐瑰嚑浣曚綅缃紙鐜粫澶勮櫧鐒朵細鈥滄姌鍥炩€濓紝浣嗗 Gephi 涓嶅奖鍝嶏級
        x = (i1 + i2) / 2.0
        y = (j1 + j2) / 2.0
        edge_type = "horizontal" if i1 != i2 else "vertical"

        Lp.add_node(
            e,
            u=str(u),
            v=str(v),
            x=float(x),
            y=float(y),
            edge_type=edge_type,
        )

    # 2. 绾垮浘閮ㄥ垎锛氭瘡涓師鍥鹃《鐐圭殑 incident edges 褰㈡垚涓€涓洟
    incident = {}  # vertex -> list of incident edges (瑙勮寖鍖栧舰寮?
    for u, v in G.edges():
        e = _norm_edge(u, v)
        for w in (u, v):
            incident.setdefault(w, []).append(e)

    for edges_at_vertex in incident.values():
        k = len(edges_at_vertex)
        for i in range(k):
            for j in range(i + 1, k):
                e1 = edges_at_vertex[i]
                e2 = edges_at_vertex[j]
                if e1 != e2:
                    Lp.add_edge(e1, e2, reason="line")

    # 3. 4-鍦堝杈归儴鍒?
    #   姣忎釜鈥滄牸瀛愨€濆氨鏄竴涓?4-鍦?
    #   (i,j)-(i+1,j)-(i+1,j+1)-(i,j+1)-(i,j)
    for i in range(m):
        for j in range(n):
            v00 = (i, j)
            v10 = ((i + 1) % m, j)
            v11 = ((i + 1) % m, (j + 1) % n)
            v01 = (i, (j + 1) % n)

            e1 = _norm_edge(v00, v10)  # 涓嬫按骞?
            e2 = _norm_edge(v10, v11)  # 鍙冲瀭鐩?
            e3 = _norm_edge(v01, v11)  # 涓婃按骞?
            e4 = _norm_edge(v00, v01)  # 宸﹀瀭鐩?

            # 姘村钩瀵硅竟
            if e1 in Lp and e3 in Lp:
                Lp.add_edge(e1, e3, reason="4cycle_opposite")

            # 鍨傜洿瀵硅竟
            if e2 in Lp and e4 in Lp:
                Lp.add_edge(e2, e4, reason="4cycle_opposite")

    return Lp


# ====================== 涓荤▼搴?& 璺緞閫昏緫 ======================

def main():
    parser = argparse.ArgumentParser(
        description="鏋勯€?C_m 鈻?C_n 鐨勬墿寮犵嚎鍥?L^+(G) 骞跺鍑哄埌 results/鏃ユ湡_瀹為獙鍚? 鐩綍銆?
    )
    parser.add_argument("m", type=int, help="绗竴涓湀 C_m 鐨勯暱搴?m (鈮?)")
    parser.add_argument("n", type=int, help="绗簩涓湀 C_n 鐨勯暱搴?n (鈮?)")

    parser.add_argument(
        "--exp-name",
        type=str,
        default=None,
        help="瀹為獙鍚嶏紙鐢ㄤ簬 results/YYYY-MM-DD_瀹為獙鍚?鐩綍鍚嶏紱榛樿 Lplus_Cm_Cn锛?,
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default=None,
        help="鑷畾涔夎緭鍑虹洰褰曪紙缁濆璺緞锛屾垨鐩稿宸ョ▼鏍圭洰褰曠殑璺緞锛?,
    )
    parser.add_argument(
        "--export-base",
        action="store_true",
        help="鍚屾椂瀵煎嚭鍘熷浘 C_m 鈻?C_n锛坆ase_Cm_Cn.gexf锛?,
    )

    args = parser.parse_args()
    m, n = args.m, args.n

    # --------- 纭畾宸ョ▼鏍圭洰褰曪紙project_root锛?----------
    # 鏈枃浠跺簲浣嶄簬 project_root/scripts/ 涓嬶紝鍥犳宸ョ▼鏍圭洰褰曟槸瀹冪殑鐖剁骇鐩綍
    project_root = Path(__file__).resolve().parents[1]
    results_root = project_root / "results"

    # --------- 纭畾杈撳嚭鐩綍 out_dir ----------
    if args.out_dir:
        # 鐢ㄦ埛鎸囧畾浜嗚緭鍑虹洰褰?
        out_dir = Path(args.out_dir)
        if not out_dir.is_absolute():
            # 鐩稿璺緞鎸夊伐绋嬫牴鐩綍瑙ｆ瀽
            out_dir = project_root / out_dir
    else:
        # 鑷姩鎸夋棩鏈?+ 瀹為獙鍚嶇粍缁囩洰褰?
        exp_name = args.exp_name or f"Lplus_C{m}_C{n}"
        date_str = date.today().isoformat()
        out_dir = results_root / f"{date_str}_{exp_name}"

    out_dir.mkdir(parents=True, exist_ok=True)

    # --------- 鏋勯€犳墿寮犵嚎鍥?----------
    Lp = extended_line_graph_Cmn(m, n)

    # 鏂囦欢鍚?
    lplus_path = out_dir / f"Lplus_C{m}_C{n}.gexf"
    nx.write_gexf(Lp, lplus_path)
    print(f"[OK] 鎵╁紶绾垮浘 L^+(C_{m} 鈻?C_{n}) 宸插鍑哄埌:\n  {lplus_path}")

    if args.export_base:
        G = build_cartesian_cycles(m, n)
        base_path = out_dir / f"base_C{m}_C{n}.gexf"
        nx.write_gexf(G, base_path)
        print(f"[OK] 鍘熷浘 C_{m} 鈻?C_{n} 宸插鍑哄埌:\n  {base_path}")

    print("\n鎻愮ず锛氬彲浠ュ湪 Gephi 涓洿鎺ユ墦寮€涓婅堪 .gexf 鏂囦欢杩涜鍙鍖栥€?)


if __name__ == "__main__":
    main()


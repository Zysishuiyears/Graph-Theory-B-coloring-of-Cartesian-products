import os
import sys
from typing import List, Tuple

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "code"))
from approxmc_cmsat import estimate_log2_model_count


# ========== 你的 CMSAT 路径 ==========
CMSAT_PATH = r"E:\cryptominisat5-win\cryptominisat5.exe"

# ========== 问题参数 ==========
m, n, k = 3, 7, 5

# 是否加入 “每个4-圈四色” 约束（你如果要就 True）
USE_C4_DISTINCT = False

# 近似计数参数（越大越稳也越慢）
TRIALS_PER_FIX = 17

# XOR 随机密度
P_IN_XOR = 0.5


# ---------- 边编号与变量 ----------
def num_edges(m, n): return 2*m*n

def edge_id_H(i, j, m, n): return i*n + j
def edge_id_V(i, j, m, n): return m*n + i*n + j

def var_color(e, c, k):  # 1..E*k
    return e*k + c + 1

def var_rep(e, c, E, k):  # E*k+1..2*E*k
    return E*k + e*k + c + 1


# ---------- incident & neighbors ----------
def incident_edges_at_vertex(i, j, m, n):
    return [
        edge_id_H(i, j, m, n),
        edge_id_H(i, (j-1) % n, m, n),
        edge_id_V(i, j, m, n),
        edge_id_V((i-1) % m, j, m, n),
    ]

def build_edge_neighbors(m, n):
    E = num_edges(m, n)
    neigh = [set() for _ in range(E)]
    for i in range(m):
        for j in range(n):
            inc = incident_edges_at_vertex(i, j, m, n)
            for a in range(len(inc)):
                for b in range(a+1, len(inc)):
                    e1, e2 = inc[a], inc[b]
                    neigh[e1].add(e2)
                    neigh[e2].add(e1)
    return [list(s) for s in neigh]

def enumerate_C4_edge_sets(m, n):
    C4s = []
    for i in range(m):
        for j in range(n):
            e_left  = edge_id_V(i, j, m, n)                 # (i,j)-(i+1,j)
            e_bot   = edge_id_H(i, j, m, n)                 # (i,j)-(i,j+1)
            e_right = edge_id_V(i, (j+1) % n, m, n)         # (i,j+1)-(i+1,j+1)
            e_top   = edge_id_H((i+1) % m, j, m, n)         # (i+1,j)-(i+1,j+1)
            C4s.append([e_left, e_bot, e_right, e_top])
    return C4s


# ---------- base CNF (NO symmetry breaking!) ----------
def build_base_bedge_bcolor_cnf(m, n, k, use_c4=False) -> Tuple[int, List[List[int]]]:
    E = num_edges(m, n)
    neigh = build_edge_neighbors(m, n)
    clauses: List[List[int]] = []

    # 1) each edge exactly one color
    for e in range(E):
        clauses.append([var_color(e, c, k) for c in range(k)])
        for c1 in range(k):
            for c2 in range(c1+1, k):
                clauses.append([-var_color(e, c1, k), -var_color(e, c2, k)])

    # 2) proper at vertices
    for i in range(m):
        for j in range(n):
            inc = incident_edges_at_vertex(i, j, m, n)
            for c in range(k):
                for a in range(len(inc)):
                    for b in range(a+1, len(inc)):
                        clauses.append([-var_color(inc[a], c, k), -var_color(inc[b], c, k)])

    # 3) optional: each 4-cycle uses all distinct colors (at most one per color)
    if use_c4:
        C4s = enumerate_C4_edge_sets(m, n)
        for cyc in C4s:
            for c in range(k):
                for a in range(4):
                    for b in range(a+1, 4):
                        clauses.append([-var_color(cyc[a], c, k), -var_color(cyc[b], c, k)])

    # 4) B-coloring dominance via rep(e,c)
    # (a) each color has at least one representative
    for c in range(k):
        clauses.append([var_rep(e, c, E, k) for e in range(E)])

    # (b) rep(e,c) -> color(e,c)
    for e in range(E):
        for c in range(k):
            clauses.append([-var_rep(e, c, E, k), var_color(e, c, k)])

    # (c) rep(e,c) -> for each d!=c exists neighbor f with color d
    for e in range(E):
        N = neigh[e]
        for c in range(k):
            for d in range(k):
                if d == c:
                    continue
                clause = [-var_rep(e, c, E, k)]
                clause += [var_color(f, d, k) for f in N]
                clauses.append(clause)

    nvars = 2 * E * k
    return nvars, clauses


# ---------- automorphisms of C_m □ C_n : D_m × D_n (translations + axis flips) ----------
def map_vertex(i, j, m, n, a, b, flip_i, flip_j):
    if flip_i:
        i = (-i) % m
    if flip_j:
        j = (-j) % n
    return ((i + a) % m, (j + b) % n)

def map_edge_id(eid, m, n, a, b, flip_i, flip_j):
    E0 = m*n
    if eid < E0:
        # horizontal H(i,j): (i,j)-(i,j+1)
        i = eid // n
        j = eid % n
        u = map_vertex(i, j, m, n, a, b, flip_i, flip_j)
        v = map_vertex(i, (j+1) % n, m, n, a, b, flip_i, flip_j)
        (iu, ju), (iv, jv) = u, v
        if iu != iv:
            raise RuntimeError("H mapped to non-H under chosen aut (shouldn't happen here)")
        i2 = iu
        if (jv - ju) % n == 1:
            j2 = ju
        else:
            j2 = jv
        return edge_id_H(i2, j2, m, n)
    else:
        # vertical V(i,j): (i,j)-(i+1,j)
        t = eid - E0
        i = t // n
        j = t % n
        u = map_vertex(i, j, m, n, a, b, flip_i, flip_j)
        v = map_vertex((i+1) % m, j, m, n, a, b, flip_i, flip_j)
        (iu, ju), (iv, jv) = u, v
        if ju != jv:
            raise RuntimeError("V mapped to non-V under chosen aut (shouldn't happen here)")
        j2 = ju
        if (iv - iu) % m == 1:
            i2 = iu
        else:
            i2 = iv
        return edge_id_V(i2, j2, m, n)


# ---------- S5 conjugacy classes representatives ----------
# We represent a color permutation sigma by an array perm[c]=sigma(c), with c in [0..k-1].
def s5_conj_reps_and_sizes():
    # cycle types for S5: sizes: 1,10,15,20,20,30,24
    # reps below written on colors 0..4
    reps = []
    # 1^5
    reps.append(([0,1,2,3,4], 1))
    # 2·1^3 : swap 0 and 1
    reps.append(([1,0,2,3,4], 10))
    # 2^2·1 : (0 1)(2 3)
    reps.append(([1,0,3,2,4], 15))
    # 3·1^2 : (0 1 2)
    reps.append(([1,2,0,3,4], 20))
    # 3·2 : (0 1 2)(3 4)
    reps.append(([1,2,0,4,3], 20))
    # 4·1 : (0 1 2 3)
    reps.append(([1,2,3,0,4], 30))
    # 5 : (0 1 2 3 4)
    reps.append(([1,2,3,4,0], 24))
    return reps


# ---------- invariance constraints: x_{e,c} <-> x_{phi(e), sigma(c)} ----------
def add_invariance_clauses(base: List[List[int]], m, n, k, phi_params, sigma_perm):
    a, b, flip_i, flip_j = phi_params
    E = num_edges(m, n)
    out = list(base)

    for eid in range(E):
        eid2 = map_edge_id(eid, m, n, a, b, flip_i, flip_j)
        for c in range(k):
            c2 = sigma_perm[c]
            x1 = var_color(eid, c, k)
            x2 = var_color(eid2, c2, k)
            out.append([-x1, x2])
            out.append([-x2, x1])

    return out


def approx_fix_log2(nvars, clauses, seed):
    log2_est, _ = estimate_log2_model_count(
        nvars=nvars,
        cnf_clauses=clauses,
        cmsat_path=CMSAT_PATH,
        trials=TRIALS_PER_FIX,
        p_in_xor=P_IN_XOR,
        seed=seed
    )
    return log2_est


def main():
    nvars, base = build_base_bedge_bcolor_cnf(m, n, k, use_c4=USE_C4_DISTINCT)

    # enumerate Aut = D_m × D_n : translations + flips
    auts = []
    for a in range(m):
        for b in range(n):
            for fi in (False, True):
                for fj in (False, True):
                    auts.append((a, b, fi, fj))
    assert len(auts) == 4*m*n  # = 84 for 3x7

    reps = s5_conj_reps_and_sizes()

    # Burnside sum in log-space is tricky; we’ll do it in float space using 2**log2.
    # For your sizes it may be huge, but Python float still gives you magnitude.
    burn_sum = 0.0
    G_size = (4*m*n) * 120  # |Aut| * |S5|

    seed0 = 20251218

    for ai, phi in enumerate(auts):
        for ri, (sigma, cls_size) in enumerate(reps):
            clauses = add_invariance_clauses(base, m, n, k, phi, sigma)
            log2_fix = approx_fix_log2(nvars, clauses, seed=seed0 + 100000*ai + 1000*ri)
            fix_est = 2.0 ** log2_fix

            burn_sum += cls_size * fix_est

            print(f"aut {ai+1:3d}/{len(auts)}  rep {ri+1}/7  log2Fix≈{log2_fix:.2f}  contrib≈{cls_size}*2^{log2_fix:.2f}")

    orbit_est = burn_sum / G_size
    print("\n==== Approx orbit count (non-similar solutions) ====")
    print(f"Group size |Aut|*|S5| = {G_size}")
    print(f"Approx orbits ≈ {orbit_est:.4e}")
    print("(This is an approximate value; increase TRIALS_PER_FIX for stability.)")


if __name__ == "__main__":
    main()

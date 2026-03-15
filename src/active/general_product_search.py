import os
import random
from datetime import datetime
from itertools import product
from typing import Dict, List, Optional, Set, Tuple

try:
    from src.active.general_product_core import (
        ProductGraph,
        build_base_cnf,
        build_problem_name,
        decode_model_to_edge_colors,
        env_int,
        make_product_spec,
        run_cmsat,
        var_color,
        write_dimacs,
    )
except ModuleNotFoundError:
    from general_product_core import (
        ProductGraph,
        build_base_cnf,
        build_problem_name,
        decode_model_to_edge_colors,
        env_int,
        make_product_spec,
        run_cmsat,
        var_color,
        write_dimacs,
    )


CMSAT_PATH = os.getenv("CMSAT_PATH", "cryptominisat5")
CMSAT_THREADS = max(1, env_int("CMSAT_THREADS", min(12, os.cpu_count() or 1)))
CMSAT_VERB = max(0, env_int("CMSAT_VERB", 0))

USE_C4_DISTINCT = True

MAX_SOLVER_CALLS = 10000
MAX_NO_NEW_CLASS = 500
MAX_NO_BEST_IMPROVE = 2000

W_LAYER = 1.0
W_TRANS = 1.0
W_UNIF = 0.1

MAX_SHIFT_SAMPLES = 2000
RANDOM_SEED = 0

RESULTS_ROOT = os.getenv("BCOLOR_RESULTS_ROOT", os.path.join("results", "runs"))


def canonical_relabel_key(colors: List[int]) -> Tuple[int, ...]:
    mapping: Dict[int, int] = {}
    next_color = 0
    relabeled = []
    for color in colors:
        if color not in mapping:
            mapping[color] = next_color
            next_color += 1
        relabeled.append(mapping[color])
    return tuple(relabeled)


def transform_vertex(
    v: Tuple[int, ...],
    sizes: List[int],
    periodic: List[bool],
    shift: List[int],
    flip: List[bool],
) -> Tuple[int, ...]:
    vv = list(v)
    for d in range(len(vv)):
        L = sizes[d]
        x = vv[d]
        if flip[d]:
            if periodic[d]:
                x = (-x) % L
            else:
                x = L - 1 - x
        if periodic[d]:
            x = (x + shift[d]) % L
        vv[d] = x
    return tuple(vv)


def apply_aut_to_coloring(
    graph: ProductGraph,
    edge_colors: List[int],
    shift: List[int],
    flip: List[bool],
) -> List[int]:
    new_colors = [-1] * graph.E
    for eid, (u, v) in enumerate(graph.edges):
        u2 = transform_vertex(u, graph.sizes, graph.periodic, shift, flip)
        v2 = transform_vertex(v, graph.sizes, graph.periodic, shift, flip)
        key = (u2, v2) if u2 <= v2 else (v2, u2)
        eid2 = graph.edge_of_uv[key]
        new_colors[eid2] = edge_colors[eid]
    assert all(color >= 0 for color in new_colors)
    return new_colors


def canonical_key_under_aut_color(
    graph: ProductGraph,
    edge_colors: List[int],
) -> Tuple[int, ...]:
    random.seed(RANDOM_SEED)

    cycle_dims = [d for d, per in enumerate(graph.periodic) if per]
    total_shifts = 1
    for d in cycle_dims:
        total_shifts *= graph.sizes[d]

    def iter_shifts():
        if total_shifts <= MAX_SHIFT_SAMPLES:
            for vals in product(*[range(graph.sizes[d]) for d in cycle_dims]):
                shift = [0] * graph.D
                for idx, d in enumerate(cycle_dims):
                    shift[d] = vals[idx]
                yield shift
            return

        for _ in range(MAX_SHIFT_SAMPLES):
            shift = [0] * graph.D
            for d in cycle_dims:
                shift[d] = random.randrange(graph.sizes[d])
            yield shift

    best: Optional[Tuple[int, ...]] = None
    for shift in iter_shifts():
        for flip_mask in product([False, True], repeat=graph.D):
            flip = list(flip_mask)
            colors2 = apply_aut_to_coloring(graph, edge_colors, shift, flip)
            key = canonical_relabel_key(colors2)
            if best is None or key < best:
                best = key
    assert best is not None
    return best


def uniform_penalty(edge_colors: List[int], k: int) -> float:
    counts = [0] * k
    for color in edge_colors:
        counts[color] += 1
    ideal = len(edge_colors) / k
    return sum((count - ideal) ** 2 for count in counts)


def translation_symmetry_score(graph: ProductGraph, edge_colors: List[int]) -> float:
    cycle_dims = [d for d, per in enumerate(graph.periodic) if per]
    if not cycle_dims:
        return 0.0

    total = 0
    good = 0
    for vals in product(*[range(graph.sizes[d]) for d in cycle_dims]):
        if all(v == 0 for v in vals):
            continue
        shift = [0] * graph.D
        for idx, d in enumerate(cycle_dims):
            shift[d] = vals[idx]
        colors2 = apply_aut_to_coloring(graph, edge_colors, shift, [False] * graph.D)
        total += 1
        if colors2 == edge_colors:
            good += 1

    return good / total if total > 0 else 0.0


def layer_regular_score_generic(
    graph: ProductGraph,
    edge_colors: List[int],
    k: int,
) -> float:
    per_dim_lines: Dict[int, Dict[Tuple[int, ...], List[Tuple[int, int]]]] = {
        d: {} for d in range(graph.D)
    }

    for eid, (u, v) in enumerate(graph.edges):
        diff = [i for i in range(graph.D) if u[i] != v[i]]
        if len(diff) != 1:
            continue
        d = diff[0]
        t = min(u[d], v[d])
        line_key = tuple(u[i] for i in range(graph.D) if i != d)
        per_dim_lines[d].setdefault(line_key, []).append((t, edge_colors[eid]))

    scores = []
    for d in range(graph.D):
        lines = per_dim_lines[d].values()
        if not lines:
            continue
        hit_sum = 0
        len_sum = 0
        for line in lines:
            line_sorted = sorted(line, key=lambda x: x[0])
            ts = [x[0] for x in line_sorted]
            cs = [x[1] for x in line_sorted]
            if not cs:
                continue

            best = 0
            for a in range(k):
                for b in range(k):
                    hits = 0
                    for idx in range(len(cs)):
                        if cs[idx] == (a * ts[idx] + b) % k:
                            hits += 1
                    if hits > best:
                        best = hits
            hit_sum += best
            len_sum += len(cs)
        scores.append(hit_sum / len_sum if len_sum else 0.0)

    return sum(scores) / len(scores) if scores else 0.0


def score_coloring(
    graph: ProductGraph,
    edge_colors: List[int],
    k: int,
) -> Tuple[float, float, float, float]:
    # Frozen search objective:
    # S = W_LAYER * L + W_TRANS * T - W_UNIF * U
    L = layer_regular_score_generic(graph, edge_colors, k)
    T = translation_symmetry_score(graph, edge_colors)
    U = uniform_penalty(edge_colors, k)
    S = W_LAYER * L + W_TRANS * T - W_UNIF * U
    return S, L, T, U


def make_blocking_clause(edge_colors: List[int], k: int) -> List[int]:
    return [-var_color(e, edge_colors[e], k) for e in range(len(edge_colors))]


def _run_search_cmsat(cnf_path: str):
    result = run_cmsat(
        cnf_path,
        cmsat_path=CMSAT_PATH,
        threads=CMSAT_THREADS,
        verb=CMSAT_VERB,
        need_model=True,
    )
    if result.status == "UNSAT":
        return None
    return result.model


def search_best(cycles: List[int], paths: List[int], k: int = 5):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    name = build_problem_name(cycles, paths, k)
    out_dir = os.path.join(RESULTS_ROOT, f"{stamp}_{name}")
    os.makedirs(out_dir, exist_ok=True)

    sizes, periodic = make_product_spec(cycles, paths)
    graph = ProductGraph(sizes, periodic)
    max_deg = max((len(graph.incident[v]) for v in graph.vertices), default=0)
    print(
        f"Graph: D={graph.D}, |V|={graph.VN}, |E|={graph.E}, "
        f"maxdeg={max_deg}, squares={len(graph.C4s)}"
    )

    nvars, base_clauses = build_base_cnf(
        graph,
        k,
        use_c4_distinct=USE_C4_DISTINCT,
    )

    blocking_clauses: List[List[int]] = []
    seen_keys: Set[Tuple[int, ...]] = set()

    best = None
    no_new = 0
    no_improve = 0

    cnf_path = os.path.join(out_dir, "working.cnf")

    for call_id in range(1, MAX_SOLVER_CALLS + 1):
        all_clauses = base_clauses + blocking_clauses
        write_dimacs(nvars, all_clauses, cnf_path)

        model = _run_search_cmsat(cnf_path)
        if model is None:
            print(f"[{call_id}] UNSAT: search space exhausted under current blocking.")
            break

        edge_colors = decode_model_to_edge_colors(model, graph.E, k)
        blocking_clauses.append(make_blocking_clause(edge_colors, k))

        key = canonical_key_under_aut_color(graph, edge_colors)
        if key in seen_keys:
            no_new += 1
            no_improve += 1
        else:
            seen_keys.add(key)
            no_new = 0

            score, L, T, U = score_coloring(graph, edge_colors, k)
            if best is None or score > best[0]:
                best = (score, L, T, U, edge_colors)
                no_improve = 0
                print(
                    f"[{call_id}] NEW BEST  S={score:.4f}  "
                    f"L={L:.3f}  T={T:.3f}  U={U:.1f}  unique={len(seen_keys)}"
                )
            else:
                no_improve += 1

        if no_new >= MAX_NO_NEW_CLASS:
            print(
                f"[{call_id}] Stop: no new equivalence class for "
                f"{MAX_NO_NEW_CLASS} consecutive calls. unique={len(seen_keys)}"
            )
            break
        if no_improve >= MAX_NO_BEST_IMPROVE:
            print(
                f"[{call_id}] Stop: no best-score improvement for "
                f"{MAX_NO_BEST_IMPROVE} consecutive calls. unique={len(seen_keys)}"
            )
            break

        if call_id % 200 == 0:
            print(
                f"[{call_id}] progress: unique={len(seen_keys)} "
                f"blocking={len(blocking_clauses)} no_new={no_new} "
                f"no_improve={no_improve}"
            )

    if best is None:
        print("No feasible coloring found.")
        return None

    score, L, T, U, best_colors = best
    sol_path = os.path.join(out_dir, "best_solution.txt")
    with open(sol_path, "w", encoding="utf-8") as f:
        f.write(f"# Best S={score:.6f}, LayerReg={L:.6f}, TransSym={T:.6f}, Penalty={U:.3f}\n")
        f.write(f"# Graph sizes={graph.sizes}, periodic={graph.periodic}, E={graph.E}\n")
        f.write("# Format: eid : u -- v : color(1..k)\n")
        for eid, (u, v) in enumerate(graph.edges):
            f.write(f"{eid} : {u} -- {v} : {best_colors[eid] + 1}\n")

    print("\n==== DONE ====")
    print("Results dir:", out_dir)
    print("Best solution:", sol_path)
    return out_dir, best


if __name__ == "__main__":
    # search_best(cycles=[7], paths=[3], k=4)
    # search_best(cycles=[3, 5], paths=[4], k=6)
    pass

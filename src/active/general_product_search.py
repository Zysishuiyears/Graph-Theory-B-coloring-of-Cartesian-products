import os
import random
import time
from datetime import datetime
from itertools import product
from typing import Dict, List, Optional, Set, Tuple

try:
    from src.active.general_product_core import (
        ProductGraph,
        RunLogger,
        build_base_cnf,
        build_problem_name,
        cnf_summary_line,
        decode_model_to_edge_colors,
        env_int,
        format_elapsed,
        get_graph_stats,
        graph_summary_line,
        make_product_spec,
        run_cmsat,
        var_color,
        write_dimacs,
    )
except ModuleNotFoundError:
    from general_product_core import (
        ProductGraph,
        RunLogger,
        build_base_cnf,
        build_problem_name,
        cnf_summary_line,
        decode_model_to_edge_colors,
        env_int,
        format_elapsed,
        get_graph_stats,
        graph_summary_line,
        make_product_spec,
        run_cmsat,
        var_color,
        write_dimacs,
    )


CMSAT_PATH = os.getenv("CMSAT_PATH", "cryptominisat5")
CMSAT_THREADS = max(1, env_int("CMSAT_THREADS", min(12, os.cpu_count() or 1)))
CMSAT_VERB = max(0, env_int("CMSAT_VERB", 0))
HEARTBEAT_SEC = max(0, env_int("BCOLOR_HEARTBEAT_SEC", 15))

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
    L = layer_regular_score_generic(graph, edge_colors, k)
    T = translation_symmetry_score(graph, edge_colors)
    U = uniform_penalty(edge_colors, k)
    S = W_LAYER * L + W_TRANS * T - W_UNIF * U
    return S, L, T, U


def make_blocking_clause(edge_colors: List[int], k: int) -> List[int]:
    return [-var_color(e, edge_colors[e], k) for e in range(len(edge_colors))]


def _run_search_cmsat(
    cnf_path: str,
    logger: RunLogger,
    call_id: int,
    blocking_count: int,
    unique_count: int,
):
    def _heartbeat(elapsed_sec: float) -> None:
        logger.log(
            f"[{call_id}] solver heartbeat: elapsed={format_elapsed(elapsed_sec)} "
            f"blocking={blocking_count} unique={unique_count}; still solving, not stuck."
        )

    return run_cmsat(
        cnf_path,
        cmsat_path=CMSAT_PATH,
        threads=CMSAT_THREADS,
        verb=CMSAT_VERB,
        need_model=True,
        heartbeat_sec=HEARTBEAT_SEC,
        heartbeat_cb=_heartbeat,
    )


def _write_search_summary(
    out_dir: str,
    stop_reason: str,
    total_elapsed_sec: float,
    total_solver_calls: int,
    unique_classes: int,
    graph_stats,
    nvars: int,
    base_clause_count: int,
    best,
) -> None:
    summary_path = os.path.join(out_dir, "search_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(f"stop_reason={stop_reason}\n")
        f.write(f"total_elapsed_sec={total_elapsed_sec:.6f}\n")
        f.write(f"total_solver_calls={total_solver_calls}\n")
        f.write(f"unique_classes={unique_classes}\n")
        f.write(f"graph_D={graph_stats.D}\n")
        f.write(f"graph_V={graph_stats.vertex_count}\n")
        f.write(f"graph_E={graph_stats.edge_count}\n")
        f.write(f"graph_maxdeg={graph_stats.max_degree}\n")
        f.write(f"graph_C4={graph_stats.c4_count}\n")
        f.write(f"nvars={nvars}\n")
        f.write(f"base_clause_count={base_clause_count}\n")
        if best is None:
            f.write("best_found=False\n")
            return
        score, L, T, U, _ = best
        f.write("best_found=True\n")
        f.write(f"best_S={score:.6f}\n")
        f.write(f"best_L={L:.6f}\n")
        f.write(f"best_T={T:.6f}\n")
        f.write(f"best_U={U:.6f}\n")


def search_best(cycles: List[int], paths: List[int], k: int = 5):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    name = build_problem_name(cycles, paths, k)
    out_dir = os.path.join(RESULTS_ROOT, f"{stamp}_{name}")
    os.makedirs(out_dir, exist_ok=True)
    logger = RunLogger(out_dir)

    logger.log(f"Results dir: {out_dir}")
    logger.log(f"Progress log: {logger.path}")
    logger.log("Search objective: S = W_LAYER * L + W_TRANS * T - W_UNIF * U")
    logger.log(
        f"Search stops: MAX_SOLVER_CALLS={MAX_SOLVER_CALLS}, "
        f"MAX_NO_NEW_CLASS={MAX_NO_NEW_CLASS}, "
        f"MAX_NO_BEST_IMPROVE={MAX_NO_BEST_IMPROVE}"
    )
    logger.log(
        f"Solver config: cmsat_path={CMSAT_PATH}, threads={CMSAT_THREADS}, heartbeat={HEARTBEAT_SEC}s"
    )

    overall_started_at = time.perf_counter()
    logger.log("Constructing product graph...")
    sizes, periodic = make_product_spec(cycles, paths)
    graph = ProductGraph(sizes, periodic)
    graph_stats = get_graph_stats(graph)
    logger.log(graph_summary_line(graph_stats))

    logger.log("Building base search CNF...")
    cnf_started_at = time.perf_counter()
    nvars, base_clauses = build_base_cnf(
        graph,
        k,
        use_c4_distinct=USE_C4_DISTINCT,
    )
    logger.log(cnf_summary_line(nvars, len(base_clauses)))
    logger.log(f"Base search CNF built in {format_elapsed(time.perf_counter() - cnf_started_at)}")

    blocking_clauses: List[List[int]] = []
    seen_keys: Set[Tuple[int, ...]] = set()

    best = None
    no_new = 0
    no_improve = 0
    total_solver_calls = 0
    stop_reason = "max_solver_calls_reached"

    cnf_path = os.path.join(out_dir, "working.cnf")
    logger.log(f"Working CNF path: {cnf_path}")

    for call_id in range(1, MAX_SOLVER_CALLS + 1):
        total_solver_calls = call_id
        call_started_at = time.perf_counter()
        all_clauses = base_clauses + blocking_clauses
        logger.log(
            f"[{call_id}] start: clauses={len(all_clauses)} blocking={len(blocking_clauses)} unique={len(seen_keys)}"
        )
        write_dimacs(nvars, all_clauses, cnf_path)

        result = _run_search_cmsat(
            cnf_path,
            logger,
            call_id,
            blocking_count=len(blocking_clauses),
            unique_count=len(seen_keys),
        )
        logger.log(
            f"[{call_id}] solver finished: status={result.status} elapsed={format_elapsed(result.elapsed_sec)}"
        )

        if result.status == "UNSAT":
            stop_reason = "search_space_exhausted"
            logger.log(f"[{call_id}] UNSAT: search space exhausted under current blocking.")
            break

        if result.model is None:
            raise RuntimeError("Expected a model for SAT result in search_best, but none was returned.")

        edge_colors = decode_model_to_edge_colors(result.model, graph.E, k)
        blocking_clauses.append(make_blocking_clause(edge_colors, k))

        key = canonical_key_under_aut_color(graph, edge_colors)
        if key in seen_keys:
            no_new += 1
            no_improve += 1
            logger.log(
                f"[{call_id}] duplicate equivalence class: unique={len(seen_keys)} "
                f"no_new={no_new} no_improve={no_improve}"
            )
        else:
            seen_keys.add(key)
            no_new = 0

            score, L, T, U = score_coloring(graph, edge_colors, k)
            if best is None or score > best[0]:
                best = (score, L, T, U, edge_colors)
                no_improve = 0
                logger.log(
                    f"[{call_id}] NEW BEST  S={score:.4f}  L={L:.3f}  T={T:.3f}  U={U:.1f}  unique={len(seen_keys)}"
                )
            else:
                no_improve += 1
                logger.log(
                    f"[{call_id}] feasible new class but not better: S={score:.4f} "
                    f"unique={len(seen_keys)} no_improve={no_improve}"
                )

        logger.log(
            f"[{call_id}] call complete: total_elapsed={format_elapsed(time.perf_counter() - call_started_at)} "
            f"blocking={len(blocking_clauses)} unique={len(seen_keys)}"
        )

        if no_new >= MAX_NO_NEW_CLASS:
            stop_reason = "max_no_new_class"
            logger.log(
                f"[{call_id}] Stop: no new equivalence class for {MAX_NO_NEW_CLASS} consecutive calls. unique={len(seen_keys)}"
            )
            break
        if no_improve >= MAX_NO_BEST_IMPROVE:
            stop_reason = "max_no_best_improve"
            logger.log(
                f"[{call_id}] Stop: no best-score improvement for {MAX_NO_BEST_IMPROVE} consecutive calls. unique={len(seen_keys)}"
            )
            break
    else:
        logger.log(f"Stop: reached MAX_SOLVER_CALLS={MAX_SOLVER_CALLS}.")

    total_elapsed_sec = time.perf_counter() - overall_started_at
    _write_search_summary(
        out_dir,
        stop_reason=stop_reason,
        total_elapsed_sec=total_elapsed_sec,
        total_solver_calls=total_solver_calls,
        unique_classes=len(seen_keys),
        graph_stats=graph_stats,
        nvars=nvars,
        base_clause_count=len(base_clauses),
        best=best,
    )

    if best is None:
        logger.log("No feasible coloring found under current constraints.")
        logger.log(f"Total search time: {format_elapsed(total_elapsed_sec)}")
        logger.log(f"Results dir: {out_dir}")
        return None

    score, L, T, U, best_colors = best
    sol_path = os.path.join(out_dir, "best_solution.txt")
    with open(sol_path, "w", encoding="utf-8") as f:
        f.write(f"# Best S={score:.6f}, LayerReg={L:.6f}, TransSym={T:.6f}, Penalty={U:.3f}\n")
        f.write(f"# Graph sizes={graph.sizes}, periodic={graph.periodic}, E={graph.E}\n")
        f.write("# Format: eid : u -- v : color(1..k)\n")
        for eid, (u, v) in enumerate(graph.edges):
            f.write(f"{eid} : {u} -- {v} : {best_colors[eid] + 1}\n")

    logger.log("==== DONE ====")
    logger.log(f"Best score: S={score:.4f}, L={L:.3f}, T={T:.3f}, U={U:.1f}")
    logger.log(f"Total search time: {format_elapsed(total_elapsed_sec)}")
    logger.log(f"Results dir: {out_dir}")
    logger.log(f"Best solution: {sol_path}")
    return out_dir, best


if __name__ == "__main__":
    # search_best(cycles=[7], paths=[3], k=4)
    # search_best(cycles=[3, 5], paths=[4], k=6)
    pass

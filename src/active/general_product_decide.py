import os
import subprocess
from datetime import datetime
import time
from typing import List, Optional, Tuple

try:
    from src.active.general_product_core import (
        ProductGraph,
        build_base_cnf,
        build_problem_name,
        decode_model_to_edge_colors,
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
        make_product_spec,
        run_cmsat,
        var_color,
        write_dimacs,
    )


def _env_int(name: str, default: int) -> int:
    raw = os.getenv(name)
    if raw is None:
        return default
    try:
        return int(raw)
    except ValueError:
        return default


def _env_flag(name: str, default: bool) -> bool:
    raw = os.getenv(name)
    if raw is None:
        return default
    return raw.strip().lower() not in {"0", "false", "no", "off"}


CMSAT_PATH = os.getenv("CMSAT_PATH", "cryptominisat5")
CMSAT_THREADS = max(1, _env_int("CMSAT_THREADS", min(12, os.cpu_count() or 1)))
CMSAT_VERB = max(0, _env_int("CMSAT_VERB", 0))
CMSAT_SOLVER_MODE = os.getenv("CMSAT_SOLVER_MODE", "single").strip().lower()
CMSAT_PORTFOLIO_WORKERS = max(
    1,
    _env_int("CMSAT_PORTFOLIO_WORKERS", min(4, max(1, CMSAT_THREADS // 3))),
)
BCOLOR_SYMMETRY_BREAKING = _env_flag("BCOLOR_SYMMETRY_BREAKING", False)
BCOLOR_HEARTBEAT_SEC = max(0, _env_int("BCOLOR_HEARTBEAT_SEC", 0))

USE_C4_DISTINCT = True
RESULTS_ROOT = os.getenv("BCOLOR_RESULTS_ROOT", os.path.join("results", "runs"))


def _ordered_root_incident_edges(graph: ProductGraph) -> List[int]:
    if not graph.vertices:
        return []
    root = min(graph.vertices)
    return sorted(graph.incident[root], key=lambda eid: graph.edges[eid])


def build_color_symmetry_breaking_clauses(graph: ProductGraph, k: int) -> List[List[int]]:
    clauses: List[List[int]] = []
    incident = _ordered_root_incident_edges(graph)
    for color_idx, eid in enumerate(incident[: min(len(incident), k)]):
        clauses.append([var_color(eid, color_idx, k)])
    return clauses


def _parse_model(output: str) -> List[int]:
    model: List[int] = []
    for line in output.splitlines():
        if line.startswith("v ") or line.startswith("V "):
            for token in line.split()[1:]:
                if token == "0":
                    continue
                try:
                    model.append(int(token))
                except ValueError:
                    pass
    return model


def _model_from_output(output: str) -> Optional[List[int]]:
    if "UNSAT" in output:
        return None
    if "SAT" not in output:
        raise RuntimeError(f"Cannot parse CryptoMiniSat output.\n{output.strip()}")
    return _parse_model(output)


def _terminate_process(proc: subprocess.Popen) -> None:
    if proc.poll() is not None:
        return
    proc.terminate()
    try:
        proc.wait(timeout=1.0)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait(timeout=1.0)


def _portfolio_profiles(total_threads: int) -> List[Tuple[int, int, List[str]]]:
    worker_count = min(max(1, CMSAT_PORTFOLIO_WORKERS), max(1, total_threads))
    threads_per_worker = max(1, total_threads // worker_count)
    variants = [
        ("auto", "glue"),
        ("stable", "geom"),
        ("rnd", "luby"),
        ("false", "glue"),
    ]
    profiles = []
    for worker_id in range(worker_count):
        polar, restart = variants[worker_id % len(variants)]
        profiles.append(
            (
                worker_id,
                threads_per_worker,
                ["--random", str(worker_id), "--restart", restart, "--polar", polar],
            )
        )
    return profiles


def _solve_with_portfolio(cnf_path: str) -> Optional[List[int]]:
    profiles = _portfolio_profiles(CMSAT_THREADS)
    if len(profiles) == 1:
        _, threads, extra_args = profiles[0]
        return run_cmsat(cnf_path, CMSAT_PATH, threads, CMSAT_VERB, extra_args=extra_args)

    running = []
    started_at = time.perf_counter()
    last_heartbeat = 0.0

    try:
        for worker_id, threads, extra_args in profiles:
            command = [
                CMSAT_PATH,
                "-t",
                str(max(1, threads)),
                "--verb",
                str(max(0, CMSAT_VERB)),
            ]
            command.extend(extra_args)
            command.append(cnf_path)
            try:
                proc = subprocess.Popen(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
            except FileNotFoundError as exc:
                for item in running:
                    _terminate_process(item["proc"])
                raise RuntimeError(f"Solver not found: {CMSAT_PATH}") from exc
            running.append({"worker_id": worker_id, "proc": proc})

        while running:
            finished_indices = []
            sat_model: Optional[List[int]] = None

            for idx, item in enumerate(running):
                proc = item["proc"]
                if proc.poll() is None:
                    continue
                stdout, stderr = proc.communicate()
                model = _model_from_output(stdout + "\n" + stderr)
                finished_indices.append(idx)
                if model is not None:
                    sat_model = model
                    break

            for idx in reversed(finished_indices):
                running.pop(idx)

            if sat_model is not None:
                for item in running:
                    _terminate_process(item["proc"])
                    item["proc"].communicate()
                return sat_model

            if not running:
                break

            elapsed_sec = time.perf_counter() - started_at
            if BCOLOR_HEARTBEAT_SEC > 0 and elapsed_sec - last_heartbeat >= BCOLOR_HEARTBEAT_SEC:
                print(
                    f"[solver] portfolio heartbeat: elapsed={elapsed_sec:.1f}s, "
                    f"workers_running={len(running)}",
                    flush=True,
                )
                last_heartbeat = elapsed_sec

            time.sleep(0.02)
    finally:
        for item in running:
            _terminate_process(item["proc"])
            item["proc"].communicate()

    return None


def _solve_decision_cnf(cnf_path: str) -> Optional[List[int]]:
    if CMSAT_SOLVER_MODE == "portfolio" and CMSAT_PORTFOLIO_WORKERS > 1 and CMSAT_THREADS > 1:
        return _solve_with_portfolio(cnf_path)
    return run_cmsat(cnf_path, CMSAT_PATH, CMSAT_THREADS, CMSAT_VERB)


def _edge_dim_and_start(
    graph: ProductGraph,
    u: Tuple[int, ...],
    v: Tuple[int, ...],
) -> Tuple[int, Tuple[int, ...]]:
    diff = [d for d in range(graph.D) if u[d] != v[d]]
    if len(diff) != 1:
        raise RuntimeError(f"edge {u} -- {v} is not 1-dimensional")

    dim = diff[0]
    if not graph.periodic[dim]:
        return (dim, u if u[dim] < v[dim] else v)

    size = graph.sizes[dim]
    if (u[dim] + 1) % size == v[dim]:
        return dim, u
    if (v[dim] + 1) % size == u[dim]:
        return dim, v
    raise RuntimeError(f"cannot orient periodic edge {u} -- {v} in dimension {dim}")


def _build_uv_matrices(
    graph: ProductGraph,
    edge_colors: List[int],
) -> Tuple[List[List[int]], List[List[int]]]:
    if graph.D != 2:
        raise RuntimeError("U/V matrices are only defined for two-factor products")

    u_rows = graph.sizes[0] if graph.periodic[0] else graph.sizes[0] - 1
    u_cols = graph.sizes[1]
    v_rows = graph.sizes[0]
    v_cols = graph.sizes[1] if graph.periodic[1] else graph.sizes[1] - 1

    U = [[0] * u_cols for _ in range(u_rows)]
    V = [[0] * v_cols for _ in range(v_rows)]

    for eid, (u, v) in enumerate(graph.edges):
        dim, start = _edge_dim_and_start(graph, u, v)
        color = edge_colors[eid] + 1
        if dim == 0:
            row = start[0]
            col = start[1]
            if not (0 <= row < u_rows and 0 <= col < u_cols):
                raise RuntimeError(f"U matrix index out of range for edge {u} -- {v}")
            if U[row][col] != 0:
                raise RuntimeError(f"duplicate U entry at ({row}, {col})")
            U[row][col] = color
        elif dim == 1:
            row = start[0]
            col = start[1]
            if not (0 <= row < v_rows and 0 <= col < v_cols):
                raise RuntimeError(f"V matrix index out of range for edge {u} -- {v}")
            if V[row][col] != 0:
                raise RuntimeError(f"duplicate V entry at ({row}, {col})")
            V[row][col] = color
        else:
            raise RuntimeError(f"unexpected dimension {dim} for 2D graph")

    if any(cell == 0 for row in U for cell in row):
        raise RuntimeError("U matrix is incomplete")
    if any(cell == 0 for row in V for cell in row):
        raise RuntimeError("V matrix is incomplete")
    return U, V


def _write_solution_file(
    path: str,
    graph: ProductGraph,
    edge_colors: List[int],
    k: int,
) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"# SAT solution: sizes={graph.sizes}, periodic={graph.periodic}, k={k}\n")
        if graph.D == 2:
            U, V = _build_uv_matrices(graph, edge_colors)
            f.write("# U matrix: factor-0 direction\n")
            for row in U:
                f.write(" ".join(map(str, row)) + "\n")
            f.write("# V matrix: factor-1 direction\n")
            for row in V:
                f.write(" ".join(map(str, row)) + "\n")
            f.write("\n")

        f.write("# eid : u -- v : color(1..k)\n")
        for eid, (u, v) in enumerate(graph.edges):
            f.write(f"{eid} : {u} -- {v} : {edge_colors[eid] + 1}\n")


def decide_existence_B(
    cycles: List[int],
    paths: List[int],
    k: int = 5,
    dump_one_solution: bool = True,
):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    name = build_problem_name(cycles, paths, k)
    out_dir = os.path.join(RESULTS_ROOT, f"{stamp}_{name}_DECIDE_B")
    os.makedirs(out_dir, exist_ok=True)

    sizes, periodic = make_product_spec(cycles, paths)
    graph = ProductGraph(sizes, periodic)
    max_deg = max((len(graph.incident[v]) for v in graph.vertices), default=0)
    c4_count = len(graph.C4s)

    print(
        f"Graph: D={graph.D}, |V|={graph.VN}, |E|={graph.E}, "
        f"maxdeg={max_deg}, #C4={c4_count}"
    )
    print(
        "Decide existence for: proper edge-coloring + rainbow 4-cycles, "
        f"k={k}, mode={CMSAT_SOLVER_MODE}, symmetry_breaking={BCOLOR_SYMMETRY_BREAKING}"
    )

    if k < max_deg:
        print(f"RESULT: UNSAT (k={k} < Δ(G)={max_deg}, proper edge-coloring is impossible)")
        return False, out_dir, None
    if c4_count > 0 and k < 4:
        print(f"RESULT: UNSAT (k={k} < 4 while #C4={c4_count}, rainbow 4-cycles are impossible)")
        return False, out_dir, None

    extra_clauses = (
        build_color_symmetry_breaking_clauses(graph, k) if BCOLOR_SYMMETRY_BREAKING else None
    )
    nvars, clauses = build_base_cnf(
        graph,
        k,
        use_c4_distinct=USE_C4_DISTINCT,
        extra_clauses=extra_clauses,
    )
    print(f"CNF: vars={nvars}, clauses={len(clauses)}")

    cnf_path = os.path.join(out_dir, "instance.cnf")
    write_dimacs(nvars, clauses, cnf_path)

    model = _solve_decision_cnf(cnf_path)
    if model is None:
        print("RESULT: UNSAT (不存在满足 proper edge-coloring + rainbow 4-cycles 的 k-色方案)")
        return False, out_dir, None

    print("RESULT: SAT (存在满足 proper edge-coloring + rainbow 4-cycles 的 k-色方案)")

    if not dump_one_solution:
        return True, out_dir, None

    edge_colors = decode_model_to_edge_colors(model, graph.E, k)
    sol_path = os.path.join(out_dir, "one_solution.txt")
    _write_solution_file(sol_path, graph, edge_colors, k)
    print("One solution dumped to:", sol_path)
    return True, out_dir, sol_path


if __name__ == "__main__":
    decide_existence_B(cycles=[4, 6 ,6], paths=[], k=6)

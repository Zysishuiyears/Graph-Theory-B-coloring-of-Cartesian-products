import os
import subprocess
import time
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional, Tuple

try:
    from src.active.general_product_core import (
        CmsatResult,
        ProductGraph,
        build_base_cnf,
        build_cmsat_command,
        build_problem_name,
        decode_model_to_edge_colors,
        env_flag,
        env_int,
        make_product_spec,
        result_from_cmsat_output,
        var_color,
        write_dimacs,
    )
except ModuleNotFoundError:
    from general_product_core import (
        CmsatResult,
        ProductGraph,
        build_base_cnf,
        build_cmsat_command,
        build_problem_name,
        decode_model_to_edge_colors,
        env_flag,
        env_int,
        make_product_spec,
        result_from_cmsat_output,
        var_color,
        write_dimacs,
    )


CMSAT_PATH = os.getenv("CMSAT_PATH", "cryptominisat5")
CMSAT_THREADS = max(1, env_int("CMSAT_THREADS", min(12, os.cpu_count() or 1)))
CMSAT_VERB = max(0, env_int("CMSAT_VERB", 0))
CMSAT_SOLVER_MODE = os.getenv("CMSAT_SOLVER_MODE", "portfolio").strip().lower()
CMSAT_PORTFOLIO_WORKERS = max(
    1,
    env_int("CMSAT_PORTFOLIO_WORKERS", min(4, max(1, CMSAT_THREADS // 3))),
)
BCOLOR_SYMMETRY_BREAKING = env_flag("BCOLOR_SYMMETRY_BREAKING", True)

USE_C4_DISTINCT = True
RESULTS_ROOT = os.getenv("BCOLOR_RESULTS_ROOT", os.path.join("results", "runs"))


_DECIDE_BASE_CACHE: Dict[
    Tuple[Tuple[int, ...], Tuple[int, ...], int, bool],
    Tuple[ProductGraph, int, List[List[int]]],
] = {}


@dataclass(frozen=True)
class SolverProfile:
    worker_id: int
    threads: int
    random_seed: int
    polar: str
    restart: str


@dataclass
class SolverAttempt:
    profile: SolverProfile
    result: CmsatResult


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


def _decision_cache_key(
    cycles: List[int],
    paths: List[int],
    k: int,
    symmetry_breaking: bool,
) -> Tuple[Tuple[int, ...], Tuple[int, ...], int, bool]:
    return tuple(cycles), tuple(paths), k, symmetry_breaking


def _get_decision_base(
    cycles: List[int],
    paths: List[int],
    k: int,
    symmetry_breaking: bool,
) -> Tuple[ProductGraph, int, List[List[int]]]:
    key = _decision_cache_key(cycles, paths, k, symmetry_breaking)
    cached = _DECIDE_BASE_CACHE.get(key)
    if cached is not None:
        return cached

    sizes, periodic = make_product_spec(cycles, paths)
    graph = ProductGraph(sizes, periodic)
    extra_clauses = (
        build_color_symmetry_breaking_clauses(graph, k) if symmetry_breaking else None
    )
    nvars, clauses = build_base_cnf(
        graph,
        k,
        use_c4_distinct=USE_C4_DISTINCT,
        extra_clauses=extra_clauses,
    )
    cached = (graph, nvars, clauses)
    _DECIDE_BASE_CACHE[key] = cached
    return cached


def _portfolio_profiles(total_threads: int) -> List[SolverProfile]:
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
            SolverProfile(
                worker_id=worker_id,
                threads=threads_per_worker,
                random_seed=worker_id,
                polar=polar,
                restart=restart,
            )
        )
    return profiles


def _build_portfolio_args(profile: SolverProfile) -> List[str]:
    return [
        "--random",
        str(profile.random_seed),
        "--restart",
        profile.restart,
        "--polar",
        profile.polar,
    ]


def _terminate_process(proc: subprocess.Popen) -> None:
    if proc.poll() is not None:
        return
    proc.terminate()
    try:
        proc.wait(timeout=1.0)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait(timeout=1.0)


def _solve_cnf_single(cnf_path: str, need_model: bool) -> SolverAttempt:
    profile = SolverProfile(
        worker_id=0,
        threads=CMSAT_THREADS,
        random_seed=0,
        polar="auto",
        restart="glue",
    )
    command = build_cmsat_command(
        cnf_path,
        CMSAT_PATH,
        profile.threads,
        verb=CMSAT_VERB,
        need_model=need_model,
        extra_args=_build_portfolio_args(profile),
    )
    started_at = time.perf_counter()
    try:
        proc = subprocess.run(command, capture_output=True, text=True, check=False)
    except FileNotFoundError as exc:
        raise RuntimeError(f"Solver not found: {CMSAT_PATH}") from exc
    elapsed_sec = time.perf_counter() - started_at
    output = proc.stdout + "\n" + proc.stderr
    result = result_from_cmsat_output(output, need_model, command, elapsed_sec)
    if result.status == "ERROR":
        raise RuntimeError(f"Cannot parse CryptoMiniSat output.\n{output.strip()}")
    return SolverAttempt(profile=profile, result=result)


def _solve_cnf_portfolio(
    cnf_path: str,
    need_model: bool,
) -> Tuple[SolverAttempt, List[SolverAttempt]]:
    profiles = _portfolio_profiles(CMSAT_THREADS)
    if len(profiles) == 1:
        attempt = _solve_cnf_single(cnf_path, need_model)
        return attempt, [attempt]

    running = []
    finished: List[SolverAttempt] = []

    try:
        for profile in profiles:
            command = build_cmsat_command(
                cnf_path,
                CMSAT_PATH,
                profile.threads,
                verb=CMSAT_VERB,
                need_model=need_model,
                extra_args=_build_portfolio_args(profile),
            )
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
            running.append(
                {
                    "profile": profile,
                    "command": command,
                    "proc": proc,
                    "started_at": time.perf_counter(),
                }
            )

        while running:
            ready_indices = []
            sat_attempt: Optional[SolverAttempt] = None

            for idx, item in enumerate(running):
                proc = item["proc"]
                if proc.poll() is None:
                    continue
                stdout, stderr = proc.communicate()
                elapsed_sec = time.perf_counter() - item["started_at"]
                output = stdout + "\n" + stderr
                result = result_from_cmsat_output(
                    output,
                    need_model,
                    item["command"],
                    elapsed_sec,
                )
                if result.status == "ERROR":
                    raise RuntimeError(
                        "Cannot parse CryptoMiniSat output from worker "
                        f"{item['profile'].worker_id}.\n{output.strip()}"
                    )
                attempt = SolverAttempt(profile=item["profile"], result=result)
                finished.append(attempt)
                ready_indices.append(idx)
                if result.status == "SAT":
                    sat_attempt = attempt
                    break

            for idx in reversed(ready_indices):
                running.pop(idx)

            if sat_attempt is not None:
                for item in running:
                    proc = item["proc"]
                    if proc.poll() is None:
                        _terminate_process(proc)
                    try:
                        proc.communicate(timeout=0.2)
                    except subprocess.TimeoutExpired:
                        proc.kill()
                        proc.communicate()
                    elapsed_sec = time.perf_counter() - item["started_at"]
                    finished.append(
                        SolverAttempt(
                            profile=item["profile"],
                            result=CmsatResult(
                                status="ABORTED",
                                model=None,
                                command=item["command"],
                                elapsed_sec=elapsed_sec,
                                output="",
                            ),
                        )
                    )
                return sat_attempt, finished

            if not ready_indices:
                time.sleep(0.05)
    finally:
        for item in running:
            _terminate_process(item["proc"])
            try:
                item["proc"].communicate(timeout=0.2)
            except subprocess.TimeoutExpired:
                item["proc"].kill()
                item["proc"].communicate()

    if finished and all(attempt.result.status == "UNSAT" for attempt in finished):
        return finished[0], finished

    raise RuntimeError("Portfolio solving finished without a SAT/UNSAT conclusion.")


def solve_decision_cnf(
    cnf_path: str,
    need_model: bool,
) -> Tuple[str, SolverAttempt, List[SolverAttempt]]:
    use_portfolio = (
        CMSAT_SOLVER_MODE == "portfolio"
        and CMSAT_PORTFOLIO_WORKERS > 1
        and CMSAT_THREADS > 1
    )
    if use_portfolio:
        winner, attempts = _solve_cnf_portfolio(cnf_path, need_model)
        return "portfolio", winner, attempts

    attempt = _solve_cnf_single(cnf_path, need_model)
    return "single", attempt, [attempt]


def _write_solver_summary(
    out_dir: str,
    mode: str,
    symmetry_breaking: bool,
    winner: SolverAttempt,
    attempts: List[SolverAttempt],
) -> None:
    summary_path = os.path.join(out_dir, "solver_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(f"mode={mode}\n")
        f.write(f"symmetry_breaking={symmetry_breaking}\n")
        f.write(f"total_threads={CMSAT_THREADS}\n")
        f.write(f"portfolio_workers={CMSAT_PORTFOLIO_WORKERS}\n")
        f.write(f"winner_status={winner.result.status}\n")
        f.write(f"winner_worker={winner.profile.worker_id}\n")
        f.write(f"winner_elapsed_sec={winner.result.elapsed_sec:.6f}\n")
        for attempt in sorted(attempts, key=lambda item: item.profile.worker_id):
            f.write(
                "worker={worker_id} status={status} elapsed_sec={elapsed:.6f} "
                "threads={threads} random={seed} polar={polar} restart={restart}\n".format(
                    worker_id=attempt.profile.worker_id,
                    status=attempt.result.status,
                    elapsed=attempt.result.elapsed_sec,
                    threads=attempt.profile.threads,
                    seed=attempt.profile.random_seed,
                    polar=attempt.profile.polar,
                    restart=attempt.profile.restart,
                )
            )


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

    symmetry_breaking = BCOLOR_SYMMETRY_BREAKING
    graph, nvars, clauses = _get_decision_base(cycles, paths, k, symmetry_breaking)

    max_deg = max((len(graph.incident[v]) for v in graph.vertices), default=0)
    print(
        f"Graph: D={graph.D}, |V|={graph.VN}, |E|={graph.E}, "
        f"maxdeg={max_deg}, #C4={len(graph.C4s)}"
    )
    print(
        "Decide existence for: proper edge-coloring + rainbow 4-cycles, "
        f"k={k}, mode={CMSAT_SOLVER_MODE}, symmetry_breaking={symmetry_breaking}"
    )

    cnf_path = os.path.join(out_dir, "instance.cnf")
    write_dimacs(nvars, clauses, cnf_path)

    mode, winner, attempts = solve_decision_cnf(
        cnf_path,
        need_model=dump_one_solution,
    )
    _write_solver_summary(out_dir, mode, symmetry_breaking, winner, attempts)

    if winner.result.status == "UNSAT":
        print("RESULT: UNSAT")
        return False, out_dir, None

    if winner.result.status != "SAT":
        raise RuntimeError(f"Unexpected solver status: {winner.result.status}")

    print("RESULT: SAT")

    if not dump_one_solution:
        return True, out_dir, None

    if winner.result.model is None:
        raise RuntimeError("Expected a model for SAT result, but none was returned.")

    edge_colors = decode_model_to_edge_colors(winner.result.model, graph.E, k)
    sol_path = os.path.join(out_dir, "one_solution.txt")
    with open(sol_path, "w", encoding="utf-8") as f:
        f.write(f"# SAT solution: sizes={graph.sizes}, periodic={graph.periodic}, k={k}\n")
        f.write("# eid : u -- v : color(1..k)\n")
        for eid, (u, v) in enumerate(graph.edges):
            f.write(f"{eid} : {u} -- {v} : {edge_colors[eid] + 1}\n")

    print("One solution dumped to:", sol_path)
    return True, out_dir, sol_path


if __name__ == "__main__":
    decide_existence_B(cycles=[4, 6, 10], paths=[], k=7)
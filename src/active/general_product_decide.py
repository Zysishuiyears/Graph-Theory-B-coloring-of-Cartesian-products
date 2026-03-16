import os
import subprocess
import time
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional, Tuple

try:
    from src.active.general_product_core import (
        CmsatResult,
        GraphStats,
        ProductGraph,
        RunLogger,
        build_base_cnf,
        build_cmsat_command,
        build_problem_name,
        cnf_summary_line,
        decode_model_to_edge_colors,
        env_flag,
        env_int,
        format_elapsed,
        get_graph_stats,
        graph_summary_line,
        make_product_spec,
        result_from_cmsat_output,
        run_cmsat,
        terminate_process,
        var_color,
        write_dimacs,
    )
except ModuleNotFoundError:
    from general_product_core import (
        CmsatResult,
        GraphStats,
        ProductGraph,
        RunLogger,
        build_base_cnf,
        build_cmsat_command,
        build_problem_name,
        cnf_summary_line,
        decode_model_to_edge_colors,
        env_flag,
        env_int,
        format_elapsed,
        get_graph_stats,
        graph_summary_line,
        make_product_spec,
        result_from_cmsat_output,
        run_cmsat,
        terminate_process,
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
HEARTBEAT_SEC = max(0, env_int("BCOLOR_HEARTBEAT_SEC", 15))

USE_C4_DISTINCT = True
RESULTS_ROOT = os.getenv("BCOLOR_RESULTS_ROOT", os.path.join("results", "runs"))


_DECIDE_GRAPH_CACHE: Dict[Tuple[Tuple[int, ...], Tuple[int, ...]], ProductGraph] = {}
_DECIDE_BASE_CACHE: Dict[
    Tuple[Tuple[int, ...], Tuple[int, ...], int, bool],
    Tuple[int, List[List[int]]],
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


@dataclass
class ObviousBoundInfo:
    lower_bound: int
    explanations: List[str]
    failing_reasons: List[str]


@dataclass
class DecisionSummary:
    result_status: str
    result_explanation: str
    stop_reason: str
    solver_invoked: bool
    total_elapsed_sec: float
    mode: str
    symmetry_breaking: bool
    graph_stats: GraphStats
    k: int
    obvious_info: ObviousBoundInfo
    cnf_path: Optional[str] = None
    nvars: Optional[int] = None
    clause_count: Optional[int] = None
    winner: Optional[SolverAttempt] = None
    attempts: Optional[List[SolverAttempt]] = None
    solution_path: Optional[str] = None


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


def _graph_cache_key(
    cycles: List[int],
    paths: List[int],
) -> Tuple[Tuple[int, ...], Tuple[int, ...]]:
    return tuple(cycles), tuple(paths)


def _decision_cache_key(
    cycles: List[int],
    paths: List[int],
    k: int,
    symmetry_breaking: bool,
) -> Tuple[Tuple[int, ...], Tuple[int, ...], int, bool]:
    return tuple(cycles), tuple(paths), k, symmetry_breaking


def _get_decision_graph(cycles: List[int], paths: List[int]) -> ProductGraph:
    key = _graph_cache_key(cycles, paths)
    cached = _DECIDE_GRAPH_CACHE.get(key)
    if cached is not None:
        return cached

    sizes, periodic = make_product_spec(cycles, paths)
    graph = ProductGraph(sizes, periodic)
    _DECIDE_GRAPH_CACHE[key] = graph
    return graph


def _get_decision_base(
    cycles: List[int],
    paths: List[int],
    k: int,
    symmetry_breaking: bool,
) -> Tuple[int, List[List[int]]]:
    key = _decision_cache_key(cycles, paths, k, symmetry_breaking)
    cached = _DECIDE_BASE_CACHE.get(key)
    if cached is not None:
        return cached

    graph = _get_decision_graph(cycles, paths)
    extra_clauses = (
        build_color_symmetry_breaking_clauses(graph, k) if symmetry_breaking else None
    )
    nvars, clauses = build_base_cnf(
        graph,
        k,
        use_c4_distinct=USE_C4_DISTINCT,
        extra_clauses=extra_clauses,
    )
    cached = (nvars, clauses)
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


def _compute_obvious_bounds(stats: GraphStats, k: int) -> ObviousBoundInfo:
    lower_bound = stats.max_degree
    explanations = [
        f"proper edge-coloring requires k >= Δ(G) = {stats.max_degree}",
    ]
    failing_reasons: List[str] = []

    if stats.c4_count > 0:
        lower_bound = max(lower_bound, 4)
        explanations.append(
            f"rainbow 4-cycles require k >= 4 because #C4 = {stats.c4_count}"
        )

    if k < stats.max_degree:
        failing_reasons.append(
            f"k = {k} < Δ(G) = {stats.max_degree}, so proper edge-coloring is impossible"
        )
    if stats.c4_count > 0 and k < 4:
        failing_reasons.append(
            f"k = {k} < 4 while #C4 = {stats.c4_count}, so a rainbow 4-cycle is impossible"
        )

    return ObviousBoundInfo(
        lower_bound=lower_bound,
        explanations=explanations,
        failing_reasons=failing_reasons,
    )


def _format_solver_profile(profile: SolverProfile) -> str:
    return (
        f"worker {profile.worker_id}: threads={profile.threads}, random={profile.random_seed}, "
        f"polar={profile.polar}, restart={profile.restart}"
    )


def _format_worker_states(
    profiles: List[SolverProfile],
    running: List[Dict[str, object]],
    finished: List[SolverAttempt],
) -> str:
    running_map = {
        item["profile"].worker_id: time.perf_counter() - item["started_at"]
        for item in running
    }
    finished_map = {attempt.profile.worker_id: attempt for attempt in finished}

    parts = []
    for profile in profiles:
        if profile.worker_id in finished_map:
            attempt = finished_map[profile.worker_id]
            parts.append(
                f"w{profile.worker_id}={attempt.result.status.lower()}"
                f"({format_elapsed(attempt.result.elapsed_sec)})"
            )
        elif profile.worker_id in running_map:
            parts.append(
                f"w{profile.worker_id}=running({format_elapsed(running_map[profile.worker_id])})"
            )
        else:
            parts.append(f"w{profile.worker_id}=pending")
    return ", ".join(parts)


def _solve_cnf_single(
    cnf_path: str,
    need_model: bool,
    logger: RunLogger,
) -> SolverAttempt:
    profile = SolverProfile(
        worker_id=0,
        threads=CMSAT_THREADS,
        random_seed=0,
        polar="auto",
        restart="glue",
    )
    logger.log(f"[solver] starting single mode with {_format_solver_profile(profile)}")

    def _heartbeat(elapsed_sec: float) -> None:
        logger.log(
            "[solver] heartbeat: "
            f"elapsed={format_elapsed(elapsed_sec)}, worker=0 status=running; still solving, not stuck."
        )

    result = run_cmsat(
        cnf_path,
        cmsat_path=CMSAT_PATH,
        threads=profile.threads,
        verb=CMSAT_VERB,
        need_model=need_model,
        extra_args=_build_portfolio_args(profile),
        heartbeat_sec=HEARTBEAT_SEC,
        heartbeat_cb=_heartbeat,
    )
    logger.log(
        f"[solver] worker 0 finished with {result.status} in {format_elapsed(result.elapsed_sec)}"
    )
    return SolverAttempt(profile=profile, result=result)


def _solve_cnf_portfolio(
    cnf_path: str,
    need_model: bool,
    logger: RunLogger,
) -> Tuple[SolverAttempt, List[SolverAttempt], str]:
    profiles = _portfolio_profiles(CMSAT_THREADS)
    if len(profiles) == 1:
        attempt = _solve_cnf_single(cnf_path, need_model, logger)
        stop_reason = f"single_{attempt.result.status.lower()}"
        return attempt, [attempt], stop_reason

    logger.log("[solver] starting portfolio mode")
    for profile in profiles:
        logger.log(f"[solver] {_format_solver_profile(profile)}")

    running: List[Dict[str, object]] = []
    finished: List[SolverAttempt] = []
    started_at = time.perf_counter()
    last_heartbeat = 0.0

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
                    terminate_process(item["proc"])
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
                logger.log(
                    f"[solver] worker {attempt.profile.worker_id} finished with "
                    f"{attempt.result.status} in {format_elapsed(attempt.result.elapsed_sec)}"
                )
                if result.status == "SAT":
                    sat_attempt = attempt
                    break

            for idx in reversed(ready_indices):
                running.pop(idx)

            if sat_attempt is not None:
                for item in running:
                    proc = item["proc"]
                    if proc.poll() is None:
                        terminate_process(proc)
                    try:
                        proc.communicate(timeout=0.2)
                    except subprocess.TimeoutExpired:
                        proc.kill()
                        proc.communicate()
                    elapsed_sec = time.perf_counter() - item["started_at"]
                    aborted_attempt = SolverAttempt(
                        profile=item["profile"],
                        result=CmsatResult(
                            status="ABORTED",
                            model=None,
                            command=item["command"],
                            elapsed_sec=elapsed_sec,
                            output="",
                        ),
                    )
                    finished.append(aborted_attempt)
                    logger.log(
                        f"[solver] worker {aborted_attempt.profile.worker_id} aborted after "
                        f"{format_elapsed(aborted_attempt.result.elapsed_sec)} because worker "
                        f"{sat_attempt.profile.worker_id} found SAT"
                    )
                stop_reason = f"sat_found_by_worker_{sat_attempt.profile.worker_id}"
                return sat_attempt, finished, stop_reason

            elapsed_total = time.perf_counter() - started_at
            if HEARTBEAT_SEC > 0 and elapsed_total - last_heartbeat >= HEARTBEAT_SEC:
                logger.log(
                    "[solver] portfolio heartbeat: "
                    f"elapsed={format_elapsed(elapsed_total)}, workers: "
                    f"{_format_worker_states(profiles, running, finished)}; still solving, not stuck."
                )
                last_heartbeat = elapsed_total

            if not ready_indices:
                time.sleep(0.05)
    finally:
        for item in running:
            terminate_process(item["proc"])
            try:
                item["proc"].communicate(timeout=0.2)
            except subprocess.TimeoutExpired:
                item["proc"].kill()
                item["proc"].communicate()

    if finished and all(attempt.result.status == "UNSAT" for attempt in finished):
        logger.log(
            f"[solver] all portfolio workers finished UNSAT in {format_elapsed(time.perf_counter() - started_at)}"
        )
        return finished[0], finished, "all_workers_unsat"

    raise RuntimeError("Portfolio solving finished without a SAT/UNSAT conclusion.")


def solve_decision_cnf(
    cnf_path: str,
    need_model: bool,
    logger: RunLogger,
) -> Tuple[str, SolverAttempt, List[SolverAttempt], str]:
    use_portfolio = (
        CMSAT_SOLVER_MODE == "portfolio"
        and CMSAT_PORTFOLIO_WORKERS > 1
        and CMSAT_THREADS > 1
    )
    if use_portfolio:
        winner, attempts, stop_reason = _solve_cnf_portfolio(cnf_path, need_model, logger)
        return "portfolio", winner, attempts, stop_reason

    attempt = _solve_cnf_single(cnf_path, need_model, logger)
    stop_reason = f"single_{attempt.result.status.lower()}"
    return "single", attempt, [attempt], stop_reason


def _write_solver_summary_file(out_dir: str, summary: DecisionSummary) -> None:
    summary_path = os.path.join(out_dir, "solver_summary.txt")
    attempts = summary.attempts or []
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(f"result_status={summary.result_status}\n")
        f.write(f"result_explanation={summary.result_explanation}\n")
        f.write(f"stop_reason={summary.stop_reason}\n")
        f.write(f"solver_invoked={summary.solver_invoked}\n")
        f.write(f"mode={summary.mode}\n")
        f.write(f"symmetry_breaking={summary.symmetry_breaking}\n")
        f.write(f"k={summary.k}\n")
        f.write(f"total_elapsed_sec={summary.total_elapsed_sec:.6f}\n")
        f.write(f"graph_D={summary.graph_stats.D}\n")
        f.write(f"graph_V={summary.graph_stats.vertex_count}\n")
        f.write(f"graph_E={summary.graph_stats.edge_count}\n")
        f.write(f"graph_maxdeg={summary.graph_stats.max_degree}\n")
        f.write(f"graph_C4={summary.graph_stats.c4_count}\n")
        f.write(f"obvious_lower_bound={summary.obvious_info.lower_bound}\n")
        for idx, text in enumerate(summary.obvious_info.explanations, start=1):
            f.write(f"obvious_explanation_{idx}={text}\n")
        for idx, text in enumerate(summary.obvious_info.failing_reasons, start=1):
            f.write(f"obvious_failure_{idx}={text}\n")
        if summary.cnf_path is not None:
            f.write(f"cnf_path={summary.cnf_path}\n")
        if summary.nvars is not None:
            f.write(f"nvars={summary.nvars}\n")
        if summary.clause_count is not None:
            f.write(f"clause_count={summary.clause_count}\n")
        if summary.solution_path is not None:
            f.write(f"solution_path={summary.solution_path}\n")
        if summary.winner is not None:
            f.write(f"winner_status={summary.winner.result.status}\n")
            f.write(f"winner_worker={summary.winner.profile.worker_id}\n")
            f.write(f"winner_elapsed_sec={summary.winner.result.elapsed_sec:.6f}\n")
        f.write(f"total_threads={CMSAT_THREADS}\n")
        f.write(f"portfolio_workers={CMSAT_PORTFOLIO_WORKERS}\n")
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


def _result_message(status: str) -> str:
    if status == "SAT":
        return "RESULT: SAT (存在满足 proper edge-coloring + rainbow 4-cycles 的 k-色方案)"
    if status == "UNSAT":
        return "RESULT: UNSAT (不存在满足 proper edge-coloring + rainbow 4-cycles 的 k-色方案)"
    return f"RESULT: {status}"


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
    logger = RunLogger(out_dir)

    logger.log(f"Results dir: {out_dir}")
    logger.log(f"Progress log: {logger.path}")
    logger.log(
        "Decide existence for: proper edge-coloring + rainbow 4-cycles, "
        f"k={k}, mode={CMSAT_SOLVER_MODE}, symmetry_breaking={BCOLOR_SYMMETRY_BREAKING}"
    )

    graph_started_at = time.perf_counter()
    logger.log("Constructing product graph...")
    graph = _get_decision_graph(cycles, paths)
    graph_elapsed = time.perf_counter() - graph_started_at
    graph_stats = get_graph_stats(graph)
    logger.log(graph_summary_line(graph_stats))
    logger.log(f"Graph built in {format_elapsed(graph_elapsed)}")

    obvious_info = _compute_obvious_bounds(graph_stats, k)
    logger.log(f"Obvious lower bound: k >= {obvious_info.lower_bound}")
    for text in obvious_info.explanations:
        logger.log(f"Obvious check: {text}")

    if obvious_info.failing_reasons:
        for text in obvious_info.failing_reasons:
            logger.log(f"Early UNSAT reason: {text}")
        result_message = _result_message("UNSAT")
        logger.log(result_message)
        summary = DecisionSummary(
            result_status="UNSAT",
            result_explanation=result_message,
            stop_reason="obvious_necessary_condition_failed",
            solver_invoked=False,
            total_elapsed_sec=0.0,
            mode="skipped",
            symmetry_breaking=BCOLOR_SYMMETRY_BREAKING,
            graph_stats=graph_stats,
            k=k,
            obvious_info=obvious_info,
        )
        _write_solver_summary_file(out_dir, summary)
        logger.log("Solver skipped because an obvious necessary condition already failed.")
        logger.log(f"Results dir: {out_dir}")
        return False, out_dir, None

    logger.log("Current k passes obvious checks; SAT solver will test sufficiency.")

    cnf_started_at = time.perf_counter()
    logger.log("Building decision CNF...")
    nvars, clauses = _get_decision_base(cycles, paths, k, BCOLOR_SYMMETRY_BREAKING)
    clause_count = len(clauses)
    logger.log(cnf_summary_line(nvars, clause_count))
    logger.log(f"Decision CNF built in {format_elapsed(time.perf_counter() - cnf_started_at)}")

    cnf_path = os.path.join(out_dir, "instance.cnf")
    write_dimacs(nvars, clauses, cnf_path)
    logger.log(f"CNF path: {cnf_path}")
    logger.log(
        f"Solver config: mode={CMSAT_SOLVER_MODE}, total_threads={CMSAT_THREADS}, "
        f"portfolio_workers={CMSAT_PORTFOLIO_WORKERS}, heartbeat={HEARTBEAT_SEC}s"
    )

    solve_started_at = time.perf_counter()
    mode, winner, attempts, stop_reason = solve_decision_cnf(
        cnf_path,
        need_model=dump_one_solution,
        logger=logger,
    )
    total_elapsed_sec = time.perf_counter() - solve_started_at

    if winner.result.status == "UNSAT":
        result_message = _result_message("UNSAT")
        logger.log(result_message)
        logger.log(f"Total solve time: {format_elapsed(total_elapsed_sec)}")
        summary = DecisionSummary(
            result_status="UNSAT",
            result_explanation=result_message,
            stop_reason=stop_reason,
            solver_invoked=True,
            total_elapsed_sec=total_elapsed_sec,
            mode=mode,
            symmetry_breaking=BCOLOR_SYMMETRY_BREAKING,
            graph_stats=graph_stats,
            k=k,
            obvious_info=obvious_info,
            cnf_path=cnf_path,
            nvars=nvars,
            clause_count=clause_count,
            winner=winner,
            attempts=attempts,
        )
        _write_solver_summary_file(out_dir, summary)
        logger.log("Solver skipped because an obvious necessary condition already failed.")
        logger.log(f"Results dir: {out_dir}")
        return False, out_dir, None

    if winner.result.status != "SAT":
        raise RuntimeError(f"Unexpected solver status: {winner.result.status}")

    result_message = _result_message("SAT")
    logger.log(result_message)

    solution_path: Optional[str] = None
    if dump_one_solution:
        if winner.result.model is None:
            raise RuntimeError("Expected a model for SAT result, but none was returned.")

        edge_colors = decode_model_to_edge_colors(winner.result.model, graph.E, k)
        solution_path = os.path.join(out_dir, "one_solution.txt")
        with open(solution_path, "w", encoding="utf-8") as f:
            f.write(f"# SAT solution: sizes={graph.sizes}, periodic={graph.periodic}, k={k}\n")
            f.write("# eid : u -- v : color(1..k)\n")
            for eid, (u, v) in enumerate(graph.edges):
                f.write(f"{eid} : {u} -- {v} : {edge_colors[eid] + 1}\n")
        logger.log(f"One solution dumped to: {solution_path}")
    else:
        logger.log("dump_one_solution=False, so no solution file was written.")

    logger.log(f"Total solve time: {format_elapsed(total_elapsed_sec)}")
    logger.log(f"Results dir: {out_dir}")

    summary = DecisionSummary(
        result_status="SAT",
        result_explanation=result_message,
        stop_reason=stop_reason,
        solver_invoked=True,
        total_elapsed_sec=total_elapsed_sec,
        mode=mode,
        symmetry_breaking=BCOLOR_SYMMETRY_BREAKING,
        graph_stats=graph_stats,
        k=k,
        obvious_info=obvious_info,
        cnf_path=cnf_path,
        nvars=nvars,
        clause_count=clause_count,
        winner=winner,
        attempts=attempts,
        solution_path=solution_path,
    )
    _write_solver_summary_file(out_dir, summary)
    return True, out_dir, solution_path


if __name__ == "__main__":
    decide_existence_B(cycles=[4, 6], paths=[], k=5)

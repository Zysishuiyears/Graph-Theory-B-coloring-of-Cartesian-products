import os
import random
import subprocess
import time
from dataclasses import dataclass
from datetime import datetime
from itertools import product
from typing import Dict, List, Optional, Sequence, Set, Tuple


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


USE_SOLVER = "cmsat"
CMSAT_PATH = os.getenv("CMSAT_PATH", "cryptominisat5")
CMSAT_THREADS = max(1, _env_int("CMSAT_THREADS", min(12, os.cpu_count() or 1)))
CMSAT_VERB = max(0, _env_int("CMSAT_VERB", 0))
CMSAT_SOLVER_MODE = os.getenv("CMSAT_SOLVER_MODE", "portfolio").strip().lower()
CMSAT_PORTFOLIO_WORKERS = max(
    1,
    _env_int("CMSAT_PORTFOLIO_WORKERS", min(4, max(1, CMSAT_THREADS // 3))),
)
BCOLOR_SYMMETRY_BREAKING = _env_flag("BCOLOR_SYMMETRY_BREAKING", True)

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


_DECIDE_BASE_CACHE: Dict[
    Tuple[Tuple[int, ...], Tuple[int, ...], int, bool],
    Tuple["ProductGraph", int, List[List[int]]],
] = {}


@dataclass(frozen=True)
class SolverProfile:
    worker_id: int
    threads: int
    random_seed: int
    polar: str
    restart: str


@dataclass
class SolverResult:
    status: str
    model: Optional[List[int]]
    profile: SolverProfile
    command: List[str]
    elapsed_sec: float
    output: str = ""


def var_color(e: int, c: int, k: int) -> int:
    return e * k + c + 1


def var_rep(e: int, c: int, E: int, k: int) -> int:
    return E * k + e * k + c + 1


class ProductGraph:
    def __init__(self, sizes: List[int], periodic: List[bool]):
        assert len(sizes) == len(periodic)
        self.sizes = sizes
        self.periodic = periodic
        self.D = len(sizes)

        self.vertices = list(product(*[range(L) for L in sizes]))
        self.VN = len(self.vertices)

        self.edges: List[Tuple[Tuple[int, ...], Tuple[int, ...]]] = []
        self.edge_of_uv: Dict[Tuple[Tuple[int, ...], Tuple[int, ...]], int] = {}

        self._build_edges()
        self.E = len(self.edges)

        self.incident: Dict[Tuple[int, ...], List[int]] = {v: [] for v in self.vertices}
        for eid, (u, v) in enumerate(self.edges):
            self.incident[u].append(eid)
            self.incident[v].append(eid)

        self.neigh = [set() for _ in range(self.E)]
        for v in self.vertices:
            inc = self.incident[v]
            for i in range(len(inc)):
                for j in range(i + 1, len(inc)):
                    e1, e2 = inc[i], inc[j]
                    self.neigh[e1].add(e2)
                    self.neigh[e2].add(e1)
        self.neigh = [list(s) for s in self.neigh]

        self.C4s = self._enumerate_squares()

    def _canon_uv(self, a: Tuple[int, ...], b: Tuple[int, ...]):
        return (a, b) if a <= b else (b, a)

    def _step(
        self, v: Tuple[int, ...], dim: int, delta: int
    ) -> Optional[Tuple[int, ...]]:
        L = self.sizes[dim]
        x = v[dim]
        y = x + delta
        if self.periodic[dim]:
            y %= L
        elif not (0 <= y < L):
            return None
        vv = list(v)
        vv[dim] = y
        return tuple(vv)

    def _build_edges(self) -> None:
        for v in self.vertices:
            for d in range(self.D):
                if self.periodic[d]:
                    w = self._step(v, d, +1)
                    assert w is not None
                else:
                    if v[d] >= self.sizes[d] - 1:
                        continue
                    w = self._step(v, d, +1)
                    assert w is not None

                key = self._canon_uv(v, w)
                if key in self.edge_of_uv:
                    continue
                eid = len(self.edges)
                self.edges.append((key[0], key[1]))
                self.edge_of_uv[key] = eid

    def _enumerate_squares(self) -> List[List[int]]:
        c4s: List[List[int]] = []
        for a in range(self.D):
            for b in range(a + 1, self.D):
                for v00 in self.vertices:
                    if (not self.periodic[a]) and (v00[a] >= self.sizes[a] - 1):
                        continue
                    if (not self.periodic[b]) and (v00[b] >= self.sizes[b] - 1):
                        continue

                    v10 = self._step(v00, a, +1)
                    v01 = self._step(v00, b, +1)
                    v11 = self._step(v10, b, +1)  # type: ignore[arg-type]
                    if v10 is None or v01 is None or v11 is None:
                        continue

                    e_a0 = self.edge_of_uv[self._canon_uv(v00, v10)]
                    e_b0 = self.edge_of_uv[self._canon_uv(v00, v01)]
                    e_a1 = self.edge_of_uv[self._canon_uv(v01, v11)]
                    e_b1 = self.edge_of_uv[self._canon_uv(v10, v11)]
                    c4s.append([e_a0, e_b0, e_a1, e_b1])

        for d in range(self.D):
            if not (self.periodic[d] and self.sizes[d] == 4):
                continue
            other_ranges = [range(self.sizes[i]) for i in range(self.D) if i != d]
            for fixed in product(*other_ranges):
                base = [0] * self.D
                idx = 0
                for i in range(self.D):
                    if i == d:
                        continue
                    base[i] = fixed[idx]
                    idx += 1

                verts = []
                for t in range(4):
                    vv = list(base)
                    vv[d] = t
                    verts.append(tuple(vv))

                edges = []
                for t in range(4):
                    u = verts[t]
                    v = verts[(t + 1) % 4]
                    edges.append(self.edge_of_uv[self._canon_uv(u, v)])
                c4s.append(edges)
        return c4s


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


def build_base_cnf(
    graph: ProductGraph,
    k: int,
    use_c4_distinct: bool = True,
    symmetry_breaking: bool = False,
) -> Tuple[int, List[List[int]]]:
    E = graph.E
    clauses: List[List[int]] = []

    for e in range(E):
        clauses.append([var_color(e, c, k) for c in range(k)])
        for c1 in range(k):
            for c2 in range(c1 + 1, k):
                clauses.append([-var_color(e, c1, k), -var_color(e, c2, k)])

    for v in graph.vertices:
        inc = graph.incident[v]
        for c in range(k):
            for i in range(len(inc)):
                for j in range(i + 1, len(inc)):
                    clauses.append([-var_color(inc[i], c, k), -var_color(inc[j], c, k)])

    if use_c4_distinct:
        for cyc in graph.C4s:
            for i in range(4):
                for j in range(i + 1, 4):
                    for c in range(k):
                        clauses.append([-var_color(cyc[i], c, k), -var_color(cyc[j], c, k)])

    if symmetry_breaking:
        clauses.extend(build_color_symmetry_breaking_clauses(graph, k))

    return E * k, clauses


def build_cnf_proper_plus_rainbowC4(
    graph: ProductGraph, k: int
) -> Tuple[int, List[List[int]]]:
    return build_base_cnf(
        graph,
        k,
        use_c4_distinct=True,
        symmetry_breaking=False,
    )


def write_dimacs(nvars: int, clauses: List[List[int]], path: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"p cnf {nvars} {len(clauses)}\n")
        for cls in clauses:
            f.write(" ".join(map(str, cls)) + " 0\n")


def _decision_cache_key(
    cycles: Sequence[int],
    paths: Sequence[int],
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

    sizes = cycles + paths
    periodic = [True] * len(cycles) + [False] * len(paths)
    graph = ProductGraph(sizes, periodic)
    nvars, clauses = build_base_cnf(
        graph,
        k,
        use_c4_distinct=USE_C4_DISTINCT,
        symmetry_breaking=symmetry_breaking,
    )
    cached = (graph, nvars, clauses)
    _DECIDE_BASE_CACHE[key] = cached
    return cached


def _build_solver_command(
    cnf_path: str,
    need_model: bool,
    profile: SolverProfile,
) -> List[str]:
    cmd = [
        CMSAT_PATH,
        "--threads",
        str(profile.threads),
        "--verb",
        str(CMSAT_VERB),
        "--random",
        str(profile.random_seed),
        "--restart",
        profile.restart,
        "--polar",
        profile.polar,
    ]
    if not need_model:
        cmd.extend(["--printsol,s", "0"])
    cmd.append(cnf_path)
    return cmd


def _parse_model_from_output(output: str) -> List[int]:
    model: List[int] = []
    for line in output.splitlines():
        if line.startswith("v ") or line.startswith("V "):
            for token in line.split()[1:]:
                if token == "0":
                    continue
                try:
                    model.append(int(token))
                except ValueError:
                    continue
    return model


def _build_solver_result(
    output: str,
    need_model: bool,
    profile: SolverProfile,
    command: List[str],
    elapsed_sec: float,
) -> SolverResult:
    if "UNSAT" in output:
        return SolverResult(
            status="UNSAT",
            model=None,
            profile=profile,
            command=command,
            elapsed_sec=elapsed_sec,
            output=output,
        )

    if "SAT" in output:
        model = _parse_model_from_output(output) if need_model else None
        return SolverResult(
            status="SAT",
            model=model,
            profile=profile,
            command=command,
            elapsed_sec=elapsed_sec,
            output=output,
        )

    return SolverResult(
        status="ERROR",
        model=None,
        profile=profile,
        command=command,
        elapsed_sec=elapsed_sec,
        output=output,
    )


def _solve_cnf_single(cnf_path: str, need_model: bool) -> SolverResult:
    profile = SolverProfile(
        worker_id=0,
        threads=max(1, CMSAT_THREADS),
        random_seed=0,
        polar="auto",
        restart="glue",
    )
    command = _build_solver_command(cnf_path, need_model, profile)
    started_at = time.perf_counter()
    try:
        proc = subprocess.run(command, capture_output=True, text=True, check=False)
    except FileNotFoundError as exc:
        raise RuntimeError(f"Solver not found: {CMSAT_PATH}") from exc
    elapsed_sec = time.perf_counter() - started_at
    output = proc.stdout + "\n" + proc.stderr
    result = _build_solver_result(output, need_model, profile, command, elapsed_sec)
    if result.status == "ERROR":
        raise RuntimeError(f"Cannot parse CryptoMiniSat output.\n{output.strip()}")
    return result


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


def _terminate_process(proc: subprocess.Popen) -> None:
    if proc.poll() is not None:
        return
    proc.terminate()
    try:
        proc.wait(timeout=1.0)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait(timeout=1.0)


def _solve_cnf_portfolio(
    cnf_path: str,
    need_model: bool,
) -> Tuple[SolverResult, List[SolverResult]]:
    profiles = _portfolio_profiles(CMSAT_THREADS)
    if len(profiles) == 1:
        result = _solve_cnf_single(cnf_path, need_model)
        return result, [result]

    running = []
    finished: List[SolverResult] = []

    try:
        for profile in profiles:
            command = _build_solver_command(cnf_path, need_model, profile)
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
            sat_result: Optional[SolverResult] = None

            for idx, item in enumerate(running):
                proc = item["proc"]
                if proc.poll() is None:
                    continue
                stdout, stderr = proc.communicate()
                elapsed_sec = time.perf_counter() - item["started_at"]
                output = stdout + "\n" + stderr
                result = _build_solver_result(
                    output,
                    need_model,
                    item["profile"],
                    item["command"],
                    elapsed_sec,
                )
                if result.status == "ERROR":
                    raise RuntimeError(
                        f"Cannot parse CryptoMiniSat output from worker {result.profile.worker_id}.\n"
                        f"{output.strip()}"
                    )
                finished.append(result)
                ready_indices.append(idx)
                if result.status == "SAT":
                    sat_result = result
                    break

            for idx in reversed(ready_indices):
                running.pop(idx)

            if sat_result is not None:
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
                        SolverResult(
                            status="ABORTED",
                            model=None,
                            profile=item["profile"],
                            command=item["command"],
                            elapsed_sec=elapsed_sec,
                            output="",
                        )
                    )
                return sat_result, finished

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

    if finished and all(result.status == "UNSAT" for result in finished):
        return finished[0], finished

    raise RuntimeError("Portfolio solving finished without a SAT/UNSAT conclusion.")


def solve_decision_cnf(
    cnf_path: str,
    need_model: bool,
) -> Tuple[str, SolverResult, List[SolverResult]]:
    use_portfolio = (
        CMSAT_SOLVER_MODE == "portfolio"
        and CMSAT_PORTFOLIO_WORKERS > 1
        and CMSAT_THREADS > 1
    )
    if use_portfolio:
        result, all_results = _solve_cnf_portfolio(cnf_path, need_model)
        return "portfolio", result, all_results

    result = _solve_cnf_single(cnf_path, need_model)
    return "single", result, [result]


def run_cmsat(cnf_path: str) -> Optional[List[int]]:
    result = _solve_cnf_single(cnf_path, need_model=True)
    return None if result.status == "UNSAT" else result.model


def decode_model_to_edge_colors(model: List[int], E: int, k: int) -> List[int]:
    max_color_var = E * k
    pos = {v for v in model if 0 < v <= max_color_var}

    edge_color = [-1] * E
    for e in range(E):
        for c in range(k):
            if var_color(e, c, k) in pos:
                edge_color[e] = c
                break
        if edge_color[e] < 0:
            raise RuntimeError("decode failed: some edge has no true color.")
    return edge_color


def _write_solver_summary(
    out_dir: str,
    mode: str,
    symmetry_breaking: bool,
    result: SolverResult,
    all_results: List[SolverResult],
) -> None:
    summary_path = os.path.join(out_dir, "solver_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(f"mode={mode}\n")
        f.write(f"symmetry_breaking={symmetry_breaking}\n")
        f.write(f"total_threads={CMSAT_THREADS}\n")
        f.write(f"portfolio_workers={CMSAT_PORTFOLIO_WORKERS}\n")
        f.write(f"winner_status={result.status}\n")
        f.write(f"winner_worker={result.profile.worker_id}\n")
        f.write(f"winner_elapsed_sec={result.elapsed_sec:.6f}\n")
        for worker in sorted(all_results, key=lambda item: item.profile.worker_id):
            f.write(
                "worker={worker_id} status={status} elapsed_sec={elapsed:.6f} "
                "threads={threads} random={seed} polar={polar} restart={restart}\n".format(
                    worker_id=worker.profile.worker_id,
                    status=worker.status,
                    elapsed=worker.elapsed_sec,
                    threads=worker.profile.threads,
                    seed=worker.profile.random_seed,
                    polar=worker.profile.polar,
                    restart=worker.profile.restart,
                )
            )


def canonical_relabel_key(colors: List[int]) -> Tuple[int, ...]:
    mp: Dict[int, int] = {}
    nxt = 0
    out = []
    for x in colors:
        if x not in mp:
            mp[x] = nxt
            nxt += 1
        out.append(mp[x])
    return tuple(out)


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
    assert all(c >= 0 for c in new_colors)
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


def decide_existence_B(
    cycles: List[int],
    paths: List[int],
    k: int = 5,
    dump_one_solution: bool = True,
):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    name = "x".join([f"C{c}" for c in cycles] + [f"P{p}" for p in paths]) + f"_k{k}"
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

    mode, result, all_results = solve_decision_cnf(
        cnf_path,
        need_model=dump_one_solution,
    )
    _write_solver_summary(out_dir, mode, symmetry_breaking, result, all_results)

    if result.status == "UNSAT":
        print("RESULT: UNSAT")
        return False, out_dir, None

    if result.status != "SAT":
        raise RuntimeError(f"Unexpected solver status: {result.status}")

    print("RESULT: SAT")

    if not dump_one_solution:
        return True, out_dir, None

    if result.model is None:
        raise RuntimeError("Expected a model for SAT result, but none was returned.")

    edge_colors = decode_model_to_edge_colors(result.model, graph.E, k)
    sol_path = os.path.join(out_dir, "one_solution.txt")
    with open(sol_path, "w", encoding="utf-8") as f:
        f.write(f"# SAT solution: sizes={graph.sizes}, periodic={graph.periodic}, k={k}\n")
        f.write("# eid : u -- v : color(1..k)\n")
        for eid, (u, v) in enumerate(graph.edges):
            f.write(f"{eid} : {u} -- {v} : {edge_colors[eid] + 1}\n")

    print("One solution dumped to:", sol_path)
    return True, out_dir, sol_path


def uniform_penalty(edge_colors: List[int], k: int) -> float:
    cnt = [0] * k
    for c in edge_colors:
        cnt[c] += 1
    ideal = len(edge_colors) / k
    return sum((x - ideal) ** 2 for x in cnt)


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


def search_best(cycles: List[int], paths: List[int], k: int = 5):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    name = "x".join([f"C{c}" for c in cycles] + [f"P{p}" for p in paths]) + f"_k{k}"
    out_dir = os.path.join(RESULTS_ROOT, f"{stamp}_{name}")
    os.makedirs(out_dir, exist_ok=True)

    sizes = cycles + paths
    periodic = [True] * len(cycles) + [False] * len(paths)

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
        symmetry_breaking=False,
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

        model = run_cmsat(cnf_path)
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
    decide_existence_B(cycles=[4, 6], paths=[], k=5)
    # decide_existence_B(cycles=[3, 5], paths=[4], k=6)

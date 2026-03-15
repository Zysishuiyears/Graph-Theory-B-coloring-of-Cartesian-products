import os
import subprocess
import time
from dataclasses import dataclass
from itertools import product
from typing import Dict, List, Optional, Sequence, Tuple


def env_int(name: str, default: int) -> int:
    raw = os.getenv(name)
    if raw is None:
        return default
    try:
        return int(raw)
    except ValueError:
        return default


def env_flag(name: str, default: bool) -> bool:
    raw = os.getenv(name)
    if raw is None:
        return default
    return raw.strip().lower() not in {"0", "false", "no", "off"}


def build_problem_name(cycles: Sequence[int], paths: Sequence[int], k: int) -> str:
    return "x".join([f"C{c}" for c in cycles] + [f"P{p}" for p in paths]) + f"_k{k}"


def make_product_spec(
    cycles: Sequence[int],
    paths: Sequence[int],
) -> Tuple[List[int], List[bool]]:
    sizes = list(cycles) + list(paths)
    periodic = [True] * len(cycles) + [False] * len(paths)
    return sizes, periodic


@dataclass
class CmsatResult:
    status: str
    model: Optional[List[int]]
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


def build_base_cnf(
    graph: ProductGraph,
    k: int,
    use_c4_distinct: bool = True,
    extra_clauses: Optional[Sequence[Sequence[int]]] = None,
) -> Tuple[int, List[List[int]]]:
    clauses: List[List[int]] = []

    for e in range(graph.E):
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

    if extra_clauses:
        clauses.extend([list(clause) for clause in extra_clauses])

    return graph.E * k, clauses


def write_dimacs(nvars: int, clauses: List[List[int]], path: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"p cnf {nvars} {len(clauses)}\n")
        for cls in clauses:
            f.write(" ".join(map(str, cls)) + " 0\n")


def build_cmsat_command(
    cnf_path: str,
    cmsat_path: str,
    threads: int,
    verb: int = 0,
    need_model: bool = True,
    extra_args: Optional[Sequence[str]] = None,
) -> List[str]:
    command = [
        cmsat_path,
        "--threads",
        str(max(1, threads)),
        "--verb",
        str(max(0, verb)),
    ]
    if extra_args:
        command.extend(extra_args)
    if not need_model:
        command.extend(["--printsol,s", "0"])
    command.append(cnf_path)
    return command


def parse_model_from_output(output: str) -> List[int]:
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


def result_from_cmsat_output(
    output: str,
    need_model: bool,
    command: List[str],
    elapsed_sec: float,
) -> CmsatResult:
    if "UNSAT" in output:
        return CmsatResult(
            status="UNSAT",
            model=None,
            command=command,
            elapsed_sec=elapsed_sec,
            output=output,
        )

    if "SAT" in output:
        model = parse_model_from_output(output) if need_model else None
        return CmsatResult(
            status="SAT",
            model=model,
            command=command,
            elapsed_sec=elapsed_sec,
            output=output,
        )

    return CmsatResult(
        status="ERROR",
        model=None,
        command=command,
        elapsed_sec=elapsed_sec,
        output=output,
    )


def run_cmsat(
    cnf_path: str,
    cmsat_path: str,
    threads: int,
    verb: int = 0,
    need_model: bool = True,
    extra_args: Optional[Sequence[str]] = None,
) -> CmsatResult:
    command = build_cmsat_command(
        cnf_path,
        cmsat_path,
        threads,
        verb=verb,
        need_model=need_model,
        extra_args=extra_args,
    )
    started_at = time.perf_counter()
    try:
        proc = subprocess.run(command, capture_output=True, text=True, check=False)
    except FileNotFoundError as exc:
        raise RuntimeError(f"Solver not found: {cmsat_path}") from exc
    elapsed_sec = time.perf_counter() - started_at
    output = proc.stdout + "\n" + proc.stderr
    result = result_from_cmsat_output(output, need_model, command, elapsed_sec)
    if result.status == "ERROR":
        raise RuntimeError(f"Cannot parse CryptoMiniSat output.\n{output.strip()}")
    return result


def decode_model_to_edge_colors(model: List[int], E: int, k: int) -> List[int]:
    max_color_var = E * k
    positives = {v for v in model if 0 < v <= max_color_var}

    edge_color = [-1] * E
    for e in range(E):
        for c in range(k):
            if var_color(e, c, k) in positives:
                edge_color[e] = c
                break
        if edge_color[e] < 0:
            raise RuntimeError("decode failed: some edge has no true color.")
    return edge_color

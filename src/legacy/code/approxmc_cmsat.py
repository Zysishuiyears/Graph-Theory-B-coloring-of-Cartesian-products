import os
import random
import subprocess
import tempfile
from typing import List, Tuple, Optional


def parse_cmsat_sat(stdout: str) -> Optional[bool]:
    # 返回 True(SAT) / False(UNSAT) / None(无法判断)
    if "UNSAT" in stdout:
        return False
    if "SAT" in stdout:
        return True
    return None


def run_cmsat_on_dimacs(dimacs_path: str, cmsat_path: str) -> bool:
    proc = subprocess.run([cmsat_path, dimacs_path], capture_output=True, text=True)
    sat = parse_cmsat_sat(proc.stdout + "\n" + proc.stderr)
    if sat is None:
        raise RuntimeError("Cannot parse CryptoMiniSat output.")
    return sat


def write_dimacs_with_xors(
    nvars: int,
    cnf_clauses: List[List[int]],
    xor_clauses: List[List[int]],
    out_path: str
) -> None:
    """
    cnf_clauses: 普通 CNF 子句，每个子句以 0 结束（我们写的时候加 0）
    xor_clauses: XOR 子句，每个是“文字列表”，我们写成一行：x<lit1> <lit2> ... 0
                在 CryptoMiniSat 里，这表示这些文字 XOR = True；
                若要 XOR=False，只需把其中一个文字取反（调用方可随机做）。
    """
    ncls = len(cnf_clauses) + len(xor_clauses)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(f"p cnf {nvars} {ncls}\n")
        for cls in cnf_clauses:
            f.write(" ".join(map(str, cls)) + " 0\n")
        for xcls in xor_clauses:
            f.write("x" + " ".join(map(str, xcls)) + " 0\n")


def random_xor_clause(nvars: int, rng: random.Random, p: float = 0.5) -> List[int]:
    """
    生成一个随机 XOR 子句（非空）。
    选变量子集：每个变量以概率 p 进入；
    然后以 1/2 概率翻转其中一个 literal，从而把 RHS 从 True 变 False（等价）。
    """
    lits = []
    for v in range(1, nvars + 1):
        if rng.random() < p:
            lits.append(v)

    # 避免空 XOR（空 XOR 在不同实现里语义麻烦），强制至少一个变量
    if not lits:
        lits.append(rng.randint(1, nvars))

    # 以 1/2 概率把 RHS 变成 False：翻转一个文字即可 :contentReference[oaicite:1]{index=1}
    if rng.random() < 0.5:
        t = rng.randrange(len(lits))
        lits[t] = -lits[t]
    return lits


def sat_with_r_xors(
    nvars: int,
    cnf_clauses: List[List[int]],
    r: int,
    cmsat_path: str,
    seed: int,
    p: float = 0.5
) -> bool:
    rng = random.Random(seed)
    xors = [random_xor_clause(nvars, rng, p=p) for _ in range(r)]
    with tempfile.NamedTemporaryFile(delete=False, suffix=".cnf", mode="w", encoding="utf-8") as tf:
        tmp_path = tf.name
    try:
        write_dimacs_with_xors(nvars, cnf_clauses, xors, tmp_path)
        return run_cmsat_on_dimacs(tmp_path, cmsat_path)
    finally:
        try:
            os.remove(tmp_path)
        except OSError:
            pass


def estimate_log2_model_count(
    nvars: int,
    cnf_clauses: List[List[int]],
    cmsat_path: str,
    trials: int = 25,
    p_in_xor: float = 0.5,
    seed: int = 1,
    r_max: Optional[int] = None
) -> Tuple[float, List[int]]:
    """
    近似估计 log2(#models)：
    - 对每次 trial，用二分找“最小 r，使 CNF+ r 个随机 XOR 变 UNSAT”的 r*
    - r* 大致 ~ log2(#models)
    返回：log2_est（用 trials 的中位数），以及所有 trial 的 r* 列表

    注意：这是简化版（非完整 ApproxMC），但在小规模结构 CNF 上通常够用。
    """
    if r_max is None:
        r_max = nvars  # 最多加到 nvars 个 XOR

    rstars = []
    base_seed = seed

    for t in range(trials):
        # 先指数扩展找上界
        lo, hi = 0, 1
        while hi <= r_max:
            sat = sat_with_r_xors(nvars, cnf_clauses, hi, cmsat_path, seed=base_seed + 100000*t + hi, p=p_in_xor)
            if not sat:
                break
            lo = hi
            hi *= 2

        if hi > r_max:
            hi = r_max

        # 如果加到 r_max 仍 SAT，说明解很多；把 r* 记为 r_max
        sat_at_hi = sat_with_r_xors(nvars, cnf_clauses, hi, cmsat_path, seed=base_seed + 100000*t + hi + 7, p=p_in_xor)
        if sat_at_hi:
            rstars.append(hi)
            continue

        # 二分：找最小 UNSAT
        L, R = lo, hi
        while L + 1 < R:
            mid = (L + R) // 2
            sat = sat_with_r_xors(nvars, cnf_clauses, mid, cmsat_path, seed=base_seed + 100000*t + mid + 13, p=p_in_xor)
            if sat:
                L = mid
            else:
                R = mid
        rstars.append(R)

    rstars_sorted = sorted(rstars)
    median_r = rstars_sorted[len(rstars_sorted)//2]
    return float(median_r), rstars

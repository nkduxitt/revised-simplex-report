#!/usr/bin/env python3
"""
run_experiments.py (updated)

Usage:
    python run_experiments.py [--repeats N] [--warmup M] [--results PATH]

Behavior:
- Produces: results/summary.csv and results/<problem>_<solver>_run<N>.json
- By default, "results" is created as a sibling to the directory that contains this script.
  (So if this file lives in "code/", results will be created next to "code/".)
"""
from __future__ import annotations
import argparse
import csv
import json
import logging
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence

import numbers
import numpy as np

# -------------------------
# Utilities
# -------------------------
logger = logging.getLogger("run_experiments")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def to_jsonable(obj: Any) -> Any:
    """
    Recursively convert numpy types (ndarray, numpy scalars) to Python built-ins
    so the structure becomes JSON-serializable.
    """
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer, np.bool_)):
        return obj.item()
    if isinstance(obj, (numbers.Number, str, bool)) or obj is None:
        return obj
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    try:
        return str(obj)
    except Exception:
        return repr(obj)


def safe_json_dump(obj: Any, path: Path, indent: int = 2) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf8") as f:
        json.dump(to_jsonable(obj), f, indent=indent, ensure_ascii=False)


# -------------------------
# Determine results directory
# -------------------------
def get_default_results_dir(script_path: Optional[Path] = None) -> Path:
    """
    Place 'results' directory as a sibling to the folder containing this script.
    Example:
        script: /home/user/project/code/run_experiments.py
        results -> /home/user/project/results
    """
    if script_path is None:
        script_path = Path(__file__).resolve()
    script_dir = script_path.resolve().parent
    parent = script_dir.parent
    return parent / "results"


# -------------------------
# Solver discovery & fallback
# -------------------------
# Try to import user's RevisedSimplex implementation if available.
try:
    import importlib

    RevisedSimplex = None
    try:
        spec = importlib.util.find_spec("code.revised_simplex")
        if spec:
            modu = importlib.import_module("code.revised_simplex")
            RevisedSimplex = getattr(modu, "RevisedSimplex", None)
    except Exception:
        RevisedSimplex = None
except Exception:
    RevisedSimplex = None


# Fallback small revised-simplex skeleton (keeps behavior similar to your previous fallback)
try:
    from scipy.linalg import lu_factor, lu_solve
except Exception:
    lu_factor = lu_solve = None  # will raise if used


class FallbackRevised:
    def __init__(
        self,
        A: np.ndarray,
        b: np.ndarray,
        c: np.ndarray,
        B_init: Sequence[int],
        tol: float = 1e-12,
        sense: str = "max",
    ):
        self.A = A.copy().astype(float)
        self.b = b.copy().astype(float)
        self.c = c.copy().astype(float)
        self.m, self.n = A.shape
        self.B_idx = list(B_init)
        self.N_idx = [j for j in range(self.n) if j not in self.B_idx]
        self.tol = tol
        self.sense = sense  # 'max' or 'min'
        self.rebuilds = 0
        self.iterations = 0
        self._rebuild()

    def _rebuild(self):
        if lu_factor is None:
            raise RuntimeError("scipy.linalg.lu_factor is required by the fallback solver.")
        B = self.A[:, self.B_idx]
        self.lu, self.piv = lu_factor(B)
        self.xB = lu_solve((self.lu, self.piv), self.b)
        self.cB = self.c[self.B_idx]
        self.rebuilds += 1

    def _reduced(self):
        # compute reduced costs and dual
        y = lu_solve((self.lu, self.piv), self.cB, trans=1)
        red = self.c[self.N_idx] - self.A[:, self.N_idx].T.dot(y)
        return red, y

    def iterate_once(self):
        red, y = self._reduced()
        # For maximization: enter if max(reduced) > tol, choose argmax.
        # For minimization: enter if min(reduced) < -tol, choose argmin.
        if self.sense == "max":
            entering_local = int(np.argmax(red))
            if red[entering_local] <= self.tol:
                return "optimal", None, None
        else:  # 'min'
            entering_local = int(np.argmin(red))
            if red[entering_local] >= -self.tol:
                return "optimal", None, None

        q = self.N_idx[entering_local]
        a_q = self.A[:, q]
        d = lu_solve((self.lu, self.piv), a_q)
        pos = d > self.tol
        if not np.any(pos):
            return "unbounded", q, None
        ratios = np.full(d.shape, np.inf)
        ratios[pos] = self.xB[pos] / d[pos]
        leave_local = int(np.argmin(ratios))
        p = self.B_idx[leave_local]
        # swap indices
        self.B_idx[leave_local] = q
        self.N_idx[entering_local] = p
        self._rebuild()
        self.iterations += 1
        return "continue", p, q

    def solve(self, max_iters: int = 1000):
        for _ in range(max_iters):
            status, p, q = self.iterate_once()
            if status == "optimal":
                x = np.zeros(self.n)
                for i, bi in enumerate(self.B_idx):
                    x[bi] = self.xB[i]
                return {
                    "status": "optimal",
                    "x": x,
                    "objective": float(self.c.dot(x)),
                    "iterations": self.iterations,
                    "rebuilds": self.rebuilds,
                    "B_idx": list(self.B_idx),
                }
            if status == "unbounded":
                return {"status": "unbounded", "entering": q, "iterations": self.iterations, "rebuilds": self.rebuilds}
        return {"status": "max_iters", "iterations": self.iterations, "rebuilds": self.rebuilds}


def make_solver(A: np.ndarray, b: np.ndarray, c: np.ndarray, B_init: Sequence[int], sense: str = "max") -> Callable[[], Dict[str, Any]]:
    """
    Return a zero-arg function that runs the chosen solver and returns the result dict.
    Passes 'sense' to the fallback solver and attempts to pass it to a user RevisedSimplex if available.
    """
    if RevisedSimplex is not None:
        def run():
            # try to instantiate RevisedSimplex with a sense kwarg if supported, else without
            try:
                solver = RevisedSimplex(A, b, c, B_init, sense=sense)
            except TypeError:
                solver = RevisedSimplex(A, b, c, B_init)
            return solver.solve()
        return run
    else:
        def run():
            solver = FallbackRevised(A, b, c, B_init, sense=sense)
            return solver.solve()
        return run


# -------------------------
# Problems (put them behind a factory so users can easily expand)
# -------------------------
@dataclass
class Problem:
    name: str
    A: np.ndarray
    b: np.ndarray
    c: np.ndarray
    B_init: List[int]
    sense: str = "max"  # 'max' or 'min' to indicate objective sense for SciPy wrapper


def build_problems() -> List[Problem]:
    problems: List[Problem] = []

    A1 = np.array([[1.0, 1.0, 1.0, 0.0], [1.0, 3.0, 0.0, 1.0]])
    b1 = np.array([4.0, 6.0])
    c1 = np.array([3.0, 2.0, 0.0, 0.0])
    B1 = [2, 3]
    problems.append(Problem("example1", A1, b1, c1, B1, sense="max"))

    A2 = np.array([[1, 1, -1, 0, 1, 0], [1, 2, 0, -1, 0, 1]], dtype=float)
    b2 = np.array([4.0, 5.0])
    # Phase I objective is a minimization of artificials a1+a2 -> sense='min'
    c2 = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 1.0])
    B2 = [4, 5]
    problems.append(Problem("example2_phaseI", A2, b2, c2, B2, sense="min"))

    return problems


# -------------------------
# SciPy baseline wrapper (respects objective sense)
# -------------------------
def run_scipy_highs(A: np.ndarray, b: np.ndarray, c: np.ndarray, sense: str = "max") -> Dict[str, Any]:
    """
    Run SciPy/HiGHS and return keys:
      - status: 'success'/'fail'/'scipy_missing'
      - x: solution vector (if any)
      - objective: objective value (c^T x using input c)
      - iterations: iteration count if provided by the solver (res.nit)
      - message: solver message

    The function respects `sense`: if sense=='max' we convert the objective to
    minimization by passing -c to linprog; if sense=='min' we pass c unchanged.
    """
    try:
        from scipy.optimize import linprog
        # convert objective according to sense
        c_min = -c if sense == "max" else c
        res = linprog(c_min, A_eq=A, b_eq=b, bounds=[(0, None)] * A.shape[1], method="highs")
        success = bool(getattr(res, "success", False))
        x = res.x if hasattr(res, "x") and res.x is not None else np.zeros(A.shape[1])
        obj = float(c.dot(x)) if x is not None else None
        iterations = getattr(res, "nit", None)
        message = getattr(res, "message", None)
        return {"status": "success" if success else "fail", "x": x, "objective": obj, "iterations": iterations, "message": message}
    except Exception as e:
        return {"status": "scipy_missing", "error": str(e)}


# -------------------------
# Helpers to make CSV rows complete
# -------------------------
def normalize_for_csv(val: Any) -> Any:
    """
    Convert None / nan to a readable placeholder '---'. Leave numbers as-is.
    """
    if val is None:
        return "---"
    if isinstance(val, float) and (np.isnan(val) or np.isinf(val)):
        return "---"
    if isinstance(val, (np.floating, np.integer)):
        return val.item()
    return val


# -------------------------
# Runner
# -------------------------
def run_all(repeats: int, warmup: int, results_dir: Path) -> None:
    results_dir.mkdir(parents=True, exist_ok=True)
    problems = build_problems()
    summary_rows: List[Dict[str, Any]] = []

    for prob in problems:
        logger.info("Running problem %s", prob.name)
        solver_fn = make_solver(prob.A, prob.b, prob.c, prob.B_init, prob.sense)

        # warmup runs
        for _ in range(warmup):
            try:
                _ = solver_fn()
            except Exception:
                logger.warning("Warmup run failed for %s", prob.name)
                traceback.print_exc()

        # timed runs
        times: List[float] = []
        raws: List[Dict[str, Any]] = []
        for i in range(repeats - warmup):
            t0 = time.perf_counter()
            try:
                out = solver_fn()
            except Exception as e:
                logger.error("Solver failed on problem %s: %s", prob.name, e)
                traceback.print_exc()
                out = {"status": "error", "error": str(e)}
            t1 = time.perf_counter()
            elapsed = t1 - t0
            times.append(elapsed)

            # SAFE handling when solver did NOT return an x vector
            x = out.get("x", None)
            if x is None:
                residual = None
                obj_val = out.get("objective", None)
                obj_val = obj_val if obj_val is not None else None
                B_idx = out.get("B_idx", None)
                condB = None
            else:
                x = np.asarray(x, dtype=float)
                residual = float(np.max(np.abs(prob.A.dot(x) - prob.b)))
                obj_val = out.get("objective", None)
                B_idx = out.get("B_idx", prob.B_init)
                try:
                    condB = float(np.linalg.cond(prob.A[:, B_idx])) if (B_idx is not None and len(B_idx) > 0) else None
                except Exception:
                    condB = None

            raw = {"elapsed": elapsed, "out": out, "residual": residual, "condB": condB}
            raws.append(raw)

            # write per-run json
            run_name = f"{prob.name}_revised_run{i+1}.json"
            safe_json_dump(
                {
                    "elapsed": elapsed,
                    "out": to_jsonable(out),
                    "objective": to_jsonable(obj_val),
                    "residual": residual,
                    "condB": to_jsonable(condB),
                },
                results_dir / run_name,
            )

        # summarise revised solver
        import statistics as stats

        med = float(stats.median(times)) if times else float("nan")
        std = float(stats.pstdev(times) if len(times) > 1 else 0.0)
        first_out = raws[0]["out"] if raws else {}
        # create row and normalize missing values to '---' placeholder
        row = {
            "problem": prob.name,
            "solver": "revised",
            "obj": normalize_for_csv(first_out.get("objective")) if first_out else "---",
            "residual": normalize_for_csv(raws[0]["residual"]) if raws and raws[0]["residual"] is not None else "---",
            "iterations": normalize_for_csv(first_out.get("iterations")),
            "rebuilds": normalize_for_csv(first_out.get("rebuilds")),
            "condB": normalize_for_csv(raws[0].get("condB")) if raws and raws[0].get("condB") is not None else "---",
            "time_median": med,
            "time_std": std,
            "trials": len(times),
        }
        summary_rows.append(row)

        # SciPy baseline
        sp = run_scipy_highs(prob.A, prob.b, prob.c, prob.sense)
        if sp.get("status") in ("success", "fail"):
            safe_json_dump({"scipy": sp}, results_dir / f"{prob.name}_scipy_run.json")

            times_sp: List[float] = []
            for _ in range(repeats - warmup):
                t0 = time.perf_counter()
                _ = run_scipy_highs(prob.A, prob.b, prob.c, prob.sense)
                t1 = time.perf_counter()
                times_sp.append(t1 - t0)
            med_sp = float(np.median(times_sp)) if times_sp else float("nan")
            std_sp = float(np.std(times_sp, ddof=1)) if len(times_sp) > 1 else 0.0

            sci_row = {
                "problem": prob.name,
                "solver": "scipy_highs",
                "obj": normalize_for_csv(float(sp.get("objective", np.nan))),
                "residual": normalize_for_csv(float(np.max(np.abs(prob.A.dot(sp.get("x", np.zeros(prob.A.shape[1]))) - prob.b)))),
                "iterations": normalize_for_csv(sp.get("iterations")),
                "rebuilds": "---",
                "condB": "---",
                "time_median": med_sp,
                "time_std": std_sp,
                "trials": len(times_sp),
            }
            summary_rows.append(sci_row)

    # write summary csv (normalize missing to '---' for all rows uniformly)
    keys = ["problem", "solver", "obj", "residual", "iterations", "rebuilds", "condB", "time_median", "time_std", "trials"]
    summary_path = results_dir / "summary.csv"
    with summary_path.open("w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, keys)
        w.writeheader()
        for r in summary_rows:
            out_row = {}
            for k in keys:
                v = r.get(k, "---")
                if isinstance(v, (np.floating, np.integer)):
                    v = v.item()
                if v is None or (isinstance(v, float) and (np.isnan(v) or np.isinf(v))):
                    v = "---"
                out_row[k] = v
            w.writerow(out_row)

    logger.info("Wrote %s and per-run json files to %s", summary_path.name, results_dir)


# -------------------------
# CLI
# -------------------------
def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run revised simplex experiments.")
    p.add_argument("--repeats", "-r", type=int, default=6, help="Total repeats per problem (including warmup).")
    p.add_argument("--warmup", "-w", type=int, default=1, help="Number of warmup runs (not counted in measurements).")
    p.add_argument(
        "--results",
        "-o",
        type=Path,
        default=None,
        help="Path to results folder. By default placed as sibling to the folder containing this script.",
    )
    return p.parse_args(argv)


def main():
    args = parse_args()
    if args.results is None:
        results_dir = get_default_results_dir()
    else:
        results_dir = args.results
    logger.info("Results directory: %s", results_dir)
    run_all(repeats=args.repeats, warmup=args.warmup, results_dir=results_dir)


if __name__ == "__main__":
    main()

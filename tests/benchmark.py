"""Kruskal vs Prim on the SAME graphs at three sizes.

Run from the repository root:
    python tests/benchmark.py            # writes tests/benchmark_results.md

Protocol
  * Graphs: random spanning tree (vertex i>0 attaches to a random j<i, so every
    graph is connected) plus random extra edges, integer weights in
    [-1_000_000, 1_000_000]; random.Random(SEED + n).  Parallel edges and
    self-loops can occur among the extras.
  * Correctness first: both results must be valid forests with equal cost.
  * Warm-up: one untimed call per method per size.
  * Timed region: ONLY algorithms.kruskal(request) / algorithms.prim(request).
    Graph generation, JSON parsing and output are outside the timer.  Each
    method's own preprocessing (Kruskal's sort + DSU, Prim's adjacency lists
    + heap) is inside it, because it is part of the algorithm.
  * Timer: time.perf_counter_ns(); gc.collect() before each timed call, the
    collector stays enabled; REPEATS timed runs per method per size, calls
    interleaved (K, P, K, P, ...) so drift affects both equally.
  * No compilation exists for Python.
"""

import os
import platform
import random
import statistics
import sys
import time
import gc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src", "Python"))
sys.path.insert(0, HERE)

import algorithms  # noqa: E402
from checkers import check_forest  # noqa: E402

SEED = 20260101
REPEATS = 5
SIZES = [(1_000, 4_000), (10_000, 40_000), (100_000, 250_000)]
W = 1_000_000


def make_graph(n, m, seed):
    rng = random.Random(seed)
    edges = [[i, rng.randrange(i), rng.randint(-W, W)] for i in range(1, n)]
    while len(edges) < m:
        edges.append([rng.randrange(n), rng.randrange(n), rng.randint(-W, W)])
    rng.shuffle(edges)
    return edges


def cpu_name():
    try:
        with open("/proc/cpuinfo", encoding="utf-8") as f:
            for line in f:
                if line.startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def timed(fn, request):
    gc.collect()
    t0 = time.perf_counter_ns()
    result = fn(request)
    return (time.perf_counter_ns() - t0) / 1e6, result  # milliseconds


def main():
    rows = []
    for n, m in SIZES:
        edges = make_graph(n, m, SEED + n)
        kreq = {"task": "kruskal", "n": n, "edges": edges}
        preq = {"task": "prim", "n": n, "edges": edges}
        # correctness before timing
        rk, rp = algorithms.kruskal(kreq), algorithms.prim(preq)
        check_forest(n, edges, rk)
        check_forest(n, edges, rp)
        assert rk["cost"] == rp["cost"], "Kruskal and Prim disagree on cost"
        # warm-up (untimed), then interleaved timed runs
        algorithms.kruskal(kreq)
        algorithms.prim(preq)
        kt, pt = [], []
        for _ in range(REPEATS):
            kt.append(timed(algorithms.kruskal, kreq)[0])
            pt.append(timed(algorithms.prim, preq)[0])
        rows.append((n, len(edges), rk["cost"], rk["components"], kt, pt))

    out = []
    out.append("| n | m | MSF cost | components | Kruskal runs (ms) | Kruskal median (ms) "
               "| Prim runs (ms) | Prim median (ms) |")
    out.append("|---:|---:|---:|---:|---|---:|---|---:|")
    for n, m, cost, comps, kt, pt in rows:
        fmt = lambda xs: ", ".join(f"{x:.1f}" for x in xs)
        out.append(f"| {n} | {m} | {cost} | {comps} | {fmt(kt)} | "
                   f"{statistics.median(kt):.1f} | {fmt(pt)} | {statistics.median(pt):.1f} |")
    env = [
        f"- Python: {platform.python_version()} ({platform.python_implementation()})",
        f"- OS: {platform.platform()}",
        f"- CPU: {cpu_name()} ({os.cpu_count()} logical cores visible)",
        f"- Seed: {SEED} (graph for size n uses random.Random({SEED} + n)); repeats per method: {REPEATS}",
    ]
    text = "\n".join(env) + "\n\n" + "\n".join(out) + "\n"
    print(text)
    with open(os.path.join(HERE, "benchmark_results.md"), "w", encoding="utf-8") as f:
        f.write(text)


if __name__ == "__main__":
    main()

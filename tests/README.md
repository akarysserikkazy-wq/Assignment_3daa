# Own tests

All commands run from the repository root.

| Command | What it does |
|---|---|
| `python -m unittest discover -s tests -v` | 49 automated tests: boundaries, full limits (n=100000, m/q=250000), deep chains, repeated mixed requests, alternative valid answers, rejected invalid witnesses, route-helper calls |
| `python tools/run.py --input tests/my_cases.jsonl` | runs the stored hand-written inputs through the launcher and prints actual results |
| `python tests/benchmark.py` | Kruskal vs Prim benchmark, writes `tests/benchmark_results.md` |

## Files

- `test_algorithms.py` – the tests. Expected answers are justified in comments or derived by an independent oracle (brute force, Warshall closure, Kahn, Bellman–Ford, planted structures).
- `checkers.py` – validators for non-unique answers (forest, SCC labels and condensation, topological order, shortest-path witnesses). They do not import `algorithms.py`.
- `my_cases.jsonl` – stored inputs, one JSON request per line.
- `my_cases_actual.jsonl` – actual outputs of `my_cases.jsonl`, in the same order.
- `benchmark.py`, `benchmark_results.md` – benchmark script and its latest results.

The expected answer for every stored case and the explanation of each test group are in section 4 of `../ANALYSIS.md`.

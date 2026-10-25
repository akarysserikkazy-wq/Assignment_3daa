# Assignment 3 – Analysis

Language: **Python 3.12** (`submission.json`: `{"language":"Python"}`).
Code: `src/Python/algorithms.py` (algorithms and route helper), `src/Python/main.py` (supplied JSON Lines adapter, unchanged).
Own tests: `tests/`.

## 1. Build, run and test commands

Python needs no compilation. Run everything from the repository root.

```sh
python tools/run.py --setup                      # supplied: launch check
python tools/run.py --check                      # supplied: 36 public cases + repeated batches
python -m unittest discover -s tests -v          # own tests (49 tests, about 7 s)
python src/Python/main.py < tests/my_cases.jsonl # run stored cases by hand
python tests/benchmark.py                        # Kruskal vs Prim, writes tests/benchmark_results.md
```

## 2. Algorithms

Notation: `n` vertices, `m` edges, `q` union operations, `α` the inverse Ackermann function.

### 2.1 DSU (`dsu`)

**Steps.** `parent[i]=i`, `size[i]=1`, `components=n`, `largest=1` (0 if `n=0`). For every pair `[u,v]`: find both roots (two-pass iterative find, the second pass points every node on the path directly at the root), attach the smaller tree under the larger, add the sizes, decrement `components` if the roots differed, update `largest`, append `[components, largest]`.

**Invariant.** `parent` encodes a forest whose trees are exactly the current connected components; `size[root]` is the true vertex count of the tree; `components` = number of roots; `largest` = max of `size` over roots.

**Correctness.** Union merges exactly two trees (or does nothing), so the partition after each step equals the connectivity of the edges processed so far. Rank and size are kept separately: I use size only, and `largest` is read from `size`, never from a rank.

**Boundary cases.** `n=0` → `states=[]`; `n=1` → `[1,1]`; self-unions and repeated unions → roots equal, nothing changes; no initial state is emitted.

**Complexity.** Time `O(n + q·α(n))`, memory `O(n + q)` (the `q` states are part of the result).

### 2.2 Kruskal (`kruskal`)

**Steps.** Sort the edge *ids* by weight (stable), scan them, and accept edge `(u,v,w)` iff `union(u,v)` merged two different trees. Cost is a Python integer sum; `components` = `n − accepted edges` as tracked by the DSU counter.

**Invariant (cut/cycle).** After processing the first `k` edges in sorted order, the accepted edges are a minimum spanning forest of those `k` edges.

**Correctness.** An edge whose endpoints are already connected closes a cycle, and a heaviest-or-tied edge of a cycle is never needed; an accepted edge is a lightest edge crossing the cut between its tree and the rest, so it belongs to some MSF (cut property). Negative weights need no special treatment because only the order matters.

**Boundary cases.** Self-loops: `union` returns False, never selected. Parallel edges: first of the cheapest wins. Ties: stable sort, any optimal forest is valid. Isolates counted as components. `n=0` → `cost=0, edge_ids=[], components=0`. IDs are original positions because ids, not rows, are sorted.

**Complexity.** `O(n + m log(m+1))` time (sort dominates; DSU is `O(m α(n))`), `O(n+m)` memory.

### 2.3 Prim (`prim`)

**Steps.** Build adjacency lists of `(w, edge_id, neighbour)` (self-loops skipped). For each unvisited vertex start a new component: mark it visited, push its edges to a lazy binary heap; repeatedly pop the lightest `(w,id,v)`; if `v` is already visited, discard (stale entry); otherwise visit `v`, add `w` and `id`, push the edges of `v` that lead to unvisited vertices.

**Invariant.** The visited set of the current component is a tree that is part of an MSF of that component, and the heap holds every edge leaving it (plus stale ones).

**Correctness.** The popped non-stale entry is a lightest edge crossing the cut (visited, unvisited), so it is safe by the cut property. Restarting at every unvisited vertex covers all components, so the result is a minimum spanning *forest*; the number of restarts is the component count.

**Boundary cases.** Same as Kruskal; additionally the restart loop is what handles disconnected graphs and isolates (an isolate starts a component, finds an empty heap and ends).

**Complexity.** Lazy heap: at most `2m` pushes, so `O(n + m log(m+1))` time, `O(n+m)` memory.

### 2.4 Kosaraju SCC (`scc`)

**Steps.** (1) Build `G` and the transpose `Gᵀ`. (2) Iterative DFS on `G` over all start vertices, recording vertices in increasing finish time. (3) Iterate vertices in decreasing finish time; every unlabelled vertex starts a DFS on `Gᵀ` that labels everything it reaches with a new label. (4) For every edge with different endpoint labels add the key `a·k + b` to a set (`k` = number of labels) and decode the keys to pairs.

**Invariant.** The first pass gives a vertex of a source component of the condensation the largest finish time among the unprocessed ones; in `Gᵀ` the edges between components are reversed, so a DFS from that vertex cannot leave its own component.

**Correctness.** Standard Kosaraju argument: if there is an edge between components `C → C'`, then `max finish(C) > max finish(C')`; processing in decreasing finish order on `Gᵀ` therefore peels components in an order in which each DFS stays inside one SCC and reaches all of it. The set makes every cross pair appear exactly once.

**Boundary cases.** `n=0` → both arrays empty; isolates are separate SCCs; self-loops and duplicate edges are harmless (self-loops are internal, duplicates collapse in the set); labels are `0..k−1`, nonnegative and < `n`.

**Complexity.** `O(n+m)` expected time (hash set), `O(n+m)` memory.

### 2.5 Topological order (`topo`)

**Steps.** Three-state iterative DFS (0 unvisited, 1 active, 2 finished) from every unvisited vertex; each stack frame is just the vertex plus a per-vertex neighbour pointer. An edge to a state-1 vertex means a cycle: return `acyclic=false, order=[]` immediately. Otherwise append a vertex when it finishes and return the reversed finish list.

**Invariant.** State-1 vertices are exactly the vertices on the DFS stack, i.e. one directed path; state-2 vertices have all descendants finished.

**Correctness.** A directed cycle contains an edge into an active vertex in any DFS (back edge) and conversely a back edge plus the tree path forms a cycle. In an acyclic graph every edge `u→v` finishes `v` before `u` (an edge to a finished vertex is harmless), so reversed finish order puts `u` before `v`.

**Boundary cases.** `n=0` → `acyclic=true, order=[]`; self-loop → cycle; duplicate edges fine; a cycle in a later component is found because every component is started; an edge to a finished vertex (diamond) is correctly *not* a cycle.

**Complexity.** `O(n+m)` time and memory; the saved pointers mean each adjacency list is scanned once in total.

### 2.6 DAG shortest paths (`dag`) and route helper

**Steps.** Build `targets` and `out_ids` (original edge ids). Run the same three-state DFS on the **whole** graph; any cycle (even unreachable from `source`) returns `acyclic=false`, `distance=[]`, `parent_edge=[]`. Otherwise `distance[source]=0`, everything else `None`, `parent_edge` all `−1`; go through the topological order, skip vertices whose distance is `None` (never relax from an unreachable vertex), and for each outgoing edge `eid=(u,v,w)` set `distance[v]=distance[u]+w`, `parent_edge[v]=eid` when strictly better (or `v` was `None`).

**Invariant.** When `u` is taken from the topological order all edges into `u` come from earlier vertices, so `distance[u]` is final.

**Correctness.** By induction on the topological order each reachable vertex gets the minimum over all incoming edges from reachable predecessors, which is its shortest distance; negative weights are fine because no vertex is revisited. `parent_edge[v]` is the id of the edge that produced the final value, so it ends at `v`, is tight, and following it strictly decreases the topological position and reaches the source.

**Boundary cases.** `n=1` → `[0]`, `[−1]`; unreachable → `None` and `−1` (a finite distance of `−1` stays distinguishable because unreachable is `None`); self-loop → cyclic; parallel and tied edges keep the first strictly best edge.

**Complexity.** `O(n+m)` time and memory.

**Route helper** `reconstruct_route(source, target, parent_edge, edges)`: walk `target → source` through `edges[parent_edge[v]][0]`, reverse. `O(L)` for `L` vertices. Results: `[source]` for the source, the vertex list for reachable targets, **`[]` for an unreachable target** (documented choice) and `[]` for an empty `parent_edge` (cyclic result) or an inconsistent parent array (edge not ending at the vertex, or a loop). JSON `distance` still uses `null` for unreachable.

Actual calls (output copied from a run):

| Graph | Call | Result | Edge ids / check |
|---|---|---|---|
| `edges=[[0,1,3],[0,2,8],[1,2,-7],[2,3,4]]`, `n=5`, `source=0`; `parent_edge=[-1,0,2,3,-1]` | `reconstruct_route(0,3,…)` | `[0,1,2,3]` | edges 0,2,3: 0→1→2→3, weights 3−7+4 = 0 = `distance[3]` |
| same | `reconstruct_route(0,0,…)` | `[0]` | source |
| same | `reconstruct_route(0,4,…)` | `[]` | vertex 4 unreachable (`distance=null`) |
| `edges=[[0,1,4],[0,1,4],[0,2,1],[2,1,3],[1,3,0]]`, `n=4`, `source=0`; `parent_edge=[-1,0,2,4]` | `reconstruct_route(0,1,…)` | `[0,1]` | three tied routes cost 4 (edges 0, 1 or 2+3); the first strictly best edge, id 0, was kept |
| same | `reconstruct_route(0,3,…)` | `[0,1,3]` | edges 0 and 4 (direction 0→1, 1→3 verified), 4+0 = 4 |
| `edges=[[0,1,2],[0,1,-2],[1,2,1]]`, `n=3` | `reconstruct_route(0,2,…)` with `parent_edge=[-1,1,2]` | `[0,1,2]` | parallel edges: the cheaper edge 1 (−2) is chosen, −2+1 = −1 = `distance[2]` |

## 3. State reset and deep DFS

**State reset.** Every function builds all of its data (`parent`, adjacency lists, heap, `visited`, `distance`, …) from the request and keeps nothing at module level; there are no caches or globals. `tests/test_algorithms.py::TestProcess` sends 15 mixed requests (all six tasks, `ping`, cyclic then acyclic, `n=0`, a DAG whose source has no outgoing edges) to `main.py` in one process three times and checks that the three batches give identical results and one JSON line per request.

**Deep DFS.** No recursion is used. First-pass DFS, second-pass DFS, the three-state DFS for `topo`/`dag` and DSU `find` are loops. A DFS keeps a stack of vertices and a per-vertex pointer `pointer[u]` to the next unread neighbour, so a suspended vertex resumes where it stopped and the total work is `O(n+m)` instead of rescanning. Chains and cycles of 100 000 vertices are tested for `scc`, `topo`, `dag` and `dsu`; Python's recursion limit is never touched.

## 4. Own tests (`tests/`)

Files: `test_algorithms.py` (49 tests), `checkers.py` (independent validators; they never import `algorithms.py`), `my_cases.jsonl` (stored inputs), `my_cases_actual.jsonl` (actual outputs from `python src/Python/main.py < tests/my_cases.jsonl`), `benchmark.py`, `benchmark_results.md`.

| Area | What is checked and why the expected answer is justified |
|---|---|
| Examples | All README examples are reproduced (hand-derived, e.g. DAG: 0, 3, 3−7=−4, −4+4=0, vertex 4 unreachable). |
| Boundaries | `n=0`, `n=1`, isolates, empty edge/operation lists, self-unions, self-loops, parallel and tied edges, negative weights, cycles unreachable from the DAG source, a finite distance of −1 versus `null`. |
| 64-bit | Forest cost 3·10⁹ and −3·10⁹; chain of 100 000 vertices with weight 10⁹ has cost 99 999 000 000 000 = (n−1)·10⁹ for MST and DAG distance. |
| Full limits | DSU `n=100000, q=250000` (first `n−1` operations form a chain, so `states[i]=[n−1−i, i+2]`); Kruskal/Prim random `n=100000, m=250000, |w|≤10⁹` (valid forests, equal cost and components); SCC with planted components (cycle blocks plus forward block edges, so the partition and the exact condensation are known by construction); random DAG for `topo` and `dag`. |
| Deep chains | 100 000-vertex chains for `scc`, `topo` (order of the reversed chain is unique), `dag`, a 100 000-cycle for `scc`, a chain closed into a cycle for `topo`, a far unreachable cycle for `dag`. |
| Randomised oracles | MST vs brute force over all subsets (250 graphs); SCC vs Warshall closure (300); topo vs Kahn (400); DAG vs Bellman–Ford (500, including cyclic graphs); DSU vs naive relabelling (100). |
| Non-unique answers | Checked by validity: forest feasibility and optimum cost; SCC partition and exact condensation; topological precedence; DAG feasibility + tight parent paths. Tests show that a valid *alternative* answer is accepted (other spanning tree of a tie, README labels `[5,5,7,9]`, orders `[0,1,2,3]` and `[1,0,2,3]`) and that *invalid witnesses are rejected* (cycle in forest, non-spanning, repeated id, wrong cost; merged/split SCC labels, duplicate or internal condensation pair; backward topological edge, missing vertex; parent edge not tight, ending elsewhere, given to an unreachable vertex, −1 instead of `null`). |
| Mutation check | I also broke the solution on purpose (DAG cycle check only for reachable vertices, a finished vertex treated as a cycle, duplicate condensation pairs) and confirmed that the test suite fails in each case. |
| Process level | Repeated mixed batch in a single `main.py` process (see section 3). |

Stored cases and their justification (actual outputs are in `tests/my_cases_actual.jsonl`):

| Case | Expected (justification) | Actual |
|---|---|---|
| `dsu n=4 [[0,1],[1,2],[0,2],[3,3]]` | 3 comps/size 2; 2/3; repeat and self-union change nothing | `[[3,2],[2,3],[2,3],[2,3]]` |
| `kruskal`/`prim n=4 [[0,1,4],[1,2,2],[0,2,3]]` | edges 1 and 2, cost 5, vertex 3 isolated → 2 components | `cost 5, ids [1,2]` / `[2,1]`, components 2 |
| `kruskal n=3 [[0,0,-5],[0,1,7],[1,0,7],[1,1,-9]]` | loops never used; one of the two parallel 7-edges; vertex 2 isolated | `cost 7, ids [1], components 2` |
| `prim n=4` chain of three 10⁹ edges | 3·10⁹ exceeds 32 bits | `cost 3000000000` |
| `scc n=4 [[0,1],[1,0],[1,2],[0,2],[2,3]]` | {0,1}, {2}, {3}; pairs {0,1}→{2} and {2}→{3}, each once | `component [0,0,1,2]`, `[[0,1],[1,2]]` |
| `scc n=3 [[0,0],[0,1],[0,1]]` | loop is internal; the duplicate edge yields one pair | one pair |
| `topo` 0→2, 1→2, 2→3 | order with 0 and 1 before 2 before 3 | `[1,0,2,3]` |
| `topo` with cycle 2→3→4→2 | cycle in second component → `false`, `[]` | `acyclic false, order []` |
| `dag` README example | `[0,3,-4,0,null]`, parents `[-1,0,2,3,-1]` | identical |
| `dag` cycle 2↔3 unreachable from source 0 | whole-graph rule → `false`, `[]`, `[]` | identical |
| `dag` `[[0,1,-1],[2,1,5]]` source 0 | distance −1 is finite, vertex 2 unreachable → `null` | `[0,-1,null]`, `[-1,0,-1]` |

## 5. Benchmark: Kruskal vs Prim on identical graphs

**Setup.** Script: `tests/benchmark.py` (`python tests/benchmark.py`); results are also saved in `tests/benchmark_results.md`. Graphs: a random spanning tree (so each graph is connected) plus random extra edges (self-loops and parallel edges possible), weights uniform in [−10⁶, 10⁶], generated by `random.Random(20260101 + n)`; both methods receive the very same `edges` list. Sizes: `(n, m)` = (1 000, 4 000), (10 000, 40 000), (100 000, 250 000). Before timing, both results are validated as forests with equal cost.
**Timer and boundary.** `time.perf_counter_ns()`; one untimed warm-up call per method and size; 5 timed runs per method and size, interleaved K, P, K, P…; `gc.collect()` before each run, GC left enabled. The timed region is only the call `algorithms.kruskal(request)` / `algorithms.prim(request)`: graph generation, JSON parsing and output are excluded for both methods; each method's own preprocessing (Kruskal: sort and DSU; Prim: adjacency lists and heap) is included. There is no compilation in Python.

**Environment of the run below.** Python 3.12.3 (CPython), Linux 6.18 x86_64, Intel Xeon @ 2.30 GHz, 1 logical core visible.

| n | m | MSF cost | Kruskal runs (ms) | Kruskal median (ms) | Prim runs (ms) | Prim median (ms) |
|---:|---:|---:|---|---:|---|---:|
| 1 000 | 4 000 | −682 182 526 | 1.8, 1.7, 1.7, 1.7, 1.7 | **1.7** | 2.8, 2.8, 2.7, 2.6, 2.8 | **2.8** |
| 10 000 | 40 000 | −7 019 931 301 | 29.2, 37.4, 42.8, 36.8, 37.0 | **37.0** | 86.8, 119.1, 85.8, 85.7, 84.6 | **85.8** |
| 100 000 | 250 000 | −53 176 181 339 | 482.2, 417.0, 493.7, 347.0, 342.5 | **417.0** | 1149.9, 981.5, 1113.4, 769.8, 731.7 | **981.5** |

**Observations (measured).** On every size Kruskal was faster than Prim in this implementation: about 1.6× at n=1 000, 2.3× at n=10 000 and 2.4× at n=100 000 by median. Both grew faster than linearly in `m` (m×62.5 gave about ×245 for Kruskal and ×350 for Prim), which is compatible with the `O(m log m)` bounds but is also affected by cache and memory effects.

**Asymptotic argument (separate from the measurements).** Both are `O(n + m log(m+1))`: Kruskal sorts `m` edges, Prim performs up to `2m` heap pushes and pops. The asymptotics predict similar growth, not which one is faster.

**Limitations.** One graph per size and one graph family (sparse, random, connected); no dense graphs, so Prim with an indexed heap or Kruskal on dense graphs is not explored. The constants are specific to CPython: Kruskal's sort runs in C (Timsort) whereas Prim does heap operations on Python tuples and builds `2m` adjacency tuples, so the gap should not be read as an algorithmic one. The single visible core and a shared machine make timings noisy (the 100 000-vertex Kruskal runs vary from 342 to 494 ms), which is why medians are reported.

## 6. Sources and AI assistance

**[STUDENT TO COMPLETE before submitting – state what is true for your work.]**
Sources: lecture notes / the course textbook and the assignment README (graph algorithm definitions: DSU, Kruskal, Prim, Kosaraju, DFS topological order, DAG relaxation).
AI assistance: *Claude (Anthropic)* was used to draft the implementation in `src/Python/algorithms.py`, the tests in `tests/`, the benchmark script and this document. My own verification: *[describe what you actually did: read and understood each algorithm, ran the commands in section 1, ran the benchmark on my machine, etc.]*. I am responsible for the correctness and originality of the submission.

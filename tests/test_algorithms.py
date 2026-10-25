"""Own tests for Assignment 3 (run: python -m unittest discover -s tests -v).

Expected answers are justified in comments or derived independently
(brute force, Warshall closure, Kahn, Bellman-Ford, planted structures).
Non-unique answers are checked by validity, never by equality to one answer.
"""

import json
import os
import random
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "src", "Python")
sys.path.insert(0, SRC)
sys.path.insert(0, HERE)

import algorithms as A  # noqa: E402
import main as M  # noqa: E402
from checkers import (  # noqa: E402
    bellman_ford, brute_force_msf_cost, check_dag, check_forest, check_scc,
    check_scc_condensation, check_topo, count_components)

MAXN, MAXM, MAXQ = 100_000, 250_000, 250_000
BIG = 10**9


def req(task, n, **kw):
    r = {"task": task, "n": n}
    r.update(kw)
    return r


def random_edges(rng, n, m, weighted=False, lo=-20, hi=20):
    out = []
    for _ in range(m):
        e = [rng.randrange(n), rng.randrange(n)]
        if weighted:
            e.append(rng.randint(lo, hi))
        out.append(e)
    return out


def random_dag_edges(rng, n, m, weighted=True, lo=-20, hi=20):
    perm = list(range(n))
    rng.shuffle(perm)
    out = []
    for _ in range(m):
        a, b = rng.sample(range(n), 2) if n > 1 else (0, 0)
        if a > b:
            a, b = b, a
        e = [perm[a], perm[b]]
        if weighted:
            e.append(rng.randint(lo, hi))
        out.append(e)
    return out


# ------------------------------------------------------------------ DSU ----

class TestDSU(unittest.TestCase):
    def test_readme_example(self):
        # 4 singletons; {0,1}; {0,1,2}; 0-2 repeated; self-union 3-3: no change.
        r = M.handle(req("dsu", 4, operations=[[0, 1], [1, 2], [0, 2], [3, 3]]))
        self.assertEqual(r, {"states": [[3, 2], [2, 3], [2, 3], [2, 3]]})

    def test_empty_and_single(self):
        self.assertEqual(A.dsu(req("dsu", 0, operations=[])), {"states": []})
        # n=1: largest=1 always; a self-union changes nothing.
        self.assertEqual(A.dsu(req("dsu", 1, operations=[[0, 0]])),
                         {"states": [[1, 1]]})
        self.assertEqual(A.dsu(req("dsu", 5, operations=[])), {"states": []})

    def test_largest_tracks_size_not_rank(self):
        # Two pairs merge to size 4 while rank is only 2.
        r = A.dsu(req("dsu", 6, operations=[[0, 1], [2, 3], [0, 2], [4, 5]]))
        self.assertEqual(r["states"], [[5, 2], [4, 2], [3, 4], [2, 4]])

    def test_random_against_naive_relabelling(self):
        rng = random.Random(1)
        for _ in range(100):
            n = rng.randint(1, 40)
            ops = random_edges(rng, n, rng.randint(0, 60))
            label = list(range(n))
            expected = []
            for u, v in ops:
                a, b = label[u], label[v]
                if a != b:
                    label = [a if x == b else x for x in label]
                expected.append([len(set(label)),
                                 max(label.count(x) for x in set(label))])
            self.assertEqual(A.dsu(req("dsu", n, operations=ops))["states"],
                             expected)

    def test_full_limit_chain_then_repeats(self):
        # Ops 0..n-2 join the chain: after op i there are n-1-i components and
        # the largest has i+2 vertices.  The remaining ops are repeats.
        rng = random.Random(2)
        ops = [[i, i + 1] for i in range(MAXN - 1)]
        ops += [[rng.randrange(MAXN), rng.randrange(MAXN)]
                for _ in range(MAXQ - len(ops))]
        self.assertEqual(len(ops), MAXQ)
        states = A.dsu(req("dsu", MAXN, operations=ops))["states"]
        self.assertEqual(len(states), MAXQ)
        for i in (0, 1, 1000, MAXN - 2):
            self.assertEqual(states[i], [MAXN - 1 - i, i + 2])
        self.assertTrue(all(s == [1, MAXN] for s in states[MAXN - 2:]))

    def test_deep_union_chain_in_reverse_order(self):
        # Long find paths if compression/union-by-size were missing.
        ops = [[i + 1, i] for i in range(MAXN - 2, -1, -1)]
        states = A.dsu(req("dsu", MAXN, operations=ops))["states"]
        self.assertEqual(states[-1], [1, MAXN])


# ------------------------------------------------------------------ MST ----

class TestForest(unittest.TestCase):
    solvers = (("kruskal", A.kruskal), ("prim", A.prim))

    def both(self, n, edges, expected_cost=None):
        results = {}
        for name, fn in self.solvers:
            r = fn(req(name, n, edges=edges))
            check_forest(n, edges, r, expected_cost)
            results[name] = r
        self.assertEqual(results["kruskal"]["cost"], results["prim"]["cost"])
        return results

    def test_readme_example(self):
        # Edges 1 (w2) and 2 (w3) are cheapest; vertex 3 isolated -> 2 comps.
        e = [[0, 1, 4], [1, 2, 2], [0, 2, 3]]
        r = self.both(4, e, 5)
        self.assertEqual(r["kruskal"]["components"], 2)
        self.assertEqual(sorted(r["kruskal"]["edge_ids"]), [1, 2])

    def test_empty_and_isolates(self):
        for name, fn in self.solvers:
            self.assertEqual(fn(req(name, 0, edges=[])),
                             {"cost": 0, "edge_ids": [], "components": 0})
            self.assertEqual(fn(req(name, 3, edges=[])),
                             {"cost": 0, "edge_ids": [], "components": 3})
            self.assertEqual(fn(req(name, 1, edges=[])),
                             {"cost": 0, "edge_ids": [], "components": 1})

    def test_self_loops_never_selected(self):
        # A negative self-loop is tempting but cannot be in a forest.
        self.both(2, [[0, 0, -5], [0, 1, 7], [1, 1, -9]], 7)
        self.both(1, [[0, 0, -5]], 0)

    def test_parallel_edges_and_ties(self):
        # Parallel 0-1 edges (ids 0,1,2): cheapest is id 1 (w=-3).
        r = self.both(2, [[0, 1, 4], [1, 0, -3], [0, 1, -3]], -3)
        self.assertIn(r["kruskal"]["edge_ids"], ([1], [2]))
        # Square with all-equal weights: every 3-edge subset is optimal.
        self.both(4, [[0, 1, 1], [1, 2, 1], [2, 3, 1], [3, 0, 1]], 3)

    def test_negative_weights_all_taken_in_a_tree(self):
        self.both(4, [[0, 1, -1], [1, 2, -2], [2, 3, -3]], -6)

    def test_totals_beyond_32_bits(self):
        # 3 edges of 1e9 = 3e9 > 2^31; and the negative counterpart.
        self.both(4, [[0, 1, BIG], [1, 2, BIG], [2, 3, BIG]], 3 * BIG)
        self.both(4, [[0, 1, -BIG], [1, 2, -BIG], [2, 3, -BIG]], -3 * BIG)

    def test_disconnected_forest_restarts(self):
        # Components {0,1,2}, {3,4}, {5}: Prim must restart twice.
        e = [[0, 1, 5], [1, 2, 1], [0, 2, 2], [3, 4, -4]]
        r = self.both(6, e, 3 + -4)
        self.assertEqual(r["prim"]["components"], 3)

    def test_random_small_against_brute_force(self):
        rng = random.Random(3)
        for _ in range(250):
            n = rng.randint(1, 6)
            e = random_edges(rng, n, rng.randint(0, 9), weighted=True,
                             lo=-3, hi=3)
            self.both(n, e, brute_force_msf_cost(n, e))

    def test_full_limit_random_agree(self):
        rng = random.Random(4)
        n, m = MAXN, MAXM
        e = random_edges(rng, n, m, weighted=True, lo=-BIG, hi=BIG)
        k = A.kruskal(req("kruskal", n, edges=e))
        p = A.prim(req("prim", n, edges=e))
        check_forest(n, e, k)
        check_forest(n, e, p)
        self.assertEqual(k["cost"], p["cost"])
        self.assertEqual(k["components"], p["components"])
        self.assertEqual(k["components"], count_components(n, e))

    def test_deep_chain_big_total(self):
        # Path of 100000 vertices: every edge is forced; cost 99999 * 1e9.
        n = MAXN
        e = [[i, i + 1, BIG] for i in range(n - 1)]
        self.both(n, e, (n - 1) * BIG)
        self.assertEqual((n - 1) * BIG, 99_999_000_000_000)

    def test_checker_accepts_alternative_and_rejects_invalid(self):
        e = [[0, 1, 1], [1, 2, 1], [0, 2, 1]]  # triangle, any 2 edges optimal
        good = {"cost": 2, "edge_ids": [0, 1], "components": 1}
        alt = {"cost": 2, "edge_ids": [2, 1], "components": 1}
        check_forest(3, e, good, 2)
        check_forest(3, e, alt, 2)
        bad = [
            {"cost": 3, "edge_ids": [0, 1, 2], "components": 1},   # cycle
            {"cost": 1, "edge_ids": [0], "components": 1},         # not spanning
            {"cost": 2, "edge_ids": [0, 0], "components": 1},      # repeated id
            {"cost": 5, "edge_ids": [0, 1], "components": 1},      # wrong cost
            {"cost": 2, "edge_ids": [0, 1], "components": 2},      # wrong comps
        ]
        for b in bad:
            with self.assertRaises(AssertionError):
                check_forest(3, e, b, 2)
        with self.assertRaises(AssertionError):  # non-optimal forest
            check_forest(3, [[0, 1, 1], [1, 2, 5], [0, 2, 1]],
                         {"cost": 6, "edge_ids": [0, 1], "components": 1}, 2)
        with self.assertRaises(AssertionError):  # self-loop chosen
            check_forest(2, [[0, 0, 1], [0, 1, 1]],
                         {"cost": 2, "edge_ids": [0, 1], "components": 1}, 1)


# ------------------------------------------------------------------ SCC ----

def check_scc_planted(n, edges, result, block):
    """Large-graph check: labels induce exactly the planted partition and the
    condensation equals the planted cross-block pairs."""
    comp = result["component"]
    assert len(comp) == n
    label_to_block, block_to_label = {}, {}
    for v in range(n):
        if label_to_block.setdefault(comp[v], block[v]) != block[v]:
            raise AssertionError("one label covers several planted blocks")
        if block_to_label.setdefault(block[v], comp[v]) != comp[v]:
            raise AssertionError("one planted block got several labels")
    pairs = [tuple(p) for p in result["condensation"]]
    assert len(set(pairs)) == len(pairs), "duplicate condensation pair"
    expected = {(block_to_label[block[u]], block_to_label[block[v]])
                for u, v in edges if block[u] != block[v]}
    assert set(pairs) == expected, "condensation mismatch"


class TestSCC(unittest.TestCase):
    def test_readme_example_and_alternative_labels(self):
        e = [[0, 1], [1, 0], [1, 2], [0, 2], [2, 3]]
        check_scc(4, e, A.scc(req("scc", 4, edges=e)))
        # README's labels are another valid answer.
        check_scc(4, e, {"component": [5, 5, 7, 9],
                         "condensation": [[5, 7], [7, 9]]})

    def test_empty_isolates_loops_duplicates(self):
        self.assertEqual(A.scc(req("scc", 0, edges=[])),
                         {"component": [], "condensation": []})
        r = A.scc(req("scc", 3, edges=[]))
        self.assertEqual(len(set(r["component"])), 3)
        self.assertEqual(r["condensation"], [])
        e = [[0, 0], [0, 1], [0, 1], [0, 1], [1, 1]]  # loops + 3 duplicates
        r = A.scc(req("scc", 2, edges=e))
        check_scc(2, e, r)
        self.assertEqual(len(r["condensation"]), 1)  # exactly once

    def test_single_big_cycle_and_two_cycles_with_bridge(self):
        e = [[0, 1], [1, 2], [2, 0], [3, 4], [4, 3], [2, 3]]
        r = A.scc(req("scc", 5, edges=e))
        check_scc(5, e, r)
        self.assertEqual(len(set(r["component"])), 2)

    def test_random_small_against_closure(self):
        rng = random.Random(5)
        for _ in range(300):
            n = rng.randint(1, 9)
            e = random_edges(rng, n, rng.randint(0, 16))
            check_scc(n, e, A.scc(req("scc", n, edges=e)))

    def test_deep_chain_all_singletons(self):
        n = MAXN
        e = [[i, i + 1] for i in range(n - 1)]
        r = A.scc(req("scc", n, edges=e))
        self.assertEqual(len(set(r["component"])), n)
        self.assertEqual(len(r["condensation"]), n - 1)

    def test_deep_cycle_single_component(self):
        n = MAXN
        e = [[i, (i + 1) % n] for i in range(n)]
        r = A.scc(req("scc", n, edges=e))
        self.assertEqual(len(set(r["component"])), 1)
        self.assertEqual(r["condensation"], [])

    def test_full_limit_planted_components(self):
        rng = random.Random(6)
        n = MAXN
        sizes, left = [], n
        while left:
            s = min(left, rng.randint(1, 40))
            sizes.append(s)
            left -= s
        perm = list(range(n))
        rng.shuffle(perm)
        block, members, p = [0] * n, [], 0
        for b, s in enumerate(sizes):
            vs = perm[p:p + s]
            p += s
            members.append(vs)
            for v in vs:
                block[v] = b
        edges = []
        for vs in members:  # a directed cycle makes each block one SCC
            if len(vs) > 1:
                edges += [[vs[i], vs[(i + 1) % len(vs)]] for i in range(len(vs))]
        k = len(sizes)
        while len(edges) < MAXM:  # forward block-to-block edges only: a DAG
            a, b = sorted(rng.sample(range(k), 2))
            edges.append([rng.choice(members[a]), rng.choice(members[b])])
        self.assertEqual(len(edges), MAXM)
        check_scc_planted(n, edges, A.scc(req("scc", n, edges=edges)), block)

    def test_checker_rejects_invalid_answers(self):
        e = [[0, 1], [1, 0], [1, 2]]
        merged = {"component": [0, 0, 0], "condensation": []}
        split = {"component": [0, 1, 2], "condensation": [[0, 1], [1, 0], [1, 2]]}
        dup = {"component": [0, 0, 1], "condensation": [[0, 1], [0, 1]]}
        internal = {"component": [0, 0, 1], "condensation": [[0, 1], [0, 0]]}
        missing = {"component": [0, 0, 1], "condensation": []}
        reversed_ = {"component": [0, 0, 1], "condensation": [[1, 0]]}
        for bad in (merged, split, dup, internal, missing, reversed_):
            with self.assertRaises(AssertionError):
                check_scc(3, e, bad)
        check_scc(3, e, {"component": [4, 4, 2], "condensation": [[4, 2]]})


# ----------------------------------------------------------------- topo ----

class TestTopo(unittest.TestCase):
    def run_check(self, n, e):
        r = A.topo(req("topo", n, edges=e))
        check_topo(n, e, r)
        return r

    def test_readme_example_and_alternatives(self):
        e = [[0, 2], [1, 2], [2, 3]]
        r = self.run_check(4, e)
        self.assertTrue(r["acyclic"])
        check_topo(4, e, {"acyclic": True, "order": [1, 0, 2, 3]})
        check_topo(4, e, {"acyclic": True, "order": [0, 1, 2, 3]})

    def test_empty_single_and_isolates(self):
        self.assertEqual(A.topo(req("topo", 0, edges=[])),
                         {"acyclic": True, "order": []})
        self.run_check(1, [])
        self.run_check(5, [])

    def test_cycles(self):
        for n, e in [(1, [[0, 0]]), (2, [[0, 1], [1, 0]]),
                     (3, [[0, 1], [1, 2], [2, 0]]),
                     # acyclic first component, cycle in the second one
                     (5, [[0, 1], [2, 3], [3, 4], [4, 2]])]:
            r = self.run_check(n, e)
            self.assertEqual(r, {"acyclic": False, "order": []})

    def test_finished_vertex_is_not_a_cycle(self):
        # Diamond: DFS reaches 3 via 1, later sees 2->3 to a FINISHED vertex.
        e = [[0, 1], [0, 2], [1, 3], [2, 3]]
        self.assertTrue(self.run_check(4, e)["acyclic"])
        # Duplicate edges / forward edge in a DAG
        e = [[0, 1], [0, 1], [1, 2], [0, 2]]
        self.assertTrue(self.run_check(3, e)["acyclic"])

    def test_random_small_against_kahn(self):
        rng = random.Random(7)
        for _ in range(400):
            n = rng.randint(1, 8)
            self.run_check(n, random_edges(rng, n, rng.randint(0, 10)))

    def test_deep_chain_with_reversed_ids(self):
        n = MAXN
        e = [[i + 1, i] for i in range(n - 1)]  # edges point to smaller ids
        r = self.run_check(n, e)
        self.assertEqual(r["order"], list(range(n - 1, -1, -1)))  # unique order

    def test_deep_chain_closed_into_cycle(self):
        n = MAXN
        e = [[i, i + 1] for i in range(n - 1)] + [[n - 1, 0]]
        self.assertEqual(A.topo(req("topo", n, edges=e)),
                         {"acyclic": False, "order": []})

    def test_full_limit_random_dag(self):
        rng = random.Random(8)
        n, m = MAXN, MAXM
        e = random_dag_edges(rng, n, m, weighted=False)
        self.run_check(n, e)

    def test_checker_rejects_invalid_answers(self):
        e = [[0, 2], [1, 2], [2, 3]]
        for bad in ({"acyclic": True, "order": [2, 0, 1, 3]},     # backward edge
                    {"acyclic": True, "order": [0, 1, 2]},        # missing vertex
                    {"acyclic": True, "order": [0, 0, 2, 3]},     # repeated
                    {"acyclic": False, "order": []}):             # wrong verdict
            with self.assertRaises(AssertionError):
                check_topo(4, e, bad)
        with self.assertRaises(AssertionError):  # cyclic graph, non-empty order
            check_topo(2, [[0, 1], [1, 0]], {"acyclic": False, "order": [0, 1]})


# ------------------------------------------------------------------ DAG ----

class TestDAG(unittest.TestCase):
    def run_check(self, n, source, e, reference=True):
        r = A.dag(req("dag", n, source=source, edges=e))
        ref = None
        if reference and kahn_acyclic(n, e):
            ref = bellman_ford(n, e, source)
        check_dag(n, source, e, r, ref)
        return r

    def test_readme_example_and_route_helper(self):
        e = [[0, 1, 3], [0, 2, 8], [1, 2, -7], [2, 3, 4]]
        r = self.run_check(5, 0, e)
        # Hand-derived: d = 0, 3, 3-7=-4, -4+4=0, vertex 4 unreachable.
        self.assertEqual(r["distance"], [0, 3, -4, 0, None])
        self.assertEqual(r["parent_edge"], [-1, 0, 2, 3, -1])
        p = r["parent_edge"]
        # reachable target: route 0-1-2-3 over edges 0,2,3, cost 3-7+4 = 0
        self.assertEqual(A.reconstruct_route(0, 3, p, e), [0, 1, 2, 3])
        # source: [source]
        self.assertEqual(A.reconstruct_route(0, 0, p, e), [0])
        # unreachable target: documented result is []
        self.assertEqual(A.reconstruct_route(0, 4, p, e), [])
        # cost of the reconstructed route equals the reported distance
        route = A.reconstruct_route(0, 3, p, e)
        self.assertEqual(sum(e[p[v]][2] for v in route[1:]), r["distance"][3])

    def test_route_helper_tied_and_parallel_edges(self):
        # Two parallel 0->1 edges (ids 0,1) tie at cost 4; path 0->2->1 also
        # costs 1+3 = 4 (id 2, 3).  Any of the three shortest routes is valid.
        e = [[0, 1, 4], [0, 1, 4], [0, 2, 1], [2, 1, 3], [1, 3, 0]]
        r = self.run_check(4, 0, e)
        self.assertEqual(r["distance"], [0, 4, 1, 4])
        self.assertIn(r["parent_edge"][1], (0, 1, 3))
        route = A.reconstruct_route(0, 3, r["parent_edge"], e)
        self.assertIn(route, ([0, 1, 3], [0, 2, 1, 3]))
        # direction + weight verification along the chosen edges
        p = r["parent_edge"]
        total = 0
        for a, b in zip(route, route[1:]):
            eid = p[b]
            self.assertEqual((e[eid][0], e[eid][1]), (a, b))
            total += e[eid][2]
        self.assertEqual(total, r["distance"][3])

    def test_helper_rejects_inconsistent_parent_edges(self):
        e = [[0, 1, 1], [1, 2, 1]]
        self.assertEqual(A.reconstruct_route(0, 2, [-1, 1, 1], e), [])  # wrong end
        self.assertEqual(A.reconstruct_route(0, 2, [], e), [])          # cyclic result

    def test_single_vertex_and_isolates(self):
        r = self.run_check(1, 0, [])
        self.assertEqual(r, {"acyclic": True, "distance": [0], "parent_edge": [-1]})
        r = self.run_check(3, 1, [[0, 2, 5]])  # edge unreachable from source 1
        self.assertEqual(r["distance"], [None, 0, None])

    def test_cycle_unreachable_from_source_still_rejects(self):
        e = [[0, 1, 1], [2, 3, 1], [3, 2, 1]]  # cycle 2<->3, source is 0
        r = self.run_check(4, 0, e)
        self.assertEqual(r, {"acyclic": False, "distance": [], "parent_edge": []})

    def test_reachable_cycle_and_self_loop(self):
        for n, e in [(3, [[0, 1, 1], [1, 2, 1], [2, 1, -5]]),
                     (2, [[0, 1, 1], [1, 1, 0]]),
                     (1, [[0, 0, -1]])]:
            r = self.run_check(n, 0, e)
            self.assertFalse(r["acyclic"])

    def test_minus_one_is_a_real_distance_not_unreachable(self):
        e = [[0, 1, -1], [2, 1, 5]]  # vertex 2 unreachable, vertex 1 at -1
        r = self.run_check(3, 0, e)
        self.assertEqual(r["distance"], [0, -1, None])
        self.assertEqual(r["parent_edge"], [-1, 0, -1])

    def test_never_relaxes_from_unreachable_vertex(self):
        # 2 is unreachable and has a hugely negative edge to 1; it must not win.
        e = [[0, 1, 5], [2, 1, -BIG]]
        r = self.run_check(3, 0, e)
        self.assertEqual(r["distance"], [0, 5, None])
        self.assertEqual(r["parent_edge"][1], 0)

    def test_random_small_against_bellman_ford(self):
        rng = random.Random(9)
        for _ in range(400):
            n = rng.randint(1, 8)
            e = random_dag_edges(rng, n, rng.randint(0, 12))
            self.run_check(n, rng.randrange(n), e)
        for _ in range(100):  # possibly cyclic graphs: whole-graph rule
            n = rng.randint(1, 6)
            e = random_edges(rng, n, rng.randint(0, 8), weighted=True)
            self.run_check(n, rng.randrange(n), e)

    def test_deep_chain_totals_beyond_32_bits(self):
        n = MAXN
        e = [[i, i + 1, BIG] for i in range(n - 1)]
        r = A.dag(req("dag", n, source=0, edges=e))
        self.assertEqual(r["distance"][-1], (n - 1) * BIG)
        check_dag(n, 0, e, r)
        route = A.reconstruct_route(0, n - 1, r["parent_edge"], e)
        self.assertEqual(route, list(range(n)))
        # negative extreme as well
        e = [[i, i + 1, -BIG] for i in range(n - 1)]
        r = A.dag(req("dag", n, source=0, edges=e))
        self.assertEqual(r["distance"][-1], -(n - 1) * BIG)

    def test_deep_chain_with_a_far_cycle_unreachable(self):
        n = MAXN
        e = [[i, i + 1, 1] for i in range(n // 2)]
        e += [[n - 2, n - 1, 1], [n - 1, n - 2, 1]]  # cycle, not reachable
        r = A.dag(req("dag", n, source=0, edges=e))
        self.assertEqual(r, {"acyclic": False, "distance": [], "parent_edge": []})

    def test_full_limit_random_dag(self):
        rng = random.Random(10)
        n, m = MAXN, MAXM
        e = random_dag_edges(rng, n, m, lo=-BIG, hi=BIG)
        source = e[0][0]
        r = A.dag(req("dag", n, source=source, edges=e))
        check_dag(n, source, e, r)

    def test_checker_rejects_invalid_witnesses(self):
        e = [[0, 1, 3], [0, 2, 8], [1, 2, -7], [2, 3, 4]]
        good = {"acyclic": True, "distance": [0, 3, -4, 0, None],
                "parent_edge": [-1, 0, 2, 3, -1]}
        check_dag(5, 0, e, good)
        alt_ok = dict(good)  # same answer: tied alternatives accepted elsewhere
        check_dag(5, 0, e, alt_ok)
        mutations = {
            "parent not tight (edge 1 gives 8, not -4)":
                dict(good, parent_edge=[-1, 0, 1, 3, -1]),
            "parent edge ends at wrong vertex":
                dict(good, parent_edge=[-1, 0, 2, 2, -1]),
            "unreachable vertex given a parent":
                dict(good, parent_edge=[-1, 0, 2, 3, 3]),
            "wrong distance":
                dict(good, distance=[0, 3, -3, 1, None]),
            "unreachable reported as finite":
                dict(good, distance=[0, 3, -4, 0, -1]),
            "-1 used for unreachable instead of null":
                dict(good, distance=[0, 3, -4, 0, -1]),
            "cyclic verdict for acyclic graph":
                {"acyclic": False, "distance": [], "parent_edge": []},
        }
        for name, bad in mutations.items():
            with self.assertRaises(AssertionError, msg=name):
                check_dag(5, 0, e, bad)


def kahn_acyclic(n, e):
    from checkers import kahn_is_acyclic
    return kahn_is_acyclic(n, e)


# ------------------------------------------- process-level state reset ----

class TestProcess(unittest.TestCase):
    """Run the supplied JSON Lines entry point; repeat mixed requests."""

    requests = [
        {"task": "ping"},
        req("dsu", 4, operations=[[0, 1], [1, 2], [0, 2], [3, 3]]),
        req("kruskal", 4, edges=[[0, 1, 4], [1, 2, 2], [0, 2, 3]]),
        req("prim", 4, edges=[[0, 1, 4], [1, 2, 2], [0, 2, 3]]),
        req("scc", 4, edges=[[0, 1], [1, 0], [1, 2], [0, 2], [2, 3]]),
        req("topo", 4, edges=[[0, 2], [1, 2], [2, 3]]),
        req("topo", 2, edges=[[0, 1], [1, 0]]),              # cyclic, then...
        req("topo", 3, edges=[[0, 1], [1, 2]]),              # ...acyclic again
        req("dag", 5, source=0, edges=[[0, 1, 3], [0, 2, 8], [1, 2, -7], [2, 3, 4]]),
        req("dag", 4, source=0, edges=[[0, 1, 1], [2, 3, 1], [3, 2, 1]]),
        req("dag", 1, source=0, edges=[]),
        req("dsu", 0, operations=[]),
        req("kruskal", 0, edges=[]),
        req("scc", 0, edges=[]),
        req("dag", 3, source=2, edges=[[0, 1, 1]]),
    ]

    def run_main(self, reqs):
        data = "".join(json.dumps(r) + "\n" for r in reqs)
        p = subprocess.run([sys.executable, os.path.join(SRC, "main.py")],
                           input=data, capture_output=True, text=True,
                           encoding="utf-8", timeout=60, cwd=SRC)
        self.assertEqual(p.returncode, 0, p.stderr)
        lines = p.stdout.splitlines()
        self.assertEqual(len(lines), len(reqs), "one result line per request")
        return [json.loads(x) for x in lines]  # every line must be JSON

    def test_batch_repeated_three_times_is_identical(self):
        once = self.run_main(self.requests)
        thrice = self.run_main(self.requests * 3)
        for k in range(3):
            chunk = thrice[k * len(once):(k + 1) * len(once)]
            self.assertEqual(once, chunk)  # same input -> same output, no leaked state
        self.assertEqual(once[0], {"status": "ready"})

    def test_results_valid_in_process(self):
        out = self.run_main(self.requests)
        self.assertEqual(out[1], {"states": [[3, 2], [2, 3], [2, 3], [2, 3]]})
        check_forest(4, self.requests[2]["edges"], out[2], 5)
        check_forest(4, self.requests[3]["edges"], out[3], 5)
        check_scc(4, self.requests[4]["edges"], out[4])
        check_topo(4, self.requests[5]["edges"], out[5])
        self.assertEqual(out[6], {"acyclic": False, "order": []})
        check_topo(3, self.requests[7]["edges"], out[7])
        check_dag(5, 0, self.requests[8]["edges"], out[8])
        self.assertEqual(out[9], {"acyclic": False, "distance": [], "parent_edge": []})
        self.assertEqual(out[10], {"acyclic": True, "distance": [0], "parent_edge": [-1]})
        self.assertEqual(out[14]["distance"], [None, None, 0])


if __name__ == "__main__":
    unittest.main()

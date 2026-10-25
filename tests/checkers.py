"""Independent semantic checkers.

They never call algorithms.py.  Non-unique answers (forests, SCC labels,
topological orders, shortest-path witnesses) are accepted when valid and
rejected (AssertionError) when not.
"""

from itertools import combinations


def _fail(message):
    raise AssertionError(message)


class _UF:
    """Tiny union-find used only for checking (no compression tricks needed)."""

    def __init__(self, n):
        self.p = list(range(n))
        self.count = n

    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a == b:
            return False
        self.p[a] = b
        self.count -= 1
        return True


# ------------------------------------------------------------------ MST ----

def count_components(n, edges):
    uf = _UF(n)
    for e in edges:
        uf.union(e[0], e[1])
    return uf.count


def brute_force_msf_cost(n, edges):
    """Exhaustive minimum spanning forest cost (small graphs only)."""
    c = count_components(n, edges)
    need = n - c
    best = None
    for subset in combinations(range(len(edges)), need):
        uf = _UF(n)
        if all(uf.union(edges[i][0], edges[i][1]) for i in subset):
            total = sum(edges[i][2] for i in subset)
            if best is None or total < best:
                best = total
    return best


def check_forest(n, edges, result, expected_cost=None):
    """Valid spanning forest: unique ids, no loops, acyclic, spans every
    component, count n-c, reported cost == sum of weights (== expected)."""
    ids = result["edge_ids"]
    c = count_components(n, edges)
    if result["components"] != c:
        _fail(f"components {result['components']} != {c}")
    if len(ids) != n - c:
        _fail(f"{len(ids)} edges selected, forest needs {n - c}")
    if len(set(ids)) != len(ids):
        _fail("repeated edge id")
    uf = _UF(n)
    total = 0
    for i in ids:
        if not (0 <= i < len(edges)):
            _fail(f"edge id {i} out of range")
        u, v, w = edges[i]
        if u == v:
            _fail(f"self-loop {i} selected")
        if not uf.union(u, v):
            _fail(f"edge {i} closes a cycle")
        total += w
    if uf.count != c:
        _fail("selected edges do not span every component")
    if result["cost"] != total:
        _fail(f"reported cost {result['cost']} != sum of selected weights {total}")
    if expected_cost is not None and total != expected_cost:
        _fail(f"cost {total} is not the optimum {expected_cost}")


# ------------------------------------------------------------------ SCC ----

def reachability(n, edges):
    """Transitive closure as bit masks (Warshall); small n only."""
    reach = [1 << i for i in range(n)]
    for u, v in edges:
        reach[u] |= 1 << v
    for k in range(n):
        bit = 1 << k
        for i in range(n):
            if reach[i] & bit:
                reach[i] |= reach[k]
    return reach


def check_scc_condensation(n, edges, result, same_scc):
    """Labels are consistent with `same_scc(i, j)` and the condensation is
    exactly the set of distinct cross-component pairs, without duplicates."""
    comp = result["component"]
    if len(comp) != n:
        _fail("component array has wrong length")
    if any((not isinstance(x, int)) or x < 0 for x in comp):
        _fail("labels must be nonnegative integers")
    # Partition equality via canonical first-seen relabelling.
    for i in range(n):
        for j in range(i + 1, n):
            if (comp[i] == comp[j]) != same_scc(i, j):
                _fail(f"vertices {i},{j}: labels {comp[i]},{comp[j]} wrong")
    pairs = [tuple(p) for p in result["condensation"]]
    if len(set(pairs)) != len(pairs):
        _fail("duplicate condensation pair")
    expected = {(comp[u], comp[v]) for u, v in edges if comp[u] != comp[v]}
    if set(pairs) != expected:
        _fail(f"condensation {sorted(pairs)} != {sorted(expected)}")


def check_scc(n, edges, result):
    reach = reachability(n, edges)
    check_scc_condensation(
        n, edges, result,
        lambda i, j: bool(reach[i] >> j & 1) and bool(reach[j] >> i & 1))


# ----------------------------------------------------------------- topo ----

def kahn_is_acyclic(n, edges):
    indeg = [0] * n
    adj = [[] for _ in range(n)]
    for e in edges:
        adj[e[0]].append(e[1])
        indeg[e[1]] += 1
    queue = [v for v in range(n) if indeg[v] == 0]
    seen = 0
    while queue:
        u = queue.pop()
        seen += 1
        for v in adj[u]:
            indeg[v] -= 1
            if indeg[v] == 0:
                queue.append(v)
    return seen == n


def check_topo(n, edges, result):
    acyclic = kahn_is_acyclic(n, edges)
    if result["acyclic"] is not acyclic:
        _fail(f"acyclic should be {acyclic}")
    if not acyclic:
        if result["order"] != []:
            _fail("cyclic graph must return an empty order")
        return
    order = result["order"]
    if sorted(order) != list(range(n)):
        _fail("order is not a permutation of all vertices")
    pos = {v: i for i, v in enumerate(order)}
    for e in edges:
        if pos[e[0]] >= pos[e[1]]:
            _fail(f"edge {e[0]}->{e[1]} goes backwards")


# ------------------------------------------------------------------ DAG ----

def reachable_set(n, edges, source):
    adj = [[] for _ in range(n)]
    for e in edges:
        adj[e[0]].append(e[1])
    seen = {source}
    stack = [source]
    while stack:
        u = stack.pop()
        for v in adj[u]:
            if v not in seen:
                seen.add(v)
                stack.append(v)
    return seen


def bellman_ford(n, edges, source):
    """Reference distances (None = unreachable); valid because the graph is a DAG."""
    dist = [None] * n
    dist[source] = 0
    for _ in range(n):
        changed = False
        for u, v, w in edges:
            if dist[u] is not None and (dist[v] is None or dist[u] + w < dist[v]):
                dist[v] = dist[u] + w
                changed = True
        if not changed:
            break
    return dist


def check_dag(n, source, edges, result, reference=None):
    """Whole-graph cycle rule, exact distances and valid, tight parent paths.

    Optimality without an oracle: every distance is feasible
    (d[v] <= d[u]+w on every edge leaving a reachable u) and every reachable
    non-source vertex has a tight parent edge ending at it that leads back to
    the source, so d[v] is the length of a real path and no path is shorter.
    """
    if not kahn_is_acyclic(n, edges):
        if result != {"acyclic": False, "distance": [], "parent_edge": []}:
            _fail("cyclic graph must give acyclic=false and two empty arrays")
        return
    if result["acyclic"] is not True:
        _fail("acyclic should be true")
    dist, parent = result["distance"], result["parent_edge"]
    if len(dist) != n or len(parent) != n:
        _fail("arrays must have length n")
    reach = reachable_set(n, edges, source)
    if dist[source] != 0 or parent[source] != -1:
        _fail("source must have distance 0 and parent_edge -1")
    for v in range(n):
        if (v in reach) != (dist[v] is not None):
            _fail(f"vertex {v}: reachability mismatch (distance {dist[v]})")
        if v not in reach and parent[v] != -1:
            _fail(f"unreachable vertex {v} must have parent_edge -1")
    for u, v, w in edges:  # feasibility
        if dist[u] is not None and dist[v] > dist[u] + w:
            _fail(f"edge {u}->{v} can still relax")
    verified = {source}  # vertices whose parent chain is known to reach source
    for v in range(n):  # tight parents that reach the source (amortised O(n))
        if dist[v] is None or v in verified:
            continue
        path, x = [], v
        while x not in verified:
            if len(path) > n:  # more than n steps without reaching a verified vertex
                _fail("parent edges loop")
            eid = parent[x]
            if not (0 <= eid < len(edges)):
                _fail(f"vertex {x}: bad parent edge {eid}")
            a, b, w = edges[eid]
            if b != x:
                _fail(f"parent edge {eid} does not end at vertex {x}")
            if dist[a] is None or dist[x] != dist[a] + w:
                _fail(f"parent edge {eid} is not tight")
            path.append(x)
            x = a
        verified.update(path)
    if reference is not None and dist != reference:
        _fail(f"distances {dist} != reference {reference}")

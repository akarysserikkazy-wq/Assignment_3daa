"""Graph algorithms for Assignment 3: DSU, Kruskal, Prim, Kosaraju SCC,
DFS topological ordering and DAG shortest paths with route reconstruction.

Every function receives one request dict and returns one result dict.  No
module-level mutable state is used, so each request starts from a clean slate.
All DFS traversals are iterative (explicit stacks with per-vertex neighbour
pointers), so graphs with 100 000 vertices cannot overflow the call stack.
"""

import heapq


# --------------------------------------------------------------------------
# Disjoint-set union
# --------------------------------------------------------------------------

class _DisjointSet:
    """Union by size + path compression (iterative, two-pass)."""

    __slots__ = ("parent", "size", "components", "largest")

    def __init__(self, n):
        self.parent = list(range(n))
        self.size = [1] * n
        self.components = n
        self.largest = 1 if n > 0 else 0

    def find(self, x):
        parent = self.parent
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:  # compress the whole path
            nxt = parent[x]
            parent[x] = root
            x = nxt
        return root

    def union(self, a, b):
        """Merge the sets of a and b; return True if they were different."""
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return False
        size = self.size
        if size[ra] < size[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        size[ra] += size[rb]
        self.components -= 1
        if size[ra] > self.largest:
            self.largest = size[ra]
        return True


def dsu(request):
    """DSU with path compression and union by size/rank; return states."""
    n = request["n"]
    operations = request["operations"]

    ds = _DisjointSet(n)
    states = []
    for u, v in operations:
        ds.union(u, v)
        states.append([ds.components, ds.largest])
    return {"states": states}


# --------------------------------------------------------------------------
# Minimum spanning forest
# --------------------------------------------------------------------------

def kruskal(request):
    """Separate Kruskal implementation; return cost, edge_ids, components."""
    n = request["n"]
    edges = request["edges"]

    # Sort edge *ids* so original positions survive the sort (stable sort).
    order = sorted(range(len(edges)), key=lambda i: edges[i][2])
    ds = _DisjointSet(n)
    cost = 0
    edge_ids = []
    for eid in order:
        u, v, w = edges[eid]
        if ds.union(u, v):  # False for self-loops and cycle-closing edges
            cost += w
            edge_ids.append(eid)
            if ds.components == 1:
                break
    return {"cost": cost, "edge_ids": edge_ids, "components": ds.components}


def prim(request):
    """Separate Prim implementation; restart for every unvisited component."""
    n = request["n"]
    edges = request["edges"]

    adjacency = [[] for _ in range(n)]
    for eid, (u, v, w) in enumerate(edges):
        if u != v:  # a self-loop can never be in a forest
            adjacency[u].append((w, eid, v))
            adjacency[v].append((w, eid, u))

    visited = [False] * n
    cost = 0
    edge_ids = []
    components = 0
    heappush = heapq.heappush
    heappop = heapq.heappop

    for start in range(n):
        if visited[start]:
            continue
        components += 1  # new component: restart from an unvisited vertex
        visited[start] = True
        heap = []
        for item in adjacency[start]:
            heappush(heap, item)
        while heap:
            w, eid, v = heappop(heap)
            if visited[v]:  # lazy deletion: stale entry, endpoint already in tree
                continue
            visited[v] = True
            cost += w
            edge_ids.append(eid)
            for item in adjacency[v]:
                if not visited[item[2]]:
                    heappush(heap, item)
    return {"cost": cost, "edge_ids": edge_ids, "components": components}


# --------------------------------------------------------------------------
# Iterative DFS building blocks
# --------------------------------------------------------------------------

def _finish_order(n, adjacency):
    """Return vertices in increasing DFS finish time (iterative DFS)."""
    visited = [False] * n
    pointer = [0] * n
    finish = []
    for start in range(n):
        if visited[start]:
            continue
        visited[start] = True
        stack = [start]
        while stack:
            u = stack[-1]
            neighbours = adjacency[u]
            i = pointer[u]
            if i < len(neighbours):
                pointer[u] = i + 1  # resume here later: no rescanning
                v = neighbours[i]
                if not visited[v]:
                    visited[v] = True
                    stack.append(v)
            else:
                stack.pop()
                finish.append(u)
    return finish


def _topological_order(n, adjacency):
    """Three-state DFS (0 unvisited, 1 active, 2 finished) over all components.

    Returns a topological order, or None if any directed cycle exists.
    """
    state = [0] * n
    pointer = [0] * n
    finish = []
    for start in range(n):
        if state[start] != 0:
            continue
        state[start] = 1
        stack = [start]
        while stack:
            u = stack[-1]
            neighbours = adjacency[u]
            i = pointer[u]
            if i < len(neighbours):
                pointer[u] = i + 1
                v = neighbours[i]
                s = state[v]
                if s == 0:
                    state[v] = 1
                    stack.append(v)
                elif s == 1:  # edge to an active vertex = back edge = cycle
                    return None
                # s == 2: edge to a finished vertex is harmless
            else:
                state[u] = 2
                stack.pop()
                finish.append(u)
    finish.reverse()
    return finish


# --------------------------------------------------------------------------
# Kosaraju SCC
# --------------------------------------------------------------------------

def scc(request):
    """Kosaraju with two DFS passes; return component and condensation."""
    n = request["n"]
    edges = request["edges"]

    adjacency = [[] for _ in range(n)]
    transpose = [[] for _ in range(n)]
    for u, v in edges:
        adjacency[u].append(v)
        transpose[v].append(u)

    # Pass 1 on G: finish order.  Pass 2 on G^T in decreasing finish order.
    finish = _finish_order(n, adjacency)
    component = [-1] * n
    pointer = [0] * n
    count = 0
    for start in reversed(finish):
        if component[start] != -1:
            continue
        component[start] = count
        stack = [start]
        while stack:
            u = stack[-1]
            neighbours = transpose[u]
            i = pointer[u]
            if i < len(neighbours):
                pointer[u] = i + 1
                v = neighbours[i]
                if component[v] == -1:
                    component[v] = count
                    stack.append(v)
            else:
                stack.pop()
        count += 1

    # Condensation: each distinct cross-component pair exactly once.
    seen = set()
    for u, v in edges:
        a = component[u]
        b = component[v]
        if a != b:
            seen.add(a * count + b)
    condensation = [[key // count, key % count] for key in seen]
    return {"component": component, "condensation": condensation}


# --------------------------------------------------------------------------
# Topological order and DAG shortest paths
# --------------------------------------------------------------------------

def topo(request):
    """Three-state DFS and reverse finish order; return acyclic and order."""
    n = request["n"]
    edges = request["edges"]

    adjacency = [[] for _ in range(n)]
    for u, v in edges:
        adjacency[u].append(v)
    order = _topological_order(n, adjacency)
    if order is None:
        return {"acyclic": False, "order": []}
    return {"acyclic": True, "order": order}


def dag(request):
    """Whole-graph cycle check and DAG relaxation; add your own route helper."""
    n = request["n"]
    source = request["source"]
    edges = request["edges"]

    targets = [[] for _ in range(n)]   # targets[u]  -> successor vertices
    out_ids = [[] for _ in range(n)]   # out_ids[u]  -> original edge ids
    for eid, (u, v, _w) in enumerate(edges):
        targets[u].append(v)
        out_ids[u].append(eid)

    order = _topological_order(n, targets)  # checks the WHOLE graph
    if order is None:
        return {"acyclic": False, "distance": [], "parent_edge": []}

    distance = [None] * n
    parent_edge = [-1] * n
    distance[source] = 0
    for u in order:
        du = distance[u]
        if du is None:  # never relax from an unreachable vertex
            continue
        for eid in out_ids[u]:
            _u, v, w = edges[eid]
            candidate = du + w
            dv = distance[v]
            if dv is None or candidate < dv:
                distance[v] = candidate
                parent_edge[v] = eid
    return {"acyclic": True, "distance": distance, "parent_edge": parent_edge}


def reconstruct_route(source, target, parent_edge, edges):
    """Rebuild the vertex sequence source -> ... -> target from parent edge IDs.

    * target == source        -> [source]
    * reachable target        -> [source, ..., target] following parent_edge
    * unreachable target      -> []   (documented no-route result)
    * empty parent_edge (cyclic graph result) -> []

    Runs in O(L) for a route of L vertices.  Each parent edge must END at the
    vertex it is stored for; a violation (or a loop) yields [] instead of a
    wrong route.
    """
    if not parent_edge:
        return []
    if target == source:
        return [source]
    route = [target]
    v = target
    limit = len(parent_edge)
    while v != source:
        eid = parent_edge[v]
        if eid < 0 or edges[eid][1] != v or len(route) > limit:
            return []
        v = edges[eid][0]
        route.append(v)
    route.reverse()
    return route

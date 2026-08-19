// ============================================================================
//  DIJKSTRA  --  single source shortest paths, NON-NEGATIVE edge weights
// ----------------------------------------------------------------------------
//  O(E log V) with a binary heap. THE default shortest path algorithm.
//
//  HARD RULE: any negative edge and Dijkstra is WRONG (not slow -- wrong).
//  Use bellmanford.cpp instead. A "cost you can save" is a negative edge in
//  disguise; so is "profit". If all edge weights are 1, use plain BFS (O(V+E)).
//
//  API  (nodes 0-based)
//    Dij g; g.init(N);
//    g.ae(u, v, w);            // directed edge u -> v, weight w >= 0
//    g.ae(u, v, w, true);      // undirected edge
//    g.run(s);                 // fills g.dist[] and g.par[]
//    g.dist[v]                 // INF if unreachable
//    g.path(v)                 // vector of nodes s..v, empty if unreachable
//
//  CONTENTS
//    1.  Dij            heap Dijkstra + path reconstruction
//    2.  denseDijkstra  O(V^2), better when E ~ V^2 (e.g. V <= 2000 complete)
//    3.  zeroOneBFS     weights in {0,1} only, O(V+E) with a deque
//    4.  countShortest  number of shortest paths (mod p) + shortest-path DAG
//    5.  secondShortest second shortest WALK from s to t
//    6.  minimaxPath    minimise the LARGEST edge on the path (bottleneck)
//    7.  MODELLING NOTES -- layered / state graphs, the real exam content
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

const ll INF = (ll)4e18;

// ============================================================================
//  1. HEAP DIJKSTRA
// ----------------------------------------------------------------------------
//  Lazy deletion: we never decrease-key, we just push duplicates and skip a
//  node the second time we pop it. This is why the `if (d > dist[v]) continue;`
//  line is mandatory -- without it the complexity blows up.
//
//  Pick INF larger than any real distance but small enough that INF + w does
//  not overflow. With ll and 4e18 you are safe as long as you never relax out
//  of an unvisited node (we don't -- we only pop settled ones).
// ============================================================================
struct Dij {
	int N;
	vector<vector<pair<int, ll>>> adj;   // adj[u] = {(v, w), ...}
	vector<ll> dist;
	vector<int> par;

	void init(int _N) {
		N = _N;
		adj.assign(N, {});
		dist.assign(N, INF);
		par.assign(N, -1);
	}
	void ae(int u, int v, ll w, bool undirected = false) {
		adj[u].push_back({v, w});
		if (undirected) adj[v].push_back({u, w});
	}

	void run(int s) {
		dist.assign(N, INF);
		par.assign(N, -1);
		// min-heap of (distance, node)
		priority_queue<pair<ll, int>, vector<pair<ll, int>>,
		               greater<pair<ll, int>>> pq;
		dist[s] = 0;
		pq.push({0, s});
		while (!pq.empty()) {
			auto [d, u] = pq.top();
			pq.pop();
			if (d > dist[u]) continue;         // stale entry -- already settled
			for (auto [v, w] : adj[u]) {
				if (d + w < dist[v]) {          // relax
					dist[v] = d + w;
					par[v] = u;
					pq.push({dist[v], v});
				}
			}
		}
	}

	// MULTI-SOURCE: shortest distance to the NEAREST of several sources.
	// Same cost as one Dijkstra -- do NOT run it once per source.
	void runMulti(const vector<int>& sources) {
		dist.assign(N, INF);
		par.assign(N, -1);
		priority_queue<pair<ll, int>, vector<pair<ll, int>>,
		               greater<pair<ll, int>>> pq;
		for (int s : sources) { dist[s] = 0; pq.push({0, s}); }
		while (!pq.empty()) {
			auto [d, u] = pq.top();
			pq.pop();
			if (d > dist[u]) continue;
			for (auto [v, w] : adj[u])
				if (d + w < dist[v]) { dist[v] = d + w; par[v] = u; pq.push({dist[v], v}); }
		}
	}

	// Nodes from the source to v inclusive. Empty if v is unreachable.
	vector<int> path(int v) {
		if (dist[v] == INF) return {};
		vector<int> p;
		for (int x = v; x != -1; x = par[x]) p.push_back(x);
		reverse(p.begin(), p.end());
		return p;
	}
};

// ============================================================================
//  2. DENSE DIJKSTRA -- O(V^2), no heap
// ----------------------------------------------------------------------------
//  Beats the heap version when E is close to V^2 (complete graphs, geometry
//  problems where every pair of points is an edge). V <= ~5000 is comfortable.
//  Takes an adjacency MATRIX; w[u][v] = INF for "no edge".
// ============================================================================
vector<ll> denseDijkstra(const vector<vector<ll>>& w, int s) {
	int n = (int)w.size();
	vector<ll> dist(n, INF);
	vector<char> done(n, 0);
	dist[s] = 0;
	for (int it = 0; it < n; it++) {
		int u = -1;
		for (int i = 0; i < n; i++)
			if (!done[i] && (u == -1 || dist[i] < dist[u])) u = i;
		if (u == -1 || dist[u] == INF) break;
		done[u] = 1;
		for (int v = 0; v < n; v++)
			if (w[u][v] != INF && dist[u] + w[u][v] < dist[v])
				dist[v] = dist[u] + w[u][v];
	}
	return dist;
}

// ============================================================================
//  3. 0-1 BFS -- all weights are 0 or 1, O(V + E)
// ----------------------------------------------------------------------------
//  Dijkstra with a deque instead of a heap: a 0-edge goes on the FRONT (same
//  distance layer), a 1-edge on the BACK. The deque stays sorted automatically.
//
//  Shows up constantly as "moving in the current direction is free, turning
//  costs 1", "you may delete up to k walls", "some roads are free".
//  Generalises to weights in {0..k} with k+1 buckets (Dial's algorithm).
// ============================================================================
vector<ll> zeroOneBFS(const vector<vector<pair<int, int>>>& adj, int s) {
	int n = (int)adj.size();
	vector<ll> dist(n, INF);
	deque<int> dq;
	dist[s] = 0;
	dq.push_back(s);
	while (!dq.empty()) {
		int u = dq.front();
		dq.pop_front();
		for (auto [v, w] : adj[u]) {          // w must be 0 or 1
			if (dist[u] + w < dist[v]) {
				dist[v] = dist[u] + w;
				if (w == 0) dq.push_front(v);
				else        dq.push_back(v);
			}
		}
	}
	return dist;
}

// ============================================================================
//  4. COUNTING SHORTEST PATHS  +  THE SHORTEST-PATH DAG
// ----------------------------------------------------------------------------
//  ways[v] = number of distinct shortest s->v paths, mod MOD.
//  An edge (u,v,w) lies on SOME shortest path iff dist[u] + w == dist[v]; the
//  set of such edges forms a DAG. Any "best shortest path under a tiebreak"
//  question (fewest edges, max total scenery, lexicographically smallest) is a
//  DP over that DAG -- process nodes in increasing dist order.
// ============================================================================
const ll MOD = 1000000007;

vector<ll> countShortestPaths(Dij& g, int s) {
	g.run(s);
	int n = g.N;
	vector<int> ord(n);
	iota(ord.begin(), ord.end(), 0);
	// increasing distance order is a topological order of the shortest-path DAG
	sort(ord.begin(), ord.end(),
	     [&](int a, int b) { return g.dist[a] < g.dist[b]; });
	vector<ll> ways(n, 0);
	ways[s] = 1;
	for (int u : ord) {
		if (g.dist[u] == INF) continue;
		for (auto [v, w] : g.adj[u])
			if (g.dist[u] + w == g.dist[v]) ways[v] = (ways[v] + ways[u]) % MOD;
	}
	return ways;
}

// ============================================================================
//  5. SECOND SHORTEST WALK  s -> t
// ----------------------------------------------------------------------------
//  "Walk", not "path": it may revisit vertices, and by the usual definition it
//  EQUALS the shortest distance when two distinct shortest paths exist.
//
//  Keep the two best distances per node and run one Dijkstra over them. A node
//  is settled only after it has been popped twice, so push a node whenever
//  either of its two values improves.
//
//  The `else if (nd < d2[v])` branch deliberately admits nd == d1[v]: that is
//  what makes a second distinct shortest path count. If instead you want the
//  second shortest STRICTLY LONGER walk, change it to
//        else if (nd > d1[v] && nd < d2[v])
//
//  (There is also an edge-based formulation -- Dijkstra from s and from t on
//  the reversed graph, then minimise d1[u] + w + d2[v] over edges. It is easy
//  to get wrong: EVERY edge of the shortest path attains d1[t], so you cannot
//  infer "two shortest paths" just from seeing that value twice. Prefer this.)
// ============================================================================
ll secondShortest(int n, const vector<array<ll, 3>>& edges, int s, int t) {
	vector<vector<pair<int, ll>>> adj(n);
	for (auto& e : edges) adj[(int)e[0]].push_back({(int)e[1], e[2]});

	vector<ll> d1(n, INF), d2(n, INF);
	priority_queue<pair<ll, int>, vector<pair<ll, int>>,
	               greater<pair<ll, int>>> pq;
	d1[s] = 0;
	pq.push({0, s});
	while (!pq.empty()) {
		auto [d, u] = pq.top();
		pq.pop();
		if (d > d2[u]) continue;                 // worse than both known values
		for (auto [v, w] : adj[u]) {
			ll nd = d + w;
			if (nd < d1[v]) {                     // new best; old best demoted
				d2[v] = d1[v];
				d1[v] = nd;
				pq.push({nd, v});
				if (d2[v] != INF) pq.push({d2[v], v});
			} else if (nd < d2[v]) {              // new second best (ties allowed)
				d2[v] = nd;
				pq.push({nd, v});
			}
		}
	}
	return d2[t];
}

// ============================================================================
//  6. MINIMAX (BOTTLENECK) PATH -- minimise the LARGEST edge used
// ----------------------------------------------------------------------------
//  "What is the smallest weight limit that still lets a truck get from s to t?"
//  Same as Dijkstra but the path cost is max(edges) rather than sum(edges);
//  max is monotone so the greedy argument still holds.
//  Equivalent answer: the largest edge on the s-t path of the MINIMUM SPANNING
//  TREE (see mst.cpp) -- useful when you need it for all pairs at once.
// ============================================================================
vector<ll> minimaxPath(const vector<vector<pair<int, ll>>>& adj, int s) {
	int n = (int)adj.size();
	vector<ll> best(n, INF);
	priority_queue<pair<ll, int>, vector<pair<ll, int>>,
	               greater<pair<ll, int>>> pq;
	best[s] = 0;
	pq.push({0, s});
	while (!pq.empty()) {
		auto [d, u] = pq.top();
		pq.pop();
		if (d > best[u]) continue;
		for (auto [v, w] : adj[u]) {
			ll cand = max(d, w);                 // <-- the only change
			if (cand < best[v]) { best[v] = cand; pq.push({cand, v}); }
		}
	}
	return best;
}

// ============================================================================
//  7. MODELLING NOTES -- where the marks actually are
// ----------------------------------------------------------------------------
//  Nearly every hard shortest-path problem is easy Dijkstra on a graph you had
//  to invent. Build the right graph and the code is 20 lines.
//
//  LAYERED / STATE GRAPHS.  Make the node a (position, state) pair and index it
//  as  id = pos * S + state.  Standard states:
//     * "you may use up to k free edges"        -> state = free edges used
//       edges: (u,j) -> (v,j)   cost w        (pay normally)
//               (u,j) -> (v,j+1) cost 0        (use a freebie)
//     * "you may take at most k flights/refuels"-> state = count used
//     * "the cost depends on how you arrived"   -> state = incoming direction
//       (turning penalties, trains that must not reverse)
//     * "alternate between two edge types"      -> state = parity
//     * "fuel tank of capacity C"               -> state = current fuel
//     Node count multiplies by S; keep V*S*log within budget.
//
//  SUPER-NODES for "everyone in group X is mutually connected at cost c".
//  Naively that is |X|^2 edges. Instead add ONE node g per group:
//       every member -> g_out (cost 0),  g_in -> every member (cost 0),
//       g_out -> g_in (cost c).
//  This is exactly the lecture's Intercountry problem: two nodes per country,
//  one for departures and one for arrivals, so a country is "entered" once.
//  O(N + C) nodes instead of O(N^2) edges. If you use a single node per group
//  instead of an in/out pair, you accidentally allow free teleport between two
//  members -- check whether that is legal in your problem.
//
//  IMPLICIT GRAPHS. Do not build the adjacency list; generate neighbours when
//  you pop a node. Numbers-on-a-display, puzzle states, grids. Bound the state
//  space first with an argument (lecture's Two Buttons: never exceed m+1), then
//  use a flat array over that bound -- a hash map is usually too slow.
//
//  BINARY SEARCH + REACHABILITY. "Maximise the minimum clearance", "minimise
//  the maximum jump": binary search the answer X, then the check is a plain
//  BFS/DFS on the graph restricted to edges/cells allowed at X. Monotone
//  because a feasible X implies every smaller X is feasible. See gridbfs.cpp.
//  (Often a minimax Dijkstra as in (6) removes the binary search entirely.)
//
//  UNIT WEIGHTS -> BFS. Do not pay log for nothing. Weights {0,1} -> 0-1 BFS.
//  Small integer weights {0..k} -> Dial's, or split an edge of weight w into w
//  unit edges when sum(w) is small.
//
//  DAG -> no Dijkstra needed: relax edges in topological order, O(V+E), and it
//  works with negative weights too (see toposort.cpp).
// ============================================================================

int main() {
	// 0 --4-- 1
	//  \      |
	//   1     2
	//    \    |
	//     2 --3-- 3
	Dij g;
	g.init(4);
	g.ae(0, 1, 4, true);
	g.ae(0, 2, 1, true);
	g.ae(2, 3, 3, true);
	g.ae(1, 3, 2, true);

	g.run(0);
	for (int v = 0; v < 4; v++) cout << "dist[" << v << "] = " << g.dist[v] << "\n";

	cout << "path to 1:";
	for (int x : g.path(1)) cout << " " << x;
	cout << "\n";

	vector<ll> ways = countShortestPaths(g, 0);
	cout << "shortest paths to 3: " << ways[3] << "\n";

	vector<array<ll, 3>> es = {{0,1,4},{1,0,4},{0,2,1},{2,0,1},
	                           {2,3,3},{3,2,3},{1,3,2},{3,1,2}};
	cout << "second shortest 0->3: " << secondShortest(4, es, 0, 3) << "\n";
	return 0;
}

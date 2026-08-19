// ============================================================================
//  BELLMAN-FORD  --  single source shortest paths WITH NEGATIVE EDGES
// ----------------------------------------------------------------------------
//  O(V * E). Slower than Dijkstra but it is the only single-source option when
//  edges can be negative, and it DETECTS negative cycles.
//
//  Key fact: if there is no negative cycle, every shortest path uses at most
//  V-1 edges. So V-1 rounds of "relax every edge" is enough. If a V-th round
//  still improves something, a negative cycle is reachable.
//
//  CONTENTS
//    1.  bellmanFord      distances + negative cycle flag
//    2.  spfa             queue-based Bellman-Ford, much faster in practice
//    3.  ruinedVertices   which nodes have distance -infinity
//    4.  findNegCycle     extract an actual negative cycle
//    5.  diffConstraints  solve a system of x_i - x_j <= c
//    6.  NOTES
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

const ll INF = (ll)4e18;

struct Edge { int u, v; ll w; };

// ============================================================================
//  1. PLAIN BELLMAN-FORD
// ----------------------------------------------------------------------------
//  Returns true iff a negative cycle is REACHABLE FROM s. If it returns true
//  the dist[] array is meaningless -- do not use it.
//
//  Note the guard `dist[e.u] != INF`. Without it you relax out of unreachable
//  vertices and INF + (negative w) becomes a finite bogus distance.
//
//  Early exit: if a whole round changes nothing we are done, so the loop is
//  fast on easy inputs. To detect a negative cycle ANYWHERE (not just one
//  reachable from s), initialise every dist[] to 0 instead of dist[s] = 0.
// ============================================================================
bool bellmanFord(int n, const vector<Edge>& edges, int s, vector<ll>& dist,
                 vector<int>& par) {
	dist.assign(n, INF);
	par.assign(n, -1);
	dist[s] = 0;
	for (int round = 0; round < n - 1; round++) {
		bool changed = false;
		for (const Edge& e : edges) {
			if (dist[e.u] == INF) continue;
			if (dist[e.u] + e.w < dist[e.v]) {
				dist[e.v] = dist[e.u] + e.w;
				par[e.v] = e.u;
				changed = true;
			}
		}
		if (!changed) return false;              // converged, no negative cycle
	}
	// one extra round: any further improvement means a negative cycle
	for (const Edge& e : edges)
		if (dist[e.u] != INF && dist[e.u] + e.w < dist[e.v]) return true;
	return false;
}

// ============================================================================
//  2. SPFA  (Shortest Path Faster Algorithm)
// ----------------------------------------------------------------------------
//  Bellman-Ford, but only re-relax edges out of vertices whose distance
//  actually changed. Same O(VE) worst case (and adversarial tests DO exist),
//  but typically near-linear. Detects negative cycles: if any vertex enters the
//  queue V or more times, it is on/behind a negative cycle.
// ============================================================================
bool spfa(int n, const vector<vector<pair<int, ll>>>& adj, int s,
          vector<ll>& dist) {
	dist.assign(n, INF);
	vector<int> cnt(n, 0);
	vector<char> inq(n, 0);
	deque<int> q;
	dist[s] = 0;
	q.push_back(s);
	inq[s] = 1;
	while (!q.empty()) {
		int u = q.front();
		q.pop_front();
		inq[u] = 0;
		for (auto [v, w] : adj[u]) {
			if (dist[u] + w < dist[v]) {
				dist[v] = dist[u] + w;
				if (!inq[v]) {
					if (++cnt[v] >= n) return true;   // negative cycle
					inq[v] = 1;
					// SLF heuristic: cheap nodes to the front
					if (!q.empty() && dist[v] < dist[q.front()]) q.push_front(v);
					else q.push_back(v);
				}
			}
		}
	}
	return false;
}

// ============================================================================
//  3. RUINED VERTICES -- true shortest distance is -infinity
// ----------------------------------------------------------------------------
//  A vertex is ruined iff it is reachable from some negative cycle that is
//  itself reachable from s. Run the V-th relaxation round, mark every vertex it
//  improves, then flood out from those marks along the graph.
//  Typical use: "print -1 if you can earn unbounded profit, else the answer".
// ============================================================================
vector<char> ruinedVertices(int n, const vector<Edge>& edges, int s) {
	vector<ll> dist;
	vector<int> par;
	bellmanFord(n, edges, s, dist, par);
	vector<char> bad(n, 0);
	vector<vector<int>> adj(n);
	for (const Edge& e : edges) adj[e.u].push_back(e.v);
	// one more round: anything still improving sits on a negative cycle
	for (const Edge& e : edges)
		if (dist[e.u] != INF && dist[e.u] + e.w < dist[e.v]) bad[e.v] = 1;
	// everything reachable from those is also ruined
	vector<int> st;
	for (int i = 0; i < n; i++) if (bad[i]) st.push_back(i);
	while (!st.empty()) {
		int u = st.back();
		st.pop_back();
		for (int v : adj[u]) if (!bad[v]) { bad[v] = 1; st.push_back(v); }
	}
	return bad;
}

// ============================================================================
//  4. EXTRACT AN ACTUAL NEGATIVE CYCLE
// ----------------------------------------------------------------------------
//  Run V rounds. If round V relaxes vertex x, then walking par[] back from x
//  V times lands you INSIDE a cycle (not merely on a path to it); from there
//  follow par[] until you return to the same node.
//  Returns the cycle as a node list in order, or empty if none exists.
// ============================================================================
vector<int> findNegCycle(int n, const vector<Edge>& edges) {
	vector<ll> dist(n, 0);          // all-zero start = detect a cycle ANYWHERE
	vector<int> par(n, -1);
	int x = -1;
	for (int round = 0; round < n; round++) {
		x = -1;
		for (const Edge& e : edges)
			if (dist[e.u] + e.w < dist[e.v]) {
				dist[e.v] = dist[e.u] + e.w;
				par[e.v] = e.u;
				x = e.v;
			}
		if (x == -1) return {};                 // converged: no negative cycle
	}
	for (int i = 0; i < n; i++) x = par[x];    // guaranteed to be on the cycle
	vector<int> cyc;
	for (int v = x;; v = par[v]) {
		cyc.push_back(v);
		if (v == x && cyc.size() > 1) break;
	}
	reverse(cyc.begin(), cyc.end());
	return cyc;
}

// ============================================================================
//  5. SYSTEM OF DIFFERENCE CONSTRAINTS
// ----------------------------------------------------------------------------
//  Given m constraints of the form   x_i - x_j <= c,   find real x_1..x_n or
//  report that none exist. O(n*m) via Bellman-Ford.
//
//  Reduction: constraint (x_i - x_j <= c) becomes an EDGE j -> i of weight c,
//  because a shortest path obeys d_i <= d_j + c, which is exactly the
//  constraint. Add a virtual source with a 0-weight edge to every vertex so
//  all distances are finite (any constant works -- shifting every x_i by the
//  same amount preserves all differences).
//
//  No solution  <=>  negative cycle. (Sum the constraints around the cycle:
//  every variable cancels, leaving 0 <= (negative), a contradiction.)
//
//  REWRITING OTHER CONSTRAINTS into this form:
//     x_i - x_j >= c    ->    x_j - x_i <= -c
//     x_i - x_j =  c    ->    both <= c and >= c
//     x_i <= c          ->    x_i - x_0 <= c   with x_0 pinned to 0
//     x_i < c (integers)->    x_i <= c - 1
//  Classic uses: scheduling with "task i starts at least d after task j",
//  and "difference between prefix sums in a window is at most k".
// ============================================================================
bool diffConstraints(int n, const vector<Edge>& cons, vector<ll>& x) {
	// nodes 0..n-1 are the variables, node n is the virtual source
	vector<Edge> edges = cons;                  // each is (u=j, v=i, w=c)
	for (int i = 0; i < n; i++) edges.push_back({n, i, 0});
	vector<int> par;
	bool neg = bellmanFord(n + 1, edges, n, x, par);
	if (neg) return false;                      // infeasible
	x.resize(n);
	return true;
}

// ============================================================================
//  6. NOTES
// ----------------------------------------------------------------------------
//  WHICH ALGORITHM?
//     all weights >= 0            -> Dijkstra, O(E log V)
//     all weights == 1            -> BFS, O(V+E)
//     weights in {0,1}            -> 0-1 BFS, O(V+E)
//     negative edges, single src  -> Bellman-Ford / SPFA, O(VE)
//     negative edges, DAG         -> relax in topological order, O(V+E)
//     negative edges, all pairs   -> Floyd-Warshall, O(V^3)
//     negative edges, all pairs, sparse -> Johnson (reweight with Bellman-Ford,
//        then run Dijkstra from each node), O(VE log V)
//
//  MINIMUM MEAN CYCLE (Karp): dp[k][v] = min weight walk of exactly k edges
//  ending at v; answer = min over v of max over k of (dp[n][v]-dp[k][v])/(n-k).
//  O(VE). Used for "is there a cycle with average weight < X" -- but if you
//  only need that decision, subtract X from every edge and look for a negative
//  cycle with Bellman-Ford, then binary search X.
//
//  TRAPS
//    * Forgetting the `dist[u] != INF` guard -- silently wrong answers.
//    * Using Dijkstra with negative edges because "the test data looks fine".
//    * A negative cycle that is NOT reachable from s does not make the answer
//      -infinity. Detect reachable-only unless the problem says otherwise.
//    * SPFA is worst-case O(VE) and there are anti-SPFA tests in the wild; if
//      the constraints already allow plain Bellman-Ford, just write that.
//    * Undirected graph + one negative edge = an immediate negative 2-cycle.
// ============================================================================

int main() {
	// 0 -> 1 (1), 1 -> 2 (-2), 2 -> 3 (1), 0 -> 3 (5)
	int n = 4;
	vector<Edge> edges = {{0,1,1},{1,2,-2},{2,3,1},{0,3,5}};
	vector<ll> dist;
	vector<int> par;
	bool neg = bellmanFord(n, edges, 0, dist, par);
	cout << "negative cycle: " << (neg ? "yes" : "no") << "\n";
	for (int v = 0; v < n; v++) cout << "dist[" << v << "] = " << dist[v] << "\n";

	// now add 3 -> 1 (-1), making 1 -> 2 -> 3 -> 1 a negative cycle
	edges.push_back({3, 1, -1});
	cout << "with 3->1(-1), negative cycle: "
	     << (bellmanFord(n, edges, 0, dist, par) ? "yes" : "no") << "\n";
	vector<int> cyc = findNegCycle(n, edges);
	cout << "cycle:";
	for (int v : cyc) cout << " " << v;
	cout << "\n";

	// x1 - x0 <= 3, x2 - x1 <= 2, x0 - x2 <= -1
	vector<Edge> cons = {{0,1,3},{1,2,2},{2,0,-1}};
	vector<ll> x;
	cout << "difference constraints feasible: "
	     << (diffConstraints(3, cons, x) ? "yes" : "no") << "\n";
	return 0;
}

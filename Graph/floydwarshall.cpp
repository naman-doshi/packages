// ============================================================================
//  FLOYD-WARSHALL  --  ALL PAIRS shortest paths, O(V^3) time, O(V^2) memory
// ----------------------------------------------------------------------------
//  Handles negative edges. Detects negative cycles. Six lines of code.
//  V <= ~500 comfortably, ~1000 if the time limit is generous (1e9 ops but the
//  inner loop is a single min over contiguous memory, so it is fast).
//
//  The DP: dist(u,v,k) = shortest u->v path using only vertices 0..k-1 as
//  INTERMEDIATES.
//     dist(u,v,k) = min( dist(u,v,k-1),  dist(u,k,k-1) + dist(k,v,k-1) )
//  The k layer can be overwritten in place: allowing k as an intermediate on a
//  path to or from k itself never helps unless there is a negative cycle.
//
//  THE LOOP ORDER IS NOT NEGOTIABLE: k OUTERMOST. Getting it wrong compiles,
//  runs, and gives wrong answers on about half the tests.
//
//  CONTENTS
//    1.  floydWarshall     distances (+ negative cycle handling)
//    2.  withPaths         distances + path reconstruction
//    3.  transitiveClosure reachability, bitset-accelerated O(V^3/64)
//    4.  minimaxAllPairs   minimise the largest edge, for every pair
//    5.  NOTES
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

// Use an INF you can add to itself without overflowing. 1e18 + 1e18 overflows
// ll; 4e17 does not. Alternatively guard with `if (d[u][k]==INF) continue;`.
const ll INF = (ll)4e17;

// ============================================================================
//  1. CORE
// ----------------------------------------------------------------------------
//  d is the adjacency matrix: d[u][v] = weight of edge u->v, INF if no edge,
//  d[u][u] = 0. With parallel edges keep the minimum. Modified in place.
// ============================================================================
void floydWarshall(vector<vector<ll>>& d) {
	int n = (int)d.size();
	for (int k = 0; k < n; k++)                    // <-- k OUTERMOST
		for (int u = 0; u < n; u++) {
			if (d[u][k] == INF) continue;           // small but real speedup
			for (int v = 0; v < n; v++)
				if (d[k][v] != INF && d[u][k] + d[k][v] < d[u][v])
					d[u][v] = d[u][k] + d[k][v];
		}
}

// d[u][u] < 0 means u lies on a negative cycle. Distances involving such a
// vertex are meaningless; mark every pair whose path can be routed through one.
// Call AFTER floydWarshall. Sets those entries to -INF.
void markNegativeCycles(vector<vector<ll>>& d) {
	int n = (int)d.size();
	for (int u = 0; u < n; u++)
		for (int v = 0; v < n; v++)
			for (int k = 0; k < n; k++)
				if (d[u][k] < INF && d[k][k] < 0 && d[k][v] < INF)
					d[u][v] = -INF;
}

// ============================================================================
//  2. WITH PATH RECONSTRUCTION
// ----------------------------------------------------------------------------
//  nxt[u][v] = the next vertex after u on a shortest u->v path.
//  Initialise nxt[u][v] = v for every real edge, and nxt[u][u] = u.
// ============================================================================
void floydWithPaths(vector<vector<ll>>& d, vector<vector<int>>& nxt) {
	int n = (int)d.size();
	for (int k = 0; k < n; k++)
		for (int u = 0; u < n; u++) {
			if (d[u][k] == INF) continue;
			for (int v = 0; v < n; v++)
				if (d[k][v] != INF && d[u][k] + d[k][v] < d[u][v]) {
					d[u][v] = d[u][k] + d[k][v];
					nxt[u][v] = nxt[u][k];           // reroute through k
				}
		}
}

vector<int> recoverPath(const vector<vector<ll>>& d, const vector<vector<int>>& nxt,
                        int u, int v) {
	if (d[u][v] >= INF) return {};
	vector<int> p = {u};
	while (u != v) { u = nxt[u][v]; p.push_back(u); }
	return p;
}

// ============================================================================
//  3. TRANSITIVE CLOSURE -- "can u reach v?"
// ----------------------------------------------------------------------------
//  Same triple loop with OR instead of min. With bitsets the inner loop is one
//  word-parallel OR, giving O(V^3 / 64): V = 2000 is fine, V = 5000 is doable.
//  Needed for minimum path cover / Dilworth style reductions (see
//  flow/modeling.cpp) and for "is this partial order total".
// ============================================================================
template <int MAXN>
void transitiveClosure(vector<bitset<MAXN>>& reach) {
	int n = (int)reach.size();
	for (int k = 0; k < n; k++)
		for (int u = 0; u < n; u++)
			if (reach[u][k]) reach[u] |= reach[k];
}

// ============================================================================
//  4. ALL-PAIRS MINIMAX (bottleneck) -- minimise the LARGEST edge on the path
// ----------------------------------------------------------------------------
//  Replace (+, min) with (max, min). "For every pair, what is the lightest
//  bridge on the best route?" Also: all-pairs MAXIMIN (widest path, maximise
//  the smallest edge) is (min, max) -- swap both operators.
// ============================================================================
void floydMinimax(vector<vector<ll>>& d) {
	int n = (int)d.size();
	for (int k = 0; k < n; k++)
		for (int u = 0; u < n; u++)
			for (int v = 0; v < n; v++)
				d[u][v] = min(d[u][v], max(d[u][k], d[k][v]));
}

// ============================================================================
//  5. NOTES
// ----------------------------------------------------------------------------
//  WHEN TO REACH FOR IT
//    * V small (<= 500) and you need many pairs, or ANY pair on demand.
//    * Negative edges present and you need all-pairs.
//    * The "relation composition" trick: any (op, combine) pair where combine
//      distributes over op works -- reachability, bottleneck, counting paths
//      (sum, product) mod p, "is there a path of even length" (parity states).
//    * As a subroutine after contracting SCCs, or on a small set of important
//      vertices while the rest is handled some other way.
//
//  DO NOT USE IT when V is 1e5 and you only need one source -- that is
//  Dijkstra. Do not use it when V is 2000+ unless the operation is a bitset.
//
//  INITIALISATION CHECKLIST (this is where the bugs are)
//    * d[u][u] = 0 BEFORE the loop, else self-distances stay INF.
//    * Parallel edges: keep the min, do not overwrite blindly.
//    * Undirected: set both d[u][v] and d[v][u].
//    * A self-loop of negative weight is already a negative cycle.
//    * INF must survive INF + INF without overflow, or guard explicitly.
//
//  NEGATIVE CYCLES: after the main loop, d[u][u] < 0 for exactly the vertices
//  on a negative cycle. If a problem asks for a shortest path that may pass
//  through such a vertex, the answer is -infinity; use markNegativeCycles.
//
//  JOHNSON'S ALGORITHM (all pairs, negative edges, SPARSE graph): add a virtual
//  source with 0-weight edges to everything, run Bellman-Ford to get h[v], then
//  reweight w'(u,v) = w(u,v) + h[u] - h[v] >= 0 and run Dijkstra from every
//  node. O(VE log V). Recover true distances by subtracting h[u] - h[v].
//  Rarely needed in contests but it is the right answer for V = 1e4, E = 1e5.
// ============================================================================

int main() {
	int n = 4;
	vector<vector<ll>> d(n, vector<ll>(n, INF));
	vector<vector<int>> nxt(n, vector<int>(n, -1));
	for (int i = 0; i < n; i++) { d[i][i] = 0; nxt[i][i] = i; }

	auto add = [&](int u, int v, ll w) {
		if (w < d[u][v]) { d[u][v] = w; nxt[u][v] = v; }
	};
	add(0, 1, 5); add(1, 2, 3); add(2, 3, 1); add(0, 3, 100); add(3, 0, 2);

	floydWithPaths(d, nxt);
	cout << "dist 0->3 = " << d[0][3] << "\n";
	cout << "path:";
	for (int x : recoverPath(d, nxt, 0, 3)) cout << " " << x;
	cout << "\n";

	cout << "negative cycle present: "
	     << ([&] { for (int i = 0; i < n; i++) if (d[i][i] < 0) return true;
	               return false; }() ? "yes" : "no") << "\n";
	return 0;
}

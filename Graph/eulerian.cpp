// ============================================================================
//  EULERIAN PATHS AND CIRCUITS  (Hierholzer's algorithm, O(V + E))
// ----------------------------------------------------------------------------
//  An Eulerian CIRCUIT uses every EDGE exactly once and returns to the start.
//  An Eulerian PATH uses every edge exactly once, ending somewhere else.
//  (Contrast: a Hamiltonian path visits every VERTEX once and is NP-hard.
//   If a problem says "every road exactly once" it is Euler and it is easy;
//   "every city exactly once" and it is Hamiltonian and it is not.)
//
//  EXISTENCE CONDITIONS -- check these before running anything.
//    UNDIRECTED
//      circuit: every vertex has EVEN degree, and all edges are in one
//               connected component.
//      path:    exactly 0 or 2 vertices of ODD degree (start and end at those
//               two), all edges in one component.
//    DIRECTED
//      circuit: indeg(v) == outdeg(v) for every v, and all edges lie in one
//               strongly connected component (ignoring isolated vertices).
//      path:    one vertex with outdeg - indeg == +1 (the start), one with
//               indeg - outdeg == +1 (the end), all others equal.
//    Isolated vertices with no edges never matter -- only check connectivity
//    over vertices that actually have an edge.
//
//  CONTENTS
//    1.  EulerUndirected
//    2.  EulerDirected
//    3.  NOTES
// ============================================================================
#include <bits/stdc++.h>
using namespace std;

// ============================================================================
//  1. UNDIRECTED EULERIAN PATH / CIRCUIT
// ----------------------------------------------------------------------------
//  Hierholzer: walk greedily, never reusing an edge, until stuck; whatever you
//  are stuck on goes onto the answer, and you back up and keep walking. The
//  reversed order of "stuck" vertices is the tour.
//
//  Iterative, so it survives 1e5+ edges. `used[]` is indexed by EDGE ID so a
//  pair of parallel edges is correctly treated as two distinct edges, and the
//  `it[]` pointer means every incident edge is examined only once overall --
//  that is what keeps it O(V+E) rather than O(V*E).
//
//  API  (nodes 0-based)
//    EulerUndirected e; e.init(N);
//    e.ae(u, v);
//    vector<int> tour = e.solve();      // node sequence, or empty if none
//    // tour.size() == E + 1 when it succeeds
// ============================================================================
struct EulerUndirected {
	int N, M = 0;
	vector<vector<pair<int, int>>> adj;    // (neighbour, edge id)
	vector<char> used;

	void init(int _N) { N = _N; M = 0; adj.assign(N, {}); used.clear(); }
	void ae(int u, int v) {
		adj[u].push_back({v, M});
		adj[v].push_back({u, M});
		M++;
	}

	vector<int> solve() {
		used.assign(M, 0);
		if (M == 0) return {};

		// existence check: 0 or 2 odd-degree vertices
		vector<int> odd;
		for (int i = 0; i < N; i++)
			if (adj[i].size() % 2) odd.push_back(i);
		if (odd.size() != 0 && odd.size() != 2) return {};

		int start = odd.empty() ? -1 : odd[0];
		if (start == -1)
			for (int i = 0; i < N; i++) if (!adj[i].empty()) { start = i; break; }

		vector<int> it(N, 0), st = {start}, tour;
		while (!st.empty()) {
			int u = st.back();
			while (it[u] < (int)adj[u].size() && used[adj[u][it[u]].second]) it[u]++;
			if (it[u] == (int)adj[u].size()) {       // stuck: retreat
				tour.push_back(u);
				st.pop_back();
			} else {
				auto [v, id] = adj[u][it[u]];
				used[id] = 1;
				st.push_back(v);
			}
		}
		if ((int)tour.size() != M + 1) return {};   // edges were disconnected
		reverse(tour.begin(), tour.end());
		return tour;
	}
};

// ============================================================================
//  2. DIRECTED EULERIAN PATH / CIRCUIT
// ----------------------------------------------------------------------------
//  Same algorithm; the existence test is the in/out degree one instead.
//  Note the tour must be REVERSED at the end (we build it backwards).
// ============================================================================
struct EulerDirected {
	int N, M = 0;
	vector<vector<pair<int, int>>> adj;    // (destination, edge id)
	vector<int> indeg;

	void init(int _N) { N = _N; M = 0; adj.assign(N, {}); indeg.assign(N, 0); }
	void ae(int u, int v) { adj[u].push_back({v, M++}); indeg[v]++; }

	vector<int> solve() {
		if (M == 0) return {};
		int start = -1, plus = 0, minus = 0;
		for (int i = 0; i < N; i++) {
			int d = (int)adj[i].size() - indeg[i];
			if (d == 1)       { plus++;  start = i; }
			else if (d == -1) { minus++; }
			else if (d != 0)  { return {}; }
		}
		if (plus > 1 || minus > 1) return {};
		if ((plus == 0) != (minus == 0)) return {};   // both zero, or one each
		if (start == -1)
			for (int i = 0; i < N; i++) if (!adj[i].empty()) { start = i; break; }

		vector<int> it(N, 0), st = {start}, tour;
		while (!st.empty()) {
			int u = st.back();
			if (it[u] == (int)adj[u].size()) {
				tour.push_back(u);
				st.pop_back();
			} else {
				int v = adj[u][it[u]++].first;
				st.push_back(v);
			}
		}
		if ((int)tour.size() != M + 1) return {};
		reverse(tour.begin(), tour.end());
		return tour;
	}
};

// ============================================================================
//  3. NOTES
// ----------------------------------------------------------------------------
//  DO NOT WRITE THE NAIVE VERSION. Erasing edges from the adjacency list, or
//  rescanning it from index 0 each time, turns this into O(V*E) and it will
//  time out. The `it[]` cursor is the whole trick (same idea as the "uptochild"
//  pointer in Dinic's algorithm).
//
//  RECURSIVE HIERHOLZER blows the stack at ~1e5 edges. Use the iterative form
//  above.
//
//  CHINESE POSTMAN (route inspection): traverse every edge at least once,
//  minimising total length. If all degrees are even, that is exactly the Euler
//  circuit. Otherwise, pair up the odd-degree vertices at minimum total
//  shortest-path cost (a minimum-weight perfect matching on the odd vertices --
//  with few odd vertices, a bitmask DP over them; the number of odd-degree
//  vertices is always even) and duplicate those paths' edges.
//
//  DE BRUIJN SEQUENCES are an Euler circuit on the graph whose vertices are
//  length-(n-1) strings and whose edges are length-n strings. Any "shortest
//  string containing every k-mer" problem is this.
//
//  RECONSTRUCTING A SEQUENCE FROM OVERLAPS (dominoes end-to-end, words where
//  the last letter matches the next word's first letter) is a directed Euler
//  path where letters are vertices and words are edges. This is the single most
//  common disguise.
//
//  COUNTING: BEST theorem counts Eulerian circuits in a directed graph via the
//  Matrix-Tree theorem. Almost never needed.
// ============================================================================

int main() {
	// square with a diagonal: degrees 3,3,2,2 -> Eulerian PATH, not a circuit
	EulerUndirected e;
	e.init(4);
	e.ae(0, 1); e.ae(1, 2); e.ae(2, 3); e.ae(3, 0); e.ae(0, 2);
	vector<int> t = e.solve();
	cout << "undirected euler tour (" << t.size() << " nodes):";
	for (int v : t) cout << " " << v;
	cout << "\n";

	// directed 3-cycle -> Eulerian circuit
	EulerDirected d;
	d.init(3);
	d.ae(0, 1); d.ae(1, 2); d.ae(2, 0);
	vector<int> t2 = d.solve();
	cout << "directed euler circuit:";
	for (int v : t2) cout << " " << v;
	cout << "\n";
	return 0;
}

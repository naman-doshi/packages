// ============================================================================
//  HOPCROFT-KARP  --  MAXIMUM BIPARTITE MATCHING in O(E * sqrt(V))
// ----------------------------------------------------------------------------
//  Far faster than Kuhn's O(V*E) augmenting-path matching on big graphs. Left
//  vertices 0..n-1, right vertices 0..m-1 (two separate index spaces). Add an
//  edge for every left-right pair that MAY be matched, then call match().
//
//  API
//    HopcroftKarp h; h.init(n, m);      // n left nodes, m right nodes
//    h.ae(u, v);                        // left u -- right v may be matched
//    int k = h.match();                 // size of maximum matching
//    // afterwards: h.ml[u] = right node matched to left u  (-1 if unmatched)
//    //             h.mr[v] = left  node matched to right v (-1 if unmatched)
//
//  KONIG'S THEOREM  (bipartite): with a maximum matching in hand,
//    min vertex cover      = k
//    max independent set   = (n + m) - k
//    min edge cover        = (n + m) - k        (if no isolated vertices)
//  To recover the actual min vertex cover / max independent set, run an
//  alternating-path search from unmatched left nodes -- see minVertexCover().
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long

struct HopcroftKarp {
	int n, m;
	vector<vector<int>> adj;
	vector<int> ml, mr, dist;
	void init(int _n, int _m) {
		n = _n; m = _m;
		adj.assign(n, {});
		ml.assign(n, -1); mr.assign(m, -1);
	}
	void ae(int u, int v) { adj[u].push_back(v); }
	bool bfs() {
		queue<int> q; dist.assign(n, -1);
		for (int u = 0; u < n; u++) if (ml[u] < 0) { dist[u] = 0; q.push(u); }
		bool found = false;
		while (!q.empty()) {
			int u = q.front(); q.pop();
			for (int v : adj[u]) {
				int w = mr[v];
				if (w < 0) found = true;
				else if (dist[w] < 0) { dist[w] = dist[u] + 1; q.push(w); }
			}
		}
		return found;
	}
	bool dfs(int u) {
		for (int v : adj[u]) {
			int w = mr[v];
			if (w < 0 || (dist[w] == dist[u] + 1 && dfs(w))) {
				ml[u] = v; mr[v] = u; return true;
			}
		}
		dist[u] = -1; return false;
	}
	int match() {
		int res = 0;
		while (bfs())
			for (int u = 0; u < n; u++)
				if (ml[u] < 0 && dfs(u)) res++;
		return res;
	}
	// Konig: returns {leftInCover, rightInCover}. Call AFTER match().
	// Min vertex cover = (left NOT reachable) + (right reachable) via
	// alternating paths from unmatched left vertices.
	pair<vector<int>, vector<int>> minVertexCover() {
		vector<char> visL(n, 0), visR(m, 0);
		// alternating BFS from every unmatched left node
		queue<int> q;
		for (int u = 0; u < n; u++) if (ml[u] < 0) { visL[u] = 1; q.push(u); }
		while (!q.empty()) {
			int u = q.front(); q.pop();
			for (int v : adj[u]) if (!visR[v]) {
				visR[v] = 1;
				int w = mr[v];               // follow the matched edge back to left
				if (w >= 0 && !visL[w]) { visL[w] = 1; q.push(w); }
			}
		}
		vector<int> L, R;
		for (int u = 0; u < n; u++) if (!visL[u]) L.push_back(u);   // unreached left
		for (int v = 0; v < m; v++) if (visR[v]) R.push_back(v);    // reached  right
		return {L, R};
	}
};

int32_t main() {
	// 3 left, 3 right.  Edges: 0-{0,1}, 1-{0}, 2-{1,2}
	HopcroftKarp h; h.init(3, 3);
	h.ae(0, 0); h.ae(0, 1);
	h.ae(1, 0);
	h.ae(2, 1); h.ae(2, 2);

	int k = h.match();
	cout << "max matching = " << k << "\n";
	for (int u = 0; u < 3; u++)
		cout << "  left " << u << " <-> right " << h.ml[u] << "\n";

	auto [L, R] = h.minVertexCover();
	cout << "min vertex cover size = " << L.size() + R.size() << "\n";
	return 0;
}

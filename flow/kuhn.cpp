// ============================================================================
//  KUHN'S ALGORITHM  --  MAXIMUM BIPARTITE MATCHING  in O(V * E)
// ----------------------------------------------------------------------------
//  The short, easy-to-remember augmenting-path matcher. Prefer this when the
//  graph is small (a few thousand edges) and you value fewer lines over speed;
//  reach for hopcroftKarp.cpp when the graph is large.
//
//  API
//    Kuhn k; k.init(n, m);       // n left nodes, m right nodes (0-based each)
//    k.ae(u, v);                 // left u may match right v
//    int sz = k.match();         // maximum matching size
//    // k.mr[v] = left node matched to right v (-1 if none)
//    // k.ml[u] = right node matched to left u (-1 if none)
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long

struct Kuhn {
	int n, m;
	vector<vector<int>> adj;
	vector<int> ml, mr;
	vector<char> used;
	void init(int _n, int _m) {
		n = _n; m = _m;
		adj.assign(n, {});
		ml.assign(n, -1); mr.assign(m, -1);
	}
	void ae(int u, int v) { adj[u].push_back(v); }
	bool tryKuhn(int u) {
		for (int v : adj[u]) if (!used[v]) {
			used[v] = 1;
			if (mr[v] < 0 || tryKuhn(mr[v])) {
				ml[u] = v; mr[v] = u; return true;
			}
		}
		return false;
	}
	int match() {
		int res = 0;
		for (int u = 0; u < n; u++) {
			used.assign(m, 0);
			if (tryKuhn(u)) res++;
		}
		return res;
	}
};

int32_t main() {
	Kuhn k; k.init(3, 3);
	k.ae(0, 0); k.ae(0, 1);
	k.ae(1, 0);
	k.ae(2, 1); k.ae(2, 2);
	cout << "max matching = " << k.match() << "\n";
	for (int u = 0; u < 3; u++)
		cout << "  left " << u << " <-> right " << k.ml[u] << "\n";
	return 0;
}

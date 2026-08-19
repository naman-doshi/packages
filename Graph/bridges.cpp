// ============================================================================
//  BRIDGES  (undirected graph, O(V + E))
// ----------------------------------------------------------------------------
//  An edge is a bridge iff removing it increases the number of connected
//  components. One DFS pass: edge (x, y) with y a child of x is a bridge iff
//  low[y] > disc[x], i.e. the subtree of y cannot reach x or above except
//  through this edge.
//
//  Parallel edges are handled by skipping the edge INDEX we came in on (not the
//  parent NODE), so a doubled edge is correctly reported as no bridge.
//  Self-loops are never bridges. Works on disconnected graphs.
//
//  API  (nodes 0-based)
//    Bridges b; b.init(N);
//    b.ae(u, v);            // add undirected edge u -- v
//    vector<pair<int,int>> br = b.gen();   // list of bridges as (u, v) pairs
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long
#define ll long long

const int INF = 9 * 1e18;

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define trav(a, x) for(auto& a : x)
#define all(x) x.begin(), x.end()
#define print(arr) for (auto i : arr) cout << i << " "; cout << endl;
typedef vector<vector<int>> vvi;
typedef vector<pair<int, int>> vpi;
typedef vector<int> vi;
typedef pair<int, int> pii;

struct Bridges {
	int N, ti = 0;
	vector<vector<pair<int, int>>> adj;  // (neighbour, edge id)
	vector<pair<int, int>> edges, bridges;
	vector<int> disc, low;
	void init(int _N) {
		N = _N;
		adj.resize(N), disc = vector<int>(N, 0), low.resize(N);
	}
	void ae(int x, int y) {
		int id = edges.size();
		edges.push_back({x, y});
		adj[x].push_back({y, id}), adj[y].push_back({x, id});
	}
	void dfs(int x, int pe) {  // pe = edge id we entered x through
		low[x] = disc[x] = ++ti;
		for (auto [y, id] : adj[x]) {
			if (id == pe) continue;  // skip the edge in, not the parent node
			if (disc[y]) low[x] = min(low[x], disc[y]);
			else {
				dfs(y, id);
				low[x] = min(low[x], low[y]);
				if (low[y] > disc[x]) bridges.push_back(edges[id]);
			}
		}
	}
	vector<pair<int, int>> gen() {
		for (int i = 0; i < N; i++)
			if (!disc[i]) dfs(i, -1);
		return bridges;
	}
};

int32_t main() {
	Bridges b;
	int numNodes = 6;
	b.init(numNodes);

	// 0-1-2 is a cycle, 2--3 and 3--4 are bridges, 5 is isolated
	b.ae(0, 1);
	b.ae(1, 2);
	b.ae(2, 0);
	b.ae(2, 3);
	b.ae(3, 4);

	vector<pair<int, int>> br = b.gen();

	cout << "Found " << br.size() << " bridges:" << endl;
	for (auto [u, v] : br) cout << u << " -- " << v << endl;

	return 0;
}

// ============================================================================
//  MINIMUM SPANNING TREE  --  Kruskal, Prim, and the DSU they need
// ----------------------------------------------------------------------------
//  A spanning tree is a subgraph that is a tree and touches every vertex; the
//  MST is the one of minimum total edge weight.
//
//  WHY GREEDY WORKS (cut property): if e is the cheapest edge crossing some
//  partition of the vertices into two non-empty sides, then some MST contains
//  e. Proof: take an MST without e, add e -- this makes a cycle, which must
//  cross the same cut on some other edge f with w(f) >= w(e). Swap f for e.
//  Every MST algorithm is just a different rule for picking which cut to look
//  at next.
//
//  CONTENTS
//    1.  DSU          union-find (needed by Kruskal, useful everywhere)
//    2.  kruskal      O(E log E). The default.
//    3.  prim         O(E log V) sparse, and O(V^2) dense
//    4.  secondBestMST
//    5.  NOTES -- variants, uniqueness, minimax, maximum spanning tree
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

// ============================================================================
//  1. DISJOINT SET UNION  (union by size + path compression)
// ----------------------------------------------------------------------------
//  Effectively O(1) amortised (inverse Ackermann).
//  find(x)      representative of x's component
//  unite(a,b)   returns false if they were already together
//  size(x)      number of elements in x's component
// ============================================================================
struct DSU {
	vector<int> p, sz;
	int comps;
	void init(int n) {
		p.resize(n);
		iota(p.begin(), p.end(), 0);
		sz.assign(n, 1);
		comps = n;
	}
	int find(int x) { return p[x] == x ? x : p[x] = find(p[x]); }
	bool same(int a, int b) { return find(a) == find(b); }
	int size(int x) { return sz[find(x)]; }
	bool unite(int a, int b) {
		a = find(a); b = find(b);
		if (a == b) return false;
		if (sz[a] < sz[b]) swap(a, b);      // attach smaller under larger
		p[b] = a;
		sz[a] += sz[b];
		comps--;
		return true;
	}
};

struct Edge {
	int u, v;
	ll w;
	bool operator<(const Edge& o) const { return w < o.w; }
};

// ============================================================================
//  2. KRUSKAL -- O(E log E), dominated by the sort
// ----------------------------------------------------------------------------
//  Sort every edge by weight, take it if it joins two different components.
//  Returns the total weight; `used` receives the chosen edges.
//  If the graph is disconnected you get a minimum spanning FOREST -- check
//  dsu.comps == 1 to detect that (a very common "output -1" case).
// ============================================================================
ll kruskal(int n, vector<Edge> edges, vector<Edge>* used = nullptr) {
	sort(edges.begin(), edges.end());
	DSU d;
	d.init(n);
	ll total = 0;
	for (const Edge& e : edges) {
		if (d.unite(e.u, e.v)) {
			total += e.w;
			if (used) used->push_back(e);
		}
	}
	return total;           // caller should check d.comps == 1 for connectivity
}

// ============================================================================
//  3. PRIM -- grow ONE component outwards from a start vertex
// ----------------------------------------------------------------------------
//  Structurally identical to Dijkstra, except the priority is the edge weight
//  itself rather than the accumulated distance. Use it when the graph is dense
//  (the O(V^2) form needs no heap and no edge list) or already an adjacency
//  list you do not want to flatten.
// ============================================================================
ll prim(const vector<vector<pair<int, ll>>>& adj) {
	int n = (int)adj.size();
	vector<char> in(n, 0);
	priority_queue<pair<ll, int>, vector<pair<ll, int>>,
	               greater<pair<ll, int>>> pq;
	ll total = 0;
	int taken = 0;
	pq.push({0, 0});
	while (!pq.empty()) {
		auto [w, u] = pq.top();
		pq.pop();
		if (in[u]) continue;               // already absorbed, at a better price
		in[u] = 1;
		total += w;
		taken++;
		for (auto [v, ww] : adj[u]) if (!in[v]) pq.push({ww, v});
	}
	return taken == n ? total : -1;        // -1 = graph is disconnected
}

// Dense version: O(V^2), no heap. w[u][v] = INF for absent edges.
ll primDense(const vector<vector<ll>>& w) {
	const ll INF = (ll)4e18;
	int n = (int)w.size();
	vector<ll> best(n, INF);
	vector<char> in(n, 0);
	best[0] = 0;
	ll total = 0;
	for (int it = 0; it < n; it++) {
		int u = -1;
		for (int i = 0; i < n; i++)
			if (!in[i] && (u == -1 || best[i] < best[u])) u = i;
		if (best[u] == INF) return -1;     // disconnected
		in[u] = 1;
		total += best[u];
		for (int v = 0; v < n; v++)
			if (!in[v] && w[u][v] < best[v]) best[v] = w[u][v];
	}
	return total;
}

// ============================================================================
//  4. SECOND BEST MST
// ----------------------------------------------------------------------------
//  The second-best spanning tree differs from the MST by exactly one edge swap.
//  So: build the MST, then for every NON-tree edge (u,v,w), find the heaviest
//  edge on the MST path u..v and consider replacing it. Take the best swap.
//
//  Below is the O(E * V) version (walk the tree path by BFS each time), which
//  is plenty for V,E <= ~2000. For larger inputs, precompute "max edge on the
//  path to the 2^k-th ancestor" with binary lifting (see Tree/lca.cpp) and each
//  query becomes O(log V).
// ============================================================================
ll secondBestMST(int n, const vector<Edge>& edges) {
	const ll INF = (ll)4e18;
	vector<Edge> tree;
	ll mstW = kruskal(n, edges, &tree);

	// adjacency of the MST, remembering each edge's weight
	vector<vector<pair<int, ll>>> adj(n);
	set<pair<int, int>> inTree;
	for (const Edge& e : tree) {
		adj[e.u].push_back({e.v, e.w});
		adj[e.v].push_back({e.u, e.w});
		inTree.insert({min(e.u, e.v), max(e.u, e.v)});
	}

	ll best = INF;
	for (const Edge& e : edges) {
		if (inTree.count({min(e.u, e.v), max(e.u, e.v)})) continue;
		// heaviest edge on the tree path e.u -> e.v, by BFS
		vector<ll> mx(n, -1);
		vector<char> seen(n, 0);
		queue<int> q;
		q.push(e.u);
		seen[e.u] = 1;
		mx[e.u] = 0;
		while (!q.empty()) {
			int u = q.front();
			q.pop();
			for (auto [v, w] : adj[u])
				if (!seen[v]) { seen[v] = 1; mx[v] = max(mx[u], w); q.push(v); }
		}
		if (mx[e.v] > 0) best = min(best, mstW - mx[e.v] + e.w);
	}
	return best == INF ? -1 : best;
}

// ============================================================================
//  5. NOTES
// ----------------------------------------------------------------------------
//  KRUSKAL OR PRIM?  Kruskal, almost always -- you usually already have an edge
//  list, and DSU is code you want anyway. Prim wins on dense graphs (V^2 form,
//  no edge list, no sort) and on implicit graphs where edges are generated on
//  demand (e.g. complete graph on points with Euclidean distance).
//
//  MAXIMUM SPANNING TREE: sort descending (or negate the weights). Same proof.
//
//  MINIMAX / BOTTLENECK: the path between u and v in the MST minimises the
//  maximum edge weight, simultaneously for every pair. So "smallest weight
//  limit letting a truck get from u to v" = heaviest edge on the MST path.
//  Build the MST, then answer path-max queries with binary lifting.
//  Equivalently: add edges in increasing order and note when u and v first
//  become connected in the DSU -- that edge's weight is the answer (this is
//  the offline "Kruskal reconstruction tree" idea).
//
//  UNIQUENESS: the MST is unique if all edge weights are distinct. With ties it
//  may not be; to test uniqueness, check whether any non-tree edge ties the
//  heaviest edge on its tree path (that is secondBestMST returning == mstW).
//
//  MST IS INVARIANT under any strictly increasing transform of the weights
//  (the algorithm only ever compares them). It is NOT the shortest path tree --
//  never use an MST to answer a distance query.
//
//  COMMON VARIANTS
//    * Some edges are forced in: unite their endpoints in the DSU first, then
//      run Kruskal on the rest.
//    * Some edges are forbidden: just drop them.
//    * "Connect all cities, but building a power plant in city i costs c_i":
//      add a virtual node 0 with an edge of weight c_i to every city, then MST.
//      This virtual-node trick handles most "or pay a fixed alternative" MSTs.
//    * Degree-constrained / Steiner tree: NP-hard, not an MST problem.
//    * Boruvka's algorithm (O(E log V), every component picks its cheapest
//      outgoing edge each round) is the one that parallelises and the one used
//      for MST on XOR distances with a trie.
//
//  DSU BEYOND MST: connectivity queries, cycle detection in undirected graphs,
//  Kruskal-style offline processing, small-to-large merging, and "DSU with
//  rollback" (skip path compression, union by size only, keep an undo stack)
//  for dynamic connectivity offline.
// ============================================================================

int main() {
	int n = 5;
	vector<Edge> edges = {{0,1,2},{0,3,6},{1,2,3},{1,3,8},{1,4,5},{2,4,7},{3,4,9}};

	vector<Edge> used;
	cout << "MST weight = " << kruskal(n, edges, &used) << "\n";
	for (const Edge& e : used) cout << "  " << e.u << " -- " << e.v << "  (" << e.w << ")\n";

	vector<vector<pair<int, ll>>> adj(n);
	for (const Edge& e : edges) { adj[e.u].push_back({e.v, e.w}); adj[e.v].push_back({e.u, e.w}); }
	cout << "prim = " << prim(adj) << "\n";
	cout << "second best MST = " << secondBestMST(n, edges) << "\n";
	return 0;
}

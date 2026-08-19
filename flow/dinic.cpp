// ============================================================================
//  DINIC'S MAX FLOW  (directed graph, O(V^2 * E) general, O(E*sqrt(V)) on unit
//  capacity / bipartite matching graphs -- the fastest general-purpose max flow
//  you will realistically ever need in a contest)
// ----------------------------------------------------------------------------
//  Computes the maximum flow from a source s to a sink t. By the max-flow /
//  min-cut theorem the returned value also equals the minimum cut. After
//  maxflow(), call minCutSide(s) to recover WHICH nodes lie on the source side
//  of a minimum cut, and cutEdges(s) to list the actual edges crossing it.
//
//  API  (nodes 0-based)
//    Dinic d; d.init(N);
//    d.ae(u, v, cap);          // directed edge u -> v with capacity cap
//    d.ae(u, v, cap, cap);     // UNdirected edge (both directions cap)
//    ll f = d.maxflow(s, t);   // value of the max flow  (== min cut)
//    // residual capacities now live in d.edges; flow on forward edge id is
//    //   (original cap) - edges[id].cap, i.e. edges[id^1].cap.
//    auto side = d.minCutSide(s);   // side[v] = true iff v is on s's side
//    auto cut  = d.cutEdges(s);     // vector of {u,v} edges crossing the cut
//
//  COMMON MODELLING
//    * vertex capacity c: split v into v_in -> v_out with capacity c.
//    * many sources/sinks: super-source S -> each source (cap = supply),
//      each sink -> super-sink T (cap = demand).
//    * bipartite matching: S -> left (cap 1), left -> right (cap 1),
//      right -> T (cap 1); max matching == maxflow.  (see hopcroftKarp.cpp)
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long
#define ll long long

const ll INF = 9e18;

struct Dinic {
	struct Edge { int to; ll cap; };
	int N;
	vector<Edge> edges;              // paired: forward id even, backward id^1 odd
	vector<vector<int>> g;           // g[u] = ids of edges leaving u
	vector<int> level, it;
	void init(int n) { N = n; g.assign(n, {}); edges.clear(); }
	// add edge u->v with capacity cap; rcap>0 makes it undirected
	void ae(int u, int v, ll cap, ll rcap = 0) {
		g[u].push_back(edges.size()); edges.push_back({v, cap});
		g[v].push_back(edges.size()); edges.push_back({u, rcap});
	}
	bool bfs(int s, int t) {
		level.assign(N, -1); level[s] = 0;
		queue<int> q; q.push(s);
		while (!q.empty()) {
			int u = q.front(); q.pop();
			for (int id : g[u]) {
				auto &e = edges[id];
				if (e.cap > 0 && level[e.to] < 0) {
					level[e.to] = level[u] + 1;
					q.push(e.to);
				}
			}
		}
		return level[t] >= 0;
	}
	ll dfs(int u, int t, ll f) {
		if (u == t) return f;
		for (int &i = it[u]; i < (int)g[u].size(); i++) {
			int id = g[u][i]; auto &e = edges[id];
			if (e.cap > 0 && level[e.to] == level[u] + 1) {
				ll d = dfs(e.to, t, min(f, e.cap));
				if (d > 0) { e.cap -= d; edges[id ^ 1].cap += d; return d; }
			}
		}
		return 0;
	}
	ll maxflow(int s, int t) {
		ll flow = 0;
		while (bfs(s, t)) {
			it.assign(N, 0);
			for (ll f; (f = dfs(s, t, INF)) > 0; ) flow += f;
		}
		return flow;
	}
	// call AFTER maxflow(): nodes reachable from s in the residual graph.
	// side[v] == true  <=>  v is on the source side of a minimum cut.
	vector<char> minCutSide(int s) {
		vector<char> side(N, 0);
		queue<int> q; q.push(s); side[s] = 1;
		while (!q.empty()) {
			int u = q.front(); q.pop();
			for (int id : g[u]) if (edges[id].cap > 0 && !side[edges[id].to]) {
				side[edges[id].to] = 1; q.push(edges[id].to);
			}
		}
		return side;
	}
	// call AFTER maxflow(): the saturated edges crossing the min cut.
	vector<pair<int,int>> cutEdges(int s) {
		auto side = minCutSide(s);
		vector<pair<int,int>> res;
		for (int id = 0; id < (int)edges.size(); id += 2) {
			int u = edges[id ^ 1].to, v = edges[id].to;   // forward edge u->v
			if (side[u] && !side[v]) res.push_back({u, v});
		}
		return res;
	}
};

int32_t main() {
	//        1 --3-- 3
	//      / |       | \
	//    s=0 2       4  t=5      (a small network)
	//      \ |       | /
	//        2 --1-- 4
	Dinic d; d.init(6);
	d.ae(0, 1, 3);
	d.ae(0, 2, 2);
	d.ae(1, 3, 3);
	d.ae(1, 2, 2);
	d.ae(2, 4, 1);
	d.ae(3, 5, 2);
	d.ae(3, 4, 4);
	d.ae(4, 5, 3);

	ll f = d.maxflow(0, 5);
	cout << "max flow = min cut = " << f << "\n";

	cout << "min-cut edges (source side -> sink side):\n";
	for (auto [u, v] : d.cutEdges(0)) cout << "  " << u << " -> " << v << "\n";
	return 0;
}

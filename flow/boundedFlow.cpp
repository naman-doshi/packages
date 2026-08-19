// ============================================================================
//  MAX FLOW WITH LOWER BOUNDS  (feasible flow / bounded circulation)
//  built on top of Dinic. O(V^2 * E).
// ----------------------------------------------------------------------------
//  Each edge carries a LOWER bound L and an UPPER bound C on its flow: every
//  unit L..C. The classic super-source/super-sink construction reduces this to
//  an ordinary max flow:
//    * for edge u->v with bounds [L, C], add a residual edge u->v of cap C-L
//      and record excess[v] += L, excess[u] -= L.
//    * add super-source SS and super-sink TT; for each node x:
//        excess[x] > 0 : SS -> x  (cap  excess[x])
//        excess[x] < 0 : x  -> TT (cap -excess[x])
//    * a feasible flow exists  <=>  the SS->TT max flow saturates every SS edge.
//
//  This struct handles the three standard queries:
//    feasibleCirculation()      -- is there a valid circulation? (add t->s inf
//                                   edge yourself for s/t flow, or use below)
//    maxFlowWithLowerBounds(s,t)-- max s->t flow respecting all lower bounds
//    minFlowWithLowerBounds(s,t)-- MIN feasible s->t flow
//
//  API  (nodes 0-based, 0..N-1 are the real nodes)
//    BoundedFlow bf; bf.init(N);
//    bf.ae(u, v, lo, hi);                 // edge u->v with flow in [lo, hi]
//    ll f; bool ok = bf.maxFlowWithLowerBounds(s, t, f);   // ok=false if infeasible
//    // real flow on the i-th added edge = bf.flowOnEdge(i);
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long
#define ll long long

const ll INF = 9e18;

struct Dinic {
	struct Edge { int to; ll cap; };
	int N;
	vector<Edge> edges;
	vector<vector<int>> g;
	vector<int> level, it;
	void init(int n) { N = n; g.assign(n, {}); edges.clear(); }
	int ae(int u, int v, ll cap) {
		int id = edges.size();
		g[u].push_back(edges.size()); edges.push_back({v, cap});
		g[v].push_back(edges.size()); edges.push_back({u, 0});
		return id;
	}
	bool bfs(int s, int t) {
		level.assign(N, -1); level[s] = 0;
		queue<int> q; q.push(s);
		while (!q.empty()) {
			int u = q.front(); q.pop();
			for (int id : g[u]) if (edges[id].cap > 0 && level[edges[id].to] < 0) {
				level[edges[id].to] = level[u] + 1; q.push(edges[id].to);
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
		while (bfs(s, t)) { it.assign(N, 0); for (ll f; (f = dfs(s, t, INF)) > 0; ) flow += f; }
		return flow;
	}
};

struct BoundedFlow {
	int N, SS, TT;
	Dinic d;
	vector<ll> excess, low;
	vector<int> edgeId;          // Dinic edge id of each user edge (the C-L edge)
	void init(int n) {
		N = n; SS = n; TT = n + 1;
		d.init(n + 2);
		excess.assign(n + 2, 0);
		low.clear(); edgeId.clear();
	}
	// edge u->v with flow in [lo, hi]; returns its index (for flowOnEdge)
	int ae(int u, int v, ll lo, ll hi) {
		excess[u] -= lo; excess[v] += lo;
		int id = d.ae(u, v, hi - lo);
		low.push_back(lo); edgeId.push_back(id);
		return (int)low.size() - 1;
	}
	// wire the super-source/sink from the excess vector; returns total demand
	ll wireExcess() {
		ll need = 0;
		for (int x = 0; x < N; x++) {
			if (excess[x] > 0) { d.ae(SS, x, excess[x]); need += excess[x]; }
			else if (excess[x] < 0) d.ae(x, TT, -excess[x]);
		}
		return need;
	}
	// feasible circulation (no distinguished s/t). true iff one exists.
	bool feasibleCirculation() {
		ll need = wireExcess();
		return d.maxflow(SS, TT) == need;
	}
	// max feasible s->t flow. writes it to outFlow; returns false if infeasible.
	bool maxFlowWithLowerBounds(int s, int t, ll &outFlow) {
		int back = d.ae(t, s, INF);         // circulation edge t->s
		ll need = wireExcess();
		if (d.maxflow(SS, TT) != need) return false;   // infeasible
		// flow already on t->s edge is a feasible base; free it and push more s->t
		ll base = d.edges[back ^ 1].cap;    // flow pushed around the cycle
		d.edges[back].cap = d.edges[back ^ 1].cap = 0;  // cut the t->s edge
		outFlow = base + d.maxflow(s, t);
		return true;
	}
	// min feasible s->t flow.
	bool minFlowWithLowerBounds(int s, int t, ll &outFlow) {
		int back = d.ae(t, s, INF);
		ll need = wireExcess();
		ll f1 = d.maxflow(SS, TT);
		ll base = d.edges[back ^ 1].cap;
		// remove t->s edge, then saturate any remaining SS->TT demand
		d.edges[back].cap = d.edges[back ^ 1].cap = 0;
		f1 += d.maxflow(SS, TT);
		if (f1 != need) return false;
		outFlow = base;                     // whatever had to flow is the minimum
		return true;
	}
	// actual flow on the i-th user edge = lower bound + used residual
	ll flowOnEdge(int i) { return low[i] + d.edges[edgeId[i] ^ 1].cap; }
};

int32_t main() {
	// simple feasibility: 0->1 must carry between 2 and 4, 1->2 between 1 and 3,
	// plus a return path so a circulation can exist.
	BoundedFlow bf; bf.init(3);
	bf.ae(0, 1, 2, 4);
	bf.ae(1, 2, 1, 3);
	bf.ae(2, 0, 0, 5);
	cout << "feasible circulation? " << (bf.feasibleCirculation() ? "yes" : "no") << "\n";

	BoundedFlow bf2; bf2.init(3);
	int e0 = bf2.ae(0, 1, 1, 3);
	int e1 = bf2.ae(1, 2, 1, 3);
	ll f;
	bool ok = bf2.maxFlowWithLowerBounds(0, 2, f);
	cout << "max s->t flow with lower bounds = " << (ok ? f : -1) << "\n";
	if (ok) cout << "  flow on 0->1 = " << bf2.flowOnEdge(e0)
	             << ", 1->2 = " << bf2.flowOnEdge(e1) << "\n";
	return 0;
}

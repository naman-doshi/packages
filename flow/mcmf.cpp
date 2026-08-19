// ============================================================================
//  MIN-COST MAX-FLOW  (SPFA / Bellman-Ford shortest augmenting path)
//  O(F * V * E) worst case; fast in practice. Handles NEGATIVE edge costs.
// ----------------------------------------------------------------------------
//  Pushes the maximum possible flow from s to t and, among all max flows,
//  achieves the minimum total cost. Each augmentation follows a cheapest
//  s->t path in the residual graph, so partial results are also optimal:
//    - want max-flow-min-cost   -> read {flow, cost} from mcmf(s,t).
//    - want min cost of exactly K units -> stop once flow == K (see note).
//    - want min cost ignoring flow amount (edges may have negative cost and
//      you only send flow while it decreases cost) -> break when dist[t] >= 0.
//
//  API  (nodes 0-based)
//    MCMF g; g.init(N);
//    g.ae(u, v, cap, cost);          // directed edge, cost is PER UNIT of flow
//    auto [flow, cost] = g.mcmf(s, t);
//
//  NOTES
//    * SPFA tolerates negative-cost edges but NOT negative-cost cycles.
//    * for a flow cap of exactly K: cap the source out-edges at K, or clamp
//      the pushed amount so total flow never exceeds K.
//    * assignment problem (min-cost perfect matching) is a special case; for
//      dense n<=~500 the Hungarian algorithm is faster (see hungarian.cpp).
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long
#define ll long long

const ll INF = 9e18;

struct MCMF {
	struct Edge { int to; ll cap, cost; };
	int N;
	vector<Edge> edges;                // forward id even, backward id^1 odd
	vector<vector<int>> g;
	vector<ll> dist;
	vector<int> pe;                    // pe[v] = edge id used to reach v
	vector<char> inq;
	void init(int n) { N = n; g.assign(n, {}); edges.clear(); }
	void ae(int u, int v, ll cap, ll cost) {
		g[u].push_back(edges.size()); edges.push_back({v, cap, cost});
		g[v].push_back(edges.size()); edges.push_back({u, 0, -cost});
	}
	bool spfa(int s, int t) {
		dist.assign(N, INF); pe.assign(N, -1); inq.assign(N, 0);
		dist[s] = 0; deque<int> q; q.push_back(s);
		while (!q.empty()) {
			int u = q.front(); q.pop_front(); inq[u] = 0;
			for (int id : g[u]) {
				auto &e = edges[id];
				if (e.cap > 0 && dist[u] + e.cost < dist[e.to]) {
					dist[e.to] = dist[u] + e.cost;
					pe[e.to] = id;
					if (!inq[e.to]) {
						inq[e.to] = 1;
						// SLF: cheaper nodes to the front for a smaller constant
						if (!q.empty() && dist[e.to] < dist[q.front()]) q.push_front(e.to);
						else q.push_back(e.to);
					}
				}
			}
		}
		return dist[t] < INF;
	}
	// returns {max flow, min cost of that flow}
	pair<ll,ll> mcmf(int s, int t) {
		ll flow = 0, cost = 0;
		while (spfa(s, t)) {
			// bottleneck along the found shortest path
			ll f = INF;
			for (int v = t; v != s; v = edges[pe[v] ^ 1].to)
				f = min(f, edges[pe[v]].cap);
			for (int v = t; v != s; v = edges[pe[v] ^ 1].to) {
				edges[pe[v]].cap -= f;
				edges[pe[v] ^ 1].cap += f;
			}
			flow += f;
			cost += f * dist[t];
		}
		return {flow, cost};
	}
};

int32_t main() {
	// source 0, sink 3. Two paths of different cost.
	MCMF g; g.init(4);
	g.ae(0, 1, 2, 1);   // 0->1 cap 2 cost 1
	g.ae(0, 2, 2, 3);   // 0->2 cap 2 cost 3
	g.ae(1, 3, 2, 1);   // 1->3 cap 2 cost 1
	g.ae(2, 3, 2, 1);   // 2->3 cap 2 cost 1
	g.ae(1, 2, 1, 1);   // 1->2 cap 1 cost 1

	auto [flow, cost] = g.mcmf(0, 3);
	cout << "max flow  = " << flow << "\n";
	cout << "min cost  = " << cost << "\n";
	return 0;
}

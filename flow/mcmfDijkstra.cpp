// ============================================================================
//  MIN-COST MAX-FLOW, DIJKSTRA + JOHNSON POTENTIALS   O(F * E log V)
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Same answer and the same augmenting scheme as mcmf.cpp -- successive
//    shortest paths -- differing only in how each shortest path is FOUND.
//
//    *** F IS THE SAME IN BOTH. *** Both push the bottleneck along a cheapest
//    augmenting path, so both perform the same sequence of augmentations; this
//    file does NOT reduce how many there are. The whole win is per path:
//        mcmf.cpp    SPFA / Bellman-Ford   O(V*E)      per augmentation
//        this file   Dijkstra + heap       O(E log V)  per augmentation
//    a factor of about V/log V, WHATEVER F happens to be. So the rule is
//    "prefer this once the graph is big enough for the per-path cost to
//    matter", not "prefer this when there is a lot of flow".
//
//  WHAT NEITHER OF THEM FIXES
//    F is bounded only by the FLOW VALUE -- with integer capacities each
//    augmentation moves at least 1 unit -- so BOTH are pseudo-polynomial.
//    Capacities of 1e9 can mean 1e9 augmentations, and F * E log V is just as
//    dead as F * V * E. The fix for that is a model with fewer units, an
//    explicit cap via mcmf(s, t, K), or capacity scaling -- not this file.
//    Where the per-path saving is the entire game is when F is naturally
//    small: unit capacities (F <= V), assignment problems (F = n), or a given
//    K. That is the common case, which is why this is the better default.
//
//  A PRACTICAL CAVEAT  SPFA's O(V*E) is a worst case it rarely reaches; on
//    ordinary graphs it behaves closer to O(E) and the two are neck and neck.
//    The gap opens on inputs built to punish SPFA -- which you cannot spot
//    from the statement, which is the real argument for defaulting to this.
//
//  THE TRICK (Johnson)  Give every node a potential h[v] and use the REDUCED
//    cost  w'(u,v) = w(u,v) + h[u] - h[v].  Along any path the h terms
//    telescope, so shortest paths are unchanged; and if h is itself a shortest
//    distance then w' >= 0 everywhere, which is what lets Dijkstra run. After
//    each augmentation add the new distances into h to keep that true --
//    residual (reverse) edges get reduced cost exactly 0, so they stay legal.
//
//  API  (identical to mcmf.cpp, plus a flow cap)
//    MCMFD g; g.init(N);
//    g.ae(u, v, cap, cost);              // directed, cost is PER UNIT
//    auto [flow, cost] = g.mcmf(s, t);         // max flow, cheapest such
//    auto [flow, cost] = g.mcmf(s, t, K);      // send AT MOST K units
//
//  NEGATIVE COSTS are allowed: init() runs one Bellman-Ford pass to seed the
//    potentials. Negative COST CYCLES are not (no algorithm here handles them).
//    If you know every cost is >= 0 you can delete that pass.
//
//  COMPLEXITY  O(F * E log V), F = number of augmentations. Each pushes a full
//    bottleneck, so in practice F sits well below the flow value -- but the
//    only guarantee is F <= flow value, which is the pseudo-polynomial catch
//    above. mcmf.cpp is O(F * V * E) with the SAME F.
//
//  PITFALLS
//    * ae() adds a reverse edge of capacity 0 and cost -cost. Do NOT add both
//      directions yourself for an undirected edge -- add two separate ae()
//      calls if you really mean an undirected edge with cost both ways.
//    * "Cheapest flow, not necessarily maximum" (when some costs are negative
//      and pushing more makes it worse): stop as soon as the path cost turns
//      positive -- see the commented line in mcmf().
//    * Potentials of nodes unreachable in the residual graph are deliberately
//      left alone. That is safe: once a node is unreachable it stays that way,
//      because augmenting only ever adds reverse edges inside the reachable
//      set. Do not "fix" it by clamping them.
//    * Costs are ll and get multiplied by the pushed amount; with 1e9 caps and
//      1e9 costs you are at 1e18. Watch it.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long
#define ll long long

struct MCMFD {
	struct Edge { int to; ll cap, cost; };
	int N;
	vector<Edge> edges;              // paired: forward id even, backward id^1
	vector<vector<int>> g;
	vector<ll> h, dist;
	vector<int> pv;                  // pv[v] = id of the edge used to reach v
	static constexpr ll BIG = (ll)4e18;   // constexpr: usable as a default arg

	void init(int n) { N = n; g.assign(n, {}); edges.clear(); }
	void ae(int u, int v, ll cap, ll cost) {
		g[u].push_back(edges.size()); edges.push_back({v, cap, cost});
		g[v].push_back(edges.size()); edges.push_back({u, 0, -cost});
	}

	// {flow pushed, total cost}. Sends at most maxf units.
	pair<ll, ll> mcmf(int s, int t, ll maxf = BIG) {
		ll flow = 0, cost = 0;
		// seed potentials with Bellman-Ford from a virtual source joined to
		// every node at cost 0 -- i.e. start them all at 0 and relax. Delete
		// this loop if every cost is already >= 0.
		h.assign(N, 0);
		for (int pass = 0; pass < N; pass++) {
			bool upd = false;
			for (int id = 0; id < (int)edges.size(); id++) {
				if (edges[id].cap <= 0) continue;
				int u = edges[id ^ 1].to, v = edges[id].to;
				if (h[u] + edges[id].cost < h[v]) { h[v] = h[u] + edges[id].cost; upd = true; }
			}
			if (!upd) break;
		}

		while (flow < maxf) {
			dist.assign(N, BIG);
			pv.assign(N, -1);
			priority_queue<pair<ll, int>, vector<pair<ll, int>>, greater<>> pq;
			dist[s] = 0; pq.push({0, s});
			while (!pq.empty()) {
				auto [d, u] = pq.top(); pq.pop();
				if (d > dist[u]) continue;
				for (int id : g[u]) {
					auto &e = edges[id];
					if (e.cap <= 0) continue;
					ll nd = d + e.cost + h[u] - h[e.to];      // reduced, >= 0
					if (nd < dist[e.to]) {
						dist[e.to] = nd; pv[e.to] = id;
						pq.push({nd, e.to});
					}
				}
			}
			if (dist[t] >= BIG) break;                        // t unreachable

			// true (unreduced) cost of this path, before h changes under us
			ll pathCost = dist[t] + h[t] - h[s];
			// for "cheapest, not maximum" flow, stop once paths stop helping:
			// if (pathCost > 0) break;

			for (int v = 0; v < N; v++) if (dist[v] < BIG) h[v] += dist[v];

			ll push = maxf - flow;
			for (int v = t; v != s; v = edges[pv[v] ^ 1].to)
				push = min(push, edges[pv[v]].cap);
			for (int v = t; v != s; v = edges[pv[v] ^ 1].to) {
				edges[pv[v]].cap -= push;
				edges[pv[v] ^ 1].cap += push;
			}
			flow += push;
			cost += push * pathCost;
		}
		return {flow, cost};
	}
};

// ---------------------------------------------------------------------------
// reference implementation for the stress test: plain Bellman-Ford augmenting.
// Deliberately a DIFFERENT algorithm, so agreement means something.
// ---------------------------------------------------------------------------
struct MCMFref {
	struct Edge { int to; ll cap, cost; };
	int N; vector<Edge> edges; vector<vector<int>> g;
	void init(int n) { N = n; g.assign(n, {}); edges.clear(); }
	void ae(int u, int v, ll cap, ll cost) {
		g[u].push_back(edges.size()); edges.push_back({v, cap, cost});
		g[v].push_back(edges.size()); edges.push_back({u, 0, -cost});
	}
	pair<ll, ll> mcmf(int s, int t, ll maxf = (ll)4e18) {
		ll flow = 0, cost = 0;
		while (flow < maxf) {
			vector<ll> d(N, (ll)4e18); vector<int> pv(N, -1);
			d[s] = 0;
			for (int pass = 0; pass < N; pass++)
				for (int id = 0; id < (int)edges.size(); id++) {
					if (edges[id].cap <= 0) continue;
					int u = edges[id ^ 1].to, v = edges[id].to;
					if (d[u] < (ll)4e18 && d[u] + edges[id].cost < d[v]) {
						d[v] = d[u] + edges[id].cost; pv[v] = id;
					}
				}
			if (d[t] >= (ll)4e18) break;
			ll push = maxf - flow;
			for (int v = t; v != s; v = edges[pv[v] ^ 1].to) push = min(push, edges[pv[v]].cap);
			for (int v = t; v != s; v = edges[pv[v] ^ 1].to) {
				edges[pv[v]].cap -= push; edges[pv[v] ^ 1].cap += push;
			}
			flow += push; cost += push * d[t];
		}
		return {flow, cost};
	}
};

int32_t main() {
	// transportation: 2 plants -> 2 shops, per-unit shipping costs
	{
		MCMFD g; g.init(6);
		int s = 4, t = 5;
		g.ae(s, 0, 3, 0); g.ae(s, 1, 4, 0);      // supplies
		g.ae(2, t, 4, 0); g.ae(3, t, 3, 0);      // demands
		g.ae(0, 2, 5, 1); g.ae(0, 3, 5, 3);
		g.ae(1, 2, 5, 4); g.ae(1, 3, 5, 2);
		auto [f, c] = g.mcmf(s, t);
		// with x units on the cheap 0->2 route the total is 25 - 4x, and plant
		// 0 caps x at 3, so the optimum is 7 units at cost 13
		cout << "transportation: flow " << f << ", cost " << c << " (expect 7, 13)\n";
	}

	// --- randomised check against Bellman-Ford augmenting -------------------
	mt19937 rng(41);
	int bad = 0, trials = 400;
	for (int trial = 0; trial < trials; trial++) {
		int n = 3 + rng() % 6;
		bool allowNeg = trial % 3 == 0;          // exercise negative costs too
		MCMFD a; MCMFref b;
		a.init(n); b.init(n);
		for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) {
			if (u == v || rng() % 10 >= 4) continue;
			if (u > v && !allowNeg) continue;    // keep it a DAG unless testing neg
			ll cap = 1 + rng() % 9;
			ll cost = (ll)(rng() % 11) - (allowNeg && u < v ? 5 : 0);
			if (u > v) continue;                 // no back edges: no negative cycles
			a.ae(u, v, cap, cost); b.ae(u, v, cap, cost);
		}
		ll cap = trial % 4 == 0 ? (ll)(1 + rng() % 15) : (ll)4e18;
		auto ra = a.mcmf(0, n - 1, cap);
		auto rb = b.mcmf(0, n - 1, cap);
		if (ra != rb) bad++;
	}
	cout << trials - bad << "/" << trials
	     << " random networks match Bellman-Ford augmenting (flow and cost)\n";
	return 0;
}

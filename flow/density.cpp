// ============================================================================
//  MAXIMUM DENSITY SUBGRAPH  --  Goldberg: binary search + max closure
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Pick a non-empty set of vertices S maximising
//          density(S) = (edges with BOTH ends in S) / |S|
//    "Find the most tightly-knit group", "the densest community", "pick a
//    subgroup where everyone knows as many others as possible". Note this is
//    NOT max clique (which is NP-hard) -- density is an average, and that
//    average is exactly what makes it polynomial.
//
//  THE REDUCTION
//    Fix a guess lambda and ask: is there an S with |E_S| - lambda*|S| > 0?
//    That is a MAX CLOSURE / project selection (modeling.cpp #9):
//        s -> (edge e)   cap 1        each edge is a "project" worth 1
//        (edge e) -> u   cap INF      taking an edge forces both endpoints
//        (edge e) -> v   cap INF
//        (vertex v) -> t cap lambda   each vertex "costs" lambda
//        best = (#edges) - mincut,  and the answer set is the source side.
//    density* > lambda  <=>  best > 0, so binary search on lambda.
//
//    Everything here is INTEGER: lambda is carried as N/Q with a fixed
//    Q = 2*n*n, so the caps above become Q per edge and N per vertex. Q is
//    chosen so the search separates any two distinct densities -- two of them
//    differ by at least 1/n^2 (denominators are at most n), and the bracket is
//    1/Q = 1/(2n^2), which is finer. So the LAST successful lambda already
//    hands back a maximum-density set, and the exact answer is then just read
//    off that set as a fraction. No floating point anywhere.
//
//  API
//    MaxDensity md; md.init(n);
//    md.ae(u, v);                      // undirected edge
//    auto r = md.solve();
//    r.num / r.den                     // the density, as an exact fraction
//    r.verts                           // a subgraph attaining it
//
//  COMPLEXITY  O(log(n^2 m)) max flows -- about 30-35 of them -- each on
//    (n + m + 2) nodes. Comfortable to n, m in the low thousands; it is the
//    edge-nodes that dominate the graph size.
//
//  PITFALLS
//    * UNDIRECTED, and every edge counts once. Call ae(u,v) once per edge, not
//      once per direction. Parallel edges are allowed and do count twice.
//    * Self loops: a loop at v contributes to |E_S| whenever v is in S, which
//      is usually NOT what "density" is meant to mean. Filter them out first.
//    * m == 0 gives density 0 and a single arbitrary vertex.
//    * Ties: one maximum-density set comes back. There can be many, and the
//      SMALLEST such set is a different (also solvable) question -- take the
//      source side of the cut at the final lambda, which is what you get.
//    * The whole graph is often the answer; that is not a bug.
//    * Weighted edges: give the s->e cap the edge's weight * Q instead of Q,
//      and raise Q if weights are not integers.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long
#define ll long long

struct Dinic {
	struct Edge { int to; ll cap; };
	int N;
	vector<Edge> edges;
	vector<vector<int>> g;
	vector<int> level, it;
	void init(int n) { N = n; g.assign(n, {}); edges.clear(); }
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
				if (e.cap > 0 && level[e.to] < 0) { level[e.to] = level[u] + 1; q.push(e.to); }
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
			for (ll f; (f = dfs(s, t, (ll)9e18)) > 0; ) flow += f;
		}
		return flow;
	}
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
};

struct MaxDensity {
	int n;
	vector<pair<int, int>> es;
	void init(int n_) { n = n_; es.clear(); }
	void ae(int u, int v) { es.push_back({u, v}); }

	struct Result { ll num, den; vector<int> verts; };   // density = num/den

	// vertices of the best S for lambda = N/Q, or {} if none beats it
	vector<int> best(ll N, ll Q) {
		int m = (int)es.size();
		int s = n + m, t = n + m + 1;
		ll INF = Q * m + N * n + 1;
		Dinic d; d.init(n + m + 2);
		for (int i = 0; i < m; i++) {
			d.ae(s, n + i, Q);                       // the edge is worth Q
			d.ae(n + i, es[i].first, INF);           // ... but forces both ends
			d.ae(n + i, es[i].second, INF);
		}
		for (int v = 0; v < n; v++) d.ae(v, t, N);   // each vertex costs N
		ll profit = Q * m - d.maxflow(s, t);
		if (profit <= 0) return {};
		auto side = d.minCutSide(s);
		vector<int> S;
		for (int v = 0; v < n; v++) if (side[v]) S.push_back(v);
		return S;
	}

	Result solve() {
		int m = (int)es.size();
		if (n == 0) return {0, 1, {}};
		if (m == 0) return {0, 1, {0}};
		ll Q = 2LL * n * n;
		// largest N with "some S has density > N/Q"; N = 0 always works (m >= 1)
		ll lo = 0, hi = Q * m;
		while (lo < hi) {
			ll mid = lo + (hi - lo + 1) / 2;
			if (!best(mid, Q).empty()) lo = mid; else hi = mid - 1;
		}
		vector<int> S = best(lo, Q);
		// read the exact density straight off the set
		vector<char> in(n, 0);
		for (int v : S) in[v] = 1;
		ll cnt = 0;
		for (auto &e : es) if (in[e.first] && in[e.second]) cnt++;
		ll g = gcd(cnt, (ll)S.size());
		if (g == 0) g = 1;
		return {cnt / g, (ll)S.size() / g, S};
	}
};

int32_t main() {
	// a 4-clique hanging off a long path: the clique (density 6/4 = 3/2) wins
	{
		MaxDensity md; md.init(8);
		md.ae(0,1); md.ae(0,2); md.ae(0,3); md.ae(1,2); md.ae(1,3); md.ae(2,3);
		md.ae(3,4); md.ae(4,5); md.ae(5,6); md.ae(6,7);
		auto r = md.solve();
		cout << "clique + tail: density " << r.num << "/" << r.den << " (expect 3/2), verts {";
		for (int v : r.verts) cout << " " << v;
		cout << " }\n";
	}
	// a plain triangle: 3 edges over 3 vertices
	{
		MaxDensity md; md.init(3);
		md.ae(0,1); md.ae(1,2); md.ae(0,2);
		auto r = md.solve();
		cout << "triangle: density " << r.num << "/" << r.den << " (expect 1/1)\n";
	}

	// --- randomised check against brute force ------------------------------
	mt19937 rng(31);
	int bad = 0, trials = 300;
	for (int trial = 0; trial < trials; trial++) {
		int n = 1 + rng() % 7;
		MaxDensity md; md.init(n);
		vector<pair<int, int>> es;
		for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++)
			if (rng() % 10 < 5) { es.push_back({u, v}); md.ae(u, v); }
		auto r = md.solve();

		// brute force: best |E_S|/|S| over all non-empty S, compared exactly
		ll bn = 0, bd = 1;
		for (int mask = 1; mask < (1 << n); mask++) {
			ll c = 0, sz = __builtin_popcount(mask);
			for (auto &e : es) if ((mask >> e.first & 1) && (mask >> e.second & 1)) c++;
			if (c * bd > bn * sz) { bn = c; bd = sz; }
		}
		// the reported fraction must equal the optimum, and the returned set
		// must actually achieve it
		vector<char> in(n, 0);
		for (int v : r.verts) in[v] = 1;
		ll c = 0, sz = (ll)r.verts.size();
		for (auto &e : es) if (in[e.first] && in[e.second]) c++;
		if (r.num * bd != bn * r.den || sz == 0 || c * bd != bn * sz) bad++;
	}
	cout << trials - bad << "/" << trials
	     << " random graphs match brute force (exact fraction and set)\n";
	return 0;
}

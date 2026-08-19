// ============================================================================
//  WEIGHTED KONIG  --  min-weight vertex cover / max-weight independent set
//  on a BIPARTITE graph, by min cut.   REQUIRES the Dinic below.
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Each vertex of a bipartite graph carries a NON-NEGATIVE weight. An edge
//    means "these two conflict". Then
//      min-weight vertex cover      = min cut of the network below
//      max-weight independent set   = (total weight) - (min-weight cover)
//    "Choose a maximum-weight set of items with no conflicting pair, where the
//    conflicts only ever run between two groups" is the phrasing to watch for.
//
//  ***  THE TRAP  ***
//    The UNWEIGHTED Konig identity (min cover = MAX MATCHING, modeling.cpp #6)
//    DOES NOT GENERALISE. With weights you must run a min cut; a maximum
//    matching answers a different question and is silently wrong. Setting every
//    weight to 1 here reproduces the unweighted answer, so this file supersedes
//    that identity whenever weights appear.
//
//  THE NETWORK   (left u, right v, an edge u--v means they conflict)
//      s -> u        cap  wl[u]
//      u -> v        cap  INF        for every conflict edge  (never cut)
//      v -> t        cap  wr[v]
//    A finite cut cannot contain an INF edge, so for every conflict it must cut
//    s->u or v->t -- i.e. the cut vertices form a vertex COVER, and any cover
//    gives a cut of the same weight. Minimum cut = minimum cover.
//    Recovering it: let R = vertices reachable from s in the residual graph.
//      cover  = { left u : u NOT in R }  union  { right v : v in R }
//      indep  = the complement of that = { left u in R } u { right v not in R }
//
//  API   (left 0..nL-1 and right 0..nR-1 are SEPARATE index spaces)
//    WeightedCover wc; wc.init(nL, nR);
//    wc.setL(u, w);  wc.setR(v, w);     // non-negative weights, default 0
//    wc.ae(u, v);                       // u and v conflict
//    auto r = wc.solve();
//    r.cover / r.indep                  // the two totals
//    r.coverL, r.coverR                 // vertices forming a minimum cover
//    r.indepL, r.indepR                 // vertices forming a max-weight set
//
//  COMPLEXITY  one max flow: O(V^2 E) worst case, in practice far less.
//
//  PITFALLS
//    * BIPARTITE ONLY. On a general graph max-weight independent set is
//      NP-hard -- there is no flow model, so check the bipartition first.
//    * Weights must be >= 0. A negative weight means "never take it": drop the
//      vertex (and its edges) before building.
//    * INF here is (total weight + 1), not 9e18, so the augmenting sums cannot
//      overflow. Do not lower it.
//    * Zero-weight vertices land in the cover (free) rather than the
//      independent set. The TOTALS are still optimal; if you need such a
//      vertex in the set, add it afterwards if it has no conflict inside.
//    * Ties: any one optimal cover/set comes back, not all of them.
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

struct WeightedCover {
	int nL, nR;
	vector<ll> wl, wr;
	vector<pair<int, int>> es;

	void init(int L, int R) {
		nL = L; nR = R;
		wl.assign(L, 0); wr.assign(R, 0); es.clear();
	}
	void setL(int u, ll w) { wl[u] = w; }
	void setR(int v, ll w) { wr[v] = w; }
	void ae(int u, int v) { es.push_back({u, v}); }

	struct Result {
		ll cover, indep;
		vector<int> coverL, coverR, indepL, indepR;
	};

	Result solve() {
		ll total = 0;
		for (ll x : wl) total += x;
		for (ll x : wr) total += x;
		ll INF = total + 1;                 // big enough to never be cut, small
		                                    // enough that sums cannot overflow
		int s = nL + nR, t = nL + nR + 1;
		Dinic d; d.init(nL + nR + 2);
		for (int u = 0; u < nL; u++) if (wl[u] > 0) d.ae(s, u, wl[u]);
		for (int v = 0; v < nR; v++) if (wr[v] > 0) d.ae(nL + v, t, wr[v]);
		for (auto &e : es) d.ae(e.first, nL + e.second, INF);

		Result r;
		r.cover = d.maxflow(s, t);
		r.indep = total - r.cover;
		auto side = d.minCutSide(s);
		for (int u = 0; u < nL; u++) (side[u] ? r.indepL : r.coverL).push_back(u);
		for (int v = 0; v < nR; v++) (side[nL + v] ? r.coverR : r.indepR).push_back(v);
		return r;
	}
};

int32_t main() {
	// tiny worked instance: left {a=5, b=1}, right {x=4, y=2}, conflicts a-x, b-x
	{
		WeightedCover wc; wc.init(2, 2);
		wc.setL(0, 5); wc.setL(1, 1); wc.setR(0, 4); wc.setR(1, 2);
		wc.ae(0, 0); wc.ae(1, 0);
		auto r = wc.solve();
		cout << "cover " << r.cover << ", independent set " << r.indep
		     << " (total 12); set = L{";
		for (int u : r.indepL) cout << " " << u;
		cout << " } R{";
		for (int v : r.indepR) cout << " " << v;
		cout << " }\n";
	}

	// --- randomised check against brute force ------------------------------
	mt19937 rng(23);
	int bad = 0, trials = 500;
	for (int trial = 0; trial < trials; trial++) {
		int nL = 1 + rng() % 4, nR = 1 + rng() % 4;
		vector<ll> wl(nL), wr(nR);
		WeightedCover wc; wc.init(nL, nR);
		for (int u = 0; u < nL; u++) { wl[u] = rng() % 10; wc.setL(u, wl[u]); }
		for (int v = 0; v < nR; v++) { wr[v] = rng() % 10; wc.setR(v, wr[v]); }
		vector<pair<int, int>> es;
		for (int u = 0; u < nL; u++) for (int v = 0; v < nR; v++)
			if (rng() % 10 < 5) { es.push_back({u, v}); wc.ae(u, v); }
		auto r = wc.solve();

		// brute force: best independent set over all subsets of both sides
		ll best = 0;
		for (int a = 0; a < (1 << nL); a++) for (int b = 0; b < (1 << nR); b++) {
			bool okSet = true;
			for (auto &e : es)
				if ((a >> e.first & 1) && (b >> e.second & 1)) { okSet = false; break; }
			if (!okSet) continue;
			ll c = 0;
			for (int u = 0; u < nL; u++) if (a >> u & 1) c += wl[u];
			for (int v = 0; v < nR; v++) if (b >> v & 1) c += wr[v];
			best = max(best, c);
		}

		// value optimal, AND the returned set really is independent with that weight
		ll got = 0;
		for (int u : r.indepL) got += wl[u];
		for (int v : r.indepR) got += wr[v];
		bool indep = true;
		for (auto &e : es) {
			bool inL = find(r.indepL.begin(), r.indepL.end(), e.first) != r.indepL.end();
			bool inR = find(r.indepR.begin(), r.indepR.end(), e.second) != r.indepR.end();
			if (inL && inR) indep = false;
		}
		// and the cover really covers every edge
		bool covers = true;
		for (auto &e : es) {
			bool inL = find(r.coverL.begin(), r.coverL.end(), e.first) != r.coverL.end();
			bool inR = find(r.coverR.begin(), r.coverR.end(), e.second) != r.coverR.end();
			if (!inL && !inR) covers = false;
		}
		if (r.indep != best || got != best || !indep || !covers) bad++;
	}
	cout << trials - bad << "/" << trials
	     << " random bipartite instances match brute force (value, set, cover)\n";
	return 0;
}

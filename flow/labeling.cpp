// ============================================================================
//  BINARY LABELING BY MIN CUT   (graph cuts / the Kolmogorov-Zabih reduction)
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Every node picks one of two labels, A or B. You pay
//      unary     costA[u] / costB[u]        for u's own choice
//      pairwise  a full 2x2 table per pair, E_XY = cost when u picks X, w picks Y
//    and you want the assignment of MINIMUM total cost. That is exactly a
//    minimum s-t cut, with source side = A and sink side = B.
//
//    Reach for this when a problem reads "each thing is one of two kinds, there
//    is a cost for each thing on its own, AND a cost when two related things
//    disagree". Project selection (modeling.cpp #9) is the special case where
//    the pairwise table is 0 / INF; this handles an arbitrary table.
//
//  THE ONE CONDITION -- SUBMODULARITY, checked per pair:
//        E_AB + E_BA >= E_AA + E_BB
//    In words: disagreeing must cost at least as much as agreeing. If it fails,
//    minimising is NP-hard in general and NO min-cut construction exists, so
//    addPairwise returns false instead of quietly handing back a wrong answer.
//    (A pure "penalty p when they differ" table, 0/p/p/0, always satisfies it.)
//
//  API   (nodes 0-based; the source and sink are added internally)
//    BinaryLabeling bl; bl.init(n);
//    bl.addUnary(u, costA, costB);                  // u's own cost
//    bl.addPairwise(u, w, E_AA, E_AB, E_BA, E_BB);  // false => not submodular
//    auto [cost, lab] = bl.solve();                 // lab[u]: 1 = A, 0 = B
//
//  HOW THE 2x2 TABLE BECOMES EDGES
//    Write x_u = 0 for "u is on the source side (A)", 1 for the sink side (B).
//    The table splits into a constant, two terms that depend on ONE node, and
//    one term that fires only on the ordered pair (A at u, B at w):
//        E = E_AA
//          + (E_BA - E_AA) * [x_u = 1]
//          + (E_BB - E_BA) * [x_w = 1]
//          + (E_AB + E_BA - E_AA - E_BB) * [x_u = 0][x_w = 1]
//    Check it against all four corners -- each reproduces its table entry. The
//    [x = 1] terms are terminal edges s->node; the last term is the edge u->w,
//    and its coefficient is exactly the submodularity slack, hence >= 0.
//
//  COMPLEXITY  one max flow on n+2 nodes with (unary edges + one edge per
//    pair). Dinic is O(V^2 E) in theory and far better here in practice.
//
//  PITFALLS
//    * Capacities are ll, i.e. INTEGER costs. If the costs are fractional,
//      scale them all to integers first -- do NOT swap in a double Dinic, the
//      residual comparisons stop being reliable.
//    * The pair is ORDERED. addPairwise(u, w, ...) reads E_AB as "u picks A and
//      w picks B". Swapping u and w means swapping E_AB with E_BA.
//    * addPairwise returning false leaves the instance UNCHANGED, so check it.
//      Submodularity is a property of each table on its own, not of the graph.
//    * Terminal capacities can go negative mid-build; solve() repairs that (see
//      the comment there) so it is not a bug -- but it does mean you must not
//      read toS / toT yourself before calling solve().
//    * Ties are broken arbitrarily: any minimum-cost labelling may come back.
//      The COST is unique, the labelling need not be.
//    * Maximising instead? Negate every cost, solve, negate the answer -- but
//      the submodularity test flips with it, so check it on the negated table.
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
};

// ----------------------------------------------------------------- labeling
struct BinaryLabeling {
	int n;
	vector<ll> toS, toT;                // capacities of s->u and u->t
	vector<array<ll, 3>> pairEdges;     // {u, w, cap}: cut iff u picks A, w picks B
	ll constant;
	bool submodular;                    // false once any table was rejected

	void init(int n_) {
		n = n_;
		toS.assign(n, 0); toT.assign(n, 0);
		pairEdges.clear(); constant = 0; submodular = true;
	}

	// u pays costA if it ends up labelled A, costB if it ends up B
	void addUnary(int u, ll costA, ll costB) {
		toT[u] += costA;                // u->t is cut exactly when u is on s's side (A)
		toS[u] += costB;                // s->u is cut exactly when u is on t's side (B)
	}

	// the full 2x2 table for the ORDERED pair (u, w). Returns false and changes
	// nothing if the table is not submodular.
	bool addPairwise(int u, int w, ll E_AA, ll E_AB, ll E_BA, ll E_BB) {
		ll rem = E_AB + E_BA - E_AA - E_BB;      // the submodularity slack
		if (rem < 0) { submodular = false; return false; }
		constant += E_AA;
		toS[u] += E_BA - E_AA;          // these two may go NEGATIVE for now;
		toS[w] += E_BB - E_BA;          // solve() shifts them back up
		if (rem > 0) pairEdges.push_back({(ll)u, (ll)w, rem});
		return true;
	}

	// {minimum total cost, labels} with lab[u] = 1 for A, 0 for B
	pair<ll, vector<char>> solve() {
		// A negative capacity is not an edge. But adding c to BOTH terminal
		// edges of one node raises EVERY cut by exactly c -- every cut severs
		// exactly one of the two, since u is on one side or the other -- so it
		// is a pure constant. Shift up until non-negative, subtract it back out.
		for (int u = 0; u < n; u++)
			for (int rep = 0; rep < 2; rep++) {
				ll &cap = rep ? toT[u] : toS[u];
				if (cap < 0) { ll sh = -cap; toS[u] += sh; toT[u] += sh; constant -= sh; }
			}

		int s = n, t = n + 1;
		Dinic d; d.init(n + 2);
		for (int u = 0; u < n; u++) {
			if (toS[u] > 0) d.ae(s, u, toS[u]);
			if (toT[u] > 0) d.ae(u, t, toT[u]);
		}
		for (auto &e : pairEdges) d.ae(e[0], e[1], e[2]);

		ll cut = d.maxflow(s, t);
		auto side = d.minCutSide(s);
		vector<char> lab(n);
		for (int u = 0; u < n; u++) lab[u] = side[u];   // source side = A
		return {cut + constant, lab};
	}
};

int32_t main() {
	// --- a worked two-node instance, full 2x2 table -------------------------
	{
		BinaryLabeling bl; bl.init(2);
		bl.addUnary(0, 4, 3);                        // u: A costs 4, B costs 3
		bl.addUnary(1, 2, 5);                        // w: A costs 2, B costs 5
		bl.addPairwise(0, 1, 1, 6, 3, 2);            // AA AB BA BB
		auto [cost, lab] = bl.solve();
		cout << "two nodes: u=" << (lab[0] ? 'A' : 'B')
		     << " w=" << (lab[1] ? 'A' : 'B') << ", cost " << cost << "\n";
	}

	// --- three tasks on two machines, penalty only when neighbours differ ---
	{
		BinaryLabeling bl; bl.init(3);
		ll onA[3] = {2, 5, 5}, onB[3] = {7, 3, 3};
		for (int i = 0; i < 3; i++) bl.addUnary(i, onA[i], onB[i]);
		bl.addPairwise(0, 1, 0, 4, 4, 0);            // differ -> 4
		bl.addPairwise(1, 2, 0, 1, 1, 0);            // differ -> 1
		auto [cost, lab] = bl.solve();
		cout << "three tasks: ";
		for (int i = 0; i < 3; i++) cout << (lab[i] ? 'A' : 'B');
		cout << ", cost " << cost << "\n";
	}

	// --- a table that is NOT submodular is refused, not fudged --------------
	{
		BinaryLabeling bl; bl.init(2);
		bool taken = bl.addPairwise(0, 1, 5, 0, 0, 5);   // 0+0 < 5+5
		cout << "non-submodular table accepted? " << (taken ? "yes" : "no") << "\n";
	}

	// --- randomised check against brute force ------------------------------
	mt19937 rng(7);
	int bad = 0, trials = 500;
	for (int trial = 0; trial < trials; trial++) {
		int n = 2 + rng() % 5;
		vector<array<ll, 2>> un(n);
		vector<array<ll, 6>> pw;
		BinaryLabeling bl; bl.init(n);
		for (int u = 0; u < n; u++) {
			un[u] = {(ll)(rng() % 10), (ll)(rng() % 10)};
			bl.addUnary(u, un[u][0], un[u][1]);
		}
		for (int u = 0; u < n; u++) for (int w = u + 1; w < n; w++)
			if (rng() % 10 < 6) {
				ll AA = rng() % 6, BB = rng() % 6, slack = rng() % 7;
				ll AB = rng() % (AA + BB + slack + 1);   // force AB+BA = AA+BB+slack
				ll BA = AA + BB + slack - AB;
				pw.push_back({(ll)u, (ll)w, AA, AB, BA, BB});
				bl.addPairwise(u, w, AA, AB, BA, BB);
			}
		auto [cost, lab] = bl.solve();

		// cost of an explicit assignment (bit set = that node is labelled A)
		auto costOf = [&](int mask) {
			ll c = 0;
			for (int u = 0; u < n; u++) c += (mask >> u & 1) ? un[u][0] : un[u][1];
			for (auto &e : pw) {
				int xu = mask >> e[0] & 1, xw = mask >> e[1] & 1;
				c += xu ? (xw ? e[2] : e[3]) : (xw ? e[4] : e[5]);
			}
			return c;
		};
		ll best = LLONG_MAX;
		for (int m = 0; m < (1 << n); m++) best = min(best, costOf(m));

		// the VALUE must be optimal, and the returned LABELS must achieve it
		int got = 0;
		for (int u = 0; u < n; u++) if (lab[u]) got |= 1 << u;
		if (cost != best || costOf(got) != cost) bad++;
	}
	cout << trials - bad << "/" << trials
	     << " random instances match brute force (value and labels)\n";
	return 0;
}

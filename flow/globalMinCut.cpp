// ============================================================================
//  GLOBAL MIN CUT  --  STOER-WAGNER,  O(V^3),  UNDIRECTED weighted graph
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    The cheapest way to split the vertices into TWO non-empty groups. Unlike
//    dinic.cpp there is no source and no sink -- you are not told which side
//    anything belongs on, you are asking for the weakest place in the whole
//    graph. "Cut the graph in two as cheaply as possible", "minimum cost to
//    disconnect the network", "how robust is this network".
//
//    Do NOT reach for this when you DO have a fixed s and t -- that is an
//    ordinary max flow (dinic.cpp), and it is faster.
//
//  API
//    vector<vector<ll>> w(n, vector<ll>(n, 0));   // symmetric, w[i][i] unused
//    w[u][v] += c; w[v][u] += c;                  // undirected edge, weight c
//    auto [cut, side] = globalMinCut(w);
//    // cut  = total weight crossing the cheapest 2-way split
//    // side = the vertices on ONE side of it (the other side is the rest)
//
//  HOW IT WORKS  Repeat n-1 times: grow a set from an arbitrary vertex, always
//    absorbing the vertex most strongly attached to what you have so far (a
//    "maximum adjacency ordering"). The KEY LEMMA is that for the LAST two
//    vertices s, t added this way, {t} alone is a minimum s-t cut. So record
//    the weight of {t} vs the rest ("cut of the phase"), then MERGE s and t and
//    repeat: any global min cut either separates s from t -- and we just
//    measured the best such cut -- or keeps them together, in which case
//    merging them loses nothing. The smallest cut of a phase is the answer.
//
//  COMPLEXITY  O(V^3) with the dense array scan below; the graph is stored as
//    a matrix, so this wants n up to ~500, not 1e5. Weights must be
//    NON-NEGATIVE (the lemma fails otherwise).
//
//  PITFALLS
//    * UNDIRECTED only. Fill BOTH w[u][v] and w[v][u]. A directed global min
//      cut is a different (harder) problem.
//    * Parallel edges: add the weights (`+=`, as above). Self loops are
//      ignored, which is correct -- they never cross a cut.
//    * A DISCONNECTED graph has a global min cut of 0, and that is what you
//      get back, with `side` one of the components. If you wanted connectivity
//      of a connected graph, check for 0 first.
//    * n < 2 has no 2-way split at all; you get {0, {}}.
//    * Unweighted edge connectivity = run this with every weight 1.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long
#define ll long long

// {weight of the cheapest 2-way split, the vertices on one side of it}
pair<ll, vector<int>> globalMinCut(vector<vector<ll>> w) {
	int n = (int)w.size();
	if (n < 2) return {0, {}};

	vector<vector<int>> grp(n);              // grp[i] = original vertices merged into i
	for (int i = 0; i < n; i++) grp[i] = {i};
	vector<char> alive(n, 1);
	ll best = LLONG_MAX;
	vector<int> bestSide;

	for (int phase = 0; phase + 1 < n; phase++) {
		// maximum adjacency ordering over the vertices still alive
		vector<ll> wt(n, 0);
		vector<char> added(n, 0);
		int prev = -1, last = -1, cnt = 0;
		for (int i = 0; i < n; i++) if (alive[i]) cnt++;

		for (int it = 0; it < cnt; it++) {
			int sel = -1;                    // most strongly attached so far
			for (int i = 0; i < n; i++)
				if (alive[i] && !added[i] && (sel < 0 || wt[i] > wt[sel])) sel = i;
			added[sel] = 1;
			prev = last; last = sel;
			for (int i = 0; i < n; i++)
				if (alive[i] && !added[i]) wt[i] += w[sel][i];
		}

		// cut of the phase: {last} against everything else. wt[last] was frozen
		// when it was picked, and that is exactly its weight to the rest.
		if (wt[last] < best) { best = wt[last]; bestSide = grp[last]; }

		// merge `last` into `prev`
		grp[prev].insert(grp[prev].end(), grp[last].begin(), grp[last].end());
		for (int i = 0; i < n; i++) if (alive[i] && i != prev && i != last) {
			w[prev][i] += w[last][i];
			w[i][prev] = w[prev][i];
		}
		alive[last] = 0;
	}
	return {best, bestSide};
}

int32_t main() {
	// the classic Stoer-Wagner paper example (8 vertices, min cut 4)
	{
		int n = 8;
		vector<vector<ll>> w(n, vector<ll>(n, 0));
		auto ae = [&](int a, int b, ll c) { w[a][b] += c; w[b][a] += c; };
		ae(0,1,2); ae(0,4,3); ae(1,2,3); ae(1,4,2); ae(1,5,2); ae(2,3,4);
		ae(2,6,2); ae(3,6,2); ae(3,7,2); ae(4,5,3); ae(5,6,1); ae(6,7,3);
		auto [cut, side] = globalMinCut(w);
		cout << "paper example: min cut " << cut << " (expect 4), side {";
		for (int v : side) cout << " " << v;
		cout << " }\n";
	}

	// disconnected -> 0
	{
		vector<vector<ll>> w(4, vector<ll>(4, 0));
		w[0][1] = w[1][0] = 5; w[2][3] = w[3][2] = 7;
		cout << "disconnected: min cut " << globalMinCut(w).first << " (expect 0)\n";
	}

	// --- randomised check against brute force over all 2-way splits ---------
	mt19937 rng(11);
	int bad = 0, trials = 400;
	for (int trial = 0; trial < trials; trial++) {
		int n = 2 + rng() % 7;
		vector<vector<ll>> w(n, vector<ll>(n, 0));
		for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++)
			if (rng() % 10 < 6) { ll c = 1 + rng() % 9; w[i][j] = w[j][i] = c; }

		auto [cut, side] = globalMinCut(w);

		ll best = LLONG_MAX;
		for (int m = 1; m < (1 << n) - 1; m++) {      // every non-trivial split
			ll c = 0;
			for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++)
				if ((m >> i & 1) != (m >> j & 1)) c += w[i][j];
			best = min(best, c);
		}
		// the value must be optimal AND `side` must actually realise it
		int mask = 0;
		for (int v : side) mask |= 1 << v;
		ll got = 0;
		for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++)
			if ((mask >> i & 1) != (mask >> j & 1)) got += w[i][j];
		if (cut != best || got != cut || side.empty() || (int)side.size() == n) bad++;
	}
	cout << trials - bad << "/" << trials
	     << " random graphs match brute force (value and side)\n";
	return 0;
}

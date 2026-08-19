// ============================================================================
//  CLOSEST PAIR, MANHATTAN DISTANCE, MANHATTAN MST   -- REQUIRES point.cpp
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    "Nearest two of n points" and everything L1 (taxicab). The L1 <-> L_inf
//    rotation at the bottom is the single most useful trick for grid-distance
//    problems -- it turns "sum of |dx|+|dy|" into two independent 1-D problems.
//
//  API
//    closestPair(pts)        -> {i, j}, the two nearest points. O(n log n)
//    closestDist2(pts)       -> the squared distance (exact integer)
//    farthestPair(pts)       -> in convex.cpp (hull + rotating calipers)
//    toChebyshev(p)          (x+y, x-y): L1 distance becomes L_inf distance
//    toManhattan(p)          the inverse (halves are implicit -- see PITFALLS)
//    manhattanFarthest(pts)  largest L1 distance between any two points, O(n)
//    manhattanMSTEdges(pts)  the O(n) candidate edges of an L1 MST
//    manhattanMST(pts)       -> {total weight, edges}, O(n log n)
//
//  COMPLEXITY  closestPair O(n log n) sweep. manhattanMST O(n log n):
//    only 4n candidate edges survive, then Kruskal.
//
//  PITFALLS
//    * closestPair needs at least 2 points; duplicates give distance 0 and are
//      reported correctly.
//    * closestDist2 is EXACT (long long). Do not sqrt it until you print.
//      Coordinates up to 1e9 give a squared distance up to 8e18 -- that fits,
//      but only just.
//    * toChebyshev doubles the scale: L1(p,q) == L_inf(T(p),T(q)) exactly, but
//      going back divides by 2. Stay in Chebyshev space and convert the final
//      ANSWER, don't round intermediate points.
//    * manhattanMST returns the MST of the COMPLETE graph under L1 -- you do
//      not need to build n^2 edges, that is the whole point.
//    * The minimum edge of the Manhattan MST is the closest pair under L1.
//    * For a EUCLIDEAN MST use delaunay.cpp (the EMST is a subgraph of the
//      Delaunay triangulation).
// ============================================================================

// indices of the two closest points (Euclidean). O(n log n) sweep.
pair<int, int> closestPair(const vP &pts) {
	int n = (int)pts.size();
	if (n < 2) return {0, 0};
	vector<pair<P, int>> v(n);
	for (int i = 0; i < n; i++) v[i] = {pts[i], i};
	sort(v.begin(), v.end(), [](const pair<P, int> &a, const pair<P, int> &b) {
		return a.first.y != b.first.y ? a.first.y < b.first.y : a.first.x < b.first.x;
	});
	set<pair<P, int>> S;                                   // ordered by x, then y
	ll best = LLONG_MAX;
	pair<int, int> res{v[0].second, v[1].second};
	int j = 0;
	for (int i = 0; i < n; i++) {
		P p = v[i].first;
		ll d = best == LLONG_MAX ? LLONG_MAX / 4 : 1 + (ll)sqrtl((long double)best);
		while (j < i && v[j].first.y <= p.y - d) S.erase(v[j]), j++;
		auto lo = S.lower_bound({P(p.x - d, LLONG_MIN), INT_MIN});
		auto hi = S.upper_bound({P(p.x + d, LLONG_MAX), INT_MAX});
		for (auto it = lo; it != hi; ++it) {
			ll dd = dist2(it->first, p);
			if (dd < best) best = dd, res = {it->second, v[i].second};
		}
		S.insert(v[i]);
	}
	return res;
}
ll closestDist2(const vP &pts) {
	auto r = closestPair(pts);
	return dist2(pts[r.first], pts[r.second]);
}

// ---------------------------------------------------------------- L1 tricks
// L1 distance in the original space == L_inf distance in Chebyshev space
P toChebyshev(P p) { return P(p.x + p.y, p.x - p.y); }
P toManhattan(P p) { return P((p.x + p.y) / 2, (p.x - p.y) / 2); }  // exact iff even

// largest |dx| + |dy| over all pairs, in O(n): max over the 4 sign patterns of
// (max of s1*x + s2*y) - (min of s1*x + s2*y)
ll manhattanFarthest(const vP &p) {
	ll best = 0;
	for (int s = 0; s < 2; s++) {
		ll mn = LLONG_MAX, mx = LLONG_MIN;
		for (auto &q : p) {
			ll v = s ? q.x - q.y : q.x + q.y;
			mn = min(mn, v);
			mx = max(mx, v);
		}
		best = max(best, mx - mn);
	}
	return best;
}

// the O(4n) candidate edges of a Manhattan MST, as {weight, i, j}
vector<array<ll, 3>> manhattanMSTEdges(vP ps) {
	int n = (int)ps.size();
	vector<int> id(n);
	for (int i = 0; i < n; i++) id[i] = i;
	vector<array<ll, 3>> edges;
	for (int k = 0; k < 4; k++) {
		sort(id.begin(), id.end(), [&](int i, int j) {
			return ps[i].x + ps[i].y < ps[j].x + ps[j].y;
		});
		map<ll, int> sweep;                       // key = -y, value = point index
		for (int i : id) {
			for (auto it = sweep.lower_bound(-ps[i].y); it != sweep.end();
			     sweep.erase(it++)) {
				int j = it->second;
				P d = ps[i] - ps[j];
				if (d.y > d.x) break;
				edges.push_back({d.x + d.y, (ll)i, (ll)j});
			}
			sweep[-ps[i].y] = i;
		}
		for (auto &p : ps) {                      // rotate/reflect into the next octant
			if (k & 1) p.x = -p.x;
			else swap(p.x, p.y);
		}
	}
	return edges;
}

// {total weight, tree edges as {w, i, j}} of the L1 minimum spanning tree
pair<ll, vector<array<ll, 3>>> manhattanMST(const vP &ps) {
	int n = (int)ps.size();
	auto edges = manhattanMSTEdges(ps);
	sort(edges.begin(), edges.end());
	vector<int> par(n);
	for (int i = 0; i < n; i++) par[i] = i;
	function<int(int)> find = [&](int x) { return par[x] == x ? x : par[x] = find(par[x]); };
	ll total = 0;
	vector<array<ll, 3>> used;
	for (auto &e : edges) {
		int a = find((int)e[1]), b = find((int)e[2]);
		if (a == b) continue;
		par[a] = b;
		total += e[0];
		used.push_back(e);
	}
	return {total, used};
}

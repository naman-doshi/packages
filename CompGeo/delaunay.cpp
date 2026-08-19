// ============================================================================
//  DELAUNAY / VORONOI / EUCLIDEAN MST   -- REQUIRES point.cpp, circles.cpp
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    The Delaunay triangulation is the triangulation in which no point lies
//    inside the circumcircle of any triangle. Reach for it when you need:
//      * the EUCLIDEAN minimum spanning tree (it is a subgraph of Delaunay,
//        so you go from n^2 candidate edges down to O(n))
//      * "nearest other point for every point" (also a Delaunay edge)
//      * the Voronoi diagram (its exact dual)
//      * a triangulation that avoids skinny triangles
//
//    Bowyer-Watson incremental insertion, driven entirely by the EXACT integer
//    inCircleDet predicate from circles.cpp. No floating point anywhere in the
//    construction, so co-circular and collinear input is handled exactly rather
//    than by luck.
//
//  API
//    vector<array<int,3>> tri = delaunay(pts);        // index triples, CCW
//    vector<pair<int,int>> e  = delaunayEdges(pts);   // deduplicated, i < j
//    auto [w, edges] = euclideanMST(pts);             // w is a double
//    vPd vc = voronoiVertices(pts, tri);              // circumcentres
//    int j = nearestNeighborOf(pts, i, tri);          // via Delaunay edges
//
//  COMPLEXITY  O(n^2) worst case, and that is what you should budget for:
//    every insertion scans the current triangle list. Comfortable to n ~ 3000.
//    For an MST at larger n, check whether the problem is really Manhattan --
//    closestPair.cpp's manhattanMST is O(n log n).
//
//  PITFALLS
//    * COORDINATE LIMIT: inCircleDet is a degree-4 determinant in __int128, so
//      keep |x|, |y| <= ~1e9. Beyond that, shift and scale first. (This is a
//      hard limit -- past it the predicate silently returns nonsense.)
//    * COLLINEAR input has no triangulation: delaunay returns {}. euclideanMST
//      handles that case on its own, but check it if you use the triangles.
//    * Duplicate points are removed internally. The returned indices refer to
//      the ORIGINAL vector, so a duplicate simply never appears in a triangle.
//    * Four exactly co-circular points admit two valid triangulations; you get
//      one of them, deterministically. Any of them is fine for the MST.
//    * euclideanMST returns a double weight. If you need it exactly, sum the
//      squared lengths of the returned edges instead.
// ============================================================================

// triangles of the Delaunay triangulation, as indices into pts (CCW)
vector<array<int, 3>> delaunay(const vP &pts) {
	int n = (int)pts.size();
	vector<array<int, 3>> res;
	if (n < 3) return res;
	// drop duplicates but remember where each survivor came from
	vector<int> id;
	{
		vector<int> ord(n);
		for (int i = 0; i < n; i++) ord[i] = i;
		sort(ord.begin(), ord.end(), [&](int a, int b) { return pts[a] < pts[b]; });
		for (int i = 0; i < n; i++)
			if (i == 0 || !(pts[ord[i]] == pts[ord[i - 1]])) id.push_back(ord[i]);
	}
	int m = (int)id.size();
	if (m < 3) return res;
	bool flat = true;                                  // all collinear -> no triangles
	for (int i = 2; i < m && flat; i++)
		if (orient(pts[id[0]], pts[id[1]], pts[id[i]]) != 0) flat = false;
	if (flat) return res;

	vP q;
	for (int i : id) q.push_back(pts[i]);

	// One symbolic vertex "at infinity" stands in for the whole outer region:
	// each hull edge carries a GHOST triangle (x, y, FAR) covering the outside.
	// That is what makes this exact -- a finite bounding triangle would have to
	// be astronomically large to avoid distorting the triangles near the hull.
	const int FAR = m;   // the symbolic vertex "at infinity"
	int a = 0, b = 1, c = -1;
	for (int i = 2; i < m && c < 0; i++)
		if (orient(q[a], q[b], q[i]) != 0) c = i;
	if (orient(q[a], q[b], q[c]) < 0) swap(a, b);
	vector<array<int, 3>> tri{{a, b, c}, {b, a, FAR}, {c, b, FAR}, {a, c, FAR}};

	// does inserting p destroy triangle t?
	auto dead = [&](const array<int, 3> &t, int p) {
		int k = -1;
		for (int e = 0; e < 3; e++)
			if (t[e] == FAR) k = e;
		if (k < 0)                                     // ordinary triangle (CCW)
			return inCircleDet(q[t[0]], q[t[1]], q[t[2]], q[p]) > 0;
		// ghost: its "circumcircle" degenerates to the half-plane left of x->y
		P x = q[t[(k + 1) % 3]], y = q[t[(k + 2) % 3]];
		int o = orient(x, y, q[p]);
		if (o != 0) return o > 0;
		// p on the line: only conflicts when it lies ON the hull edge itself,
		// which is what stops us from emitting a zero-area triangle for it
		return sgn(dot(q[p] - x, q[p] - y)) < 0;
	};

	for (int i = 0; i < m; i++) {
		if (i == a || i == b || i == c) continue;
		map<pair<int, int>, int> border;
		vector<array<int, 3>> keep;
		for (auto &t : tri) {
			if (dead(t, i))
				for (int e = 0; e < 3; e++) border[{t[e], t[(e + 1) % 3]}]++;
			else
				keep.push_back(t);
		}
		// the cavity boundary is every directed edge whose reverse did NOT die
		for (auto &e : border)
			if (!border.count({e.first.second, e.first.first}))
				keep.push_back({e.first.first, e.first.second, i});
		tri = keep;
	}
	for (auto &t : tri)
		if (t[0] != FAR && t[1] != FAR && t[2] != FAR)
			res.push_back({id[t[0]], id[t[1]], id[t[2]]});
	return res;
}

// every Delaunay edge, once, as i < j
vector<pair<int, int>> delaunayEdges(const vP &pts) {
	set<pair<int, int>> s;
	for (auto &t : delaunay(pts))
		for (int e = 0; e < 3; e++) {
			int a = t[e], b = t[(e + 1) % 3];
			s.insert({min(a, b), max(a, b)});
		}
	return vector<pair<int, int>>(s.begin(), s.end());
}

// {total length, tree edges}. The EMST is a subgraph of the Delaunay graph.
pair<ld, vector<pair<int, int>>> euclideanMST(const vP &pts) {
	int n = (int)pts.size();
	vector<pair<int, int>> cand = delaunayEdges(pts);
	// Always add the x-sorted chain as well. It costs O(n) extra candidates and
	// it is what keeps the degenerate cases connected: duplicate points never
	// receive a Delaunay edge, and fully collinear input has no triangles at
	// all. Extra real edges can never make Kruskal's answer too small.
	{
		vector<int> ord(n);
		for (int i = 0; i < n; i++) ord[i] = i;
		sort(ord.begin(), ord.end(), [&](int a, int b) { return pts[a] < pts[b]; });
		for (int i = 0; i + 1 < n; i++) cand.push_back({ord[i], ord[i + 1]});
	}
	vector<pair<ld, pair<int, int>>> es;
	for (auto &e : cand) es.push_back({dist(pts[e.first], pts[e.second]), e});
	sort(es.begin(), es.end());
	vector<int> par(n);
	for (int i = 0; i < n; i++) par[i] = i;
	function<int(int)> find = [&](int x) { return par[x] == x ? x : par[x] = find(par[x]); };
	ld total = 0;
	vector<pair<int, int>> used;
	for (auto &e : es) {
		int a = find(e.second.first), b = find(e.second.second);
		if (a == b) continue;
		par[a] = b;
		total += e.first;
		used.push_back(e.second);
	}
	return {total, used};
}

// the Voronoi vertices: one circumcentre per Delaunay triangle. The Voronoi
// cell of point i is bounded by the circumcentres of the triangles touching i,
// in angular order around it.
vPd voronoiVertices(const vP &pts, const vector<array<int, 3>> &tri) {
	vPd v;
	for (auto &t : tri) v.push_back(circumcircle(pts[t[0]], pts[t[1]], pts[t[2]]).c);
	return v;
}

// nearest other point to pts[i]; the nearest neighbour is always a Delaunay edge
int nearestNeighborOf(const vP &pts, int i, const vector<array<int, 3>> &tri) {
	int best = -1;
	ll bd = LLONG_MAX;
	for (auto &t : tri)
		for (int e = 0; e < 3; e++)
			if (t[e] == i)
				for (int o = 0; o < 3; o++)
					if (o != e && dist2(pts[i], pts[t[o]]) < bd)
						bd = dist2(pts[i], pts[t[o]]), best = t[o];
	return best;
}

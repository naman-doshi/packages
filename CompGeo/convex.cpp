// ============================================================================
//  CONVEX HULL + CONVEX POLYGON TOOLS -- REQUIRES point.cpp, lines.cpp, polygon.cpp
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    convexHull is the single most useful routine in geometry. Anything of the
//    form "extreme / farthest / maximum linear objective over a point set" only
//    ever cares about hull vertices. Everything below assumes the hull comes
//    out of convexHull: CCW, starting at the lexicographically smallest point.
//
//  API
//    convexHull(pts, keepCollinear=false)   CCW hull. EXACT for P.
//    hullIdx(pts, keepCollinear=false)      the same, as indices into pts
//    hullDiameter(h)      -> {i, j}, the FARTHEST pair of the point set
//    farthestPair(pts)    -> {i, j} indices into pts (hull + calipers)
//    hullWidth(h)         min width over all directions (a double)
//    minAreaRect(h)       -> {area, 4 corners}      rotating calipers
//    minPerimeterRect(h)  -> {perimeter, 4 corners}
//    maxDot(h, dir)       max of dot(dir, v) over vertices, O(log n)
//    extreme(h, dir)      the index attaining it,        O(log n)
//    inConvex(h, p)       0 out / 1 boundary / 2 inside, O(log n)
//    tangents(h, p)       {right, left} tangent indices from an OUTSIDE p, O(n)
//    minkowski(A, B)      Minkowski sum of two convex polygons, O(n+m)
//    convexInterArea(A,B) area of the intersection of two convex polygons
//    polyDist(A, B)       min distance between two polygons (0 if they touch)
//    maxTriangle(h)       2 * max area of a triangle with vertices on h
//    convexLayers(pts)    onion peeling: hull, then hull of the rest, ...
//
//  COMPLEXITY  hull O(n log n) (the sort); calipers O(n); extreme / inConvex
//    O(log n); minkowski O(n+m); maxTriangle O(n^2); convexLayers O(n^2 log n).
//
//  PITFALLS
//    * keepCollinear=false (the default) drops points lying in the middle of a
//      hull edge. Almost every algorithm here NEEDS that -- calipers, extreme
//      and inConvex assume no three consecutive vertices are collinear. Only
//      pass true when the problem literally asks for "points on the boundary".
//    * With keepCollinear=true and ALL points collinear you get the segment
//      traversed out and back (2k-2 vertices), which is correct but surprising.
//    * The hull of 0/1/2 distinct points is returned as 0/1/2 points, not a
//      polygon. Guard h.size() < 3 before calling the calipers routines.
//    * hullDiameter returns squared distance if you take dist2 -- keep it
//      squared and integral, don't sqrt until output.
//    * inConvex needs the hull CCW. If you built the polygon yourself, run
//      makeCCW first.
//    * minkowski's output can contain collinear vertices; re-run convexHull if
//      you need it strict.
// ============================================================================

// CCW convex hull, starting at the lexicographically smallest point.
template <class T> vector<Pt<T>> convexHull(vector<Pt<T>> p, bool keepCollinear = false) {
	sort(p.begin(), p.end());
	p.erase(unique(p.begin(), p.end()), p.end());
	int n = (int)p.size(), k = 0;
	if (n <= 2) return p;
	auto bad = [&](Pt<T> o, Pt<T> a, Pt<T> b) {
		int c = orient(o, a, b);
		return keepCollinear ? c < 0 : c <= 0;
	};
	vector<Pt<T>> h(2 * n);
	for (int i = 0; i < n; i++) {                       // lower hull
		while (k >= 2 && bad(h[k - 2], h[k - 1], p[i])) k--;
		h[k++] = p[i];
	}
	for (int i = n - 2, t = k + 1; i >= 0; i--) {       // upper hull
		while (k >= t && bad(h[k - 2], h[k - 1], p[i])) k--;
		h[k++] = p[i];
	}
	h.resize(k - 1);
	return h;
}

// same hull but as indices into the ORIGINAL vector (duplicates: keeps one)
template <class T> vector<int> hullIdx(const vector<Pt<T>> &p, bool keepCollinear = false) {
	int n = (int)p.size();
	vector<int> id(n);
	for (int i = 0; i < n; i++) id[i] = i;
	sort(id.begin(), id.end(), [&](int a, int b) { return p[a] < p[b]; });
	id.erase(unique(id.begin(), id.end(), [&](int a, int b) { return p[a] == p[b]; }), id.end());
	int m = (int)id.size(), k = 0;
	if (m <= 2) return id;
	auto bad = [&](int o, int a, int b) {
		int c = orient(p[o], p[a], p[b]);
		return keepCollinear ? c < 0 : c <= 0;
	};
	vector<int> h(2 * m);
	for (int i = 0; i < m; i++) {
		while (k >= 2 && bad(h[k - 2], h[k - 1], id[i])) k--;
		h[k++] = id[i];
	}
	for (int i = m - 2, t = k + 1; i >= 0; i--) {
		while (k >= t && bad(h[k - 2], h[k - 1], id[i])) k--;
		h[k++] = id[i];
	}
	h.resize(k - 1);
	return h;
}

// ------------------------------------------------------- rotating calipers
// indices of the two farthest-apart vertices of a strictly convex CCW hull
template <class T> pair<int, int> hullDiameter(const vector<Pt<T>> &h) {
	int n = (int)h.size();
	if (n < 2) return {0, 0};
	if (n == 2) return {0, 1};
	pair<T, pair<int, int>> best{0, {0, 0}};
	for (int i = 0, j = 1; i < n; i++) {
		int ni = (i + 1) % n;
		while (abs(cross(h[i], h[ni], h[(j + 1) % n])) > abs(cross(h[i], h[ni], h[j])))
			j = (j + 1) % n;
		best = max(best, {dist2(h[i], h[j]), {i, j}});
		best = max(best, {dist2(h[ni], h[j]), {ni, j}});
	}
	return best.second;
}
// farthest pair over an arbitrary point set, as indices into pts
template <class T> pair<int, int> farthestPair(const vector<Pt<T>> &pts) {
	vector<int> id = hullIdx(pts);
	vector<Pt<T>> h;
	for (int i : id) h.push_back(pts[i]);
	auto d = hullDiameter(h);
	return {id[d.first], id[d.second]};
}

// min width: the smallest distance between two parallel supporting lines
template <class T> ld hullWidth(const vector<Pt<T>> &h) {
	int n = (int)h.size();
	if (n < 3) return 0;
	ld best = 1e18;
	for (int i = 0, j = 1; i < n; i++) {
		int ni = (i + 1) % n;
		while (abs(cross(h[i], h[ni], h[(j + 1) % n])) > abs(cross(h[i], h[ni], h[j])))
			j = (j + 1) % n;
		best = min(best, (ld)abs((ld)cross(h[i], h[ni], h[j])) / dist(h[i], h[ni]));
	}
	return best;
}

// The minimum-area (or minimum-perimeter) enclosing rectangle always has one
// side flush with a hull edge -- so try every edge, with three calipers.
template <class T>
pair<ld, array<Pd, 4>> minRect(const vector<Pt<T>> &h, bool byPerimeter) {
	int n = (int)h.size();
	array<Pd, 4> corners{};
	if (n < 3) {
		for (int i = 0; i < 4 && n > 0; i++) corners[i] = toD(h[i % n]);
		return {0, corners};
	}
	ld best = 1e18;
	int up = 1, rt = 1, lf = 1;
	for (int i = 0; i < n; i++) {
		Pt<T> e = h[(i + 1) % n] - h[i];
		while (cross(e, h[(up + 1) % n] - h[i]) > cross(e, h[up] - h[i])) up = (up + 1) % n;
		while (dot(e, h[(rt + 1) % n]) > dot(e, h[rt])) rt = (rt + 1) % n;
		if (i == 0) lf = rt;
		// <= (not <) so a tie at the maximum does not stall the walk immediately
		while (dot(e, h[(lf + 1) % n]) <= dot(e, h[lf])) lf = (lf + 1) % n;
		ld L = len(e);
		ld w = (ld)(dot(e, h[rt]) - dot(e, h[lf])) / L;
		ld ht = (ld)cross(e, h[up] - h[i]) / L;
		ld val = byPerimeter ? 2 * (w + ht) : w * ht;
		if (val < best) {
			best = val;
			Pd A = toD(h[i]), u = toD(e) / L, v = perp(u);
			ld dl = dot(u, toD(h[lf]) - A), dr = dot(u, toD(h[rt]) - A);
			corners[0] = A + u * dl;
			corners[1] = A + u * dr;
			corners[2] = corners[1] + v * ht;
			corners[3] = corners[0] + v * ht;
		}
	}
	return {best, corners};
}
template <class T> pair<ld, array<Pd, 4>> minAreaRect(const vector<Pt<T>> &h) {
	return minRect(h, false);
}
template <class T> pair<ld, array<Pd, 4>> minPerimeterRect(const vector<Pt<T>> &h) {
	return minRect(h, true);
}

// ------------------------------------------------------ O(log n) queries
// index of a vertex maximising dot(dir, h[i]). h must be STRICTLY convex, CCW.
template <class T> int extreme(const vector<Pt<T>> &h, Pt<T> dir) {
	int n = (int)h.size();
	if (n < 3) {
		int b = 0;
		for (int i = 1; i < n; i++)
			if (dot(dir, h[i]) > dot(dir, h[b])) b = i;
		return b;
	}
	auto cmp = [&](int i, int j) { return sgn(dot(dir, h[j % n] - h[i % n])); };
	auto extr = [&](int i) { return cmp(i + 1, i) >= 0 && cmp(i, i + n - 1) < 0; };
	int lo = 0, hi = n;
	if (extr(0)) return 0;
	while (lo + 1 < hi) {
		int m = (lo + hi) / 2;
		if (extr(m)) return m;
		int ls = cmp(lo + 1, lo), ms = cmp(m + 1, m);
		((ls < ms || (ls == ms && ls == cmp(lo, m))) ? hi : lo) = m;
	}
	return lo;
}
template <class T> T maxDot(const vector<Pt<T>> &h, Pt<T> dir) {
	return dot(dir, h[extreme(h, dir)]);
}

// 0 outside, 1 on the boundary, 2 strictly inside. h CCW, O(log n).
template <class T> int inConvex(const vector<Pt<T>> &h, Pt<T> q) {
	int n = (int)h.size();
	if (n == 0) return 0;
	if (n == 1) return q == h[0] ? 1 : 0;
	if (n == 2) return onSegment(q, h[0], h[1]) ? 1 : 0;
	if (orient(h[0], h[1], q) < 0 || orient(h[0], h[n - 1], q) > 0) return 0;
	int lo = 1, hi = n - 1;
	while (hi - lo > 1) {                                 // wedge from h[0]
		int m = (lo + hi) / 2;
		(orient(h[0], h[m], q) >= 0 ? lo : hi) = m;
	}
	int s = orient(h[lo], h[lo + 1], q);
	if (s < 0) return 0;
	if (s == 0) return 1;
	if (orient(h[0], h[1], q) == 0 && onSegment(q, h[0], h[1])) return 1;
	if (orient(h[0], h[n - 1], q) == 0 && onSegment(q, h[0], h[n - 1])) return 1;
	return 2;
}

// tangent vertices from a point OUTSIDE the hull: everything lies in the wedge
// p->h[res.first] .. p->h[res.second], turning CCW. O(n).
template <class T> pair<int, int> tangents(const vector<Pt<T>> &h, Pt<T> p) {
	int n = (int)h.size(), r = 0, l = 0;
	for (int i = 1; i < n; i++) {
		if (orient(p, h[r], h[i]) < 0) r = i;             // clockwise-most
		if (orient(p, h[l], h[i]) > 0) l = i;             // counter-clockwise-most
	}
	return {r, l};
}

// ------------------------------------------------------------- combination
// Minkowski sum {a + b}. A and B must be convex and CCW. O(n + m).
template <class T> vector<Pt<T>> minkowski(vector<Pt<T>> A, vector<Pt<T>> B) {
	auto reorder = [](vector<Pt<T>> &v) {                 // lowest, then leftmost
		int k = 0;
		for (int i = 1; i < (int)v.size(); i++)
			if (v[i].y < v[k].y || (v[i].y == v[k].y && v[i].x < v[k].x)) k = i;
		rotate(v.begin(), v.begin() + k, v.end());
	};
	if (A.empty() || B.empty()) return {};
	if (A.size() == 1 || B.size() == 1) {
		vector<Pt<T>> r;
		for (auto &a : A) for (auto &b : B) r.push_back(a + b);
		return convexHull(r);
	}
	reorder(A); reorder(B);
	int n = (int)A.size(), m = (int)B.size();
	A.push_back(A[0]); A.push_back(A[1]);
	B.push_back(B[0]); B.push_back(B[1]);
	vector<Pt<T>> res;
	int i = 0, j = 0;
	while (i < n || j < m) {
		res.push_back(A[i] + B[j]);
		T c = cross(A[i + 1] - A[i], B[j + 1] - B[j]);
		if (sgn(c) >= 0 && i < n) i++;
		if (sgn(c) <= 0 && j < m) j++;
	}
	return res;
}

// area of the intersection of two CONVEX CCW polygons
template <class T> ld convexInterArea(const vector<Pt<T>> &A, const vector<Pt<T>> &B) {
	vPd a, b;
	for (auto &p : A) a.push_back(toD(p));
	for (auto &p : B) b.push_back(toD(p));
	vPd c = convexClip(a, b);
	return c.size() < 3 ? 0 : polyArea(c);
}

// minimum distance between the boundaries of two polygons (0 if they cross or
// touch). O(n*m) -- plenty for the usual constraints.
template <class T> ld polyDist(const vector<Pt<T>> &A, const vector<Pt<T>> &B) {
	int n = (int)A.size(), m = (int)B.size();
	ld best = 1e18;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			best = min(best, segDist(A[i], A[(i + 1) % n], B[j], B[(j + 1) % m]));
	return best;
}

// 2 * the largest triangle area with all three vertices on a convex CCW hull
template <class T> T maxTriangle(const vector<Pt<T>> &h) {
	int n = (int)h.size();
	T best = 0;
	if (n < 3) return 0;
	for (int i = 0; i < n; i++) {
		int k = (i + 2) % n;
		for (int j = (i + 1) % n; j != i && (j + 1) % n != i; j = (j + 1) % n) {
			if (k == j) k = (k + 1) % n;
			while ((k + 1) % n != i &&
			       abs(cross(h[i], h[j], h[(k + 1) % n])) >= abs(cross(h[i], h[j], h[k])))
				k = (k + 1) % n;
			best = max(best, abs(cross(h[i], h[j], h[k])));
		}
	}
	return best;
}

// onion peeling: layer 0 is the hull, layer 1 the hull of what is left, ...
template <class T> vector<vector<Pt<T>>> convexLayers(vector<Pt<T>> p) {
	vector<vector<Pt<T>>> res;
	sort(p.begin(), p.end());
	p.erase(unique(p.begin(), p.end()), p.end());
	while (!p.empty()) {
		vector<Pt<T>> h = convexHull(p);
		res.push_back(h);
		vector<Pt<T>> rest;
		for (auto &q : p)
			if (find(h.begin(), h.end(), q) == h.end()) rest.push_back(q);
		p = rest;
	}
	return res;
}

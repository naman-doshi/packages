// ============================================================================
//  ANGULAR SWEEP  -- the "rotate a ray and count" archetypes
//                    REQUIRES point.cpp
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    A family of problems that all look different and are all the same trick:
//    fix one point, sort everything else by ANGLE around it, then sweep.
//    Recognise it when the statement mentions "a line through", "a circle of
//    radius r placed anywhere", "contains the origin", "sees", or "collinear".
//
//  API
//    maxPointsOnLine(pts)         most points on one straight line
//    maxCollinearThrough(pts, i)  most points on a line through pts[i]
//    maxPointsInCircle(pts, r)    most points a disk of radius r can cover
//    countTrianglesContaining(pts, q)   triangles from pts that contain q
//    countPointsInHalf(pts, o, d) points strictly left of the ray o -> o+d
//    angularOrder(pts, o)         indices of pts sorted CCW around o
//
//  COMPLEXITY  maxPointsOnLine and maxPointsInCircle are O(n^2 log n) -- that
//    is the intended complexity for these, n is typically <= 2000.
//    countTrianglesContaining is O(n log n).
//
//  PITFALLS
//    * maxPointsOnLine normalises each direction by its gcd and forces a
//      canonical sign, so it is EXACT -- never compare slopes as doubles.
//    * Duplicate points: maxPointsOnLine counts them all as collinear (they
//      are), which is usually what you want. maxPointsInCircle also counts
//      them. Deduplicate first if the problem says "distinct".
//    * maxPointsInCircle counts points ON the boundary as covered, and uses a
//      1e-9 tolerance to do it. If the problem says STRICTLY inside, shrink r
//      by an epsilon.
//    * countTrianglesContaining assumes q is not one of the points and that no
//      two points are collinear WITH q (general position). Degenerate input
//      makes "contains" ambiguous anyway -- read the statement.
// ============================================================================

// indices of pts sorted counter-clockwise around o, starting from the +x axis
template <class T> vector<int> angularOrder(const vector<Pt<T>> &pts, Pt<T> o) {
	int n = (int)pts.size();
	vector<int> id(n);
	for (int i = 0; i < n; i++) id[i] = i;
	sort(id.begin(), id.end(), [&](int a, int b) { return polarCmp(pts[a] - o, pts[b] - o); });
	return id;
}

// most points on a single line through pts[i] (counting pts[i] itself)
int maxCollinearThrough(const vP &pts, int i) {
	int n = (int)pts.size(), same = 0;
	map<pair<ll, ll>, int> cnt;
	for (int j = 0; j < n; j++) {
		if (j == i) continue;
		P d = pts[j] - pts[i];
		if (d.x == 0 && d.y == 0) { same++; continue; }
		ll g = gcd(llabs(d.x), llabs(d.y));
		d.x /= g; d.y /= g;
		if (d.x < 0 || (d.x == 0 && d.y < 0)) d.x = -d.x, d.y = -d.y;  // canonical
		cnt[{d.x, d.y}]++;
	}
	int best = 0;
	for (auto &e : cnt) best = max(best, e.second);
	return best + same + 1;
}
int maxPointsOnLine(const vP &pts) {
	int n = (int)pts.size(), best = n ? 1 : 0;
	for (int i = 0; i < n; i++) best = max(best, maxCollinearThrough(pts, i));
	return best;
}

// the largest number of points a disk of radius r can cover, placed anywhere.
// An optimal disk can always be slid until a point sits on its boundary, so fix
// each point on the boundary and sweep the centre around it.
int maxPointsInCircle(const vP &pts, ld r) {
	int n = (int)pts.size();
	if (n == 0) return 0;
	int best = 1;
	ld R2 = 4 * r * r;
	for (int i = 0; i < n; i++) {
		vector<pair<ld, int>> ev;                 // (angle, +1 enter / -1 leave)
		for (int j = 0; j < n; j++) {
			if (i == j) continue;
			ld d2 = (ld)dist2(pts[i], pts[j]);
			if (d2 > R2 + EPS) continue;
			ld d = sqrt(d2);
			ld A = atan2((ld)(pts[j].y - pts[i].y), (ld)(pts[j].x - pts[i].x));
			ld B = acos(max((ld)-1, min((ld)1, d / (2 * r))));
			ld lo = fmod(A - B, 2 * PI);
			if (lo < 0) lo += 2 * PI;
			ld hi = lo + 2 * B;
			if (hi > 2 * PI) {
				ev.push_back({lo, 1}); ev.push_back({2 * PI, -1});
				ev.push_back({0, 1}); ev.push_back({hi - 2 * PI, -1});
			} else {
				ev.push_back({lo, 1}); ev.push_back({hi, -1});
			}
		}
		sort(ev.begin(), ev.end(), [](const pair<ld, int> &a, const pair<ld, int> &b) {
			return a.first != b.first ? a.first < b.first : a.second > b.second;
		});
		int cur = 1;                              // pts[i] is always covered
		for (auto &e : ev) best = max(best, cur += e.second);
	}
	return best;
}

// points strictly to the LEFT of the ray o -> o+d
template <class T> int countPointsInHalf(const vector<Pt<T>> &pts, Pt<T> o, Pt<T> d) {
	int c = 0;
	for (auto &p : pts) c += sgn(cross(d, p - o)) > 0;
	return c;
}

// number of triples of pts whose triangle strictly contains q. Complement
// counting: a triangle FAILS to contain q exactly when all three points lie in
// some half-plane through q, i.e. inside a 180 degree wedge.
ll countTrianglesContaining(const vP &pts, P q) {
	int n = (int)pts.size();
	if (n < 3) return 0;
	vector<ld> a;
	for (auto &p : pts) a.push_back(atan2((ld)(p.y - q.y), (ld)(p.x - q.x)));
	sort(a.begin(), a.end());
	ll total = (ll)n * (n - 1) * (n - 2) / 6, bad = 0;
	int j = 0;
	for (int i = 0; i < n; i++) {                 // count points in (a[i], a[i]+pi)
		if (j < i + 1) j = i + 1;
		while (j < i + n && a[j % n] + (j >= n ? 2 * PI : 0) - a[i] < PI - EPS) j++;
		ll k = j - i - 1;                         // how many fit in the half-plane
		bad += k * (k - 1) / 2;
	}
	return total - bad;
}

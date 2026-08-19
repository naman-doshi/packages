// ============================================================================
//  POLYGONS (any simple polygon)   -- REQUIRES point.cpp, lines.cpp
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Area, containment, cutting, triangulation, lattice counting. Works for any
//    SIMPLE polygon (no self-intersections); convex-only tricks that are faster
//    (O(log n) containment, tangents, Minkowski) live in convex.cpp.
//
//    Convention: vertices in a vector, in order, NO repeated last vertex.
//    Most routines want CCW order -- call makeCCW(poly) once and stop worrying.
//
//  API
//    polyArea2(poly)        2 * signed area, EXACT integer. >0 iff CCW.
//    polyArea(poly)         |area| as a double
//    makeCCW(poly)          reverse in place if it was clockwise
//    perimeter(poly)
//    isConvex(poly)         true for convex (collinear vertices allowed)
//    isSimple(poly)         O(n^2) self-intersection check
//    inPolygon(poly, p)     0 = outside, 1 = ON the boundary, 2 = strictly in
//    windingNumber(poly, p) turns around p; !=0 means inside (works for
//                           self-intersecting polygons too)
//    centroid(poly)         area-weighted centre of mass
//    distToPoly(poly, p)    distance to the boundary (0 if p is on it)
//    polygonCut(poly, a, b) the part of poly LEFT of the directed line a->b
//    convexClip(poly, clip) intersection of a polygon with a CONVEX polygon
//    lineInterPoly(poly, a, b)  EVERY point where the infinite line a->b meets
//                           the boundary, in order along a->b, deduplicated
//    lineInsidePoly(poly, a, b) the pieces of that line lying inside poly, as
//                           {entry, exit} pairs -- sum their lengths for "how
//                           much of the line the polygon covers"
//    triangulate(poly)      ear clipping -> list of vertex-index triples
//    pointInTriangle(p,a,b,c)
//    boundaryLattice(poly)  B: lattice points on the boundary
//    interiorLattice(poly)  I: lattice points strictly inside  (Pick's theorem)
//    normalizePoly(poly)    rotate so the smallest vertex is first (for ==)
//    minRotation(v)         start of the least cyclic rotation, O(n) (Booth)
//    canonicalPoly(poly)    canonical rotation/direction of the vertex cycle
//    polyHash(poly)         64-bit hash of that -- same polygon, same hash
//    congruenceHash(poly)   hash of the SHAPE: equal for congruent polygons
//                           (translated / rotated / mirrored / relisted)
//    forEachConvexPolygon(pts, f)  every convex polygon over a POINT SET, one
//                           at a time (callback); allConvexPolygons(pts) to
//                           collect them, countConvexPolygons(pts) to count
//    PolyCount(verts, pts)  how many of pts are inside a query polygon whose
//                           vertices come from verts, O(k log m) per query
//
//  COMPLEXITY  O(n) except isSimple / triangulate (O(n^2)) and convexClip
//    (O(n*m)). inPolygon is O(n) PER QUERY -- if the polygon is convex and you
//    have many queries use convex.cpp's inConvex, which is O(log n).
//
//  PITFALLS
//    * polyArea2 returns TWICE the area so it stays an exact integer. Divide by
//      2.0 (not 2) at the very end. It is SIGNED -- take abs() unless you are
//      using the sign to detect orientation.
//    * inPolygon's boundary case is the one that gets you. It returns 1 there;
//      decide deliberately whether the problem counts the boundary as inside.
//    * Pick's theorem (A = I + B/2 - 1) needs INTEGER vertices. It is the fast
//      way to count lattice points -- never loop over the bounding box.
//    * triangulate assumes the polygon is simple. It will silently produce
//      garbage on a self-intersecting one; run isSimple first if unsure.
//    * lineInterPoly counts a vertex the line merely TOUCHES as an
//      intersection, and gives two points for an edge lying along the line, so
//      its size is not "number of times the line crosses in or out" and can be
//      odd. Use lineInsidePoly if what you mean is inside.
//    * polygonCut/convexClip return doubles. If you only need the AREA of a
//      convex clip and inputs are integral, prefer halfplane.cpp or keep the
//      cut symbolic.
// ============================================================================

// 2 * signed area (shoelace). Positive iff the vertices are CCW. EXACT for P.
template <class T> T polyArea2(const vector<Pt<T>> &v) {
	T a = 0;
	int n = (int)v.size();
	for (int i = 0; i < n; i++) a += cross(v[i], v[(i + 1) % n]);
	return a;
}
template <class T> ld polyArea(const vector<Pt<T>> &v) {
	return (ld)abs((ld)polyArea2(v)) / 2.0;
}
template <class T> void makeCCW(vector<Pt<T>> &v) {
	if (polyArea2(v) < 0) reverse(v.begin(), v.end());
}
template <class T> ld perimeter(const vector<Pt<T>> &v) {
	ld s = 0;
	int n = (int)v.size();
	for (int i = 0; i < n; i++) s += dist(v[i], v[(i + 1) % n]);
	return s;
}

// convex in either orientation; collinear vertices are tolerated
template <class T> bool isConvex(const vector<Pt<T>> &v) {
	int n = (int)v.size(), pos = 0, neg = 0;
	if (n < 3) return true;
	for (int i = 0; i < n; i++) {
		int o = orient(v[i], v[(i + 1) % n], v[(i + 2) % n]);
		if (o > 0) pos++;
		if (o < 0) neg++;
	}
	return !pos || !neg;
}

// no two non-adjacent edges touch; adjacent edges meet only at their vertex
template <class T> bool isSimple(const vector<Pt<T>> &v) {
	int n = (int)v.size();
	if (n < 3) return false;
	for (int i = 0; i < n; i++) {
		if (v[i] == v[(i + 1) % n]) return false;            // zero-length edge
		// consecutive edges must not fold back over each other
		Pt<T> a = v[i], s = v[(i + 1) % n], c = v[(i + 2) % n];
		if (orient(a, s, c) == 0 && sgn(dot(a - s, c - s)) > 0) return false;
	}
	for (int i = 0; i < n; i++)
		for (int j = i + 2; j < n; j++) {
			if (i == 0 && j == n - 1) continue;              // adjacent, done above
			if (segInter(v[i], v[(i + 1) % n], v[j], v[(j + 1) % n])) return false;
		}
	return true;
}

// 0 outside, 1 on the boundary, 2 strictly inside. EXACT for P.
template <class T> int inPolygon(const vector<Pt<T>> &v, Pt<T> p) {
	int n = (int)v.size(), cnt = 0;
	for (int i = 0; i < n; i++) {
		Pt<T> a = v[i], b = v[(i + 1) % n];
		if (onSegment(p, a, b)) return 1;
		bool up = a.y < b.y;                 // orient the edge upward
		if (!up) swap(a, b);
		// half-open [a.y, b.y): counts each crossing of the +x ray exactly once
		if (sgn(a.y - p.y) <= 0 && sgn(b.y - p.y) > 0 && orient(a, b, p) > 0) cnt++;
	}
	return (cnt & 1) ? 2 : 0;
}

// how many times the boundary winds around p (0 = outside). Handles polygons
// that self-intersect, where the even-odd rule above is ambiguous.
template <class T> int windingNumber(const vector<Pt<T>> &v, Pt<T> p) {
	int n = (int)v.size(), w = 0;
	for (int i = 0; i < n; i++) {
		Pt<T> a = v[i], b = v[(i + 1) % n];
		if (sgn(a.y - p.y) <= 0) {
			if (sgn(b.y - p.y) > 0 && orient(a, b, p) > 0) w++;
		} else if (sgn(b.y - p.y) <= 0 && orient(a, b, p) < 0)
			w--;
	}
	return w;
}

// area-weighted centre of mass of the (filled) polygon
template <class T> Pd centroid(const vector<Pt<T>> &v) {
	int n = (int)v.size();
	ld a = 0;
	Pd c(0, 0);
	for (int i = 0; i < n; i++) {
		Pd p = toD(v[i]), q = toD(v[(i + 1) % n]);
		ld cr = cross(p, q);
		a += cr;
		c += (p + q) * cr;
	}
	return c / (3 * a);
}

template <class T> ld distToPoly(const vector<Pt<T>> &v, Pd p) {
	int n = (int)v.size();
	ld best = 1e18;
	for (int i = 0; i < n; i++) best = min(best, distToSeg(p, v[i], v[(i + 1) % n]));
	return best;
}

// inclusive of the boundary; a, b, c in any orientation
template <class T> bool pointInTriangle(Pt<T> p, Pt<T> a, Pt<T> b, Pt<T> c) {
	int d1 = orient(a, b, p), d2 = orient(b, c, p), d3 = orient(c, a, p);
	bool neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
	bool pos = (d1 > 0) || (d2 > 0) || (d3 > 0);
	return !(neg && pos);
}

// ------------------------------------------------------------------ cutting
// keep the part of poly on the LEFT of the directed line a->b (the half-plane
// cross(b-a, x-a) >= 0). poly may be any simple polygon; the result of one cut
// is a polygon only if poly is convex, otherwise use it as a clipping step.
vPd polygonCut(const vPd &poly, Pd a, Pd b) {
	vPd res;
	int n = (int)poly.size();
	for (int i = 0; i < n; i++) {
		Pd cur = poly[i], nxt = poly[(i + 1) % n];
		ld s1 = cross(a, b, cur), s2 = cross(a, b, nxt);
		if (sgn(s1) >= 0) res.push_back(cur);
		if (sgn(s1) * sgn(s2) < 0) res.push_back(cur + (nxt - cur) * (s1 / (s1 - s2)));
	}
	return res;
}
// Sutherland-Hodgman: poly (any) clipped to the inside of the CONVEX polygon
// `clip` (which must be CCW). Result may be empty.
vPd convexClip(vPd poly, const vPd &clip) {
	int m = (int)clip.size();
	for (int i = 0; i < m && !poly.empty(); i++)
		poly = polygonCut(poly, clip[i], clip[(i + 1) % m]);
	return poly;
}

// ------------------------------------------- a LINE against the polygon
// Every point where the INFINITE line a->b meets the BOUNDARY of poly, sorted
// along a->b and deduplicated. Any simple polygon, either orientation. The
// side tests are EXACT on P; only a crossing in the middle of an edge has to
// be computed, so a vertex hit comes back with its coordinates untouched.
//
//   crossing        the ordinary case: one point per edge the line passes
//                   through
//   vertex hit      reported ONCE, not once per incident edge, and reported
//                   whether the line crosses the boundary there or merely
//                   touches it (a tangent vertex is still an intersection)
//   collinear edge  the line contains a whole edge -> BOTH its endpoints are
//                   reported, so the overlap is the piece between them
//
// Do not read the result as entry/exit pairs -- with touches and collinear
// edges in play the count can be odd, and consecutive points need not bracket
// the interior. lineInsidePoly below is the one that answers "inside".
//
// For a SEGMENT rather than a line, keep the points with onSegment(p, a, b).
//
// TIME O(n log n), the sort.  a must differ from b.
template <class T> vPd lineInterPoly(const vector<Pt<T>> &poly, Pt<T> a, Pt<T> b) {
	int n = (int)poly.size();
	vPd res;
	if (n < 3 || a == b) return res;
	Pt<T> d = b - a;
	for (int i = 0; i < n; i++) {
		Pt<T> u = poly[i], v = poly[(i + 1) % n];
		T cu = cross(d, u - a), cv = cross(d, v - a);   // signed side of each end
		int su = sgn(cu), sv = sgn(cv);
		if (!su && !sv) { res.push_back(toD(u)); res.push_back(toD(v)); }  // edge on it
		else if (!su) res.push_back(toD(u));            // touches at this end only
		else if (!sv) res.push_back(toD(v));
		else if (su * sv < 0) {                         // strictly straddles: solve
			// the subtraction is done in ld: cu - cv can overflow T even when
			// each cross() still fits
			ld t = (ld)cu / ((ld)cu - (ld)cv);
			res.push_back(toD(u) + (toD(v) - toD(u)) * t);
		}
	}
	Pd A = toD(a), D = toD(d);
	// every point is on the line, so the projection onto d orders them
	sort(res.begin(), res.end(), [&](Pd p, Pd q) { return dot(p - A, D) < dot(q - A, D); });
	res.erase(unique(res.begin(), res.end(),
	                 [](Pd p, Pd q) { return near(p, q); }), res.end());
	return res;
}

// The maximal pieces of the line a->b that lie INSIDE poly (the boundary
// counts as inside), as {entry, exit} pairs in order along a->b.
//   ld covered = 0; for (auto &s : lineInsidePoly(poly, a, b)) covered += dist(s.first, s.second);
// Consecutive boundary hits bracket a stretch that is entirely in or entirely
// out -- nothing can change in between -- so one containment test per gap
// settles it, which is why reflex vertices, tangent vertices and edges lying
// along the line all fall out correctly instead of needing special cases.
// Pieces touching end to end are merged, so the output is disjoint and maximal.
//
// Zero-length pieces are NOT reported: a line that only grazes one vertex, or
// runs along an edge of a polygon it never enters, gives {} here even though
// lineInterPoly lists those points. Ask lineInterPoly if you need them.
//
// TIME O(n^2) -- O(n) containment tests of O(n) each. For a CONVEX polygon
// there are at most two hits, so it is O(n).
template <class T> vector<pair<Pd, Pd>> lineInsidePoly(const vector<Pt<T>> &poly, Pt<T> a, Pt<T> b) {
	vPd cut = lineInterPoly(poly, a, b), pd;
	vector<pair<Pd, Pd>> res;
	for (auto &p : poly) pd.push_back(toD(p));
	for (int i = 0; i + 1 < (int)cut.size(); i++) {
		Pd mid = (cut[i] + cut[i + 1]) * (ld)0.5;
		if (inPolygon(pd, mid) == 0) continue;
		if (!res.empty() && near(res.back().second, cut[i])) res.back().second = cut[i + 1];
		else res.push_back({cut[i], cut[i + 1]});
	}
	return res;
}

// ------------------------------------------------------------ triangulation
// ear clipping on a SIMPLE polygon; returns n-2 triangles as index triples.
template <class T> vector<array<int, 3>> triangulate(vector<Pt<T>> v) {
	int n = (int)v.size();
	vector<array<int, 3>> res;
	if (n < 3) return res;
	vector<int> id(n);
	for (int i = 0; i < n; i++) id[i] = i;
	if (polyArea2(v) < 0) reverse(id.begin(), id.end());   // work CCW
	while ((int)id.size() > 3) {
		int m = (int)id.size();
		bool cut = false;
		for (int i = 0; i < m; i++) {
			int a = id[(i + m - 1) % m], b = id[i], c = id[(i + 1) % m];
			if (orient(v[a], v[b], v[c]) <= 0) continue;   // reflex / straight
			bool ok = true;
			for (int j : id) {
				if (j == a || j == b || j == c) continue;
				if (pointInTriangle(v[j], v[a], v[b], v[c])) { ok = false; break; }
			}
			if (!ok) continue;
			res.push_back({a, b, c});
			id.erase(id.begin() + i);
			cut = true;
			break;
		}
		if (!cut) break;                                   // degenerate input
	}
	if (id.size() == 3) res.push_back({id[0], id[1], id[2]});
	return res;
}

// ------------------------------------------------------- lattice / Pick's
// B = lattice points ON the boundary
ll boundaryLattice(const vP &v) {
	ll b = 0;
	int n = (int)v.size();
	for (int i = 0; i < n; i++) b += gcd(llabs(v[i].x - v[(i + 1) % n].x),
	                                     llabs(v[i].y - v[(i + 1) % n].y));
	return b;
}
// I = lattice points strictly INSIDE.  Pick: A = I + B/2 - 1  ->  I = A - B/2 + 1
ll interiorLattice(const vP &v) {
	ll a2 = llabs(polyArea2(v)), b = boundaryLattice(v);
	return (a2 - b + 2) / 2;
}

// rotate so the lexicographically smallest vertex comes first -- lets you
// compare two polygons of the same orientation with ==
template <class T> void normalizePoly(vector<Pt<T>> &v) {
	rotate(v.begin(), min_element(v.begin(), v.end()), v.end());
}

// -------------------------------------------------------------- hashing
// For "have I seen this polygon before?" -- bucket polygons in a set / map /
// unordered_map instead of comparing them pairwise. Two flavours:
//
//   polyHash(poly)        the SAME polygon: same vertices in the same cyclic
//                         order, whatever vertex the list happens to start at
//                         and (by default) whichever way round it is listed.
//                         Position matters -- a translate hashes differently.
//   congruenceHash(poly)  the same SHAPE: equal for any translate, rotate
//                         (by ANY angle, not just multiples of 90) and, unless
//                         you pass false, mirror image of the polygon.
//
// Both need INTEGER coordinates -- there is no such thing as a tolerant hash,
// so a Pd whose vertices are off by 1e-12 will land in a different bucket.
// Neither cares whether the polygon is simple or convex.
//
// A hash is a filter, not a proof: on a match compare canonicalPoly(a) ==
// canonicalPoly(b) if a false positive would be fatal (2^-64 per pair).
typedef unsigned long long ull;

// start index of the lexicographically least cyclic rotation of s. O(n),
// and unlike min_element it is well defined when values repeat.  (Booth)
template <class T> int minRotation(vector<T> s) {
	int n = (int)s.size(), a = 0;
	if (n == 0) return 0;
	s.insert(s.end(), s.begin(), s.end());
	for (int b = 0; b < n; b++)
		for (int k = 0; k < n; k++) {
			if (a + k == b || s[a + k] < s[b + k]) { b += max(0, k - 1); break; }
			if (s[a + k] > s[b + k]) { a = b; break; }
		}
	return a;
}

ull mixHash(ull x) {                                            // splitmix64
	x += 0x9e3779b97f4a7c15ULL;
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
	x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
	return x ^ (x >> 31);
}
void hashAdd(ull &h, ll x) { h = mixHash(h ^ mixHash((ull)x)); }   // order matters

// the vertex cycle rotated (and flipped) into its unique smallest listing, so
// that two spellings of one polygon become the same vector -- exact ==
template <class T> vector<Pt<T>> canonicalPoly(vector<Pt<T>> v, bool anyDirection = true) {
	if (v.empty()) return v;
	rotate(v.begin(), v.begin() + minRotation(v), v.end());
	if (anyDirection) {
		vector<Pt<T>> r(v.rbegin(), v.rend());
		rotate(r.begin(), r.begin() + minRotation(r), r.end());
		if (r < v) v.swap(r);
	}
	return v;
}
template <class T> ull polyHash(const vector<Pt<T>> &v, bool anyDirection = true) {
	ull h = mixHash(v.size() * 2 + 1);
	for (auto &p : canonicalPoly(v, anyDirection)) { hashAdd(h, (ll)p.x); hashAdd(h, (ll)p.y); }
	return h;
}

// Shape only. Each vertex becomes (outgoing |edge|^2, cross, dot of the two
// edges at it): the edge lengths and turn angles, which pin the polygon down
// up to a rigid motion and are untouched by translation and rotation. A mirror
// negates every cross, listing it backwards reverses the walk -- so canonicise
// over all four and take the smallest.
template <class T> ull congruenceHash(vector<Pt<T>> v, bool allowReflection = true) {
	int n = (int)v.size();
	ull h = mixHash(n * 2 + 1);
	if (!n) return h;
	auto sig = [&](const vector<Pt<T>> &p) {
		vector<array<ll, 3>> s(n);
		for (int i = 0; i < n; i++) {
			Pt<T> in = p[i] - p[(i + n - 1) % n], out = p[(i + 1) % n] - p[i];
			s[i] = {(ll)len2(out), (ll)cross(in, out), (ll)dot(in, out)};
		}
		rotate(s.begin(), s.begin() + minRotation(s), s.end());
		return s;
	};
	vector<Pt<T>> r(v.rbegin(), v.rend());
	auto best = min(sig(v), sig(r));
	if (allowReflection) {
		for (auto &p : v) p.x = -p.x;
		for (auto &p : r) p.x = -p.x;
		best = min(best, min(sig(v), sig(r)));
	}
	for (auto &t : best) for (ll x : t) hashAdd(h, x);
	return h;
}

// ------------------------------------- every convex polygon over a POINT SET
// Reports each subset of >= 3 of the points that is in STRICTLY convex
// position (no vertex sitting on the line through its two neighbours) exactly
// once, as a CCW list of INDICES into pts starting at the polygon's
// bottom-most (then leftmost) vertex. Other points of the set are allowed to
// lie inside -- counting only the EMPTY convex polygons is a different DP.
// Points must be DISTINCT: run dedup(pts) first, or a repeated point reports
// the same polygon twice.
//
//   forEachConvexPolygon(pts, [&](const vector<int> &poly) { ... });
//
// HOW  Fix the bottom-most vertex b. The other vertices of any convex polygon
// through b appear, in CCW order, sorted by angle around b -- and every one of
// those angles is in [0, pi), so there is no wraparound and a cross product
// sorts them. Grow a chain in that order keeping every turn a left turn; the
// turn at the last vertex and the turn back at b are then automatically left,
// so EVERY chain closes into a polygon and the search never hits a dead end.
//
// TIME  O(n^2 log n) + O(n) per polygon reported -- output sensitive, and the
// output is the bottleneck: n points in convex position give 2^n - O(n^2)
// polygons (60 random points -> 6e7 polygons in ~0.4s). Use the callback and
// consume them as they come; allConvexPolygons stores them all.
template <class T, class F> void forEachConvexPolygon(const vector<Pt<T>> &pts, F f) {
	int n = (int)pts.size();
	vector<int> ord(n);
	for (int i = 0; i < n; i++) ord[i] = i;
	sort(ord.begin(), ord.end(), [&](int a, int b) {           // bottom-most first
		return pts[a].y != pts[b].y ? pts[a].y < pts[b].y : pts[a].x < pts[b].x;
	});
	vector<int> c, nxt, chain;
	for (int t = 0; t + 2 < n; t++) {
		int b = ord[t];
		c.assign(ord.begin() + t + 1, ord.end());              // everything above b
		int m = (int)c.size();
		sort(c.begin(), c.end(), [&](int u, int v) {           // by angle around b
			int s = sgn(cross(pts[b], pts[u], pts[v]));
			return s ? s > 0 : len2(pts[u] - pts[b]) < len2(pts[v] - pts[b]);
		});
		nxt.assign(m, m);                    // first index at a STRICTLY larger angle
		for (int i = m - 2; i >= 0; i--)
			nxt[i] = sgn(cross(pts[b], pts[c[i]], pts[c[i + 1]])) ? i + 1 : nxt[i + 1];
		chain.assign(1, b);
		auto dfs = [&](auto &&self, int pi, int i) -> void {
			chain.push_back(c[i]);
			if (chain.size() >= 3) f((const vector<int> &)chain);
			for (int j = nxt[i]; j < m; j++)   // pi < 0: the turn at c[i] is free
				if (pi < 0 || orient(pts[c[pi]], pts[c[i]], pts[c[j]]) > 0) self(self, i, j);
			chain.pop_back();
		};
		for (int i = 0; i < m; i++) dfs(dfs, -1, i);
	}
}
template <class T> vector<vector<int>> allConvexPolygons(const vector<Pt<T>> &pts) {
	vector<vector<int>> res;
	forEachConvexPolygon(pts, [&](const vector<int> &poly) { res.push_back(poly); });
	return res;
}

// Just how MANY there are, without listing them: same chains, but counted with
// a DP over (previous, current) instead of walked. O(n^4 / 4), n = 150 in
// ~0.02s. The answer is around 2^n, so take it mod something if n > 60.
template <class T> ll countConvexPolygons(const vector<Pt<T>> &pts) {
	int n = (int)pts.size();
	vector<int> ord(n);
	for (int i = 0; i < n; i++) ord[i] = i;
	sort(ord.begin(), ord.end(), [&](int a, int b) {
		return pts[a].y != pts[b].y ? pts[a].y < pts[b].y : pts[a].x < pts[b].x;
	});
	ll tot = 0;
	vector<int> c, nxt;
	for (int t = 0; t + 2 < n; t++) {
		int b = ord[t];
		c.assign(ord.begin() + t + 1, ord.end());
		int m = (int)c.size();
		sort(c.begin(), c.end(), [&](int u, int v) {
			int s = sgn(cross(pts[b], pts[u], pts[v]));
			return s ? s > 0 : len2(pts[u] - pts[b]) < len2(pts[v] - pts[b]);
		});
		nxt.assign(m, m);
		for (int i = m - 2; i >= 0; i--)
			nxt[i] = sgn(cross(pts[b], pts[c[i]], pts[c[i + 1]])) ? i + 1 : nxt[i + 1];
		vector<vector<ll>> dp(m, vector<ll>(m, 0));    // dp[i][j]: chains b..c[i]->c[j]
		for (int j = 0; j < m; j++)
			for (int k = nxt[j]; k < m; k++) {
				ll s = 1;                                  // the triangle b, c[j], c[k]
				for (int i = 0; i < j; i++)
					if (dp[i][j] && orient(pts[c[i]], pts[c[j]], pts[c[k]]) > 0) s += dp[i][j];
				dp[j][k] = s;
				tot += s;                                  // every chain closes
			}
	}
	return tot;
}

// ------------------------ points inside a QUERY polygon, O(k log m) per query
// Static set of m points, plus a static set of n candidate vertices. Each query
// names a polygon as indices into the candidates, and you get back how many of
// the m points are inside it. Any SIMPLE polygon, either orientation -- convex
// is not required. Build O(n m log m) time and O(n m) ints of memory, so this
// wants n, m of a few thousand, not 1e5.
//
//   PolyCount pc(red, black);          // vertices, points to count
//   pc.inside({3, 1, 0, 2});           // strictly inside      O(k log m)
//   pc.onBoundary({3, 1, 0, 2});       // on an edge or vertex O(k log m)
//
// WHY IT IS FAST  Points inside = sum over edges of the signed trapezoid under
// that edge (the shoelace decomposition, counting points instead of area). One
// trapezoid is a slab-clipped half plane, which normally needs heavy machinery
// -- but it splits into a difference of two WEDGES, one apexed at each end of
// the edge:
//     under(a->b) = #{q.x > a.x, q right of ab} - #{q.x > b.x, q right of ab}
// Both conditions in each term are half planes THROUGH THAT VERTEX, so each
// wedge is an angular interval about a candidate vertex. Sort the m points by
// angle around each of the n candidates once, and every wedge is a binary
// search. That is the whole trick: only vertices from a fixed set can be used,
// which is exactly what makes the preprocessing possible.
//
// DEGENERACIES  All handled, exactly: points on an edge, points sitting on a
// vertex, duplicate points, vertical edges, collinear vertices, reflex
// vertices. inside() counts STRICTLY inside; add onBoundary() for the closed
// count. Coordinates must be integers (cross/dot are used exactly).
struct PolyCount {
	vector<P> vs, ps;                  // candidate vertices / the points to count
	vector<vector<int>> ang;           // ps indices, by angle around vs[r] from DOWN
	vector<int> down, same;            // size of the straight-down block / ps == vs[r]

	// angles are measured CCW from straight down, so a wedge "right of the edge
	// and to the right of the apex" is always a prefix-to-prefix range
	static bool hf(P d) { return d.x < 0 || (d.x == 0 && d.y > 0); }
	static bool angBefore(P d, P e) {
		if (hf(d) != hf(e)) return hf(d) < hf(e);
		return cross(d, e) > 0;
	}

	PolyCount(vector<P> verts, vector<P> pts) : vs(verts), ps(pts) {
		int n = (int)vs.size(), m = (int)ps.size();
		ang.resize(n); down.assign(n, 0); same.assign(n, 0);
		for (int r = 0; r < n; r++) {
			P o = vs[r];
			for (int i = 0; i < m; i++) {
				if (ps[i] == o) same[r]++;                  // no angle -- count apart
				else ang[r].push_back(i);
			}
			sort(ang[r].begin(), ang[r].end(), [&](int i, int j) {
				P a = ps[i] - o, b = ps[j] - o;
				if (hf(a) != hf(b)) return hf(a) < hf(b);
				ll c = cross(a, b);
				return c ? c > 0 : len2(a) < len2(b);       // same ray: nearest first
			});
			for (int i : ang[r]) { P d = ps[i] - o; if (d.x == 0 && d.y < 0) down[r]++; }
		}
	}
	int lb(int r, P e) const {                              // first at angle >= e
		const vector<int> &A = ang[r]; P o = vs[r];
		int lo = 0, hi = (int)A.size();
		while (lo < hi) {
			int mid = (lo + hi) / 2;
			if (angBefore(ps[A[mid]] - o, e)) lo = mid + 1; else hi = mid;
		}
		return lo;
	}
	int ub(int r, P e) const {                              // first at angle > e
		const vector<int> &A = ang[r]; P o = vs[r];
		int lo = 0, hi = (int)A.size();
		while (lo < hi) {
			int mid = (lo + hi) / 2;
			if (angBefore(e, ps[A[mid]] - o)) hi = mid; else lo = mid + 1;
		}
		return lo;
	}
	// #points q with (q - vs[r]).x > 0 and q strictly right of the line vs[r]+t*e
	ll wedge(int r, P e) const { return lb(r, e) - down[r]; }
	// #points on the ray vs[r] + t*e with 0 < t <= 1, or 0 < t < 1 if !closed
	ll onRay(int r, P e, bool closed) const {
		const vector<int> &A = ang[r]; P o = vs[r];
		int i = lb(r, e), lo = i, hi = ub(r, e);            // [i,hi): exactly this ray
		ll L = len2(e);
		while (lo < hi) {                                   // sorted by length inside
			int mid = (lo + hi) / 2;
			ll d = len2(ps[A[mid]] - o);
			if (closed ? d <= L : d < L) lo = mid + 1; else hi = mid;
		}
		return lo - i;
	}
	// signed count in the trapezoid under the directed edge a->b, down to
	// y = -inf, over the half-open x range (left, right]
	ll under(int a, int b) const {
		if (vs[a].x == vs[b].x) return 0;                   // vertical: nothing under
		int lo = a, hi = b, sg = 1;
		if (vs[a].x < vs[b].x) { swap(lo, hi); sg = -1; }   // vs[lo].x > vs[hi].x
		P e = vs[lo] - vs[hi];                              // e.x > 0
		return sg * (wedge(hi, e) - wedge(lo, e));
	}

	// The sum counts a point q exactly when q + (-d, eps) is strictly inside, for
	// 0 < d << eps << 1: the half-open x range nudges q left, "strictly below the
	// edge" nudges it up, and up wins. So a point ON the boundary is counted iff
	// that direction u leads into the interior there -- which is what inside()
	// subtracts back off. sideU is the sign of cross(a, u).
	static int sideU(P a) { return a.x ? sgn(a.x) : sgn(a.y); }

	ll onBoundary(const vector<int> &poly) const {          // each point once
		int k = (int)poly.size();
		ll r = 0;
		for (int i = 0; i < k; i++) {
			int a = poly[i], b = poly[(i + 1) % k];
			r += same[a] + onRay(a, vs[b] - vs[a], false);  // vertex a + edge interior
		}
		return r;
	}
	ll inside(vector<int> poly) const {                     // STRICTLY inside
		int k = (int)poly.size();
		if (k < 3) return 0;
		vector<P> pv;
		for (int i : poly) pv.push_back(vs[i]);
		if (polyArea2(pv) < 0) reverse(poly.begin(), poly.end());   // work CCW
		ll s = 0;
		for (int i = 0; i < k; i++) {
			int pr = poly[(i + k - 1) % k], a = poly[i], b = poly[(i + 1) % k];
			s += under(a, b);
			P e = vs[b] - vs[a];
			if (sideU(e) > 0) s -= onRay(a, e, false);      // u enters through this edge
			P A = e, B = vs[pr] - vs[a];                    // interior wedge at vs[a],
			ll c = cross(A, B);                             // CCW from A round to B
			bool in = c > 0 ? (sideU(A) > 0 && sideU(B) < 0)          // convex corner
			        : c < 0 ? (sideU(A) > 0 || sideU(B) < 0)          // reflex corner
			                : sideU(A) > 0;                           // straight
			if (in) s -= same[a];
		}
		return s;
	}
};

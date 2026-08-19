// ============================================================================
//  CIRCLES   -- REQUIRES point.cpp, lines.cpp
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Everything circular. All of it is floating point -- circles leave the
//    lattice the moment you intersect anything -- so read the EPS notes.
//
//  API
//    Circle C{Pd(x,y), r};
//    inCircle(C, p)            0 outside / 1 on / 2 inside
//    circleLine(C, a, b)       0,1,2 points where the INFINITE line ab cuts C
//    circleSeg(C, a, b)        the same, restricted to the segment
//    circleCircle(C1, C2)      0,1,2 points (empty if identical -- check first)
//    tangentPoints(C, p)       the 2 points of C where a line from p is tangent
//    tangentLen(C, p)          length of that tangent
//    commonTangents(C1,C2,inner)  -> list of {touch on C1, touch on C2}
//    circumcircle(a, b, c)     circle through 3 points
//    incircle(a, b, c)         inscribed circle of a triangle
//    minEnclosingCircle(pts)   smallest circle covering every point, O(n) exp.
//    circleInterArea(C1, C2)   area of the lens where two circles overlap
//    circlePolyArea(C, poly)   area of (circle AND polygon)
//    circleUnionArea(cs)       area covered by at least one circle
//    segmentArea(r, th)        area of a circular segment cut by angle th
//    circleFrom2R(a, b, r)     the (up to 2) circles of radius r through a, b
//
//  COMPLEXITY  O(1) each; minEnclosingCircle O(n) expected (shuffled);
//    circlePolyArea O(n); circleUnionArea O(n^2 log n).
//
//  PITFALLS
//    * minEnclosingCircle MUST see a shuffled input -- it shuffles internally,
//      do not remove that, the O(n) bound is entirely due to randomisation.
//    * circleCircle returns {} both for "no intersection" and for "identical
//      circles" (infinitely many). Test for identical yourself if it matters.
//    * Tangent/intersection routines clamp their acos argument, so a point
//      exactly on the circle gives one point rather than NaN. Don't remove the
//      clamps.
//    * circleUnionArea deletes circles contained in another and de-duplicates
//      identical ones; zero-radius circles contribute nothing.
//    * EPS is absolute (1e-9). With coordinates ~1e9 that is far below the
//      representable precision of a double -- scale your input down, or switch
//      ld to long double, before trusting boundary cases.
//    * For "is this lattice point inside the circumcircle of 3 lattice points"
//      use the EXACT inCircleDet predicate at the bottom, not circumcircle().
// ============================================================================

struct Circle {
	Pd c;
	ld r;
	Circle(Pd c = Pd(), ld r = 0) : c(c), r(r) {}
};

ld arcLength(ld r, ld th) { return r * th; }
ld sectorArea(ld r, ld th) { return r * r * th / 2; }
// area cut off by a chord subtending angle th at the centre
ld segmentArea(ld r, ld th) { return r * r * (th - sin(th)) / 2; }

template <class T> int inCircle(Circle C, Pt<T> p) {
	ld d = dist(C.c, toD(p));
	return sgn(C.r - d) < 0 ? 0 : (sgn(C.r - d) == 0 ? 1 : 2);
}

// where the INFINITE line ab meets C, ordered along a -> b
template <class T> vPd circleLine(Circle C, Pt<T> a, Pt<T> b) {
	Pd A = toD(a), B = toD(b);
	Pd f = proj(C.c, A, B);
	ld h2 = C.r * C.r - dist2(C.c, f);
	if (sgn(h2) < 0) return {};
	if (sgn(h2) == 0) return {f};
	Pd v = unit(B - A) * sqrt(h2);
	return {f - v, f + v};
}
// ...restricted to the segment ab, still ordered along a -> b
template <class T> vPd circleSeg(Circle C, Pt<T> a, Pt<T> b) {
	vPd r, hits = circleLine(C, a, b);
	for (auto &p : hits)
		if (sgn(dot(toD(a) - p, toD(b) - p)) <= 0) r.push_back(p);
	return r;
}

// 0, 1 or 2 points. IDENTICAL circles also return {} -- check for that first.
vPd circleCircle(Circle A, Circle B) {
	ld d = dist(A.c, B.c);
	if (sgn(d) == 0) return {};                       // concentric (or identical)
	if (sgn(d - (A.r + B.r)) > 0) return {};          // too far apart
	if (sgn(d - fabs(A.r - B.r)) < 0) return {};      // one strictly inside
	ld x = (d * d + A.r * A.r - B.r * B.r) / (2 * d);
	ld h2 = A.r * A.r - x * x;
	Pd u = (B.c - A.c) / d, mid = A.c + u * x;
	if (sgn(h2) <= 0) return {mid};
	Pd v = perp(u) * sqrt(h2);
	return {mid - v, mid + v};
}

ld tangentLen(Circle C, Pd p) { return sqrt(max((ld)0, dist2(C.c, p) - C.r * C.r)); }
// the two points of C at which a line through p is tangent (p must be outside)
vPd tangentPoints(Circle C, Pd p) {
	ld L = dist2(C.c, p) - C.r * C.r;
	if (sgn(L) < 0) return {};
	return circleCircle(C, Circle(p, sqrt(max((ld)0, L))));
}

// common tangents. inner=false -> the 2 outer ones, inner=true -> the 2 inner
// ones. Each entry is {tangency point on C1, tangency point on C2}; the tangent
// line is the line through those two points.
vector<pair<Pd, Pd>> commonTangents(Circle c1, Circle c2, bool inner) {
	ld r2 = inner ? -c2.r : c2.r;
	Pd d = c2.c - c1.c;
	ld dr = c1.r - r2, d2 = len2(d), h2 = d2 - dr * dr;
	if (sgn(d2) == 0 || sgn(h2) < 0) return {};
	vector<pair<Pd, Pd>> out;
	for (ld s : {(ld)-1, (ld)1}) {
		Pd v = (d * dr + perp(d) * sqrt(max((ld)0, h2)) * s) / d2;
		out.push_back({c1.c + v * c1.r, c2.c + v * r2});
	}
	if (sgn(h2) == 0) out.pop_back();                  // the circles touch
	return out;
}

// circle through three non-collinear points
template <class T> Circle circumcircle(Pt<T> a, Pt<T> b, Pt<T> c) {
	Pd A = toD(a), B = toD(b) - A, C = toD(c) - A;
	ld d = 2 * cross(B, C);
	Pd o(C.y * len2(B) - B.y * len2(C), B.x * len2(C) - C.x * len2(B));
	o = A + o / d;
	return Circle(o, dist(o, A));
}
template <class T> Circle incircle(Pt<T> a, Pt<T> b, Pt<T> c) {
	ld A = dist(b, c), B = dist(a, c), C = dist(a, b), p = A + B + C;
	Pd o = (toD(a) * A + toD(b) * B + toD(c) * C) / p;
	return Circle(o, fabs((ld)cross(a, b, c)) / p);
}

// smallest circle containing every point. O(n) expected -- keep the shuffle.
Circle minEnclosingCircle(vPd p) {
	if (p.empty()) return Circle();
	static mt19937 mecRng(0x9e3779b9);
	shuffle(p.begin(), p.end(), mecRng);
	auto has = [&](Circle c, Pd q) { return dist(c.c, q) <= c.r + 1e-7 * (1 + c.r); };
	auto from2 = [&](Pd a, Pd b) { return Circle((a + b) / 2, dist(a, b) / 2); };
	Circle c(p[0], 0);
	for (int i = 1; i < (int)p.size(); i++)
		if (!has(c, p[i])) {
			c = Circle(p[i], 0);
			for (int j = 0; j < i; j++)
				if (!has(c, p[j])) {
					c = from2(p[i], p[j]);
					for (int k = 0; k < j; k++)
						if (!has(c, p[k])) c = circumcircle(p[i], p[j], p[k]);
				}
		}
	return c;
}

// the (0, 1 or 2) circles of radius r whose boundary passes through a and b
vPd circleFrom2R(Pd a, Pd b, ld r) {          // returns the CENTRES
	ld d2 = dist2(a, b), h2 = r * r - d2 / 4;
	if (sgn(h2) < 0) return {};
	Pd mid = (a + b) / 2;
	if (sgn(h2) == 0) return {mid};
	Pd v = perp(b - a) * sqrt(h2 / d2);
	return {mid - v, mid + v};
}

// ------------------------------------------------------------------- areas
// area of the lens where two circles overlap
ld circleInterArea(Circle A, Circle B) {
	ld d = dist(A.c, B.c);
	if (sgn(d - (A.r + B.r)) >= 0) return 0;
	if (sgn(d - fabs(A.r - B.r)) <= 0) { ld r = min(A.r, B.r); return PI * r * r; }
	ld a1 = acos(max((ld)-1, min((ld)1, (d * d + A.r * A.r - B.r * B.r) / (2 * d * A.r))));
	ld a2 = acos(max((ld)-1, min((ld)1, (d * d + B.r * B.r - A.r * A.r) / (2 * d * B.r))));
	return A.r * A.r * (a1 - sin(2 * a1) / 2) + B.r * B.r * (a2 - sin(2 * a2) / 2);
}

// signed area of (disk of radius r centred at the ORIGIN) AND (triangle O,a,b)
ld _circTri(Pd a, Pd b, ld r) {
	auto arg = [](Pd p, Pd q) { return atan2(cross(p, q), dot(p, q)); };
	bool ina = sgn(len2(a) - r * r) <= 0, inb = sgn(len2(b) - r * r) <= 0;
	if (ina && inb) return cross(a, b) / 2;
	vPd ip = circleSeg(Circle(Pd(0, 0), r), a, b);
	if (ina && !ip.empty()) return cross(a, ip.back()) / 2 + r * r * arg(ip.back(), b) / 2;
	if (inb && !ip.empty()) return r * r * arg(a, ip[0]) / 2 + cross(ip[0], b) / 2;
	if (ip.size() == 2)
		return r * r * arg(a, ip[0]) / 2 + cross(ip[0], ip[1]) / 2 + r * r * arg(ip[1], b) / 2;
	return r * r * arg(a, b) / 2;
}
// area of (circle C) AND (polygon poly, any simple polygon)
template <class T> ld circlePolyArea(Circle C, const vector<Pt<T>> &poly) {
	int n = (int)poly.size();
	ld s = 0;
	for (int i = 0; i < n; i++)
		s += _circTri(toD(poly[i]) - C.c, toD(poly[(i + 1) % n]) - C.c, C.r);
	return fabs(s);
}

// area covered by at least one circle (Green's theorem over the uncovered arcs)
ld circleUnionArea(vector<Circle> cs) {
	int n = (int)cs.size();
	vector<bool> dead(n, false);
	for (int i = 0; i < n; i++) dead[i] = sgn(cs[i].r) <= 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n && !dead[i]; j++) {
			if (i == j || dead[j]) continue;
			ld d = dist(cs[i].c, cs[j].c);
			bool same = sgn(d) == 0 && sgn(cs[i].r - cs[j].r) == 0;
			if (same ? j < i : sgn(d + cs[i].r - cs[j].r) <= 0) dead[i] = true;
		}
	ld total = 0;
	for (int i = 0; i < n; i++) {
		if (dead[i]) continue;
		ld R = cs[i].r;
		vector<pair<ld, ld>> iv;
		for (int j = 0; j < n; j++) {
			if (i == j || dead[j]) continue;
			ld d = dist(cs[i].c, cs[j].c);
			if (sgn(d - (R + cs[j].r)) >= 0) continue;        // disjoint
			if (sgn(d + cs[j].r - R) <= 0) continue;          // j inside i
			ld A = atan2(cs[j].c.y - cs[i].c.y, cs[j].c.x - cs[i].c.x);
			ld co = (d * d + R * R - cs[j].r * cs[j].r) / (2 * d * R);
			ld B = acos(max((ld)-1, min((ld)1, co)));
			ld lo = fmod(A - B, 2 * PI);
			if (lo < 0) lo += 2 * PI;
			ld hi = lo + 2 * B;
			if (hi > 2 * PI) { iv.push_back({lo, 2 * PI}); iv.push_back({0, hi - 2 * PI}); }
			else iv.push_back({lo, hi});
		}
		sort(iv.begin(), iv.end());
		// integrate 1/2 (x dy - y dx) over every arc NOT covered by another circle
		auto arc = [&](ld t1, ld t2) {
			if (t2 <= t1) return;
			total += 0.5 * (cs[i].c.x * R * (sin(t2) - sin(t1)) -
			                cs[i].c.y * R * (cos(t2) - cos(t1)) + R * R * (t2 - t1));
		};
		ld cur = 0;
		for (auto &e : iv) {
			if (e.first > cur) arc(cur, e.first);
			cur = max(cur, e.second);
		}
		arc(cur, 2 * PI);
	}
	return total;
}

// ------------------------------------------------------------------- exact
// EXACT "is d strictly inside the circle through a, b, c" for LATTICE points.
// a, b, c must be counter-clockwise. Returns >0 inside, 0 on, <0 outside.
// This is the Delaunay predicate; it needs __int128 because it is a 4x4-ish
// determinant of squared coordinates (values up to ~C^4).
__int128 inCircleDet(P a, P b, P c, P d) {
	auto f = [&](P p) { return (__int128)(p.x - d.x); };
	auto g = [&](P p) { return (__int128)(p.y - d.y); };
	auto h = [&](P p) { return f(p) * f(p) + g(p) * g(p); };
	return f(a) * (g(b) * h(c) - g(c) * h(b)) - g(a) * (f(b) * h(c) - f(c) * h(b)) +
	       h(a) * (f(b) * g(c) - f(c) * g(b));
}

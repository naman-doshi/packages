// ============================================================================
//  HALF-PLANE INTERSECTION   -- REQUIRES point.cpp
//                               (add polygon.cpp if you want polyArea of the result)
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Intersect any number of half-planes ("the region left of this directed
//    line") in O(n log n). This is the general tool for:
//      * feasible region of a system of linear inequalities a*x + b*y <= c
//      * the visible / lit region of a convex arrangement
//      * "the area a robot can reach given these walls"
//      * kernel of a polygon (the set of points that see the whole polygon)
//      * intersection of many convex polygons -- feed in all their edges
//
//  A half-plane is written as a directed line: HP(a, b) means "the closed side
//  to the LEFT of a -> b". If your input is `a*x + b*y <= c`, the inward normal
//  is (-a,-b), so a direction vector along the boundary is (b, -a):
//      pick any point p on the line, then HP(p, p + Pd(b, -a)).
//
//  API
//    vector<HP> h = { HP(Pd(0,0), Pd(1,0)), ... };
//    vPd region = halfplaneInter(h);          // CCW polygon, empty if infeasible
//    vPd region = halfplaneInter(h, 1e9);     // custom bounding box half-width
//    vPd k = polygonKernel(poly);             // kernel of a simple polygon
//
//  COMPLEXITY  O(n log n) -- one sort by angle, then a deque sweep.
//
//  PITFALLS
//    * An UNBOUNDED feasible region is clipped to the bounding box (default
//      +-1e9). If your answer touches the box, the true region was unbounded.
//      Set BOX just above your coordinate range so that is easy to detect.
//    * The result is EMPTY for an infeasible system AND for a region that is a
//      single point or a segment (zero area). Check size() < 3.
//    * All floating point. Two nearly-parallel half-planes produce a far-away
//      intersection point; EPS handles the usual cases but do not push it.
//    * Duplicate / parallel half-planes are fine, they get de-duplicated.
//    * The half-planes are CLOSED (the boundary line is included).
// ============================================================================

struct HP {
	Pd p, pq;            // the boundary is the line p -> p+pq; keep the LEFT side
	ld ang;
	HP() {}
	HP(Pd a, Pd b) : p(a), pq(b - a) { ang = atan2(pq.y, pq.x); }
	bool out(Pd r) const { return cross(pq, r - p) < -EPS; }   // r strictly outside
	bool operator<(const HP &o) const { return ang < o.ang; }
};
Pd hpInter(const HP &s, const HP &t) {
	ld a = cross(t.p - s.p, t.pq) / cross(s.pq, t.pq);
	return s.p + s.pq * a;
}

// intersection of all the half-planes, as a CCW polygon. Empty = infeasible or
// degenerate (a point / a segment).
vPd halfplaneInter(vector<HP> h, ld BOX = 1e9) {
	Pd box[4] = {Pd(BOX, BOX), Pd(-BOX, BOX), Pd(-BOX, -BOX), Pd(BOX, -BOX)};
	for (int i = 0; i < 4; i++) h.push_back(HP(box[i], box[(i + 1) % 4]));
	sort(h.begin(), h.end());
	deque<HP> dq;
	for (int i = 0; i < (int)h.size(); i++) {
		while (dq.size() > 1 && h[i].out(hpInter(dq[dq.size() - 1], dq[dq.size() - 2])))
			dq.pop_back();
		while (dq.size() > 1 && h[i].out(hpInter(dq[0], dq[1]))) dq.pop_front();
		if (!dq.empty() && fabs(cross(h[i].pq, dq.back().pq)) < EPS) {
			if (dot(h[i].pq, dq.back().pq) < 0) return {};      // opposite: empty
			if (h[i].out(dq.back().p)) dq.pop_back();           // keep the tighter one
			else continue;
		}
		dq.push_back(h[i]);
	}
	while (dq.size() > 2 && dq[0].out(hpInter(dq[dq.size() - 1], dq[dq.size() - 2])))
		dq.pop_back();
	while (dq.size() > 2 && dq.back().out(hpInter(dq[0], dq[1]))) dq.pop_front();
	if (dq.size() < 3) return {};
	vPd res(dq.size());
	for (size_t i = 0; i + 1 < dq.size(); i++) res[i] = hpInter(dq[i], dq[i + 1]);
	res.back() = hpInter(dq.back(), dq[0]);
	return res;
}

// build a half-plane for the inequality  A*x + B*y <= C
HP hpFromIneq(ld A, ld B, ld C) {
	Pd n(A, B);                              // outward normal
	Pd p = n * (C / len2(n));                // a point on the boundary line
	// the kept side must be on the LEFT of the direction, and "left of d" is
	// perp(d) = (-d.y, d.x); we need that to point inward, i.e. to equal -n,
	// which gives d = (-B, A).
	return HP(p, p + Pd(-B, A));
}

// the kernel of a simple polygon: every point that can see the WHOLE polygon.
// Non-empty exactly when the polygon is star-shaped. poly must be CCW.
template <class T> vPd polygonKernel(const vector<Pt<T>> &poly) {
	vector<HP> h;
	int n = (int)poly.size();
	for (int i = 0; i < n; i++) h.push_back(HP(toD(poly[i]), toD(poly[(i + 1) % n])));
	return halfplaneInter(h);
}

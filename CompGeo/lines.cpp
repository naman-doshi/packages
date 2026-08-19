// ============================================================================
//  LINES & SEGMENTS   -- REQUIRES point.cpp
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Everything about infinite lines, rays and segments: do they cross, where,
//    how far apart, project / reflect a point. All the predicates (intersects?
//    which side? on the segment?) are EXACT when you use P = Pt<ll>; only the
//    functions that must return a non-lattice point go through double.
//
//  API   (a,b) is one line/segment, (c,d) the other. The predicates are
//        templated and EXACT on P; the ones that return a distance or a
//        computed point take Pd -- and a P converts to a Pd automatically, so
//        `distToLine(someCircle.c, latticeA, latticeB)` just works.
//    sideOf(a, b, p)             +1 left of a->b, -1 right, 0 on the line
//    lineInter(a,b, c,d)         -> {0 parallel, 1 unique, -1 same line}, Pd hit
//    segInter(a,b, c,d)          bool: do the two SEGMENTS share a point
//    segInterPts(a,b, c,d)       0, 1 or 2 points (2 = collinear overlap ends)
//    segInterProper(a,b, c,d)    bool: cross at a single interior point
//    proj(p, a,b)                foot of perpendicular from p onto line ab
//    refl(p, a,b)                mirror image of p across line ab
//    closestOnSeg(p, a,b)        nearest point of segment ab to p
//    distToLine(p, a,b)          distance point <-> infinite line
//    distToSeg(p, a,b)           distance point <-> segment
//    segDist(a,b, c,d)           distance segment <-> segment
//    linesParallel / linesSame / perpendicular
//    perpBisector(a,b)           -> {point, direction}
//    angleBisector(o, a, b)      direction of the interior bisector at o
//    lineFromPts / evalABC / interABC     for the ax + by = c form
//    lineDist(a,b, c,d)          distance between two PARALLEL lines
//    ptOnLineAt(a, b, t)         a + t*(b-a)
//    segRect(a, b, lo, hi)       does segment ab touch the axis-aligned box
//
//  A POINT SET against a line
//    countOnLine(pts, a, b)      how many of pts lie on the LINE ab      O(n)
//    countOnSeg(pts, a, b)       ... on the SEGMENT ab                   O(n)
//    pointsOnLine(pts, a, b)     their indices, ordered along a->b       O(n log n)
//    lineKey(a, b)               canonical {A,B,C} of Ax+By=C -- equal iff same line
//    OnLine tab(pts);            many query lines against a fixed point set
//      tab.count(a, b)           how many points on the line, O(log n)
//      tab.any(a, b)             is there at least one,       O(log n)
//    LineSweep sw(pts);          the same, with a worst-case bound, and
//      sw.count / sw.any                                      O(log^2 n)
//      sw.countLeft(a, b)        half-plane count: points strictly left of a->b
//      sw.countRight(a, b)                        ... strictly right
//
//  WHICH ONE   OnLine is smaller and answers in O(log n), but only lines
//    holding >= 2 points are tabulated; the 0-or-1 case is a binary search in a
//    per-direction index it builds on demand, O(n log n) the first time it sees
//    a direction. That is free when the queries reuse directions. LineSweep has
//    no such hole -- every query is O(log^2 n) whatever the direction -- and it
//    also does half-plane counting. Both build in O(n^2 log n), so if the
//    queries are adversarial, just use LineSweep.
//
//  COMPLEXITY  everything is O(1) except the point-set block above:
//    countOnLine / countOnSeg O(n), pointsOnLine O(n log n).
//
//  PITFALLS
//    * lineInter needs a real line: a != b and c != d. Degenerate "segments"
//      that are single points make every orientation test 0 -- guard for them.
//    * distToSeg and closestOnSeg handle a == b fine (returns dist to the point).
//    * segInter treats TOUCHING as intersecting (shared endpoint -> true). If
//      the problem wants strictly-crossing, use segInterProper.
//    * Collinear overlapping segments: segInter -> true, segInterPts -> the two
//      ends of the shared piece (which may be equal if they touch at one point).
//    * Overflow: sideOf on coordinates up to 1e9 computes up to ~4e18. Safe but
//      tight; beyond that switch Pt's T to __int128.
//    * For an EXACT intersection point of integer segments use Math/frac.cpp:
//      the coordinates are rationals with denominator cross(b-a, d-c).
//    * countOnLine / lineKey / OnLine all need a REAL line: a != b. With a == b
//      every orientation test is 0, so countOnLine happily returns n and
//      lineKey divides by zero.
//    * lineKey computes A*a.x + B*a.y with A,B up to 2*C -- ~4e18 at C = 1e9.
//      Fits, but shift the coordinates or switch to __int128 past that. Same
//      bound for the cross(dir, p) that OnLine::any binary searches.
//    * OnLine::any / count keep one array of n longs per DISTINCT query
//      direction. That is the point (repeat directions are free), but q
//      all-different directions cost O(q n) memory -- clear `proj` if that
//      bites, or use LineSweep, which has no cache at all.
//    * LineSweep costs ~24 bytes per PAIR while building, so n <= ~2000. Its
//      countLeft is STRICT (orient > 0); "on or left" is countLeft + count.
//    * Looking for the BEST line rather than a given one? That is
//      maxPointsOnLine in angular.cpp -- don't build OnLine for it.
// ============================================================================

// which side of the directed line a->b does p lie on: +1 left, -1 right, 0 on
template <class T> int sideOf(Pt<T> a, Pt<T> b, Pt<T> p) { return orient(a, b, p); }

template <class T> bool linesParallel(Pt<T> a, Pt<T> b, Pt<T> c, Pt<T> d) {
	return sgn(cross(b - a, d - c)) == 0;
}
template <class T> bool linesSame(Pt<T> a, Pt<T> b, Pt<T> c, Pt<T> d) {
	return linesParallel(a, b, c, d) && orient(a, b, c) == 0;
}
template <class T> bool perpendicular(Pt<T> a, Pt<T> b, Pt<T> c, Pt<T> d) {
	return sgn(dot(b - a, d - c)) == 0;
}

// INFINITE line ab  x  infinite line cd.
//   first =  1 : unique intersection, returned in second
//   first =  0 : parallel and distinct (second is garbage)
//   first = -1 : the same line       (second is garbage)
pair<int, Pd> lineInter(Pd a, Pd b, Pd c, Pd d) {
	ld den = cross(b - a, d - c);
	if (sgn(den) == 0) return {orient(a, b, c) == 0 ? -1 : 0, Pd()};
	return {1, a + (b - a) * (cross(c - a, d - c) / den)};
}

Pd ptOnLineAt(Pd a, Pd b, ld t) { return a + (b - a) * t; }   // a + t*(b-a)

// do segments [a,b] and [c,d] share at least one point? (touching counts)
template <class T> bool segInter(Pt<T> a, Pt<T> b, Pt<T> c, Pt<T> d) {
	int d1 = orient(a, b, c), d2 = orient(a, b, d);
	int d3 = orient(c, d, a), d4 = orient(c, d, b);
	if (d1 * d2 < 0 && d3 * d4 < 0) return true;
	if (d1 == 0 && onSegment(c, a, b)) return true;
	if (d2 == 0 && onSegment(d, a, b)) return true;
	if (d3 == 0 && onSegment(a, c, d)) return true;
	if (d4 == 0 && onSegment(b, c, d)) return true;
	return false;
}
// strictly crossing: exactly one common point, interior to both segments
template <class T> bool segInterProper(Pt<T> a, Pt<T> b, Pt<T> c, Pt<T> d) {
	int d1 = orient(a, b, c), d2 = orient(a, b, d);
	int d3 = orient(c, d, a), d4 = orient(c, d, b);
	return d1 * d2 < 0 && d3 * d4 < 0;
}

// the shared point set of two segments: {} / {p} / {p,q} (collinear overlap)
template <class T> vPd segInterPts(Pt<T> a, Pt<T> b, Pt<T> c, Pt<T> d) {
	T den = cross(b - a, d - c);
	if (sgn(den) != 0) {                      // not parallel: at most one point
		if (!segInter(a, b, c, d)) return {};
		return {lineInter(a, b, c, d).second};
	}
	if (orient(a, b, c) != 0) return {};      // parallel, distinct
	// collinear: intersect the two intervals along the common direction
	vector<Pt<T>> got;
	if (onSegment(c, a, b)) got.push_back(c);
	if (onSegment(d, a, b)) got.push_back(d);
	if (onSegment(a, c, d)) got.push_back(a);
	if (onSegment(b, c, d)) got.push_back(b);
	if (got.empty()) return {};
	Pt<T> lo = got[0], hi = got[0];
	for (auto &p : got) { lo = min(lo, p); hi = max(hi, p); }
	if (lo == hi) return {toD(lo)};
	return {toD(lo), toD(hi)};
}

// ------------------------------------------------ projections and distances
// These take Pd, and a P converts to a Pd on its own, so you can freely mix
// lattice points and computed points:  distToLine(circle.c, latticeA, latticeB)
Pd proj(Pd p, Pd a, Pd b) {                  // foot of the perpendicular onto ab
	Pd v = b - a;
	return a + v * (dot(p - a, v) / len2(v));
}
Pd refl(Pd p, Pd a, Pd b) { return proj(p, a, b) * 2.0 - p; }   // mirror across ab

// nearest point of the SEGMENT ab to p (handles a == b)
Pd closestOnSeg(Pd p, Pd a, Pd b) {
	Pd v = b - a;
	if (sgn(len2(v)) == 0) return a;
	ld t = max((ld)0, min((ld)1, dot(p - a, v) / len2(v)));
	return a + v * t;
}
ld distToLine(Pd p, Pd a, Pd b) { return fabs(cross(a, b, p)) / len(b - a); }
// signed: positive when p is left of a->b
ld signedDistToLine(Pd p, Pd a, Pd b) { return cross(a, b, p) / len(b - a); }
ld distToSeg(Pd p, Pd a, Pd b) { return dist(p, closestOnSeg(p, a, b)); }
// distance between two PARALLEL lines (assumes they really are parallel, so the
// 4th point is along for the ride and never read)
ld lineDist(Pd a, Pd b, Pd c, Pd /*d*/) { return distToLine(c, a, b); }

// segment <-> segment. Templated so the "do they cross" half stays EXACT on
// lattice input; only the fallback distance goes through double.
template <class T> ld segDist(Pt<T> a, Pt<T> b, Pt<T> c, Pt<T> d) {
	if (segInter(a, b, c, d)) return 0;
	return min({distToSeg(a, c, d), distToSeg(b, c, d),
	            distToSeg(c, a, b), distToSeg(d, a, b)});
}

// ------------------------------------------------------------------ bisectors
// perpendicular bisector of ab, as {a point on it, a direction vector}
pair<Pd, Pd> perpBisector(Pd a, Pd b) { return {(a + b) * 0.5, perp(b - a)}; }
// direction of the interior angle bisector at o in the wedge o->a, o->b
Pd angleBisector(Pd o, Pd a, Pd b) { return unit(a - o) + unit(b - o); }

// ------------------------------------------------------------- ax + by = c
// Many editorials state lines as ax + by = c. Convert both ways.
struct ABC { ld a, b, c; };
ABC lineFromPts(Pd p, Pd q) {                                   // through p and q
	return ABC{q.y - p.y, p.x - q.x, (q.y - p.y) * p.x + (p.x - q.x) * p.y};
}
ld evalABC(ABC L, Pd p) { return L.a * p.x + L.b * p.y - L.c; }  // 0 on the line
// intersection of two ax+by=c lines; .first = false when parallel
pair<bool, Pd> interABC(ABC L, ABC M) {
	ld den = L.a * M.b - L.b * M.a;
	if (sgn(den) == 0) return {false, Pd()};
	return {true, Pd((L.c * M.b - L.b * M.c) / den, (L.a * M.c - L.c * M.a) / den)};
}

// -------------------------------------------------- a point set and a line
// How many of `pts` lie on the INFINITE line through a and b? Exact for P.
template <class T> int countOnLine(const vector<Pt<T>> &pts, Pt<T> a, Pt<T> b) {
	int c = 0;
	for (auto &p : pts) c += orient(a, b, p) == 0;
	return c;
}
// ... and on the SEGMENT ab only (endpoints count)
template <class T> int countOnSeg(const vector<Pt<T>> &pts, Pt<T> a, Pt<T> b) {
	int c = 0;
	for (auto &p : pts) c += onSegment(p, a, b);
	return c;
}
// indices of the points on the line, ordered along a -> b
template <class T> vector<int> pointsOnLine(const vector<Pt<T>> &pts, Pt<T> a, Pt<T> b) {
	vector<int> id;
	for (int i = 0; i < (int)pts.size(); i++)
		if (orient(a, b, pts[i]) == 0) id.push_back(i);
	Pt<T> d = b - a;
	sort(id.begin(), id.end(), [&](int i, int j) {
		return dot(d, pts[i] - a) < dot(d, pts[j] - a);
	});
	return id;
}

// Canonical integer form of the line through a != b, as Ax + By = C reduced by
// gcd with a fixed sign: two keys are equal IFF the two lines are identical.
// (gcd(|A|,|B|) always divides C, so all three stay integers -- exact, no slopes.)
typedef array<ll, 3> LineKey;
LineKey lineKey(P a, P b) {
	ll A = b.y - a.y, B = a.x - b.x, C = A * a.x + B * a.y;
	ll g = gcd(llabs(A), llabs(B));                       // > 0 since a != b
	A /= g; B /= g; C /= g;
	if (A < 0 || (A == 0 && B < 0)) A = -A, B = -B, C = -C;
	return {A, B, C};
}

// MANY query lines against one FIXED point set. Everything here is EXACT.
//   OnLine tab(pts);          build O(n^2 log n), O(#lines) <= O(n^2) memory
//   tab.count(a, b)           how many points are on the line
//   tab.any(a, b)             is there at least one
//
// Two indexes, because "2 or more" and "1 or more" are different problems:
//   * `tab`, every line holding >= 2 of the points, keyed by lineKey.
//     O(log n), always, for any query line. A hit already answers any() = true.
//   * `proj`, built lazily PER QUERY DIRECTION: the points sorted by their
//     cross product with that direction, so a line of that direction is a
//     binary search for one value. First query of a direction costs
//     O(n log n), every later one O(log n).
// So a query is O(log n) unless it is a brand-new direction whose line holds
// <= 1 point. Reusing directions (one sweep angle, all vertical lines, lines
// through a fixed pencil...) is the case that has to be fast, and it is.
// There is no O(log n) worst case for arbitrary directions without half-plane
// range counting (dual arrangement, O(n^2) space) -- not worth it here.
// Duplicate points in `pts` are all counted, same as countOnLine.
struct OnLine {
	vP pts;                                   // kept for the per-direction build
	vector<pair<LineKey, int>> tab;           // sorted; lines with >= 2 points
	mutable map<P, vector<ll>> proj;          // direction -> sorted cross(dir, p)
	OnLine(const vP &p) : pts(p) {
		vP u = p; dedup(u);                   // u sorted & unique, w = multiplicity
		int m = (int)u.size();
		vector<int> w(m, 0);
		for (auto &q : p) w[lower_bound(u.begin(), u.end(), q) - u.begin()]++;
		// every line with >= 2 points is seen at its lowest-index point, where
		// the whole rest of it shows up as one direction bucket
		vector<pair<P, int>> dir;
		for (int i = 0; i < m; i++) {
			dir.clear();
			for (int j = i + 1; j < m; j++) {
				P d = u[j] - u[i];
				ll g = gcd(llabs(d.x), llabs(d.y));
				d.x /= g; d.y /= g;
				if (d.x < 0 || (d.x == 0 && d.y < 0)) d = -d;   // canonical
				dir.push_back({d, w[j]});
			}
			sort(dir.begin(), dir.end(),
			     [](const pair<P, int> &x, const pair<P, int> &y) { return x.first < y.first; });
			for (int s = 0, e; s < (int)dir.size(); s = e) {
				int tot = 0;
				for (e = s; e < (int)dir.size() && dir[e].first == dir[s].first; e++)
					tot += dir[e].second;
				tab.push_back({lineKey(u[i], u[i] + dir[s].first), w[i] + tot});
			}
		}
		sort(tab.begin(), tab.end());
		// same line reached from a later point: keep the largest (= complete) count
		int k = 0;
		for (int i = 0; i < (int)tab.size(); i++)
			if (k && tab[k - 1].first == tab[i].first)
				tab[k - 1].second = max(tab[k - 1].second, tab[i].second);
			else tab[k++] = tab[i];
		tab.resize(k);
	}
	// canonical primitive direction of a -> b; p is on the line iff
	// cross(dir, p) == cross(dir, a), which is what `proj` is sorted by
	static P dirOf(P a, P b) {
		P d = b - a;
		ll g = gcd(llabs(d.x), llabs(d.y));                   // > 0 since a != b
		d.x /= g; d.y /= g;
		if (d.x < 0 || (d.x == 0 && d.y < 0)) d = -d;
		return d;
	}
	const vector<ll> &projOf(P d) const {                     // memoised, O(n log n) once
		auto it = proj.find(d);
		if (it != proj.end()) return it->second;
		vector<ll> &v = proj[d];
		v.reserve(pts.size());
		for (auto &p : pts) v.push_back(cross(d, p));
		sort(v.begin(), v.end());
		return v;
	}
	// the fast half of a query: -1 = "not tabulated, 0 or 1 points", else the count
	int tabbed(P a, P b) const {
		LineKey k = lineKey(a, b);
		auto it = lower_bound(tab.begin(), tab.end(), k,
		                      [](const pair<LineKey, int> &e, const LineKey &v) { return e.first < v; });
		return it != tab.end() && it->first == k ? it->second : -1;
	}
	int count(P a, P b) const {
		int c = tabbed(a, b);
		if (c >= 0) return c;
		P d = dirOf(a, b);
		const vector<ll> &v = projOf(d);
		ll t = cross(d, a);
		return (int)(upper_bound(v.begin(), v.end(), t) - lower_bound(v.begin(), v.end(), t));
	}
	bool any(P a, P b) const {
		if (tabbed(a, b) >= 0) return true;                   // >= 2 is >= 1
		P d = dirOf(a, b);
		const vector<ll> &v = projOf(d);
		ll t = cross(d, a);
		return binary_search(v.begin(), v.end(), t);
	}
};

// SAME QUESTIONS, WORST-CASE BOUND. Every query is O(log^2 n) whatever its
// direction -- no cache to miss -- and half-plane counting comes with it.
//   LineSweep sw(pts);        build O(n^2 log n) time, O(n^2) memory
//   sw.count(a, b)            points ON the line ab
//   sw.any(a, b)              is there at least one
//   sw.countLeft(a, b)        points strictly LEFT of the directed line a->b
//   sw.countRight(a, b)       ... strictly right      (left + right + count = n)
//
// HOW. Fix a direction d and sort the points by cross(d, p): p is on the line
// iff cross(d, p) == cross(d, a), so the query is a lookup in a sorted array.
// Rotate d over a half turn and that order only changes at the O(n^2)
// directions where two points tie, each change swapping neighbours -- so
// recording "rank k holds this point from this direction on" is O(n^2) in
// total, not the O(n^3) of storing every order. A query binary searches the
// events for its own direction, then binary searches the ranks, whose
// projections are sorted by construction. Two binary searches -> O(log^2 n).
// (Fractional cascading over the ranks would make it O(log n). Not worth it.)
//
// Exact, duplicate-tolerant, and it needs no general position. Costs about
// 24 bytes per pair while building, so n <= ~2000.
struct LineSweep {
	int n;
	vP pts;
	vector<P> ev;                             // event directions, in angle order
	vector<vector<pair<int, int>>> at;        // at[k] = {(from which state, point)}
	// is dir x at a strictly smaller angle than y? both folded into [0, pi)
	static bool angLess(P x, P y) { return cross(x, y) > 0; }
	// primitive direction of a->b folded into [0, pi), plus the sign that fold
	// cost: b - a == s * g * d, so orient(a, b, p) has the sign of s*(val - t)
	static pair<P, int> canon(P a, P b) {
		P d = b - a;
		int s = 1;
		ll g = gcd(llabs(d.x), llabs(d.y));                   // > 0 since a != b
		d.x /= g; d.y /= g;
		if (d.y < 0 || (d.y == 0 && d.x < 0)) { d = -d; s = -1; }
		return {d, s};
	}
	LineSweep(const vP &p) : n((int)p.size()), pts(p) {
		vector<int> ord(n), pos(n);
		iota(ord.begin(), ord.end(), 0);
		// state 0 is the direction (1,0), where cross(d, p) is just p.y. The tie
		// break by x matters: the first event group reverses the equal-y blocks,
		// and just past (1,0) they have to come out ordered by DECREASING x.
		sort(ord.begin(), ord.end(), [&](int i, int j) {
			return pts[i].y != pts[j].y ? pts[i].y < pts[j].y : pts[i].x < pts[j].x;
		});
		at.assign(n, {});
		for (int k = 0; k < n; k++) { pos[ord[k]] = k; at[k].push_back({0, ord[k]}); }
		vector<pair<P, pair<int, int>>> evp;                  // (tie direction, pair)
		for (int i = 0; i < n; i++)
			for (int j = i + 1; j < n; j++) {
				if (pts[i] == pts[j]) continue;               // never swap, never ties
				P d = pts[i] - pts[j];
				ll g = gcd(llabs(d.x), llabs(d.y));
				d.x /= g; d.y /= g;
				if (d.y < 0 || (d.y == 0 && d.x < 0)) d = -d;
				evp.push_back({d, {i, j}});
			}
		sort(evp.begin(), evp.end(), [](const pair<P, pair<int, int>> &x,
		                                const pair<P, pair<int, int>> &y) {
			return x.first == y.first ? x.second < y.second : angLess(x.first, y.first);
		});
		vector<int> who;
		for (int s = 0, e; s < (int)evp.size(); s = e) {
			P d = evp[s].first;
			for (e = s; e < (int)evp.size() && evp[e].first == d; e++) ;
			// every point that ties with another at this direction, by rank.
			// Points sharing a cross value form ONE contiguous block of ranks
			// (anything ranked between them ties too), and rotating past d
			// reverses each such block -- that is the whole state change.
			who.clear();
			for (int q = s; q < e; q++) { who.push_back(evp[q].second.first); who.push_back(evp[q].second.second); }
			sort(who.begin(), who.end(), [&](int i, int j) { return pos[i] < pos[j]; });
			who.erase(unique(who.begin(), who.end()), who.end());
			int g = (int)ev.size();                           // state g+1 starts after d
			for (int x = 0, y; x < (int)who.size(); x = y) {
				ll v = cross(d, pts[who[x]]);
				for (y = x; y < (int)who.size() && cross(d, pts[who[y]]) == v; y++) ;
				int lo = pos[who[x]], hi = pos[who[y - 1]];
				reverse(ord.begin() + lo, ord.begin() + hi + 1);
				for (int k = lo; k <= hi; k++) {
					pos[ord[k]] = k;
					if (at[k].back().second != ord[k]) at[k].push_back({g + 1, ord[k]});
				}
			}
			ev.push_back(d);
		}
	}
	// which state is in force at direction d. Landing exactly ON an event is
	// fine: the points it swaps are tied there, so the projections agree.
	int stateAt(P d) const { return (int)(lower_bound(ev.begin(), ev.end(), d, angLess) - ev.begin()); }
	ll valAt(int st, int k, P d) const {                      // projection at rank k
		const vector<pair<int, int>> &v = at[k];
		int i = (int)(upper_bound(v.begin(), v.end(), make_pair(st, INT_MAX)) - v.begin()) - 1;
		return cross(d, pts[v[i].second]);
	}
	// first rank whose projection is >= t (strict: > t). Projections rise with k.
	int rankOf(int st, P d, ll t, bool strict) const {
		int lo = 0, hi = n;
		while (lo < hi) {
			int mid = (lo + hi) / 2;
			ll v = valAt(st, mid, d);
			if (v > t || (!strict && v == t)) hi = mid; else lo = mid + 1;
		}
		return lo;
	}
	int count(P a, P b) const {
		pair<P, int> c = canon(a, b);
		int st = stateAt(c.first);
		ll t = cross(c.first, a);
		return rankOf(st, c.first, t, true) - rankOf(st, c.first, t, false);
	}
	bool any(P a, P b) const {
		pair<P, int> c = canon(a, b);
		int st = stateAt(c.first);
		ll t = cross(c.first, a);
		int k = rankOf(st, c.first, t, false);
		return k < n && valAt(st, k, c.first) == t;
	}
	// strictly left / right of the DIRECTED line a -> b (orient(a,b,p) = +1/-1)
	int countLeft(P a, P b) const {
		pair<P, int> c = canon(a, b);
		int st = stateAt(c.first);
		ll t = cross(c.first, a);
		return c.second > 0 ? n - rankOf(st, c.first, t, true) : rankOf(st, c.first, t, false);
	}
	int countRight(P a, P b) const { return countLeft(b, a); }
};

// -------------------------------------------------------------------- boxes
// does segment ab touch the axis-aligned box [lo.x,hi.x] x [lo.y,hi.y]?
template <class T> bool segRect(Pt<T> a, Pt<T> b, Pt<T> lo, Pt<T> hi) {
	auto inBox = [&](Pt<T> p) {
		return p.x >= lo.x && p.x <= hi.x && p.y >= lo.y && p.y <= hi.y;
	};
	if (inBox(a) || inBox(b)) return true;
	Pt<T> c1(lo.x, lo.y), c2(hi.x, lo.y), c3(hi.x, hi.y), c4(lo.x, hi.y);
	return segInter(a, b, c1, c2) || segInter(a, b, c2, c3) ||
	       segInter(a, b, c3, c4) || segInter(a, b, c4, c1);
}

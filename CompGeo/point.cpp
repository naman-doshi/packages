// ============================================================================
//  POINT  -- the base file. PASTE THIS FIRST, every other CompGeo file needs it
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    2-D vector type + the four primitives everything else is built from:
//    cross (turn direction / signed area), dot (projection / angle sign),
//    orient (the ONLY correct way to ask "left or right"), and polar sort.
//
//    Two instantiations are provided:
//      P  = Pt<ll>      EXACT. Use this whenever the input coordinates are
//                       integers -- which is most of the time. No epsilons,
//                       no wrong answers.
//      Pd = Pt<double>  Only when the geometry genuinely leaves the lattice
//                       (circles, intersection points, angles, rotations).
//    Algorithms are templated, so `convexHull(v)` works on either.
//
//  API
//    P a(3, 4), b(1, 2);
//    a + b   a - b   a * 3   a / 2   -a   a == b   a < b   cin >> a   cout << a
//    dot(a,b)  cross(a,b)  cross(o,a,b)          // cross(o,a,b) = cross(a-o,b-o)
//    len2(a)   len(a)   dist2(a,b)   dist(a,b)   // len2/dist2 are EXACT for P
//    orient(a,b,c)      // +1 = c is LEFT of a->b (CCW), -1 = right, 0 = collinear
//    collinear(a,b,c)   onSegment(p,a,b)         // exact for P
//    perp(a)            // rotate 90 CCW      perpCW(a) = 90 CW
//    rot(p, th)  rotAround(p, o, th)  unit(p)  scaleTo(p, len)   // Pd only
//    angle(p)           // atan2 in (-pi, pi]
//    angleBetween(a,b)  // unsigned, [0, pi]
//    signedAngle(a,b)   // (-pi, pi], + if b is CCW from a
//    polarSort(v)  polarSortAround(o, v)         // CCW from +x axis
//    latticeOnSeg(a,b)  // # lattice points strictly between two lattice points
//    area2(a,b,c)       // 2 * signed triangle area (EXACT integer)
//
//  COMPLEXITY  everything O(1) except polarSort, O(n log n).
//
//  PITFALLS
//    * OVERFLOW is the #1 killer. cross() of coordinates up to C computes ~2*C^2.
//      |C| <= 1e9 -> 2e18, that fits in ll but ONLY just. If you subtract a point
//      first (cross(o,a,b) with |o|,|a|,|b| <= 1e9) the differences are up to 2e9
//      and the cross is up to 8e18 -- OVERFLOW. Shift coordinates to be small, or
//      use __int128 (typedef T = __int128 in Pt) when |C| > ~1e9.
//    * NEVER compare doubles with ==. `sgn()` does it with EPS; use orient(),
//      not `cross(...) > 0`, so the EPS is applied for you.
//    * == and < are EXACT tie-compares even for Pd (so sort() stays a valid
//      strict weak order). For tolerant equality on doubles use `near(a,b)`.
//    * P::operator/ on integers is INTEGER division. Convert to Pd first.
//    * atan2 sorting is slower and can misorder near-equal angles. Use
//      polarSort (cross-product based), not sort-by-atan2.
//    * These files do NOT `#define int long long`, so they behave the same
//      whether or not your template has that macro on.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef double ld;                 // -> long double if you need the precision

const ld EPS = 1e-9;
const ld PI = acos((ld)-1.0);

// sign with tolerance for floats, exact for integers. Written as one template so
// it does not break under a `#define int long long` template.
template <class T> int sgn(T x) {
	if constexpr (is_integral_v<T>) return (x > 0) - (x < 0);
	else return (x > EPS) - (x < -EPS);
}

template <class T> struct Pt {
	T x, y;
	Pt(T x = 0, T y = 0) : x(x), y(y) {}
	// P -> Pd happens implicitly (widening only). Pd -> P does NOT, so you can
	// never silently truncate a computed point back onto the lattice.
	template <class U, class = enable_if_t<is_integral_v<U> && !is_integral_v<T>>>
	Pt(const Pt<U> &o) : x((T)o.x), y((T)o.y) {}

	Pt operator+(const Pt &p) const { return Pt(x + p.x, y + p.y); }
	Pt operator-(const Pt &p) const { return Pt(x - p.x, y - p.y); }
	Pt operator*(T d) const { return Pt(x * d, y * d); }
	Pt operator/(T d) const { return Pt(x / d, y / d); }   // integer div if T=ll
	Pt operator-() const { return Pt(-x, -y); }
	Pt &operator+=(const Pt &p) { x += p.x; y += p.y; return *this; }
	Pt &operator-=(const Pt &p) { x -= p.x; y -= p.y; return *this; }
	Pt &operator*=(T d) { x *= d; y *= d; return *this; }

	// exact compares: keeps sort()/set<> a valid strict weak ordering
	bool operator==(const Pt &p) const { return x == p.x && y == p.y; }
	bool operator!=(const Pt &p) const { return !(*this == p); }
	bool operator<(const Pt &p) const { return x != p.x ? x < p.x : y < p.y; }
	bool operator>(const Pt &p) const { return p < *this; }
};

template <class T> Pt<T> operator*(T d, Pt<T> p) { return p * d; }
template <class T> istream &operator>>(istream &is, Pt<T> &p) { return is >> p.x >> p.y; }
template <class T> ostream &operator<<(ostream &os, const Pt<T> &p) {
	return os << p.x << " " << p.y;
}

typedef Pt<ll> P;                  // EXACT -- the default
typedef Pt<ld> Pd;                 // floating point
typedef vector<P> vP;
typedef vector<Pd> vPd;

template <class T> Pd toD(Pt<T> p) { return Pd((ld)p.x, (ld)p.y); }
template <class T> bool near(Pt<T> a, Pt<T> b) {                  // tolerant ==
	return sgn(a.x - b.x) == 0 && sgn(a.y - b.y) == 0;
}

// ---------------------------------------------------------------- primitives
template <class T> T dot(Pt<T> a, Pt<T> b) { return a.x * b.x + a.y * b.y; }
template <class T> T cross(Pt<T> a, Pt<T> b) { return a.x * b.y - a.y * b.x; }
template <class T> T cross(Pt<T> o, Pt<T> a, Pt<T> b) { return cross(a - o, b - o); }

template <class T> T len2(Pt<T> p) { return dot(p, p); }
template <class T> ld len(Pt<T> p) { return sqrt((ld)len2(p)); }
template <class T> T dist2(Pt<T> a, Pt<T> b) { return len2(a - b); }
template <class T> ld dist(Pt<T> a, Pt<T> b) { return len(a - b); }

// +1: c strictly LEFT of the directed line a->b (a,b,c is counter-clockwise)
// -1: c strictly right        0: a, b, c collinear
template <class T> int orient(Pt<T> a, Pt<T> b, Pt<T> c) { return sgn(cross(a, b, c)); }
template <class T> bool ccw(Pt<T> a, Pt<T> b, Pt<T> c) { return orient(a, b, c) > 0; }
template <class T> bool collinear(Pt<T> a, Pt<T> b, Pt<T> c) { return orient(a, b, c) == 0; }

// 2 * signed area of triangle abc (EXACT for P). Positive iff CCW.
template <class T> T area2(Pt<T> a, Pt<T> b, Pt<T> c) { return cross(a, b, c); }

// is p on segment [a,b]?  (endpoints count)
template <class T> bool onSegment(Pt<T> p, Pt<T> a, Pt<T> b) {
	return orient(a, b, p) == 0 && sgn(dot(a - p, b - p)) <= 0;
}
// is p strictly inside the segment (not an endpoint)?
template <class T> bool onSegmentStrict(Pt<T> p, Pt<T> a, Pt<T> b) {
	return orient(a, b, p) == 0 && sgn(dot(a - p, b - p)) < 0;
}

// ------------------------------------------------------------------ rotation
template <class T> Pt<T> perp(Pt<T> p) { return Pt<T>(-p.y, p.x); }    // +90 CCW
template <class T> Pt<T> perpCW(Pt<T> p) { return Pt<T>(p.y, -p.x); }  // -90 CW

Pd rot(Pd p, ld th) {                                    // CCW by th radians
	return Pd(p.x * cos(th) - p.y * sin(th), p.x * sin(th) + p.y * cos(th));
}
Pd rotAround(Pd p, Pd o, ld th) { return o + rot(p - o, th); }
Pd unit(Pd p) { return p / len(p); }
Pd scaleTo(Pd p, ld L) { return p * (L / len(p)); }
// reflect p through the point o
template <class T> Pt<T> reflectPoint(Pt<T> p, Pt<T> o) { return o * (T)2 - p; }

// ------------------------------------------------------------------- angles
template <class T> ld angle(Pt<T> p) { return atan2((ld)p.y, (ld)p.x); }  // (-pi,pi]
template <class T> ld angleBetween(Pt<T> a, Pt<T> b) {                    // [0, pi]
	return atan2((ld)abs(cross(a, b)), (ld)dot(a, b));
}
template <class T> ld signedAngle(Pt<T> a, Pt<T> b) {                     // (-pi,pi]
	return atan2((ld)cross(a, b), (ld)dot(a, b));
}
// angle of the corner at b in the path a-b-c, in [0, pi]
template <class T> ld cornerAngle(Pt<T> a, Pt<T> b, Pt<T> c) {
	return angleBetween(a - b, c - b);
}
// is p inside the CCW angular wedge from a to b (all three seen from origin)?
template <class T> bool inAngle(Pt<T> a, Pt<T> b, Pt<T> p) {
	if (cross(a, b) >= 0) return cross(a, p) >= 0 && cross(p, b) >= 0;
	return cross(a, p) >= 0 || cross(p, b) >= 0;
}

// --------------------------------------------------------------- polar sort
// half(p) = 0 for angle in [0, pi), 1 for [pi, 2pi). The zero vector lands in 0.
template <class T> int half(Pt<T> p) {
	int s = sgn(p.y);
	return s < 0 || (s == 0 && sgn(p.x) < 0);
}
// strict weak ordering by angle CCW starting at the +x axis; ties broken by
// length (closer first), so equal-direction points end up adjacent.
template <class T> bool polarCmp(Pt<T> a, Pt<T> b) {
	int ha = half(a), hb = half(b);
	if (ha != hb) return ha < hb;
	int c = sgn(cross(a, b));
	return c ? c > 0 : len2(a) < len2(b);
}
template <class T> void polarSort(vector<Pt<T>> &v) {
	sort(v.begin(), v.end(), polarCmp<T>);
}
template <class T> void polarSortAround(Pt<T> o, vector<Pt<T>> &v) {
	sort(v.begin(), v.end(), [&](Pt<T> a, Pt<T> b) { return polarCmp(a - o, b - o); });
}
// same but starting the sweep at a chosen direction `dir` instead of +x
template <class T> void polarSortFrom(Pt<T> o, Pt<T> dir, vector<Pt<T>> &v) {
	auto hf = [&](Pt<T> p) {
		int c = sgn(cross(dir, p));
		return c < 0 || (c == 0 && sgn(dot(dir, p)) < 0);
	};
	sort(v.begin(), v.end(), [&](Pt<T> A, Pt<T> B) {
		Pt<T> a = A - o, b = B - o;
		int ha = hf(a), hb = hf(b);
		if (ha != hb) return ha < hb;
		int c = sgn(cross(a, b));
		return c ? c > 0 : len2(a) < len2(b);
	});
}

// ------------------------------------------------------------------ lattice
// number of lattice points STRICTLY between two lattice points (a != b)
ll latticeOnSeg(P a, P b) { return gcd(llabs(a.x - b.x), llabs(a.y - b.y)) - 1; }

// ---------------------------------------------------------------- utilities
// remove duplicate points (needs no particular input order)
template <class T> void dedup(vector<Pt<T>> &v) {
	sort(v.begin(), v.end());
	v.erase(unique(v.begin(), v.end()), v.end());
}
// bounding box: returns {lowerLeft, upperRight}
template <class T> pair<Pt<T>, Pt<T>> bbox(const vector<Pt<T>> &v) {
	T x0 = v[0].x, x1 = v[0].x, y0 = v[0].y, y1 = v[0].y;
	for (auto &p : v) {
		x0 = min(x0, p.x); x1 = max(x1, p.x);
		y0 = min(y0, p.y); y1 = max(y1, p.y);
	}
	return {Pt<T>(x0, y0), Pt<T>(x1, y1)};
}

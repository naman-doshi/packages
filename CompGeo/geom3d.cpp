// ============================================================================
//  3-D GEOMETRY + 3-D CONVEX HULL   -- REQUIRES point.cpp (for ld, EPS, sgn)
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Points, planes, spheres and the 3-D convex hull. The hull is also how you
//    get a Delaunay triangulation in 2-D (see delaunay.cpp), so it earns its
//    place even in problems that look flat.
//
//  API
//    P3 a(1,2,3), b(4,5,6);
//    a + b   a - b   a * 2   a / 2   a == b
//    dot3(a,b)  cross3(a,b)  len3(a)  len3sq(a)  dist3(a,b)  unit3(a)
//    tripleProduct(a,b,c)             // = dot3(a, cross3(b,c)), 6*tetra volume
//    planeFrom3(a,b,c)                // -> {normal, d}, plane is dot(n,x) = d
//    distToPlane(p, pl)               // signed
//    projectOnPlane(p, pl)
//    lineePlaneInter(a, b, pl)        // -> {hit?, point}
//    lineSphere(a, b, centre, r)      // 0/1/2 points
//    rotate3(p, axis, th)             // Rodrigues, axis need not be a unit vec
//    sphericalDist(lat1,lon1,lat2,lon2, R)   // great-circle, degrees in
//    hull3d(pts)                      // -> triangular faces, CCW from outside
//    hullVolume(pts, faces)  hullArea(pts, faces)
//    pointInHull(pts, faces, q)       // 0 out / 1 on the surface / 2 inside
//
//  COMPLEXITY  hull3d is O(n^2) (incremental with an explicit horizon walk).
//    Good to n ~ 2000-5000. Everything else is O(1).
//
//  PITFALLS
//    * hull3d needs at least 4 points that are NOT all coplanar. It searches
//      for a valid starting tetrahedron and returns {} if there is none.
//    * Points that are coplanar with an existing face are treated as NOT
//      visible, so a point lying exactly on the hull surface is silently
//      dropped from the face list. That keeps the output a valid triangulation;
//      it does mean "is this point a hull vertex" is not answered by presence.
//    * The faces are triangles even when the true face is a polygon (a cube
//      comes out as 12 triangles). Merge coplanar neighbours if you need the
//      real faces.
//    * All floating point. Scale coordinates into a sane range first: the
//      visibility test compares a triple product, which is cubic in the input.
// ============================================================================

struct P3 {
	ld x, y, z;
	P3(ld x = 0, ld y = 0, ld z = 0) : x(x), y(y), z(z) {}
	P3 operator+(const P3 &o) const { return P3(x + o.x, y + o.y, z + o.z); }
	P3 operator-(const P3 &o) const { return P3(x - o.x, y - o.y, z - o.z); }
	P3 operator*(ld d) const { return P3(x * d, y * d, z * d); }
	P3 operator/(ld d) const { return P3(x / d, y / d, z / d); }
	bool operator==(const P3 &o) const {
		return sgn(x - o.x) == 0 && sgn(y - o.y) == 0 && sgn(z - o.z) == 0;
	}
	bool operator<(const P3 &o) const {
		return x != o.x ? x < o.x : (y != o.y ? y < o.y : z < o.z);
	}
};
typedef vector<P3> vP3;
istream &operator>>(istream &is, P3 &p) { return is >> p.x >> p.y >> p.z; }
ostream &operator<<(ostream &os, const P3 &p) { return os << p.x << " " << p.y << " " << p.z; }

ld dot3(P3 a, P3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
P3 cross3(P3 a, P3 b) {
	return P3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
ld len3sq(P3 a) { return dot3(a, a); }
ld len3(P3 a) { return sqrt(len3sq(a)); }
ld dist3(P3 a, P3 b) { return len3(a - b); }
P3 unit3(P3 a) { return a / len3(a); }
// 6 * the signed volume of the tetrahedron (0, a, b, c)
ld tripleProduct(P3 a, P3 b, P3 c) { return dot3(a, cross3(b, c)); }
ld tetraVolume(P3 a, P3 b, P3 c, P3 d) {
	return fabs(tripleProduct(b - a, c - a, d - a)) / 6;
}

// plane: dot3(n, x) == d
struct Plane { P3 n; ld d; };
Plane planeFrom3(P3 a, P3 b, P3 c) {
	P3 n = cross3(b - a, c - a);
	return Plane{n, dot3(n, a)};
}
ld distToPlane(P3 p, Plane pl) { return (dot3(pl.n, p) - pl.d) / len3(pl.n); }
P3 projectOnPlane(P3 p, Plane pl) {
	return p - pl.n * ((dot3(pl.n, p) - pl.d) / len3sq(pl.n));
}
// where the line a->b meets the plane
pair<bool, P3> linePlaneInter(P3 a, P3 b, Plane pl) {
	ld den = dot3(pl.n, b - a);
	if (sgn(den) == 0) return {false, P3()};
	return {true, a + (b - a) * ((pl.d - dot3(pl.n, a)) / den)};
}
// where the line a->b meets the sphere (c, r)
vP3 lineSphere(P3 a, P3 b, P3 c, ld r) {
	P3 dir = b - a, f = a - c;
	ld A = len3sq(dir), B = 2 * dot3(f, dir), C = len3sq(f) - r * r;
	ld disc = B * B - 4 * A * C;
	if (sgn(disc) < 0) return {};
	disc = sqrt(max((ld)0, disc));
	ld t1 = (-B - disc) / (2 * A), t2 = (-B + disc) / (2 * A);
	if (sgn(disc) == 0) return {a + dir * t1};
	return {a + dir * t1, a + dir * t2};
}
// rotate p around the axis vector `ax` by th radians (right-hand rule)
P3 rotate3(P3 p, P3 ax, ld th) {
	P3 k = unit3(ax);
	return p * cos(th) + cross3(k, p) * sin(th) + k * (dot3(k, p) * (1 - cos(th)));
}
// great-circle distance; latitudes/longitudes in DEGREES
ld sphericalDist(ld lat1, ld lon1, ld lat2, ld lon2, ld R) {
	auto rad = [](ld d) { return d * PI / 180; };
	lat1 = rad(lat1); lon1 = rad(lon1); lat2 = rad(lat2); lon2 = rad(lon2);
	ld h = sin((lat2 - lat1) / 2) * sin((lat2 - lat1) / 2) +
	       cos(lat1) * cos(lat2) * sin((lon2 - lon1) / 2) * sin((lon2 - lon1) / 2);
	return 2 * R * asin(sqrt(min((ld)1, h)));
}

// -------------------------------------------------------------- convex hull
// Triangular faces of the 3-D convex hull, each CCW when viewed from OUTSIDE.
vector<array<int, 3>> hull3d(const vP3 &p) {
	int n = (int)p.size();
	vector<array<int, 3>> faces;
	if (n < 4) return faces;
	// pick a starting tetrahedron that is not degenerate
	int i0 = 0, i1 = -1, i2 = -1, i3 = -1;
	for (int i = 1; i < n && i1 < 0; i++)
		if (!(p[i] == p[i0])) i1 = i;
	if (i1 < 0) return faces;
	for (int i = 1; i < n && i2 < 0; i++)
		if (i != i1 && sgn(len3sq(cross3(p[i1] - p[i0], p[i] - p[i0]))) > 0) i2 = i;
	if (i2 < 0) return faces;
	for (int i = 1; i < n && i3 < 0; i++)
		if (i != i1 && i != i2 &&
		    sgn(tripleProduct(p[i1] - p[i0], p[i2] - p[i0], p[i] - p[i0])) != 0)
			i3 = i;
	if (i3 < 0) return faces;                       // everything is coplanar

	auto normal = [&](array<int, 3> f) {
		return cross3(p[f[1]] - p[f[0]], p[f[2]] - p[f[0]]);
	};
	auto outward = [&](array<int, 3> f, int inside) {   // orient away from `inside`
		if (sgn(dot3(normal(f), p[inside] - p[f[0]])) > 0) swap(f[1], f[2]);
		return f;
	};
	int base[4] = {i0, i1, i2, i3};
	for (int a = 0; a < 4; a++) {
		array<int, 3> f{};
		int k = 0, opp = 0;
		for (int b = 0; b < 4; b++) {
			if (b == a) opp = base[b];
			else f[k++] = base[b];
		}
		faces.push_back(outward(f, opp));
	}
	vector<bool> done(n, false);
	for (int i : base) done[i] = true;

	for (int i = 0; i < n; i++) {
		if (done[i]) continue;
		vector<char> vis(faces.size(), 0);
		bool any = false;
		for (size_t f = 0; f < faces.size(); f++)
			if (sgn(dot3(normal(faces[f]), p[i] - p[faces[f][0]])) > 0) vis[f] = 1, any = true;
		if (!any) continue;                          // inside or on the surface
		// horizon = directed edges of visible faces whose twin is not visible
		set<pair<int, int>> seen;
		for (size_t f = 0; f < faces.size(); f++)
			if (vis[f])
				for (int e = 0; e < 3; e++)
					seen.insert({faces[f][e], faces[f][(e + 1) % 3]});
		vector<array<int, 3>> keep;
		for (size_t f = 0; f < faces.size(); f++)
			if (!vis[f]) keep.push_back(faces[f]);
		for (auto &e : seen)
			if (!seen.count({e.second, e.first})) keep.push_back({e.first, e.second, i});
		faces = keep;
	}
	return faces;
}

ld hullVolume(const vP3 &p, const vector<array<int, 3>> &f) {
	ld v = 0;
	for (auto &t : f) v += tripleProduct(p[t[0]], p[t[1]], p[t[2]]);
	return fabs(v) / 6;
}
ld hullArea(const vP3 &p, const vector<array<int, 3>> &f) {
	ld a = 0;
	for (auto &t : f) a += len3(cross3(p[t[1]] - p[t[0]], p[t[2]] - p[t[0]]));
	return a / 2;
}
// 0 outside, 1 on the surface, 2 strictly inside
int pointInHull(const vP3 &p, const vector<array<int, 3>> &f, P3 q) {
	int on = 0;
	for (auto &t : f) {
		ld s = dot3(cross3(p[t[1]] - p[t[0]], p[t[2]] - p[t[0]]), q - p[t[0]]);
		if (sgn(s) > 0) return 0;
		if (sgn(s) == 0) on = 1;
	}
	return on ? 1 : 2;
}

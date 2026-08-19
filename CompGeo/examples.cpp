// ============================================================================
//  EXAMPLES -- one worked call for every function in CompGeo/
// ----------------------------------------------------------------------------
//  This file is RUNNABLE and SELF-CHECKING. Build it with every module:
//
//     cat point.cpp lines.cpp polygon.cpp convex.cpp circles.cpp halfplane.cpp \
//         closestPair.cpp sweep.cpp kdtree.cpp angular.cpp geom3d.cpp \
//         delaunay.cpp examples.cpp > /tmp/ex.cpp
//     g++ -std=c++17 -O2 -o /tmp/ex /tmp/ex.cpp && /tmp/ex
//
//  It prints "ALL EXAMPLES OK". Every line below is of the form
//
//     ok( <a real call> , "<what it returns>" );
//
//  so the call is the example and the string is the documented answer -- and
//  the program fails loudly if the two ever stop agreeing. Ctrl-F a function
//  name to see how it is used.
// ============================================================================

int bad = 0;
void ok(bool c, const string &what) {
	if (!c) { cout << "MISMATCH: " << what << "\n"; bad++; }
}
bool close(ld a, ld b, ld eps = 1e-6) { return fabs(a - b) < eps; }

signed main() {
	// =====================================================================
	//  point.cpp
	// =====================================================================
	P a(3, 4), b(1, 2), o(0, 0);

	ok(sgn(-7LL) == -1 && sgn(0LL) == 0 && sgn(2.5) == 1, "sgn: -1 / 0 / +1");
	ok((a + b) == P(4, 6), "P(3,4) + P(1,2) == P(4,6)");
	ok((a - b) == P(2, 2), "P(3,4) - P(1,2) == P(2,2)");
	ok((a * 2) == P(6, 8), "P(3,4) * 2 == P(6,8)");
	ok((a / 2) == P(1, 2), "P(3,4) / 2 == P(1,2)   -- INTEGER division!");
	ok((-a) == P(-3, -4), "-P(3,4) == P(-3,-4)");
	ok(b < a, "operator< is lexicographic: P(1,2) < P(3,4)");
	ok(near(Pd(1, 2), Pd(1 + 1e-12, 2)), "near() compares Pd with EPS; == is exact");
	ok(toD(a).x == 3.0, "toD(P) -> Pd; the reverse never happens implicitly");

	ok(dot(a, b) == 11, "dot(P(3,4),P(1,2)) == 11   (>0 same side, 0 perpendicular)");
	ok(cross(a, b) == 2, "cross(P(3,4),P(1,2)) == 2   (>0 means b is CCW from a)");
	ok(cross(o, P(1, 0), P(0, 1)) == 1, "cross(o,a,b) == cross(a-o, b-o) == 1");
	ok(len2(a) == 25 && close(len(a), 5), "len2(P(3,4))==25 exact, len(...)==5");
	ok(dist2(a, b) == 8 && close(dist(a, b), sqrt(8.0)), "dist2 is exact, dist is a double");

	ok(orient(o, P(1, 0), P(0, 1)) == 1, "orient: +1 when c is LEFT of a->b (CCW)");
	ok(orient(o, P(1, 0), P(0, -1)) == -1, "orient: -1 when c is to the right");
	ok(orient(o, P(1, 0), P(2, 0)) == 0, "orient: 0 when collinear");
	ok(ccw(o, P(1, 0), P(0, 1)), "ccw(a,b,c) is orient(a,b,c) > 0");
	ok(collinear(o, P(1, 1), P(5, 5)), "collinear(a,b,c)");
	ok(area2(o, P(4, 0), P(0, 3)) == 12, "area2 = 2x signed triangle area = 12 (area 6)");
	ok(onSegment(P(2, 2), o, P(4, 4)), "onSegment: endpoints included");
	ok(!onSegmentStrict(P(4, 4), o, P(4, 4)), "onSegmentStrict: endpoints excluded");

	ok(perp(P(1, 0)) == P(0, 1), "perp = rotate 90 CCW");
	ok(perpCW(P(1, 0)) == P(0, -1), "perpCW = rotate 90 CW");
	ok(near(rot(Pd(1, 0), PI / 2), Pd(0, 1)), "rot(p, theta) rotates about the origin");
	ok(near(rotAround(Pd(2, 1), Pd(1, 1), PI), Pd(0, 1)), "rotAround(p, centre, theta)");
	ok(close(len(unit(Pd(3, 4))), 1), "unit(p) has length 1");
	ok(close(len(scaleTo(Pd(3, 4), 10)), 10), "scaleTo(p, L) has length L");
	ok(reflectPoint(P(5, 5), P(1, 1)) == P(-3, -3), "reflectPoint(p, o): mirror through a POINT");

	ok(close(angle(P(0, 1)), PI / 2), "angle(p) = atan2, in (-pi, pi]");
	ok(close(angleBetween(P(1, 0), P(0, 1)), PI / 2), "angleBetween is unsigned, [0,pi]");
	ok(close(signedAngle(P(1, 0), P(0, -1)), -PI / 2), "signedAngle is signed, (-pi,pi]");
	ok(close(cornerAngle(P(1, 0), o, P(0, 1)), PI / 2), "cornerAngle(a,b,c): the angle AT b");
	ok(inAngle(P(1, 0), P(0, 1), P(1, 1)), "inAngle(a,b,p): is p in the CCW wedge a..b");

	ok(half(P(1, 1)) == 0 && half(P(1, -1)) == 1, "half(): 0 for the upper half plane, 1 for the lower");
	ok(polarCmp(P(1, 0), P(0, 1)), "polarCmp orders CCW from the +x axis");
	{
		vP v{P(0, -1), P(-1, 0), P(1, 1), P(1, 0)};
		polarSort(v);
		ok(v[0] == P(1, 0) && v[1] == P(1, 1) && v[2] == P(-1, 0) && v[3] == P(0, -1),
		   "polarSort -> +x, then CCW: (1,0) (1,1) (-1,0) (0,-1)");
		vP w{P(6, 5), P(4, 5), P(5, 6)};
		polarSortAround(P(5, 5), w);
		ok(w[0] == P(6, 5), "polarSortAround(centre, v) sorts by angle seen from `centre`");
		vP u{P(1, 0), P(0, 1), P(-1, 0)};
		polarSortFrom(o, P(0, 1), u);
		ok(u[0] == P(0, 1), "polarSortFrom(o, dir, v) starts the sweep at `dir` instead of +x");
	}
	ok(latticeOnSeg(P(0, 0), P(4, 6)) == 1, "latticeOnSeg: 1 lattice point strictly inside (0,0)-(4,6)");
	{
		vP v{P(1, 1), P(0, 0), P(1, 1)};
		dedup(v);
		ok(v.size() == 2, "dedup() sorts and removes duplicate points");
		auto bb = bbox(vP{P(1, 5), P(-2, 3), P(4, -1)});
		ok(bb.first == P(-2, -1) && bb.second == P(4, 5), "bbox -> {lowerLeft, upperRight}");
	}

	// =====================================================================
	//  lines.cpp     -- segment (a,b) and segment (c,d)
	// =====================================================================
	P L0(0, 0), L1(4, 0), M0(2, -2), M1(2, 2);

	ok(sideOf(L0, L1, P(1, 1)) == 1, "sideOf(a,b,p): +1 left of a->b");
	ok(linesParallel(L0, L1, P(0, 3), P(4, 3)), "linesParallel");
	ok(linesSame(L0, L1, P(7, 0), P(9, 0)), "linesSame: parallel AND collinear");
	ok(perpendicular(L0, L1, M0, M1), "perpendicular");

	{
		auto r = lineInter(L0, L1, M0, M1);
		ok(r.first == 1 && near(r.second, Pd(2, 0)), "lineInter -> {1, (2,0)}: unique crossing");
		ok(lineInter(L0, L1, P(0, 3), P(4, 3)).first == 0, "lineInter -> {0, ...}: parallel, distinct");
		ok(lineInter(L0, L1, P(7, 0), P(9, 0)).first == -1, "lineInter -> {-1, ...}: the same line");
		ok(near(ptOnLineAt(Pd(0, 0), Pd(4, 0), 0.25), Pd(1, 0)), "ptOnLineAt(a,b,t) = a + t*(b-a)");
	}
	ok(segInter(L0, L1, M0, M1), "segInter: do the two segments share a point (touching counts)");
	ok(!segInter(L0, L1, P(5, -1), P(5, 1)), "segInter: false when they miss");
	ok(segInter(L0, L1, P(4, 0), P(9, 9)), "segInter: TOUCHING at an endpoint is true...");
	ok(!segInterProper(L0, L1, P(4, 0), P(9, 9)), "...but segInterProper wants a strict crossing");
	{
		auto v = segInterPts(L0, L1, M0, M1);
		ok(v.size() == 1 && near(v[0], Pd(2, 0)), "segInterPts -> 1 point when they cross");
		auto w = segInterPts(L0, L1, P(2, 0), P(9, 0));
		ok(w.size() == 2 && near(w[0], Pd(2, 0)) && near(w[1], Pd(4, 0)),
		   "segInterPts -> 2 points: the overlap of collinear segments is (2,0)..(4,0)");
		ok(segInterPts(L0, L1, P(0, 1), P(4, 1)).empty(), "segInterPts -> {} when disjoint");
	}
	ok(near(proj(Pd(2, 5), L0, L1), Pd(2, 0)), "proj(p,a,b): foot of the perpendicular on the LINE");
	ok(near(refl(Pd(2, 5), L0, L1), Pd(2, -5)), "refl(p,a,b): mirror image across the line");
	ok(near(closestOnSeg(Pd(9, 5), L0, L1), Pd(4, 0)), "closestOnSeg clamps to the segment's end");
	ok(close(distToLine(Pd(2, 5), L0, L1), 5), "distToLine: 5");
	ok(close(signedDistToLine(Pd(2, -5), L0, L1), -5), "signedDistToLine: negative on the right");
	ok(close(distToSeg(Pd(9, 0), L0, L1), 5), "distToSeg: measured to the endpoint, so 5");
	ok(close(lineDist(L0, L1, P(0, 3), P(4, 3)), 3), "lineDist between two PARALLEL lines");
	ok(close(segDist(L0, L1, P(0, 3), P(4, 3)), 3), "segDist: segment to segment (0 if they touch)");
	{
		auto pb = perpBisector(Pd(0, 0), Pd(4, 0));
		ok(near(pb.first, Pd(2, 0)) && close(cross(pb.second, Pd(0, 1)), 0),
		   "perpBisector -> {midpoint (2,0), direction along y}");
		ok(close(cross(angleBisector(Pd(0, 0), Pd(4, 0), Pd(0, 4)), Pd(1, 1)), 0),
		   "angleBisector(o,a,b): direction (1,1) for the right angle at the origin");
	}
	{
		ABC ln = lineFromPts(Pd(0, 0), Pd(4, 0));      // the line y = 0
		ok(close(evalABC(ln, Pd(3, 0)), 0), "evalABC == 0 exactly on the line");
		auto hit = interABC(ln, lineFromPts(Pd(2, -2), Pd(2, 2)));
		ok(hit.first && near(hit.second, Pd(2, 0)), "interABC: intersection in ax+by=c form");
	}
	ok(segRect(P(-5, 1), P(5, 1), P(0, 0), P(4, 4)), "segRect: does a segment touch an axis-aligned box");
	{
		//  (0,0) (1,0) (2,0) (3,0) all on y = 0;  (1,1) (2,2) off it
		vP pts{P(0, 0), P(1, 0), P(2, 0), P(3, 0), P(1, 1), P(2, 2)};
		ok(countOnLine(pts, P(0, 0), P(10, 0)) == 4, "countOnLine: 4 of the points are on y = 0");
		ok(countOnSeg(pts, P(0, 0), P(2, 0)) == 3, "countOnSeg: only 3 of them are between (0,0) and (2,0)");
		auto id = pointsOnLine(pts, P(10, 0), P(0, 0));
		ok(id.size() == 4 && pts[id[0]] == P(3, 0), "pointsOnLine: indices ordered along a->b, so (3,0) first");
		ok(lineKey(P(0, 0), P(2, 0)) == lineKey(P(7, 0), P(-1, 0)),
		   "lineKey: same canonical {A,B,C} for two descriptions of y = 0");
		ok(lineKey(P(0, 0), P(2, 0)) != lineKey(P(0, 1), P(2, 1)), "lineKey: y = 1 is a different key");
		OnLine tab(pts);
		ok(tab.count(P(0, 0), P(10, 0)) == 4, "OnLine::count: same 4, in O(log n) after an O(n^2 log n) build");
		ok(tab.count(P(0, 0), P(2, 2)) == 3, "OnLine::count: (0,0) (1,1) (2,2) are collinear too");
		ok(tab.count(P(0, 5), P(1, 5)) == 0, "OnLine::count: 0 when the line misses everything");
		ok(tab.count(P(1, 1), P(2, 5)) == 1, "OnLine::count: exact in the 1-point case too, via the direction index");
		ok(tab.any(P(1, 1), P(2, 5)), "OnLine::any: >= 1 point on the line, O(log n)");
		ok(!tab.any(P(0, 5), P(1, 5)), "OnLine::any: false when the line holds none of them");
		ok(tab.any(P(9, 9), P(10, 9)) == false && tab.any(P(9, 2), P(10, 2)) == true,
		   "OnLine::any: later queries of a seen direction are pure binary searches");

		LineSweep sw(pts);
		ok(sw.count(P(0, 0), P(10, 0)) == 4, "LineSweep::count: 4 points on y = 0, O(log^2 n) for ANY direction");
		ok(sw.count(P(0, 0), P(2, 2)) == 3 && sw.count(P(1, 1), P(2, 5)) == 1,
		   "LineSweep::count: exact whether the line holds 3 points or 1");
		ok(sw.any(P(1, 1), P(2, 5)) && !sw.any(P(0, 5), P(1, 5)), "LineSweep::any: >= 1 point on the line?");
		ok(sw.countLeft(P(0, 0), P(10, 0)) == 2, "LineSweep::countLeft: (1,1) and (2,2) are left of ->x");
		ok(sw.countRight(P(0, 0), P(10, 0)) == 0, "LineSweep::countRight: nothing below y = 0");
		ok(sw.countLeft(P(10, 0), P(0, 0)) == 0 && sw.countRight(P(10, 0), P(0, 0)) == 2,
		   "LineSweep: reversing the line swaps left and right");
		ok(sw.countLeft(P(0, 0), P(0, 10)) + sw.countRight(P(0, 0), P(0, 10)) + sw.count(P(0, 0), P(0, 10))
		   == (int)pts.size(), "LineSweep: left + right + on == n, always");
	}

	// =====================================================================
	//  polygon.cpp
	// =====================================================================
	vP sq{P(0, 0), P(4, 0), P(4, 4), P(0, 4)};         // a CCW 4x4 square
	vP cw{P(0, 0), P(0, 4), P(4, 4), P(4, 0)};         // the same, clockwise

	ok(polyArea2(sq) == 32, "polyArea2 = 2 x SIGNED area = 32 (positive => CCW)");
	ok(polyArea2(cw) == -32, "polyArea2 < 0 for a clockwise polygon");
	ok(close(polyArea(sq), 16), "polyArea = |area| as a double = 16");
	{
		vP t = cw;
		makeCCW(t);
		ok(polyArea2(t) > 0, "makeCCW flips a clockwise polygon in place");
	}
	ok(close(perimeter(sq), 16), "perimeter = 16");
	ok(isConvex(sq), "isConvex");
	ok(!isConvex(vP{P(0, 0), P(4, 0), P(1, 1), P(0, 4)}), "isConvex is false for a dart");
	ok(isSimple(sq) && !isSimple(vP{P(0, 0), P(4, 4), P(4, 0), P(0, 4)}),
	   "isSimple: true for the square, false for the bow-tie");

	ok(inPolygon(sq, P(2, 2)) == 2, "inPolygon -> 2: strictly inside");
	ok(inPolygon(sq, P(4, 2)) == 1, "inPolygon -> 1: ON the boundary");
	ok(inPolygon(sq, P(9, 2)) == 0, "inPolygon -> 0: outside");
	ok(windingNumber(sq, P(2, 2)) == 1, "windingNumber != 0 means inside (works on self-crossing polys)");
	ok(near(centroid(sq), Pd(2, 2)), "centroid: area-weighted centre of mass");
	ok(close(distToPoly(sq, Pd(6, 2)), 2), "distToPoly: distance to the BOUNDARY");
	ok(pointInTriangle(P(1, 1), P(0, 0), P(4, 0), P(0, 4)), "pointInTriangle (boundary counts)");

	{
		vPd sqd;
		for (auto &p : sq) sqd.push_back(toD(p));
		vPd half = polygonCut(sqd, Pd(0, 2), Pd(4, 2));      // keep the LEFT of a->b
		ok(close(polyArea(half), 8), "polygonCut keeps the half left of a->b: area 8");
		vPd clip{Pd(2, 2), Pd(6, 2), Pd(6, 6), Pd(2, 6)};
		ok(close(polyArea(convexClip(sqd, clip)), 4), "convexClip: square AND clip box = 4");
	}
	{
		auto tri = triangulate(sq);
		ok(tri.size() == 2, "triangulate: an n-gon becomes n-2 triangles");
		ll s = 0;
		for (auto &t : tri) s += llabs(cross(sq[t[0]], sq[t[1]], sq[t[2]]));
		ok(s == 32, "...and their areas add back up to the polygon's");
	}
	ok(boundaryLattice(sq) == 16, "boundaryLattice B = 16 lattice points on the border");
	ok(interiorLattice(sq) == 9, "interiorLattice I = 9, via Pick's theorem A = I + B/2 - 1");
	{
		vP r{P(4, 4), P(0, 4), P(0, 0), P(4, 0)};
		normalizePoly(r);
		ok(r[0] == P(0, 0), "normalizePoly rotates the smallest vertex to the front");
	}

	// =====================================================================
	//  convex.cpp
	// =====================================================================
	vP pts{P(0, 0), P(4, 0), P(4, 4), P(0, 4), P(2, 2), P(2, 0)};

	vP h = convexHull(pts);
	ok(h.size() == 4 && h[0] == P(0, 0), "convexHull: 4 corners, CCW, starting at the smallest point");
	ok(convexHull(pts, true).size() == 5, "convexHull(pts, true) also KEEPS (2,0), collinear on an edge");
	{
		auto id = hullIdx(pts);
		ok(id.size() == 4 && pts[id[0]] == P(0, 0), "hullIdx: the same hull, as indices into pts");
	}
	{
		auto d = hullDiameter(h);
		ok(dist2(h[d.first], h[d.second]) == 32, "hullDiameter: the farthest pair, squared distance 32");
		auto f = farthestPair(pts);
		ok(dist2(pts[f.first], pts[f.second]) == 32, "farthestPair: same, but indices into the raw points");
	}
	ok(close(hullWidth(h), 4), "hullWidth: narrowest slab that still contains the hull");
	{
		auto r = minAreaRect(h);
		ok(close(r.first, 16), "minAreaRect -> {16, the 4 corners}");
		ok(close(minPerimeterRect(h).first, 16), "minPerimeterRect -> {16, corners}");
	}
	ok(maxDot(h, P(1, 1)) == 8, "maxDot(h, dir): best dot(dir, vertex) in O(log n) = 8");
	ok(h[extreme(h, P(1, 1))] == P(4, 4), "extreme(h, dir): the vertex attaining it");
	ok(inConvex(h, P(2, 2)) == 2, "inConvex -> 2 inside   (O(log n))");
	ok(inConvex(h, P(4, 2)) == 1, "inConvex -> 1 on the boundary");
	ok(inConvex(h, P(5, 2)) == 0, "inConvex -> 0 outside");
	{
		auto t = tangents(h, P(10, 2));
		for (auto &z : h) ok(orient(P(10, 2), h[t.first], z) >= 0 && orient(P(10, 2), h[t.second], z) <= 0,
		                     "tangents(h,p): every vertex lies inside the wedge from p");
	}
	{
		vP A{P(0, 0), P(2, 0), P(2, 2), P(0, 2)}, B{P(0, 0), P(1, 0), P(1, 1), P(0, 1)};
		ok(close(polyArea(minkowski(A, B)), 9), "minkowski of a 2x2 and a 1x1 square is 3x3, area 9");
		ok(close(convexInterArea(A, B), 1), "convexInterArea(A,B) = 1");
		ok(close(polyDist(A, vP{P(5, 0), P(6, 0), P(6, 1)}), 3), "polyDist: 3 between the two polygons");
	}
	ok(maxTriangle(h) == 16, "maxTriangle: 2 x the biggest inscribed triangle area (=8)");
	{
		vP many;
		for (int x = 0; x <= 3; x++)
			for (int y = 0; y <= 3; y++) many.push_back(P(x, y));
		auto layers = convexLayers(many);
		ok(layers.size() == 3 && layers[0].size() == 4 && layers[1].size() == 8,
		   "convexLayers: a 4x4 grid peels into 3 rings of 4, 8 and 4 points");
	}

	// =====================================================================
	//  circles.cpp
	// =====================================================================
	Circle C(Pd(0, 0), 5);

	ok(close(arcLength(2, PI), 2 * PI), "arcLength(r, theta)");
	ok(close(sectorArea(2, PI), 2 * PI), "sectorArea(r, theta)");
	ok(close(segmentArea(1, PI), PI / 2), "segmentArea(r, theta): the bit cut off by a chord");
	ok(inCircle(C, P(0, 0)) == 2 && inCircle(C, P(5, 0)) == 1 && inCircle(C, P(9, 0)) == 0,
	   "inCircle -> 2 inside / 1 on / 0 outside");
	{
		auto v = circleLine(C, P(-9, 3), P(9, 3));
		ok(v.size() == 2 && close(v[0].x, -4) && close(v[1].x, 4),
		   "circleLine: the line y=3 cuts the r=5 circle at x = -4 and +4 (ordered a->b)");
		ok(circleSeg(C, P(0, 3), P(9, 3)).size() == 1, "circleSeg keeps only the hits ON the segment");
		ok(circleLine(C, P(-9, 7), P(9, 7)).empty(), "circleLine -> {} when the line misses");
	}
	{
		auto v = circleCircle(C, Circle(Pd(8, 0), 5));
		ok(v.size() == 2 && close(v[0].x, 4) && close(fabs(v[0].y), 3), "circleCircle -> 2 points (4,-3),(4,3)");
		ok(circleCircle(C, Circle(Pd(20, 0), 5)).empty(), "circleCircle -> {} when they are too far apart");
	}
	ok(close(tangentLen(C, Pd(13, 0)), 12), "tangentLen: 13-5-12 right triangle");
	{
		auto v = tangentPoints(C, Pd(13, 0));
		ok(v.size() == 2, "tangentPoints: the 2 points of C where a line from p is tangent");
		for (auto &q : v) ok(close(dist(C.c, q), 5) && close(dot(q - C.c, Pd(13, 0) - q), 0),
		                     "...each is on the circle and the radius meets the tangent at 90 degrees");
	}
	{
		auto v = commonTangents(C, Circle(Pd(20, 0), 5), false);
		ok(v.size() == 2, "commonTangents(..., false): 2 OUTER tangents");
		ok(commonTangents(C, Circle(Pd(20, 0), 5), true).size() == 2, "commonTangents(..., true): 2 INNER ones");
	}
	{
		Circle cc = circumcircle(P(0, 0), P(4, 0), P(0, 4));
		ok(near(cc.c, Pd(2, 2)) && close(cc.r, sqrt(8.0)), "circumcircle through 3 points");
		Circle ic = incircle(P(0, 0), P(4, 0), P(0, 3));
		ok(close(ic.r, 1), "incircle of the 3-4-5 triangle has r = 1");
	}
	{
		Circle m = minEnclosingCircle(vPd{Pd(0, 0), Pd(4, 0), Pd(0, 4), Pd(4, 4), Pd(2, 2)});
		ok(near(m.c, Pd(2, 2)) && close(m.r, sqrt(8.0)), "minEnclosingCircle of a square: centre (2,2)");
	}
	{
		auto v = circleFrom2R(Pd(0, 0), Pd(6, 0), 5);
		ok(v.size() == 2 && close(fabs(v[0].y), 4), "circleFrom2R -> the 2 CENTRES (3,-4) and (3,4)");
	}
	ok(close(circleInterArea(C, Circle(Pd(20, 0), 5)), 0), "circleInterArea: 0 when disjoint");
	ok(close(circleInterArea(C, Circle(Pd(0, 0), 2)), 4 * PI), "circleInterArea: the smaller disk when nested");
	ok(close(circlePolyArea(Circle(Pd(0, 0), 100), sq), 16),
	   "circlePolyArea: a huge circle over the 4x4 square gives the square's area");
	ok(close(circleUnionArea(vector<Circle>{C, C}), 25 * PI), "circleUnionArea de-duplicates: still 25pi");
	ok(close(circleUnionArea(vector<Circle>{C, Circle(Pd(100, 0), 1)}), 25 * PI + PI),
	   "circleUnionArea of two disjoint circles adds up");
	ok(inCircleDet(P(0, 0), P(4, 0), P(0, 4), P(2, 2)) > 0,
	   "inCircleDet > 0: (2,2) is strictly inside the circumcircle. EXACT, a,b,c must be CCW");

	// =====================================================================
	//  halfplane.cpp
	// =====================================================================
	{
		// the square [0,4]^2 written as four "keep the left of a->b" constraints
		vector<HP> hs{HP(Pd(0, 0), Pd(4, 0)), HP(Pd(4, 0), Pd(4, 4)),
		              HP(Pd(4, 4), Pd(0, 4)), HP(Pd(0, 4), Pd(0, 0))};
		vPd reg = halfplaneInter(hs, 1e6);
		ok(reg.size() == 4 && close(polyArea(reg), 16), "halfplaneInter rebuilds the square, area 16");

		hs.push_back(HP(Pd(0, 2), Pd(4, 2)));               // ...and cut off the bottom half
		ok(close(polyArea(halfplaneInter(hs, 1e6)), 8), "one more half-plane leaves area 8");

		vector<HP> bad{HP(Pd(0, 0), Pd(1, 0)), HP(Pd(0, -5), Pd(-1, -5))};
		ok(halfplaneInter(bad, 1e6).empty(), "halfplaneInter -> {} when the system is infeasible");

		// x <= 3 as an inequality, i.e. 1*x + 0*y <= 3
		vector<HP> ineq{HP(Pd(0, 0), Pd(4, 0)), HP(Pd(4, 0), Pd(4, 4)),
		                HP(Pd(4, 4), Pd(0, 4)), HP(Pd(0, 4), Pd(0, 0)), hpFromIneq(1, 0, 3)};
		ok(close(polyArea(halfplaneInter(ineq, 1e6)), 12), "hpFromIneq(1,0,3) is `x <= 3`: area drops to 12");
	}
	{
		// an L-shaped polygon: its kernel is the part that can see every corner
		vP L{P(0, 0), P(4, 0), P(4, 2), P(2, 2), P(2, 4), P(0, 4)};
		ok(close(polyArea(polygonKernel(L)), 4), "polygonKernel of the L is the 2x2 corner square");
		ok(close(polyArea(polygonKernel(sq)), 16), "polygonKernel of a convex polygon is the polygon");
	}

	// =====================================================================
	//  closestPair.cpp
	// =====================================================================
	{
		vP v{P(0, 0), P(10, 10), P(11, 10), P(30, 5)};
		auto cp = closestPair(v);
		ok(dist2(v[cp.first], v[cp.second]) == 1, "closestPair -> the indices of (10,10) and (11,10)");
		ok(closestDist2(v) == 1, "closestDist2 -> 1 (SQUARED, and exact)");
	}
	ok(toChebyshev(P(3, 5)) == P(8, -2), "toChebyshev(x,y) = (x+y, x-y): L1 becomes L_inf");
	ok(toManhattan(P(8, -2)) == P(3, 5), "toManhattan is the inverse");
	ok(manhattanFarthest(vP{P(0, 0), P(3, 4), P(-2, 1)}) == 8,
	   "manhattanFarthest: the largest |dx|+|dy| is 8, between (3,4) and (-2,1)");
	{
		vP v{P(0, 0), P(0, 3), P(4, 0), P(4, 3)};
		ok(manhattanMSTEdges(v).size() > 0, "manhattanMSTEdges: the O(n) candidate edges");
		auto mst = manhattanMST(v);
		ok(mst.first == 10 && mst.second.size() == 3, "manhattanMST -> {weight 10, 3 tree edges}");
	}

	// =====================================================================
	//  sweep.cpp
	// =====================================================================
	{
		vector<array<ll, 4>> rs{{0, 0, 4, 4}, {2, 2, 6, 6}};      // {x1,y1,x2,y2}
		ok(rectUnionArea(rs) == 28, "rectUnionArea: 16 + 16 - 4 overlap = 28");
		ok(rectUnionPerimeter(rs) == 24, "rectUnionPerimeter of that L-shaped union = 24");
		ok(rectUnionArea(vector<array<ll, 4>>{{0, 0, 4, 4}, {0, 0, 4, 4}}) == 16,
		   "rectUnionArea counts overlapping copies once");
	}
	{
		vector<array<ll, 3>> bs{{0, 4, 3}, {2, 6, 5}};            // {left, right, height}
		auto sk = skyline(bs);
		ok(sk.size() == 3 && sk[0] == make_pair(0LL, 3LL) && sk[1] == make_pair(2LL, 5LL) &&
		       sk[2] == make_pair(6LL, 0LL),
		   "skyline -> (0,3) (2,5) (6,0): the outline's corner points");
	}
	{
		vector<pair<P, P>> segs{{P(0, 0), P(4, 0)}, {P(2, -2), P(2, 2)}, {P(9, 9), P(9, 10)}};
		auto hit = anySegIntersect(segs);
		ok(hit.first >= 0, "anySegIntersect -> a crossing pair (here segments 0 and 1)");
		ok(anySegIntersect(vector<pair<P, P>>{{P(0, 0), P(1, 0)}, {P(0, 5), P(1, 5)}}).first == -1,
		   "anySegIntersect -> {-1,-1} when nothing crosses");
		ok(countHVIntersect(segs) == 1, "countHVIntersect: 1 horizontal/vertical crossing");
	}
	ok(maxOverlap(vector<pair<ll, ll>>{{1, 5}, {2, 6}, {6, 9}}) == 2, "maxOverlap: at most 2 intervals overlap");

	// =====================================================================
	//  kdtree.cpp
	// =====================================================================
	{
		KDTree T(vP{P(0, 0), P(5, 5), P(1, 1), P(9, 0)});
		ok(T.nearest(P(2, 2)).first == 2, "nearest(p) -> {squared distance 2, the point (1,1)}");
		ok(T.nearest(P(0, 0)).first == 0, "nearest returns p itself when p is in the tree...");
		ok(T.nearest(P(0, 0), true).first == 2, "...unless you pass excludeSelf = true");
		ok(T.kNearest(P(0, 0), 2).size() == 2, "kNearest(p, k) -> k pairs, closest first");
		ok(T.kNearest(P(0, 0), 2)[1].first == 2, "...so [1] is the second closest, at distance^2 2");
		ok(T.countInRect(P(0, 0), P(5, 5)) == 3, "countInRect: 3 points in the inclusive box");
		ok(T.reportInRect(P(0, 0), P(5, 5)).size() == 3, "reportInRect returns those points");
		ok(T.countInCircle(P(0, 0), 2) == 2, "countInCircle(centre, RADIUS) -> 2");
	}

	// =====================================================================
	//  angular.cpp
	// =====================================================================
	{
		vP v{P(1, 0), P(0, 1), P(-1, 0)};
		ok(angularOrder(v, P(0, 0))[0] == 0, "angularOrder -> indices sorted CCW around a centre");
	}
	{
		vP line{P(0, 0), P(1, 1), P(2, 2), P(5, 0)};
		ok(maxCollinearThrough(line, 0) == 3, "maxCollinearThrough(pts, 0): 3 points on a line through pts[0]");
		ok(maxPointsOnLine(line) == 3, "maxPointsOnLine: 3");
	}
	ok(maxPointsInCircle(vP{P(0, 0), P(1, 0), P(0, 1), P(9, 9)}, 1.0) == 3,
	   "maxPointsInCircle: a radius-1 disk can cover 3 of these 4 points");
	ok(countPointsInHalf(vP{P(1, 1), P(1, -1), P(2, 3)}, P(0, 0), P(1, 0)) == 2,
	   "countPointsInHalf: 2 points strictly left of the ray (0,0)->(1,0)");
	ok(countTrianglesContaining(vP{P(3, 0), P(-3, 2), P(-3, -2), P(5, 5)}, P(0, 0)) == 2,
	   "countTrianglesContaining: 2 of the 4 possible triangles contain the origin");

	// =====================================================================
	//  geom3d.cpp
	// =====================================================================
	P3 u(1, 0, 0), v(0, 1, 0);
	ok(close(dot3(u, v), 0), "dot3");
	ok(cross3(u, v) == P3(0, 0, 1), "cross3: x cross y = z");
	ok(close(len3sq(P3(1, 2, 2)), 9) && close(len3(P3(1, 2, 2)), 3), "len3sq / len3");
	ok(close(dist3(P3(0, 0, 0), P3(1, 2, 2)), 3), "dist3");
	ok(close(len3(unit3(P3(1, 2, 2))), 1), "unit3");
	ok(close(tripleProduct(u, v, P3(0, 0, 1)), 1), "tripleProduct = 6 x tetra volume = 1");
	ok(close(tetraVolume(P3(0, 0, 0), u, v, P3(0, 0, 1)), 1.0 / 6), "tetraVolume = 1/6");
	{
		Plane pl = planeFrom3(P3(0, 0, 0), u, v);              // the z = 0 plane
		ok(close(distToPlane(P3(0, 0, 5), pl), 5), "distToPlane is SIGNED");
		ok(close(projectOnPlane(P3(2, 3, 7), pl).z, 0), "projectOnPlane drops it onto the plane");
		auto hit = linePlaneInter(P3(0, 0, 5), P3(0, 0, -5), pl);
		ok(hit.first && close(hit.second.z, 0), "linePlaneInter -> {true, (0,0,0)}");
	}
	ok(lineSphere(P3(-9, 0, 0), P3(9, 0, 0), P3(0, 0, 0), 2).size() == 2, "lineSphere -> 2 points");
	ok(near(Pd(rotate3(u, P3(0, 0, 1), PI / 2).x, rotate3(u, P3(0, 0, 1), PI / 2).y), Pd(0, 1)),
	   "rotate3(p, axis, theta): Rodrigues rotation");
	ok(close(sphericalDist(0, 0, 0, 180, 1), PI), "sphericalDist(lat1,lon1,lat2,lon2,R), degrees in");
	{
		vP3 cube;
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				for (int k = 0; k < 2; k++) cube.push_back(P3(3 * i, 3 * j, 3 * k));
		auto f = hull3d(cube);
		ok(f.size() == 12, "hull3d: a cube comes back as 12 TRIANGLES (faces are not merged)");
		ok(close(hullVolume(cube, f), 27), "hullVolume = 27");
		ok(close(hullArea(cube, f), 54), "hullArea = 54");
		ok(pointInHull(cube, f, P3(1.5, 1.5, 1.5)) == 2, "pointInHull -> 2 inside");
		ok(pointInHull(cube, f, P3(0, 1.5, 1.5)) == 1, "pointInHull -> 1 on the surface");
		ok(pointInHull(cube, f, P3(-1, 1.5, 1.5)) == 0, "pointInHull -> 0 outside");
	}

	// =====================================================================
	//  delaunay.cpp
	// =====================================================================
	{
		vP v{P(0, 0), P(4, 0), P(4, 4), P(0, 4), P(2, 1)};
		auto tri = delaunay(v);
		ld area = 0;
		for (auto &t : tri) area += fabs((ld)cross(v[t[0]], v[t[1]], v[t[2]])) / 2;
		ok(close(area, 16), "delaunay: the triangles exactly tile the convex hull (area 16)");
		ok(delaunayEdges(v).size() >= 5, "delaunayEdges: each edge once, as i < j");

		auto mst = euclideanMST(vP{P(0, 0), P(0, 3), P(4, 3), P(4, 0)});
		ok(close(mst.first, 10) && mst.second.size() == 3, "euclideanMST -> {10, 3 edges}");

		ok(voronoiVertices(v, tri).size() == tri.size(), "voronoiVertices: one circumcentre per triangle");
		ok(nearestNeighborOf(v, 4, tri) == 0, "nearestNeighborOf: (2,1)'s nearest neighbour is (0,0)");
	}

	cout << (bad ? "MISMATCHES: " : "ALL EXAMPLES OK ") << bad << "\n";
	return bad != 0;
}

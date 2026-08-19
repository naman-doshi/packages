# CompGeo

Computational geometry. Every file opens with a `WHAT / WHEN` → `API` →
`COMPLEXITY` → `PITFALLS` banner; read the pitfalls, that's where the WA lives.

## How to use it

**`point.cpp` is the base. Paste it first — every other file needs it.** Then
paste whichever module the problem wants. Nothing else has hidden dependencies
beyond what its banner lists.

```
point.cpp                                          -> vectors, cross/dot, polar sort
point.cpp + lines.cpp                              -> intersections, distances
point.cpp + lines.cpp + polygon.cpp                -> areas, containment, Pick's
point.cpp + lines.cpp + polygon.cpp + convex.cpp   -> hulls, calipers, tangents
point.cpp + lines.cpp + circles.cpp                -> circles (+ delaunay.cpp on top)
point.cpp alone is enough for  halfplane / closestPair / kdtree / angular / geom3d
```

Two point types, both from `point.cpp`:

| | when | exact? |
|---|---|---|
| `P`  = `Pt<ll>`     | the input coordinates are integers — **the default** | yes, no epsilons anywhere |
| `Pd` = `Pt<double>` | circles, intersections, rotations, anything computed | no, `EPS = 1e-9` |

A `P` converts to a `Pd` on its own; a `Pd` never converts back to a `P`. So
`distToLine(circle.c, latticeA, latticeB)` just compiles, and you can't
accidentally truncate a computed point onto the lattice.

`examples.cpp` has a worked, checked call for **every** function in this folder.
It is runnable and self-verifying — build it with all the modules and it prints
`ALL EXAMPLES OK`:

```sh
cat point.cpp lines.cpp polygon.cpp convex.cpp circles.cpp halfplane.cpp \
    closestPair.cpp sweep.cpp kdtree.cpp angular.cpp geom3d.cpp \
    delaunay.cpp examples.cpp > /tmp/ex.cpp
g++ -std=c++17 -O2 -o /tmp/ex /tmp/ex.cpp && /tmp/ex
```

## Index

| File | Contents |
|---|---|
| `point.cpp` | **base.** `Pt<T>`, `cross`/`dot`/`orient`, rotation, angles, polar sort, lattice helpers |
| `lines.cpp` | line & segment intersection, projection, reflection, all the distances, `ax+by=c`, counting the points of a set on a line |
| `polygon.cpp` | area, centroid, point-in-polygon, winding, cutting, clipping, ear clipping, Pick's theorem |
| `convex.cpp` | convex hull, rotating calipers (diameter / width / min rectangle), O(log n) queries, Minkowski, onion layers |
| `circles.cpp` | intersections, tangents, circumcircle, min enclosing circle, circle∩polygon and circle-union areas, exact `inCircleDet` |
| `halfplane.cpp` | half-plane intersection, linear feasible region, polygon kernel |
| `closestPair.cpp` | closest pair, Chebyshev/Manhattan rotation, Manhattan MST |
| `sweep.cpp` | rectangle union area & perimeter, skyline, segment-intersection detection, H×V crossing count |
| `kdtree.cpp` | nearest / k-nearest neighbour, points in a box or disk |
| `angular.cpp` | max points on a line, max points in a radius-r disk, triangles containing a point |
| `geom3d.cpp` | 3-D points, planes, spheres, 3-D convex hull, volume |
| `delaunay.cpp` | Delaunay (exact), Voronoi vertices, **Euclidean MST**, nearest-neighbour graph |
| `examples.cpp` | a checked usage example for every function |
| `formulas.md` | formula sheet: identities, areas, angles, the standard tricks |

## Which file do I want?

**"is this point inside / which side"**
- one polygon, few queries → `inPolygon` (`polygon.cpp`), O(n), tells you
  outside / boundary / inside
- convex polygon, many queries → `inConvex` (`convex.cpp`), O(log n)
- self-intersecting polygon → `windingNumber`
- three points, just the turn direction → `orient`. Never compare slopes.

**"largest / farthest / extreme"** → hull it first (`convexHull`), the answer is
always on the hull.
- farthest pair → `farthestPair` / `hullDiameter`
- max of `dot(dir, p)` → `maxDot`, O(log n)
- biggest inscribed triangle → `maxTriangle`
- smallest enclosing rectangle → `minAreaRect` / `minPerimeterRect`
- smallest enclosing circle → `minEnclosingCircle` (`circles.cpp`)

**"closest"**
- closest pair of a point set → `closestPair` (`closestPair.cpp`), O(n log n)
- nearest neighbour, online queries → `kdtree.cpp`
- nearest neighbour for every point → `delaunay.cpp`
- point to a segment / polygon → `distToSeg`, `distToPoly`

**"area of the union / overlap"**
- axis-aligned rectangles → `rectUnionArea` (`sweep.cpp`)
- two convex polygons → `convexInterArea`
- a circle and a polygon → `circlePolyArea`
- many circles → `circleUnionArea`

**"minimum spanning tree of points"**
- Manhattan / L1 → `manhattanMST`, O(n log n)
- Euclidean → `euclideanMST` (`delaunay.cpp`), O(n²) here — check the problem
  isn't secretly Manhattan first

**"count lattice points"** → Pick's theorem: `interiorLattice` /
`boundaryLattice`. Never loop over the bounding box.

**"a line / circle placed anywhere, maximise coverage"** → angular sweep:
`maxPointsOnLine`, `maxPointsInCircle` (`angular.cpp`).

**"how many of these points are on THIS line"** (`lines.cpp`)
- one line, or a handful → `countOnLine`, O(n). `pointsOnLine` if you want them
- many lines → `OnLine tab(pts)`, O(log n) per query, plus `tab.any` for "is
  there at least one". Its 0-or-1 case indexes per direction, so it is at its
  best when the query lines reuse directions
- many lines, adversarial directions, or you also want **half-plane counts** →
  `LineSweep sw(pts)`: `count` / `any` / `countLeft` / `countRight`, O(log² n)
  each with no cache to miss, after an O(n² log n) build (n ≤ ~2000)

**"these linear constraints"** → `halfplaneInter` (`halfplane.cpp`); build each
one with `hpFromIneq(A, B, C)` for `Ax + By <= C`.

**"do any two segments cross"** → `anySegIntersect` (`sweep.cpp`), O(n log n).
If they're all horizontal/vertical and you need a *count*, `countHVIntersect`.

## Conventions

- Polygons are a `vector<P>` in order, **with no repeated last vertex**. Most
  routines want CCW — call `makeCCW(poly)` once and stop worrying.
- Anything returning an area returns **twice** the area when it can stay an
  exact integer (`polyArea2`, `area2`, `maxTriangle`). Divide by `2.0` at the end.
- Anything returning a distance returns the **squared** distance when it can
  stay exact (`dist2`, `len2`, `closestDist2`). Don't `sqrt` until you print.
- Containment predicates return **0 = outside, 1 = on the boundary, 2 = inside**,
  so `if (inPolygon(...))` means "inside or on". Decide deliberately which one
  the problem wants.
- These files do **not** `#define int long long`, so they behave identically
  whether or not your template has that macro on. They define no macros at all,
  so nothing collides with `rep` / `all` / `f` / `s`.
- `std::gcd` is used, never `__gcd` (Apple's libc++ static-asserts on it).

## Overflow — read this once

`cross` is the whole library, and it computes a product of coordinate
*differences*. With `|x|, |y| <= C` the differences reach `2C` and the cross
product reaches `8C²`:

| C | max `cross` | verdict |
|---|---|---|
| 1e4 | 8e8 | fine |
| 1e6 | 8e12 | fine |
| 1e9 | 8e18 | **fits in `ll`, but with almost nothing to spare** |
| >1e9 | — | change `Pt<ll>` to `Pt<__int128>` |

`dist2` reaches `8C²` too, so at `C = 1e9` you must not add two of them
together. `inCircleDet` is degree 4 and already uses `__int128`; it is good to
`|coord| <= ~1e9` and no further.

## What was checked

Everything here was validated against brute force before being committed —
convex hull and all the calipers against O(n³) enumeration, `inConvex` against
`inPolygon` on every lattice point, Pick's theorem against counting, segment
intersection against exact rational solving, Delaunay against the exact
in-circle predicate on dense grids (the degenerate case), Manhattan and
Euclidean MST against Prim, the rectangle sweeps against a unit-cell grid, and
the circle/polygon areas against Monte Carlo. `examples.cpp` re-checks every
documented return value on every run.

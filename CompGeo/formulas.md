# Geometry formula sheet

Everything you'd otherwise re-derive at 2am. Notation: `a × b = ax·by − ay·bx`
(cross), `a · b = ax·bx + ay·by` (dot).

## The two primitives

| | meaning |
|---|---|
| `a × b` | signed area of the parallelogram. `> 0` ⟺ **b is counter-clockwise from a** |
| `a · b` | `\|a\|\|b\|cos θ`. `> 0` same general direction, `0` perpendicular, `< 0` opposite |
| `(b−a) × (c−a)` | `> 0` ⟺ **a→b→c turns left** (CCW). This is `orient(a,b,c)`. |

Almost every predicate is one of these. **If you are writing a slope, a
division, or an `atan2` inside a comparison, stop** — there is a cross product
that answers it exactly.

- `a × b = 0` ⟺ parallel (collinear through the origin)
- `a · b = 0` ⟺ perpendicular
- `θ` between a and b: `atan2(|a × b|, a · b)` ∈ [0, π] — better conditioned
  than `acos(dot/(|a||b|))`, which loses all precision near 0 and π
- signed θ from a to b: `atan2(a × b, a · b)` ∈ (−π, π]
- rotate 90° CCW: `(x, y) → (−y, x)`
- rotate by θ: `(x cos θ − y sin θ, x sin θ + y cos θ)`

## Areas

- **Shoelace**: `2·Area = Σ (x_i·y_{i+1} − x_{i+1}·y_i)` over the cycle.
  Positive ⟺ CCW. Keep the factor 2 and stay in integers.
- Triangle: `2·Area = (b−a) × (c−a)`
- Heron: `Area = √(s(s−a)(s−b)(s−c))`, `s = (a+b+c)/2` — numerically bad for
  thin triangles; prefer the cross product.
- Centroid of a polygon (**not** the average of the vertices):
  `C = (1/(6A)) Σ (p_i + p_{i+1})·(x_i·y_{i+1} − x_{i+1}·y_i)`
- Trapezoid form: the area between a polygon edge and the x-axis is
  `(x2 − x1)(y1 + y2)/2` — the basis of every sweep.

## Lattice points — Pick's theorem

For a simple polygon with **integer** vertices:

```
A = I + B/2 − 1          I = interior points, B = boundary points
```

- `B = Σ gcd(|dx|, |dy|)` over the edges
- `I = A − B/2 + 1 = (2A − B + 2)/2`
- lattice points strictly inside the segment `(x1,y1)–(x2,y2)`:
  `gcd(|dx|, |dy|) − 1`

This is *the* way to count lattice points. Never iterate the bounding box.

## Lines

- Through `p, q`: `(y_q − y_p)x + (x_p − x_q)y = (y_q − y_p)x_p + (x_p − x_q)y_p`
- Distance from `p` to the line `ab`: `|(b−a) × (p−a)| / |b−a|`  (drop the `| |`
  for a signed distance — positive on the left)
- Projection of `p` onto line `ab`: `a + (b−a)·((p−a)·(b−a) / |b−a|²)`
- Reflection: `2·proj(p) − p`
- Two lines `a+t·u` and `c+s·v` meet at `t = ((c−a) × v) / (u × v)`.
  `u × v = 0` ⟺ parallel.
- Segments `[a,b]`, `[c,d]` **properly** cross ⟺
  `orient(a,b,c)·orient(a,b,d) < 0` **and** `orient(c,d,a)·orient(c,d,b) < 0`.
  Touching and collinear overlap need the extra `onSegment` checks.

## Circles

- Through 3 points — circumcentre is where the perpendicular bisectors meet;
  `R = abc / (4·Area)`
- Incentre `= (a·A + b·B + c·C)/(a+b+c)` with `a = |BC|` etc.;  `r = Area / s`
- Chord of length `L` in radius `r`: half-angle `θ/2 = asin(L/(2r))`
- Circular **sector** area (angle θ): `r²θ/2`
- Circular **segment** area (cut off by the chord): `r²(θ − sin θ)/2`
- Two circles at distance `d` intersect ⟺ `|r1 − r2| ≤ d ≤ r1 + r2`.
  The radical line is at `x = (d² + r1² − r2²)/(2d)` from centre 1, and the
  half-chord is `√(r1² − x²)`.
- Tangent length from an external `p`: `√(|p − c|² − r²)`
- **In-circle test** (exact, for lattice points, `a,b,c` CCW): `d` is strictly
  inside the circle through `a,b,c` ⟺

```
| ax−dx  ay−dy  (ax−dx)²+(ay−dy)² |
| bx−dx  by−dy  (bx−dx)²+(by−dy)² |  > 0
| cx−dx  cy−dy  (cx−dx)²+(cy−dy)² |
```

  Degree 4 — needs `__int128`. This is the Delaunay predicate.

## Convexity

- A polygon is convex ⟺ all `orient(p_i, p_{i+1}, p_{i+2})` share one sign.
- A point set's convex hull has ≤ n vertices, and **every** "maximise a linear
  function / find the farthest thing" answer lies on it.
- **Rotating calipers**: as an edge direction rotates through 360°, the
  supporting vertex on the far side rotates monotonically too — so the
  two-pointer never backtracks, giving O(n) for diameter, width, minimum
  enclosing rectangle, and the distance between two convex polygons.
- The minimum-area enclosing rectangle always has one side **flush with a hull
  edge**. Same for minimum perimeter. So: try every edge.
- **Minkowski sum** `A ⊕ B = {a + b}`: for convex polygons, merge their edge
  vectors in angular order. Size `|A| + |B|`. Uses:
  - `A ∩ B ≠ ∅` ⟺ `0 ∈ A ⊕ (−B)`
  - distance between convex `A` and `B` = distance from the origin to `A ⊕ (−B)`
  - "can this shape pass through" / collision detection
- The centroid, the hull, and the diameter are unchanged by translation —
  shift coordinates to be small before crossing anything.

## Angular sweep

The pattern: **fix one point, sort the rest by angle around it, sweep.**
Recognise it from "a line through", "a circle of radius r anywhere", "sees",
"contains the origin", "collinear".

- Sort by angle with `half()` + cross product, **never** by `atan2` — exact,
  and twice as fast.
- Points in a half-plane through `q`: sort by angle, then a two-pointer over a
  window of π. Complement-count with it: the number of triangles **not**
  containing `q` is `Σ C(k_i, 2)` where `k_i` is the number of points inside the
  half-plane starting at point `i`; so the answer is `C(n,3) − Σ C(k_i, 2)`.
- A disk of radius `r` covering the most points: an optimal disk can always be
  slid until a point is on its boundary. Fix that point, sweep the centre around
  a circle of radius `r`; point `j` is covered over an arc of half-width
  `acos(d/2r)` centred on the direction to `j`.
- Most collinear points: for each `i`, bucket the directions `(dx, dy)/gcd`
  with a canonical sign. O(n² log n), exact.

## L1 / L∞ (Manhattan and Chebyshev)

```
(x, y) → (x + y, x − y)      turns L1 distance into L∞ distance
(u, v) → ((u+v)/2, (u−v)/2)  and back
```

- `L1(p,q) = max(|(x_p+y_p) − (x_q+y_q)|, |(x_p−y_p) − (x_q−y_q)|)` after the
  rotation — so the **L1 diameter is O(n)**: take the max minus the min of
  `x+y` and of `x−y`.
- L∞ distance = `max(|dx|, |dy|)`; L1 = `|dx| + |dy|`.
- Manhattan MST: only **4n** candidate edges matter (each point's nearest
  neighbour in each of 8 octants; by symmetry 4 sweeps suffice). O(n log n).
- The smallest edge of the Manhattan MST is the L1 closest pair.

## Sweep line

- **Rectangle union area**: sweep x, keep the covered y-length in a segment tree
  with a "covered count" per node. `area += covered · Δx`.
- **Rectangle union perimeter**: same sweep, both directions, summing
  `|Δ covered length|` — but process events **one at a time, additions before
  removals**. Grouping a whole x together loses the edges of two boxes touching
  at a corner; removing before adding double-counts the shared edge of two
  abutting boxes.
- **Do any two segments cross** (Shamos–Hoey): sweep x with a `set` ordered by
  the y at the current sweep position; a crossing pair must become adjacent in
  that order at some point, so only check new neighbours. O(n log n).
- **Counting** crossings needs Bentley–Ottmann (O((n+k) log n)) and is almost
  never intended — if a count is wanted, the segments are usually axis-parallel,
  which is a BIT sweep.

## Delaunay / Voronoi

- Delaunay ⟺ no point strictly inside any triangle's circumcircle.
- Lift `(x, y) → (x, y, x² + y²)`; the Delaunay triangulation is the **lower**
  convex hull in 3-D. (Elegant, but it makes co-circular points degenerate —
  `delaunay.cpp` uses exact incremental insertion instead.)
- Voronoi is its **dual**: one Voronoi vertex per Delaunay triangle (its
  circumcentre); one Voronoi edge per Delaunay edge.
- Delaunay contains: the Euclidean MST, the nearest-neighbour graph, and the
  convex hull (as its outer boundary).
- Sizes: a Delaunay triangulation of `n` points has ≤ `2n − 5` triangles and
  ≤ `3n − 6` edges.

## Euler's formula

For a connected planar graph: `V − E + F = 2` (`F` counts the outer face).
Handy for "how many regions do these lines/segments cut the plane into":
`n` lines in general position make `1 + n + C(n,2)` regions.

## Precision, in order of preference

1. **Stay in integers.** Use `P` = `Pt<ll>`, compare cross products, keep areas
   doubled and distances squared. Most geometry problems are solvable this way.
2. `__int128` when the degree gets to 4 (in-circle) or the coordinates exceed 1e9.
3. Rationals (`Math/frac.cpp`) for exact intersection coordinates — an
   intersection of integer segments has denominator `(b−a) × (d−c)`.
4. `double` with a **relative** epsilon, and only at the end.

Danger signs: `acos` near ±1, subtracting nearly equal large numbers, comparing
`atan2` values, and normalising a near-zero vector. Scale and translate the
input to be small and centred before doing any of it.

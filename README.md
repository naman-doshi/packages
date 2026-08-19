# Packages

Contest library. Every file opens with a `WHAT / WHEN` → `API` → `COMPLEXITY` →
`PITFALLS` banner; read the pitfalls before pasting.

```
Bitwise/      bit tricks and identities
CompGeo/      all computational geometry          (see CompGeo/README.md)
DP/           classic DP shapes, incl. bitmask/TSP
flow/         max flow, min cost flow, matching
Geometry/     superseded by CompGeo/ -- hull.cpp is CompGeo/convex.cpp now
Graph/        traversal, shortest paths, MST, SCC, 2-SAT, ...
Math/         number theory, combinatorics, matrices, FFT   (see Math/README.md)
Segtree/      every segment tree variant
Strings/      hashing
Structures/   fenwick, sparse table, ordered set, Li Chao
Techniques/   binary search, monotonic stack
Tree/         LCA, HLD, centroid, Euler tour
```

## Segtree/

| File | Use it for |
|---|---|
| `seg.cpp` | point update, range query. The default. Start here. |
| `lazyseg.cpp` | **range** update + range query |
| `dseg.cpp` | huge/unbounded index space, nodes created on demand |
| `logspaceseg.cpp` | block leaves; saves a log factor of memory when you need many trees |
| `seg2d.cpp` | rectangle query with updates, non-invertible ops (min/max/gcd) |
| `pseg.cpp` | persistent — every past version stays queryable |
| `pseg2d.cpp` | static 2-D points: count in rectangle, k-th y |
| `parr.cpp` | persistent array (a pseg over indices) |
| `mergeSortTree.cpp` | static "count values ≤ k in a range", simpler than pseg |
| `pointLocation.cpp` | tree descent: leftmost/rightmost index beating a threshold |
| `customSegtrees.cpp` | **guide**: how to change *what* any of them compute |

Rule of thumb: invertible op (sum/xor) and point updates → use `Structures/fenwick.cpp`
instead, it's shorter and ~3x faster. Static array + idempotent op (min/max/gcd)
→ `Structures/sparse.cpp` for O(1) queries.

## Structures/

| File | Use it for |
|---|---|
| `fenwick.cpp` | point add, prefix/range sum. The workhorse. |
| `fenwick2d.cpp` | point add, rectangle sum |
| `sparse.cpp` | static array, O(1) idempotent range query |
| `orderedSet.cpp` | pbds tree: k-th smallest, rank of a key |
| `liChao.cpp` | min over a set of lines (CHT with unsorted slopes) |

## CompGeo/

Paste `point.cpp` first — it defines `P` (exact integer points) and `Pd`
(floating) and everything else builds on it. Then the module you need:
`polygon.cpp`, `convex.cpp` (hull + calipers), `circles.cpp`, `halfplane.cpp`,
`closestPair.cpp` (incl. Manhattan MST), `sweep.cpp` (rectangle union, skyline),
`kdtree.cpp`, `angular.cpp`, `geom3d.cpp`, `delaunay.cpp` (incl. Euclidean MST).
`examples.cpp` has a checked usage example for every function, `formulas.md` is
the formula sheet. Full index and an archetype→file guide in `CompGeo/README.md`.

## Techniques/

| File | Use it for |
|---|---|
| `binarySearch.cpp` | boundary of a monotone predicate; "binary search the answer" |
| `histogram.cpp` | largest rectangle under a skyline (monotonic stack) |

## Notes

- Most files assume the template's `#define int long long`, plus `vi`, `pii`,
  `INF`. A few (`seg.cpp`, `dseg.cpp`, `logspaceseg.cpp`, `parr.cpp`) are pure
  snippets that won't compile alone — that's deliberate, they inherit those
  names from the file you paste them into.
- `orderedSet.cpp` needs libstdc++ (`ext/pb_ds`), so it compiles on the judge
  but **not** with Apple clang locally.
- Use `std::gcd` / `std::lcm`, never `__gcd` — Apple's libc++ static-asserts on
  signed types and won't compile it.

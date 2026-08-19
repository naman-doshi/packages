// ============================================================================
//  2-D FENWICK TREE (BIT)  -- point add + rectangle sum, O(log^2)
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    A grid version of the Fenwick tree: add to a cell, sum any axis-aligned
//    rectangle, both in O(log R * log C). The default choice for DYNAMIC 2-D
//    "point update, rectangle sum" (counting points in a region while inserting
//    them, 2-D prefix sums with updates, etc). Like 1-D BIT it only does
//    INVERTIBLE ops (sum, xor) -- for rectangle min/max use seg2d.cpp.
//
//  API  (0-based coords)
//    Fenwick2D f(R, C);
//    f.add(r, c, delta);                 // grid[r][c] += delta
//    f.rect(r1, c1, r2, c2);             // sum over the inclusive rectangle
//
//  COMPLEXITY  O(log R * log C) per op; O(R*C) memory.
//
//  PITFALLS
//    * Memory is O(R*C) -- if coordinates are large/sparse, COORDINATE-COMPRESS
//      first (map the distinct x's and y's to 0..k), or go offline: sort queries
//      by one axis and sweep with a 1-D BIT over the other.
//    * OVERFLOW: use long long.
//    * rect() is inclusion-exclusion of four prefixes; sumPrefix(-1, .) = 0 is
//      handled by the loop (r+1 == 0 -> no iterations).
// ============================================================================
#include <bits/stdc++.h>
using namespace std;

struct Fenwick2D {
    int R, C;
    vector<vector<long long>> bit;   // 1-indexed internally
    Fenwick2D(int R, int C) : R(R), C(C), bit(R + 1, vector<long long>(C + 1, 0)) {}

    void add(int r, int c, long long delta) {           // 0-based cell
        for (int i = r + 1; i <= R; i += i & -i)
            for (int j = c + 1; j <= C; j += j & -j)
                bit[i][j] += delta;
    }
    long long sumPrefix(int r, int c) {                 // sum over [0..r] x [0..c]
        long long s = 0;
        for (int i = r + 1; i > 0; i -= i & -i)
            for (int j = c + 1; j > 0; j -= j & -j)
                s += bit[i][j];
        return s;
    }
    long long rect(int r1, int c1, int r2, int c2) {    // inclusive rectangle sum
        if (r1 > r2 || c1 > c2) return 0;
        return sumPrefix(r2, c2) - sumPrefix(r1 - 1, c2)
             - sumPrefix(r2, c1 - 1) + sumPrefix(r1 - 1, c1 - 1);
    }
};

/*  USAGE
      Fenwick2D f(1000, 1000);
      f.add(3, 5, 1);            // put a point / weight at (3,5)
      f.add(7, 2, 4);
      f.rect(0, 0, 5, 5);        // sum (count) in rectangle [0..5]x[0..5] -> 1
    Counting points in a rectangle while inserting them = add(y,x,1) then rect(...).
*/

// ============================================================================
//  2-D SEGMENT TREE  -- point update + rectangle query, O(log^2)
//  (segment tree of segment trees; this copy does rectangle MAX)
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Point-assign a cell and query the aggregate over any rectangle in
//    O(log R * log C). Use it over a 2-D Fenwick when the op is NOT invertible
//    -- rectangle MIN / MAX / gcd with updates (a BIT can't do those). For
//    rectangle SUM with updates, prefer fenwick2d.cpp (lighter/faster).
//
//  API  (0-based coords; this copy = rectangle MAX)
//    Seg2D t(R, C);                      // all cells start at LLONG_MIN
//    t.modify(x, y, v);                  // set cell (x,y) = v  (point assign)
//    t.query(x1, x2, y1, y2);            // max over the inclusive rectangle
//
//  CUSTOMISE  it's a monoid in 2-D. To switch the op, change the merge in BOTH
//    modify() and query() and the fill value:
//      max -> LLONG_MIN, max(a,b)     min -> LLONG_MAX, min(a,b)     sum -> 0, a+b
//
//  COMPLEXITY  build/modify/query all O(log R * log C); memory O(R * C) (4x a grid).
//
//  PITFALLS
//    * Coordinates must be dense (0..R, 0..C). Large/sparse -> compress first.
//    * modify is a POINT ASSIGN (overwrite), not add. For "add", read-modify or
//      switch to a sum tree and add deltas.
//    * Empty cells read as the fill value (LLONG_MIN for max) -- guard if you
//      query rectangles that may contain never-set cells.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;

struct Seg2D {
    int R, C;
    vector<vector<long long>> t;                 // 2R x 2C, iterative
    Seg2D(int R, int C) : R(R), C(C), t(2 * R, vector<long long>(2 * C, LLONG_MIN)) {}

    void modify(int X, int Y, long long v) {     // point assign at (X,Y)
        int x = X + R, y = Y + C;
        t[x][y] = v;
        for (int j = y >> 1; j >= 1; j >>= 1) t[x][j] = max(t[x][2 * j], t[x][2 * j + 1]);
        for (int i = x >> 1; i >= 1; i >>= 1) {
            t[i][y] = max(t[2 * i][y], t[2 * i + 1][y]);
            for (int j = y >> 1; j >= 1; j >>= 1) t[i][j] = max(t[i][2 * j], t[i][2 * j + 1]);
        }
    }
    long long queryRow(int row, int y1, int y2) {         // max over columns [y1,y2] of a row-node
        long long res = LLONG_MIN;
        for (y1 += C, y2 += C + 1; y1 < y2; y1 >>= 1, y2 >>= 1) {
            if (y1 & 1) res = max(res, t[row][y1++]);
            if (y2 & 1) res = max(res, t[row][--y2]);
        }
        return res;
    }
    long long query(int x1, int x2, int y1, int y2) {     // max over inclusive rectangle
        long long res = LLONG_MIN;
        for (x1 += R, x2 += R + 1; x1 < x2; x1 >>= 1, x2 >>= 1) {
            if (x1 & 1) res = max(res, queryRow(x1++, y1, y2));
            if (x2 & 1) res = max(res, queryRow(--x2, y1, y2));
        }
        return res;
    }
};

/*  USAGE
      Seg2D t(100, 100);
      t.modify(3, 5, 42);
      t.modify(10, 20, 7);
      t.query(0, 50, 0, 50);      // max over rectangle [0..50]x[0..50] -> 42
*/

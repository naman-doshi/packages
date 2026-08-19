// ============================================================================
//  POINT LOCATION over a MAX-SEGTREE  -- "leftmost column beating a threshold"
//  (the segment-tree-descent trick behind lab3/D_Points)
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Keep one value per column x (e.g. the MAX y of the points currently in that
//    column). Then answer, in O(log n):
//      * update(i, v)              : set column i's stored value
//      * query(l, r)              : max stored value over columns [l,r]
//      * firstGreater(l, r, x)    : LEFTMOST column in [l,r] whose value > x (-1 if none)
//      * lastGreater(l, r, x)     : RIGHTMOST such column
//    firstGreater walks the tree, pruning subtrees whose max <= x, so it finds
//    the boundary column in O(log n) instead of a binary search that costs
//    O(log^2 n). This is how you do dominance / "next point above and to the
//    right" queries: store max-y per column, then firstGreater(qx+1, N-1, qy).
//
//  API  (0-based columns)
//    PointColMax s(n);          // n columns, all -INF
//    s.update(i, y);            // column i now has value y
//    s.firstGreater(l, r, x);   // leftmost column in [l,r] with value > x
//
//  COMPLEXITY  update / query / firstGreater all O(log n).
//
//  PITFALLS
//    * The stored value is ONE number per column; if a column holds several
//      points, keep the extreme you query on (e.g. max y) and, after locating
//      the column, look up the exact point in a per-column sorted set/BIT.
//    * Compress x-coordinates to 0..n-1 first if they are large.
//    * firstGreater is strict (> x); use >= by testing `>= x` in the prune.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;

struct PointColMax {
    int sz;
    vector<long long> seg;
    PointColMax(int n) { sz = 1; while (sz < max(1, n)) sz <<= 1; seg.assign(2 * sz, LLONG_MIN); }

    void update(int i, long long v) {
        for (seg[i += sz] = v, i >>= 1; i; i >>= 1)
            seg[i] = max(seg[2 * i], seg[2 * i + 1]);
    }
    long long query(int lo, int hi, int i = 1, int l = 0, int r = -1) {
        if (r < 0) r = sz - 1;
        if (lo > r || hi < l) return LLONG_MIN;
        if (lo <= l && r <= hi) return seg[i];
        int m = (l + r) / 2;
        return max(query(lo, hi, 2 * i, l, m), query(lo, hi, 2 * i + 1, m + 1, r));
    }
    // leftmost index in [lo,hi] whose value > x, or -1
    int firstGreater(int lo, int hi, long long x, int i = 1, int l = 0, int r = -1) {
        if (r < 0) r = sz - 1;
        if (lo > r || hi < l || seg[i] <= x) return -1;     // prune whole subtree
        if (l == r) return l;
        int m = (l + r) / 2;
        int res = firstGreater(lo, hi, x, 2 * i, l, m);
        if (res != -1) return res;
        return firstGreater(lo, hi, x, 2 * i + 1, m + 1, r);
    }
    // rightmost index in [lo,hi] whose value > x, or -1
    int lastGreater(int lo, int hi, long long x, int i = 1, int l = 0, int r = -1) {
        if (r < 0) r = sz - 1;
        if (lo > r || hi < l || seg[i] <= x) return -1;
        if (l == r) return l;
        int m = (l + r) / 2;
        int res = lastGreater(lo, hi, x, 2 * i + 1, m + 1, r);
        if (res != -1) return res;
        return lastGreater(lo, hi, x, 2 * i, l, m);
    }
};

/*  USAGE  (like lab3/D_Points: find a point strictly right of and above (qx,qy))
      PointColMax s(n);                 // columns = compressed x, value = max y in column
      s.update(colOfX, y);              // when a point is added
      int col = s.firstGreater(qxIndex + 1, n - 1, qy);   // leftmost column right of qx with some y>qy
      // then within that column's sorted y-set, upper_bound(qy) is the point.
*/

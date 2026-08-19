// ============================================================================
//  MERGE SORT TREE  -- static "count points in a rectangle" via sorted lists
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    A segment tree over array indices where each node stores the SORTED list of
//    the values in its index range. Answers, on a STATIC array, in O(log^2 n):
//      * countLE(l, r, k)         : how many i in [l,r] have a[i] <= k
//      * countRange(l,r,lo,hi)    : how many i in [l,r] have a[i] in [lo,hi]
//    Treat a[i] as the y of a point at x = i: countRange is exactly "points in
//    the rectangle x in [l,r], y in [lo,hi]". Simpler to write than a persistent
//    tree; slightly slower (extra log from the binary search per node).
//
//  IDEA  A query range decomposes into O(log n) tree nodes; in each fully-covered
//    node the count of values <= k is one binary search on its sorted list.
//
//  COMPLEXITY  build O(n log n) time & memory; query O(log^2 n).
//
//  PITFALLS
//    * STATIC only (values fixed). For updates use a BIT of sorted structures or
//      a wavelet/merge-sort tree with fractional cascading -- much heavier.
//    * countRange uses countLE(hi) - countLE(lo-1); mind inclusive bounds.
//    * O(n log n) memory: ~n*log2(n) longs. Fine to a few million.
//
//  ALTERNATIVE  For k-th smallest in a range, persistent segtree (pseg2d.cpp)
//    is the cleaner O(log n) tool; merge-sort tree does it in O(log^3 n) by
//    binary-searching the answer, so prefer pseg2d for k-th.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;

struct MergeSortTree {
    int n;
    vector<vector<long long>> tr;

    MergeSortTree(const vector<long long>& a) {
        n = a.size();
        tr.assign(4 * max(1, n), {});
        if (n) build(1, 0, n - 1, a);
    }
    void build(int node, int l, int r, const vector<long long>& a) {
        if (l == r) { tr[node] = {a[l]}; return; }
        int m = (l + r) / 2;
        build(2 * node, l, m, a);
        build(2 * node + 1, m + 1, r, a);
        merge(tr[2 * node].begin(), tr[2 * node].end(),
              tr[2 * node + 1].begin(), tr[2 * node + 1].end(),
              back_inserter(tr[node]));
    }
    int countLE(int node, int l, int r, int ql, int qr, long long k) {
        if (qr < l || r < ql) return 0;
        if (ql <= l && r <= qr)
            return upper_bound(tr[node].begin(), tr[node].end(), k) - tr[node].begin();
        int m = (l + r) / 2;
        return countLE(2 * node, l, m, ql, qr, k) + countLE(2 * node + 1, m + 1, r, ql, qr, k);
    }
    // count of i in [ql,qr] (0-based, inclusive) with a[i] <= k
    int countLE(int ql, int qr, long long k) {
        if (ql > qr || n == 0) return 0;
        return countLE(1, 0, n - 1, ql, qr, k);
    }
    // count of i in [ql,qr] with a[i] in [lo,hi]  (== points in rectangle)
    int countRange(int ql, int qr, long long lo, long long hi) {
        if (lo > hi) return 0;
        return countLE(ql, qr, hi) - countLE(ql, qr, lo - 1);
    }
};

/*  USAGE
      MergeSortTree mst({5, 2, 8, 4, 7});
      mst.countLE(0, 4, 5);         // values <= 5 among all           -> 3  (5,2,4)
      mst.countRange(1, 3, 3, 8);   // i in [1,3] with a[i] in [3,8]   -> 2  (8,4)
    Points view: a[i] = y at x = i;  countRange(x1,x2,y1,y2) = points in the rectangle.
*/

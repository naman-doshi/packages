// ============================================================================
//  PERSISTENT SEGMENT TREE for 2-D POINT QUERIES
//  (store a static set of points; count / k-th over any x,y rectangle)
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    You have N fixed 2-D points and want, in O(log N) each:
//      * countRect(x1,y1,x2,y2) : how many points lie in that rectangle
//      * kthY(x1,x2,k)          : the k-th smallest y among points with x in [x1,x2]
//    This is THE tool for offline/online 2-D dominance counting, "how many
//    points below-left of (x,y)", and the classic "k-th smallest value in a
//    subarray" (make points (index, value) and query x in [l,r]).
//
//  THE IDEA (why persistence)
//    Sort points by x. version[i] = a segment tree over (compressed) y that has
//    inserted the first i points. Persistence keeps EVERY prefix version alive
//    at only O(log N) extra nodes per insert. Points with x in [x1,x2] are a
//    contiguous band of versions, so any y-range count over that band is just
//        query(version[hi], yRange) - query(version[lo], yRange).
//    That subtraction of two versions is the whole trick.
//
//  COMPLEXITY  build O(N log N) time & memory; each query O(log N).
//
//  PITFALLS
//    * y is coordinate-compressed; queries use the real y values (handled).
//    * countRect x-band uses x<x1 (lower_bound) and x<=x2 (upper_bound); keep
//      those half-open or you double count / miss the boundary column.
//    * node 0 is the shared null (cnt 0) -- never allocate over it.
//    * Points are STATIC. For insert/delete over time, use an offline sweep +
//      Fenwick, or a merge-sort tree for static count-only queries.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;

struct PSeg2D {
    struct Node { int l = 0, r = 0, cnt = 0; };
    vector<Node> t;
    vector<int> root;              // root[i] = version after first i points (sorted by x)
    vector<long long> xs, ys;      // xs[i] = x of i-th inserted point (sorted); ys = sorted-unique y
    int m = 0;                     // number of distinct y

    int update(int prev, int lo, int hi, int pos) {
        int cur = t.size(); t.push_back(t[prev]); t[cur].cnt++;
        if (lo == hi) return cur;
        int mid = (lo + hi) / 2;
        if (pos <= mid) t[cur].l = update(t[prev].l, lo, mid, pos);
        else            t[cur].r = update(t[prev].r, mid + 1, hi, pos);
        return cur;
    }
    // count in one version with y-index in [ql, qr]
    int qrange(int node, int lo, int hi, int ql, int qr) {
        if (!node || qr < lo || hi < ql) return 0;
        if (ql <= lo && hi <= qr) return t[node].cnt;
        int mid = (lo + hi) / 2;
        return qrange(t[node].l, lo, mid, ql, qr) + qrange(t[node].r, mid + 1, hi, ql, qr);
    }

    PSeg2D(vector<pair<long long,long long>> pts) {
        t.push_back(Node());                        // node 0 = null
        sort(pts.begin(), pts.end());               // by x, then y
        for (auto& p : pts) ys.push_back(p.second);
        sort(ys.begin(), ys.end());
        ys.erase(unique(ys.begin(), ys.end()), ys.end());
        m = ys.size();
        root.push_back(0);                          // version 0 = empty
        for (auto& p : pts) {
            xs.push_back(p.first);
            int yi = lower_bound(ys.begin(), ys.end(), p.second) - ys.begin() + 1;  // 1-based
            root.push_back(update(root.back(), 1, max(1, m), yi));
        }
    }

    int loIdx(long long x1) { return lower_bound(xs.begin(), xs.end(), x1) - xs.begin(); } // #(x < x1)
    int hiIdx(long long x2) { return upper_bound(xs.begin(), xs.end(), x2) - xs.begin(); } // #(x <= x2)

    // number of points with x in [x1,x2] and y in [y1,y2]
    long long countRect(long long x1, long long y1, long long x2, long long y2) {
        if (x1 > x2 || y1 > y2 || m == 0) return 0;
        int lo = loIdx(x1), hi = hiIdx(x2);
        if (hi <= lo) return 0;
        int yl = lower_bound(ys.begin(), ys.end(), y1) - ys.begin() + 1;  // 1-based first y >= y1
        int yr = upper_bound(ys.begin(), ys.end(), y2) - ys.begin();      // #(y <= y2), 1-based last
        if (yl > yr) return 0;
        return qrange(root[hi], 1, m, yl, yr) - qrange(root[lo], 1, m, yl, yr);
    }

    // k-th smallest y (1-indexed) among points with x in [x1,x2]; LLONG_MIN if <k such points
    long long kthY(long long x1, long long x2, int k) {
        int lo = loIdx(x1), hi = hiIdx(x2);
        if (k < 1 || k > hi - lo) return LLONG_MIN;
        int a = root[lo], b = root[hi], l = 1, r = m;
        while (l < r) {
            int mid = (l + r) / 2;
            int leftCnt = t[t[b].l].cnt - t[t[a].l].cnt;
            if (k <= leftCnt) { a = t[a].l; b = t[b].l; r = mid; }
            else { k -= leftCnt; a = t[a].r; b = t[b].r; l = mid + 1; }
        }
        return ys[l - 1];
    }
};

/*  USAGE
      PSeg2D pst({{1,5},{2,2},{2,8},{4,4}});     // points (x,y)
      pst.countRect(2, 0, 4, 5);                 // points with 2<=x<=4, 0<=y<=5  -> 2  ({2,2},{4,4})
      pst.kthY(1, 4, 2);                         // 2nd smallest y among all 4    -> 4

    "k-th smallest in subarray a[l..r]" : build points (i, a[i]); answer = kthY(l, r, k).
*/

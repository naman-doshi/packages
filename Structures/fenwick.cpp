// ============================================================================
//  FENWICK TREE / BINARY INDEXED TREE (BIT)
// ----------------------------------------------------------------------------
//  WHAT IT IS
//    A tiny array-backed structure for  point update + prefix/range SUM  in
//    O(log n), using O(n) memory and a much smaller constant factor than a
//    segment tree. The magic: index i is responsible for the block of
//    (i & -i) elements ending at i, so walking i += i&-i / i -= i&-i visits
//    only O(log n) blocks to cover any prefix.
//
//  WHEN TO USE IT  (vs a segment tree)
//    USE Fenwick when your operation is INVERTIBLE (has subtraction) so a range
//    is prefix(r) - prefix(l-1): sums, counts, XOR (its own inverse). It's the
//    default for "point add, range sum" -- shorter, ~2-3x faster, less memory.
//    USE a segment tree instead when:
//      * the operation is NOT invertible (min / max / gcd) -- you can't subtract
//        a prefix, so range-min needs a segtree/sparse table, not a BIT;
//      * you need lazy range updates of a non-additive kind, or to descend the
//        tree (though BIT can do range-add/range-sum and k-th, see below).
//
//  THE FOUR STANDARD PATTERNS  (pick per problem)
//    1. point add, range sum          -> `Fenwick` below (the common case)
//    2. range add, point query        -> one BIT on the DIFFERENCE array
//    3. range add, range sum          -> `FenwickRange` below (two BITs)
//    4. k-th element / order-statistic-> `Fenwick::kth` (binary lifting on BIT),
//                                         when values are non-negative counts.
//
//  CANONICAL USES
//    * running prefix sums / range sums under updates
//    * counting inversions: sweep left->right, inv += (i) - query(rank a[i]); add
//    * "how many earlier values <= x": a frequency BIT over compressed values
//    * 2-D version (nested BITs) for point-add / rectangle-sum on a grid
//
//  PITFALLS
//    * OVERFLOW: sums grow fast -- use long long (this file already does).
//    * Values can exceed array size: COORDINATE-COMPRESS first, then BIT over
//      ranks. A frequency BIT is sized to the number of distinct values.
//    * Off-by-one: this API is 0-based externally but 1-indexed internally
//      (the +i / ++i). sumPrefix(-1) == 0 by design (empty prefix).
//    * A BIT can't do range min/max -- that need is a segment tree.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;

// ---- Pattern 1: point add, prefix/range sum -------------------------------
struct Fenwick {
    int n;
    vector<long long> bit;                 // 1-indexed internally

    Fenwick(int n = 0) : n(n), bit(n + 1, 0) {}

    // Build from a 0-based array in O(n) (faster than n add()s).
    Fenwick(const vector<long long>& a) : n((int)a.size()), bit(n + 1, 0) {
        for (int i = 1; i <= n; i++) {
            bit[i] += a[i - 1];
            int j = i + (i & -i);
            if (j <= n) bit[j] += bit[i];
        }
    }

    // a[i] += delta   (i is 0-based)
    void add(int i, long long delta) {
        for (++i; i <= n; i += i & -i) bit[i] += delta;
    }

    // sum of a[0..i]  (0-based, inclusive). sumPrefix(-1) == 0.
    long long sumPrefix(int i) {
        long long s = 0;
        for (++i; i > 0; i -= i & -i) s += bit[i];
        return s;
    }

    // sum of a[l..r]  (0-based, inclusive)
    long long sum(int l, int r) {
        if (r < l) return 0;
        return sumPrefix(r) - sumPrefix(l - 1);
    }

    // Order-statistic: smallest 0-based index p with sum of a[0..p] >= k.
    // Requires all values >= 0 (e.g. a frequency BIT). Returns n if none.
    // Great for "find the k-th present element" in O(log n).
    int kth(long long k) {
        int pos = 0;
        long long cur = 0;
        int LOG = (n > 0) ? (31 - __builtin_clz(n)) : 0;
        for (int pw = 1 << LOG; pw > 0; pw >>= 1)
            if (pos + pw <= n && cur + bit[pos + pw] < k) {
                pos += pw;
                cur += bit[pos];
            }
        return pos;                        // 0-based index of the k-th element
    }
};

// ---- Pattern 3: range add, range sum  (two BITs) --------------------------
//  Maintains sum over a[0..i] while supporting "add delta to every a[l..r]".
//  Trick: sum_{0..i} = (i+1)*B1.prefix(i) - B2.prefix(i), with a paired update.
struct FenwickRange {
    Fenwick b1, b2;
    FenwickRange(int n = 0) : b1(n), b2(n) {}

    void _add(Fenwick& b, int i, long long v) { if (i >= 0) b.add(i, v); }

    // add delta to every a[l..r]  (0-based, inclusive).
    // Pairing rule (0-based): each b1.add(pos,v) is mirrored by b2.add(pos,v*pos).
    void rangeAdd(int l, int r, long long delta) {
        _add(b1, l, delta);
        _add(b1, r + 1, -delta);
        _add(b2, l, delta * l);
        _add(b2, r + 1, -delta * (r + 1));
    }

    long long prefix(int i) {              // sum of a[0..i], 0-based inclusive
        return b1.sumPrefix(i) * (long long)(i + 1) - b2.sumPrefix(i);
    }

    long long sum(int l, int r) {          // sum of a[l..r]
        if (r < l) return 0;
        return prefix(r) - prefix(l - 1);
    }
};

/*  USAGE

    // Pattern 1: point add, range sum
    Fenwick f(n);                 // all zeros
    f.add(i, v);                  // a[i] += v      (0-based)
    long long s = f.sum(l, r);    // sum of a[l..r] (inclusive)
    // or build from an array:
    Fenwick g(vector<long long>{3, 1, 4, 1, 5});

    // Frequency BIT + order statistic:
    Fenwick freq(MAXV);
    freq.add(x, 1);               // insert value x
    int median = freq.kth((freq.sumPrefix(MAXV - 1) + 1) / 2);

    // Pattern 3: range add, range sum
    FenwickRange fr(n);
    fr.rangeAdd(l, r, delta);     // a[l..r] += delta
    long long s2 = fr.sum(l, r);

    // Counting inversions (values already compressed to 1..n ranks):
    Fenwick bit(n); long long inv = 0;
    for (int i = 0; i < n; i++) {
        inv += i - bit.sumPrefix(rank[i]); // how many earlier are > rank[i]
        bit.add(rank[i], 1);
    }
*/

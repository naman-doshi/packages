// ============================================================================
//  CUSTOMISING SEGMENT TREES, LAZY SEGMENT TREES & SPARSE TABLES
//  A guide to changing WHAT they compute -- with worked, compilable examples.
// ----------------------------------------------------------------------------
//  THE ONE IDEA: A MONOID
//    Every one of these structures computes a fold of a range with some binary
//    operation. To customise it you only ever supply two things:
//        (1) MERGE  : an ASSOCIATIVE operation  combine(a, b)
//        (2) IDENTITY: a value `e` with  combine(e, x) = combine(x, e) = x
//    That pair (identity, merge) is a "monoid". Examples:
//        sum: (0, +)   min: (+INF, min)   max: (-INF, max)   gcd: (0, gcd)
//        xor: (0, ^)   product: (1, *)    and: (all-ones, &)  or: (0, |)
//    If your operation is associative and has an identity, it drops straight in.
//
//    WHICH STRUCTURE?
//      * point update, range fold, IDEMPOTENT op (min/max/gcd/and/or), STATIC
//                                              -> Sparse Table (O(1) query)
//      * point update, range fold, ANY monoid, with updates
//                                              -> Segment Tree (Section 1)
//      * RANGE update + range fold             -> Lazy Segment Tree (Section 2)
//
//    Non-commutative merges (e.g. matrix product, "max subarray") are fine for
//    segment trees and lazy trees, but you must keep left/right order in the
//    merge -- the templates below already do.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;

// ===========================================================================
//  SECTION 1 -- PLAIN SEGMENT TREE (point set, range fold)
// ---------------------------------------------------------------------------
//  Customising = writing ONE monoid struct S with:
//      static S e();               // the identity
//      S operator+(const S&) const // the (associative) merge
//  Then SegTree<S> just works. Nothing else changes.
// ===========================================================================
template <class S>
struct SegTree {
    int n; vector<S> t;
    SegTree(int n) : n(n), t(2 * n, S::e()) {}
    SegTree(const vector<S>& a) : n(a.size()), t(2 * n, S::e()) {
        for (int i = 0; i < n; i++) t[n + i] = a[i];
        for (int i = n - 1; i > 0; i--) t[i] = t[2 * i] + t[2 * i + 1];
    }
    void set(int i, S v) { for (t[i += n] = v, i >>= 1; i; i >>= 1) t[i] = t[2 * i] + t[2 * i + 1]; }
    S get(int i) { return t[i + n]; }
    S query(int l, int r) {                 // half-open [l, r)
        S rl = S::e(), rr = S::e();
        for (l += n, r += n; l < r; l >>= 1, r >>= 1) {
            if (l & 1) rl = rl + t[l++];    // left pieces fold on the left...
            if (r & 1) rr = t[--r] + rr;    // ...right pieces on the right (keeps order)
        }
        return rl + rr;
    }
};

// ---- Example monoids: change ONLY these to change the operation ------------
struct Sum { long long v; static Sum e() { return {0}; } };
Sum operator+(const Sum& a, const Sum& b) { return {a.v + b.v}; }

struct Min { long long v; static Min e() { return {LLONG_MAX}; } };
Min operator+(const Min& a, const Min& b) { return {min(a.v, b.v)}; }

struct Max { long long v; static Max e() { return {LLONG_MIN}; } };
Max operator+(const Max& a, const Max& b) { return {max(a.v, b.v)}; }

struct Gcd { long long v; static Gcd e() { return {0}; } };   // gcd(0,x)=x, so 0 is identity
Gcd operator+(const Gcd& a, const Gcd& b) { return {gcd(a.v, b.v)}; }

// A NON-OBVIOUS monoid: maximum subarray sum (contiguous, non-empty).
// Store, for a range: best prefix / best suffix / best overall / total.
struct Sub {
    long long pre, suf, best, tot;
    static Sub e() { return {LLONG_MIN, LLONG_MIN, LLONG_MIN, 0}; } // identity
};
Sub operator+(const Sub& a, const Sub& b) {
    if (a.best == LLONG_MIN) return b;      // a is identity
    if (b.best == LLONG_MIN) return a;      // b is identity
    return { max(a.pre, a.tot + b.pre),
             max(b.suf, b.tot + a.suf),
             max({a.best, b.best, a.suf + b.pre}),
             a.tot + b.tot };
}
//  USAGE (Section 1)
//    SegTree<Sum> st(n);          st.set(i, {x});   st.query(l, r).v;   // sum [l,r)
//    SegTree<Min> st(vector<Min>{{3},{1},{4}});     st.query(l, r).v;   // min [l,r)
//    Max subarray:  leaf i = Sub{a[i],a[i],a[i],a[i]};  st.query(l,r).best

// ===========================================================================
//  SECTION 2 -- LAZY SEGMENT TREE (range update + range fold)
// ---------------------------------------------------------------------------
//  Here you fill in FOUR things (this matches the interface in lazyseg.cpp):
//      Node:  data + operator+(Node)         -- the MERGE (a monoid, as above)
//             + upd(Lazy, l, r)              -- APPLY a lazy op over cells [l,r]
//      Lazy:  data + operator+=(Lazy)        -- COMPOSE a new op onto a pending one
//      NID (node identity for +) and UID (a no-op lazy)  -- passed to init().
//  The three laws that must hold (or the tree lies):
//      * merge is associative, NID is its identity;
//      * applying a lazy distributes over merge (apply to parent == apply to both
//        children then merge) -- this is why sum uses the range length (r-l+1);
//      * composing lazies is associative and UID composes as a no-op.
//  Below are drop-in Node/Lazy pairs for LazySeg<Node, Lazy, SZ> from lazyseg.cpp.
// ===========================================================================

// ---- (A) range ADD, range SUM ---------------------------------------------
struct AddSumLazy { long long add; void operator+=(const AddSumLazy& b) { add += b.add; } };
struct AddSumNode {
    long long sum;
    AddSumNode operator+(const AddSumNode& b) const { return {sum + b.sum}; }
    void upd(const AddSumLazy& u, int l, int r) { sum += u.add * (r - l + 1); } // *length!
};
//  init({0}, {0});   leaf i = {a[i]};   tree.upd(l,r,{d});   tree.query(l,r).sum

// ---- (B) range ADD, range MIN ---------------------------------------------
struct AddMinLazy { long long add; void operator+=(const AddMinLazy& b) { add += b.add; } };
struct AddMinNode {
    long long mn;
    AddMinNode operator+(const AddMinNode& b) const { return {min(mn, b.mn)}; }
    void upd(const AddMinLazy& u, int, int) { mn += u.add; }   // min shifts by add; no length
};
//  init({(long long)4e18}, {0});   leaf {a[i]};   upd(l,r,{d});   query(l,r).mn

// ---- (C) range ASSIGN, range MAX ------------------------------------------
struct AssignLazy {
    long long val; bool set;                 // set=false => "no pending assign"
    void operator+=(const AssignLazy& b) { if (b.set) { val = b.val; set = true; } }
};
struct AssignMaxNode {
    long long mx;
    AssignMaxNode operator+(const AssignMaxNode& b) const { return {max(mx, b.mx)}; }
    void upd(const AssignLazy& u, int, int) { if (u.set) mx = u.val; }
};
//  init({(long long)-4e18}, {0,false});  leaf {a[i]};  upd(l,r,{v,true});  query.mx
//  NOTE UID here is {0,false} = "no assign pending" (the no-op), NOT {0,true}.

// ---- (D) range ADD *and* ASSIGN, range MAX & SUM together -----------------
//  That's exactly the struct pair already in lazyseg.cpp (Lazy has an `inc`
//  flag: inc=true adds, inc=false assigns; assign wins when composing). Read it
//  there as the canonical "two update kinds at once" example.

// ===========================================================================
//  SECTION 3 -- SPARSE TABLE (static, idempotent range fold)
// ---------------------------------------------------------------------------
//  A sparse table needs LESS than a monoid: just an associative + IDEMPOTENT
//  merge (op(x,x) = x), and NO identity (it never folds an empty range). That
//  idempotence is what lets the two overlapping half-blocks double-count safely.
//  In sparse.cpp you change ONLY the merge -- replace both `min(...)` calls:
//      min : jmp[j][i] = min(jmp[j-1][i], jmp[j-1][i+(1<<(j-1))]);
//            return     min(jmp[d][l],   jmp[d][r-(1<<d)+1]);
//      max : swap min -> max in both lines.
//      gcd : gcd(a, b)           (idempotent: gcd(x,x)=x)
//      and : (a & b)             or : (a | b)
//  DO NOT use a sparse table for sum/product/xor -- those are not idempotent, so
//  overlapping blocks corrupt the answer. Use a Fenwick or segment tree instead.
// ===========================================================================

// ===========================================================================
//  DESIGNING YOUR OWN MONOID (when the answer isn't a plain min/max/sum)
//    1. Ask: "to merge a LEFT range with a RIGHT range, what do I need to store
//       about each so the combined answer is computable in O(1)?" Store exactly
//       that. (Max-subarray needs pre/suf/best/tot -- best alone can't merge.)
//    2. Write merge(left, right). Check it's associative. Non-commutative is OK
//       (keep order), but the merge must not depend on absolute positions.
//    3. Pick the identity: the value that merges as a no-op (often "empty range").
//    4. For LAZY: also define how a bulk op transforms the stored data (using the
//       range length if the data is size-dependent, like sum), and how two bulk
//       ops compose. If a bulk op can't be expressed on your aggregate, lazy
//       won't work -- rethink the state.
//    Other useful monoids: (count of the min in a range), (max & 2nd max),
//    (assignment count / hash of a range), (matrix product for linear recurrences).
// ===========================================================================

// ============================================================================
//  LAZY SEGMENT TREE  (range update + range query, O(log n) each)
// ----------------------------------------------------------------------------
//  WHAT IT DOES
//    Maintains an array and supports, all in O(log SZ):
//      * apply an operation to EVERY element of a range [l, r]        (upd)
//      * aggregate a range [l, r]   (this file: max AND sum together) (query)
//    "Lazy" = a pending range-update is parked at a node and only pushed down
//    to its children when you actually descend past it. That deferral is the
//    whole trick that makes RANGE updates O(log n) instead of O(n).
//
//  MENTAL MODEL  (read this if it feels like a black box)
//    Two parallel arrays over the same 2*SZ node space:
//      seg[i]  = aggregate for node i's range, ASSUMING lazy[i] is applied to
//                it but NOT yet to its children.
//      lazy[i] = an update owed to node i's whole subtree, not yet pushed down.
//    push(i,l,r): "cash in" lazy[i] -- apply it to seg[i], copy it into both
//                 children's lazy, then clear lazy[i].  Called on the way DOWN.
//    pull(i):     recompute seg[i] from its (already-correct) children. Way UP.
//    Every upd/query push()es a node on entry (so what it reads is current) and
//    pull()s internal nodes on exit (so parents stay correct). That's the loop.
//
//  THE THREE THINGS YOU CUSTOMISE  (structs Node, Lazy below)
//    Node (T): the per-range aggregate you store.
//        operator+(Node)   -> MERGE two child aggregates (must be associative).
//        upd(Lazy,l,r)     -> APPLY a lazy op to this node covering cells [l,r]
//                             (needs length r-l+1 for sum-like data).
//    Lazy (U): a pending update.
//        operator+=(Lazy b)-> COMPOSE b on top of the pending op. Order matters:
//                             here an "assign" (inc=false) wipes a pending "add".
//    Identities:
//        NID = identity of operator+   (max -> -INF, sum -> 0, min -> +INF).
//        UID = a NO-OP lazy            (here {0, true} = "add 0").
//
//  DEEP DIVE: what NID / UID are, and exactly what tree.init(NID, UID) does
//    tree.init(NID, UID) does two simple things:
//      1. remembers your two "identity" values (stores them in the struct), and
//      2. fills the whole tree with them -- every seg[] cell = NID, every lazy[]
//         cell = UID. So NID and UID are just the "empty / nothing-here-yet"
//         values the tree starts from and falls back to.
//
//    NID  ("Node IDentity") = the aggregate of an EMPTY range. Requirement:
//         merging it with anything must change nothing, i.e. NID + X == X for
//         your operator+. Choose it per aggregate:
//            max  -> mx = -INF   because max(-INF, x) = x
//            min  -> mn = +INF   because min(+INF, x) = x
//            sum  -> sum = 0     because 0 + x = x
//         This tree merges max AND sum at once, so NID = {mx:-INF, sum:0}.
//         WHY IT MATTERS: query() returns NID for the parts of the tree outside
//         your [l,r]. If NID were wrong (say mx:0), then querying the max of an
//         all-negative array would wrongly return 0. NID is the safe "this piece
//         contributes nothing" value.
//
//    UID  ("Update IDentity") = the "do nothing" lazy. Requirement: applying it
//         to a node changes nothing (X.upd(UID) == X) and composing it with any
//         update leaves that update unchanged. Here a Lazy is {v, inc}; with
//         inc=true it means "add v", so {0, true} = "add 0" = a true no-op.
//         NOTE you could NOT use {0, false} as the identity: inc=false means
//         "ASSIGN 0", which is a real change, not nothing. UID is what a node's
//         lazy holds when nothing is pending, and it's the value push() resets
//         lazy[i] to after it has cashed a pending update in.
//
//    So the standard init for THIS (max+sum) tree reads:
//         tree.init({-INF, 0}, {0, true});
//                     ~~~~~~~   ~~~~~~~~~
//                     NID       UID
//         {-INF, 0} : an empty range has max -INF and sum 0.
//         {0, true} : the neutral "add 0" lazy.
//    If you re-purpose the tree (e.g. min instead of max), change NID's field to
//    +INF to match your merge. UID stays "add 0" as long as your lazy has an
//    additive mode; if it doesn't, invent a sentinel value that upd() treats as
//    "skip". (For concrete re-purposing examples see customSegtrees.cpp.)
//
//  API  (indices 0-based over [0, SZ); you use [0, n))
//    LazySeg<Node, Lazy, 1<<20> tree;          // SZ MUST be a power of two >= n
//    tree.init({-INF, 0}, {0, true});          // NID, UID
//    for (i in 0..n-1) tree[i] = {a[i], a[i]}; // seed leaves (mx=val, sum=val)
//    tree.build();                             // build internal nodes -- REQUIRED
//    tree.upd(l, r, {v, true});                // add v to every cell in [l,r]
//    tree.upd(l, r, {v, false});               // assign v to every cell in [l,r]
//    tree.query(l, r).mx;                      // max on [l,r]
//    tree.query(l, r).sum;                     // sum on [l,r]
//
//  ADAPTING TO OTHER PROBLEMS  (swap Node/Lazy, keep the skeleton)
//    range add    + range sum : Node{sum};      upd: sum += v*(r-l+1)
//    range assign + range min : Node{mn}; NID{+INF}; upd assign: mn = v
//    range add    + range min : Node{mn};       upd: mn += v
//    range add/assign + max&sum: exactly this file (both tags via `inc`)
//    Rule of thumb: a Lazy works only if it (a) composes with itself
//    (operator+=) and (b) can be applied to an aggregate given the range length.
//
//  COMPLEXITY   build O(SZ);  upd/query O(log SZ).
//    Memory ~ 2*SZ*(sizeof(Node)+sizeof(Lazy)); SZ=1<<20 => ~64MB, mind limits.
//
//  COMMON PITFALLS
//    * SZ must be a power of two and >= n, or leaf indexing (i + SZ) breaks.
//    * You MUST call build() after seeding leaves; skipping it leaves internal
//      aggregates wrong.
//    * Identities: -INF for max, +INF for min, 0 for sum -- a wrong identity
//      silently corrupts negative arrays and empty ranges.
//    * Overflow: sum of many big values needs long long (this file #defines int
//      as long long -- keep that if you adapt it).
//    * upd/query ranges are INCLUSIVE [l,r]; disjoint ranges return NID.
//    * Tag composition lives in Lazy::operator+= -- if you add tag types, get
//      their composition order right or the lazies will fight.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long
#define ll long long
 
const int INF = 1e18;
 
#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define trav(a, x) for(auto& a : x)
#define all(x) x.begin(), x.end()
#define print(arr) for (auto i : arr) cout << i << " "; cout << endl;
typedef vector<vector<int>> vvi;
typedef vector<pair<int, int>> vpi;
typedef vector<int> vi;
typedef pair<int, int> pii;
 
void setFile(string name) {
  freopen((name + "in.txt").c_str(), "r", stdin);
  freopen((name + "out.txt").c_str(), "w", stdout);
}
 
int MOD = 1e9 + 7;

struct Lazy {
    int v;
    bool inc;
    void operator+=(const Lazy &b) {
        if (b.inc) v += b.v;
        else v = b.v, inc = false;
    }
};
 
struct Node {
    int mx, sum;
    Node operator+(const Node &b) {
        return {max(mx, b.mx), sum + b.sum};
    }
    void upd(const Lazy &u, int l, int r) {
        if (u.inc) mx += u.v, sum += u.v * (r - l + 1);
        else mx = u.v, sum = u.v * (r - l + 1);
    }
};
 
template<class T, class U, int SZ> struct LazySeg {
    T NID;
    U UID;
    vector<T> seg;
    vector<U> lazy;
    void init(T _NID, U _UID) {
        NID = _NID;
        UID = _UID;
        seg.resize(2 * SZ, NID);
        lazy.resize(2 * SZ, UID);
    }
    void pull(int i) {
        seg[i] = seg[2 * i] + seg[2 * i + 1];
    }
    void push(int i, int l, int r) {
        seg[i].upd(lazy[i], l, r);
        if (l != r)  rep(j, 0, 2) lazy[2 * i + j] += lazy[i];
        lazy[i] = UID;
    }
    void build() {
        for (int i = SZ - 1; i > 0; i--) pull(i);
    }
    void upd(int lo, int hi, U val, int i = 1, int l = 0, int r = SZ - 1) {
        push(i, l, r);
        if (r < lo || l > hi) return;
        if (lo <= l && r <= hi) {
            lazy[i] += val;
            push(i, l, r);
            return;
        }
        int m = (l + r) / 2;
        upd(lo, hi, val, 2 * i, l, m);
        upd(lo, hi, val, 2 * i + 1, m + 1, r);
        pull(i);
    }
    T query(int lo = 0, int hi = SZ - 1, int i = 1, int l = 0, int r = SZ - 1) {
        push(i, l, r);
        if (r < lo || l > hi) return NID;
        if (lo <= l && r <= hi) return seg[i];
        int m = (l + r) / 2;
        return query(lo, hi, 2 * i, l, m) + query(lo, hi, 2 * i + 1, m + 1, r);
    }
    T& operator[](int i) {
        return seg[i + SZ];
    }
};

// LazySeg<Node, Lazy, 1 << 20> tree;
// tree.init({-INF, 0}, {0, true});
// F0R (i, n) tree[i] = {0, 0};
// tree.build();

// tree.query(l, r).mx; // maximum of [l, r]
// tree.upd(l, r, {1, true}); // increments [l, r] by 1
// tree.upd(l, r, {0, false}); // sets [l, r] to 0
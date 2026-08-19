// ============================================================================
//  SPARSE TABLE (RMQ)  -- O(1) range query on a STATIC array
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Precomputes the answer for every power-of-two length, so any range query
//    is answered in O(1) by overlapping two blocks. Use it when the array never
//    changes and the operation is IDEMPOTENT (safe to double-count on overlap):
//    min, max, gcd, bitwise and/or. NOT for sum/product/xor (overlap corrupts
//    them -- use a Fenwick or segment tree for those).
//
//  API  (0-based, INCLUSIVE range)
//    RMQ<int> r;
//    r.init(v);                 // v is a vector<int>, build is O(n log n)
//    int m = r.query(l, r);     // min over v[l..r], O(1)
//
//  CUSTOMISING IT  (see customSegtrees.cpp for the full guide)
//    A sparse table needs only an idempotent, associative merge -- no identity.
//    To change the operation, replace BOTH `min(...)` calls (in init and query):
//      max:  jmp[j][i] = max(...);      return max(jmp[d][l], jmp[d][r-(1<<d)+1]);
//      gcd:  __gcd(a, b)                (works because gcd is idempotent)
//      or/and: (a | b) / (a & b)
//    That's the entire change -- the index math is operation-agnostic.
//
//  COMPLEXITY  build O(n log n) time & memory;  query O(1).
//
//  PITFALLS
//    * Static only -- no updates. If values change, use a segment tree.
//    * query(l, r) is inclusive and needs l <= r; empty ranges aren't handled.
//    * v must be non-empty (init reads v.size()).
//    * Idempotent ops only -- never use this for sum.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long
#define ll long long

const int INF = 9 * 1e18;
 
#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define trav(a, x) for(auto& a : x)
#define all(x) x.begin(), x.end()
#define print(arr) for (auto i : arr) cout << i << " "; cout << endl;
typedef vector<vector<int>> vvi;
typedef vector<pair<int, int>> vpi;
typedef vector<int> vi;
typedef pair<int, int> pii; 
 
template<class T> struct RMQ {
  vector<vector<T>> jmp;
  void init(const vector<T>& v) { // assert(v.size() <= sz);
    jmp.resize(32 - __builtin_clz(v.size()), vector<T>(v.size()));
    copy(all(v), begin(jmp[0]));
    for (int j = 1; 1 << j <= v.size(); ++j) {
      rep (i, 0, v.size() - (1 << j) + 1) jmp[j][i] = min(jmp[j - 1][i], jmp[j - 1][i + (1 << (j - 1))]);
    }
  }
  T query(int l, int r) {
    int d = 31 - __builtin_clz(r - l + 1);
    return min(jmp[d][l], jmp[d][r - (1 << d) + 1]);
  }
};
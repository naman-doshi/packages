// ============================================================================
//  XOR BASIS (linear basis over GF(2))  -- the span of a set under XOR
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Treat each number as a bit-vector over GF(2). The basis is a maximal
//    independent subset: it has at most 60 elements no matter how many numbers
//    you insert, and it spans EXACTLY the same set of XOR values as the
//    original array. That instantly answers:
//      - maximum / minimum XOR of any subset
//      - is value v achievable as a XOR of some subset?
//      - how many distinct XOR values are achievable  (2^size)
//      - the k-th smallest achievable value
//      - how many subsets give a particular XOR  (2^(n - size) each, if achievable)
//
//  API
//    XorBasis b;
//    b.insert(x)        // true if x grew the basis (i.e. x was independent)
//    b.size             // number of independent vectors
//    b.canMake(v)       // is v in the span?
//    b.maxXor()         // largest achievable value  (maxXor(start) to seed it)
//    b.minXor()         // smallest NON-ZERO achievable value (0 is always in)
//    b.kth(k)           // k-th smallest achievable value, 1-indexed, 0 counts
//    b.countDistinct()  // 2^size
//    b.merge(other)     // union of two spans
//
//  COMPLEXITY  every operation is O(B) = O(60). Building over n numbers is
//              O(n * 60), which is effectively linear.
//
//  PITFALLS
//    * B = 60 covers values < 2^60. Raise it (and only it) for bigger inputs;
//      for values up to 1e18 you need B = 60, for unsigned 64-bit use 64 and
//      switch the type to unsigned long long.
//    * kth() and minXor() need the REDUCED basis -- they call reduce()
//      themselves, which mutates internal state but not the span. Don't cache
//      basis[] across a reduce().
//    * The empty basis spans only {0}: maxXor() = 0, canMake(0) = true.
//    * For "max XOR of a SUBARRAY" you need prefix XORs + a persistent /
//      offline structure, not this alone. For "max XOR pair" a binary trie is
//      simpler and faster.
//
//  RECIPE -- max XOR path in a graph (classic):
//    Run any DFS, record dist[v] = XOR of edges from the root. Every cycle's
//    XOR goes into the basis. Then the answer for u->v is
//    maxXor(dist[u] ^ dist[v]) -- seed the greedy with the tree path value.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define all(x) x.begin(), x.end()
typedef vector<int> vi;

struct XorBasis {
  static const int B = 60;            // bits; values must be < 2^B
  vi basis;                           // basis[i] has highest set bit i, or 0
  int size = 0;
  bool reduced = false;

  XorBasis() : basis(B, 0) {}

  // Insert x. Returns true if it was independent of everything already in.
  bool insert(int x) {
    for (int i = B - 1; i >= 0; i--) {
      if (!((x >> i) & 1)) continue;
      if (!basis[i]) { basis[i] = x; size++; reduced = false; return true; }
      x ^= basis[i];
    }
    return false;                     // x reduced to 0 -> already in the span
  }

  bool canMake(int x) const {
    for (int i = B - 1; i >= 0; i--) {
      if (!((x >> i) & 1)) continue;
      if (!basis[i]) return false;
      x ^= basis[i];
    }
    return true;
  }

  // Largest achievable XOR. Pass a starting value to compute
  // max over subsets S of (start ^ xor(S)).
  int maxXor(int start = 0) const {
    int r = start;
    for (int i = B - 1; i >= 0; i--)
      if (basis[i] && ((r >> i) & 1) == 0) r ^= basis[i];
    return r;
  }
  // Smallest achievable value strictly greater than 0 (0 itself is always in
  // the span via the empty subset). Returns 0 if the basis is empty.
  int minXor() {
    reduce();
    for (int i = 0; i < B; i++) if (basis[i]) return basis[i];
    return 0;
  }
  // Reduce to canonical form: each pivot bit appears in exactly one vector.
  void reduce() {
    if (reduced) return;
    for (int i = 0; i < B; i++)
      if (basis[i])
        for (int j = i + 1; j < B; j++)
          if (basis[j] && ((basis[j] >> i) & 1)) basis[j] ^= basis[i];
    reduced = true;
  }
  // k-th smallest achievable value, 1-indexed, counting 0 as the 1st.
  // Returns -1 if k exceeds 2^size.
  int kth(int k) {
    reduce();
    vi rows;
    rep(i, 0, B) if (basis[i]) rows.push_back(basis[i]);
    int cnt = rows.size();
    if (cnt < 62 && k > (1LL << cnt)) return -1;
    k--;                               // 0-indexed: bit i picks rows[i]
    int res = 0;
    rep(i, 0, cnt) if ((k >> i) & 1) res ^= rows[i];
    return res;
  }
  // 2^size distinct XOR values are reachable (may overflow past size 62).
  int countDistinct() const { return 1LL << size; }

  void merge(const XorBasis &o) { for (int v : o.basis) if (v) insert(v); }
};

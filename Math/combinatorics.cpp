// ============================================================================
//  COMBINATORICS  -- nCr, Catalan, stars & bars, derangements, Lucas
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Precompute factorials and inverse factorials once, then every binomial is
//    O(1). This is the workhorse for counting problems mod a prime.
//
//  API
//    Comb c(N);            // N = the LARGEST factorial you will ever ask for.
//                          // For catalan(k) / central binomials you need 2*k!
//    c.C(n, k)  c.P(n, k)  c.f[n]  c.fi[n]  c.inv[n]
//    c.catalan(k)  c.derange(k)  c.stars(k, boxes)  c.starsPos(k, boxes)
//    c.multinomial({a, b, c})    c.C2(n)   // n choose 2, no precompute needed
//
//    Lucas L(p);  L.C(n, k);     // C(n,k) mod SMALL prime p, n up to 1e18
//    Cbig(n, k, p);              // n huge, k small (<= 1e6), p prime
//    pascal(n, m);               // C table mod ANY m (composite included)
//
//  COMPLEXITY  Comb build O(N);  every query O(1).  Lucas build O(p), query
//              O(log_p n).  pascal O(n^2) time and memory.
//
//  PITFALLS
//    * MOD MUST BE PRIME and > N. Composite modulus => use pascal(), or CRT
//      over the prime power factorisation (Andrew Granville / generalised Lucas).
//    * C(n,k) returns 0 for k < 0 or k > n, which is usually what you want in
//      inclusion-exclusion -- but it will also silently hide an out-of-range
//      index bug, so size N generously.
//    * catalan(k) needs the table built to 2k, NOT k.
//    * Everything here is `int` = long long, and products are reduced after
//      every multiply, so no overflow as long as MOD < 3e9.
//
//  IDENTITIES worth having in your head
//    C(n,k) = C(n-1,k-1) + C(n-1,k)            Pascal
//    sum_k C(n,k) = 2^n        sum_k (-1)^k C(n,k) = [n == 0]
//    sum_k C(n,k)^2 = C(2n,n)                  Vandermonde special case
//    sum_{i<=n} C(i,k) = C(n+1,k+1)            hockey stick
//    sum_i C(a,i) C(b,k-i) = C(a+b,k)          Vandermonde
//    C(n,k) * C(k,j) = C(n,j) * C(n-j,k-j)     subset-of-a-subset swap
//    # of lattice paths (0,0)->(a,b) with unit right/up steps = C(a+b, a)
//    # paths that never cross the diagonal = Catalan / ballot numbers
//    # labelled trees on n nodes = n^(n-2) (Cayley); with fixed degrees d_i it
//      is the multinomial (n-2)! / prod (d_i - 1)!  (Prufer sequences)
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define all(x) x.begin(), x.end()
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int, int> pii;

const int MOD = 1000000007;

int powmod(int a, int e, int m) {
  int r = 1 % m; a %= m; if (a < 0) a += m;
  while (e > 0) { if (e & 1) r = r * a % m; a = a * a % m; e >>= 1; }
  return r;
}

struct Comb {
  int n;
  vi f;     // f[i]   = i!            mod MOD
  vi fi;    // fi[i]  = (i!)^-1       mod MOD
  vi inv;   // inv[i] = i^-1          mod MOD
  vi der;   // der[i] = derangements of i items

  Comb(int N) : n(max(N, 2LL)) {
    f.assign(n + 1, 1); fi.assign(n + 1, 1);
    inv.assign(n + 1, 1); der.assign(n + 1, 0);
    // inv[i] from p = (p/i)*i + p%i  =>  i^-1 = -(p/i) * (p%i)^-1
    rep(i, 2, n + 1) inv[i] = (MOD - (MOD / i) * inv[MOD % i] % MOD) % MOD;
    rep(i, 1, n + 1) f[i] = f[i - 1] * i % MOD;
    rep(i, 1, n + 1) fi[i] = fi[i - 1] * inv[i] % MOD;
    der[0] = 1; der[1] = 0;
    rep(i, 2, n + 1) der[i] = (i - 1) % MOD * ((der[i - 1] + der[i - 2]) % MOD) % MOD;
  }

  // ---- the two basics ----
  int C(int a, int b) const {                    // choose: order doesn't matter
    if (b < 0 || a < 0 || b > a) return 0;
    return f[a] * fi[b] % MOD * fi[a - b] % MOD;
  }
  int P(int a, int b) const {                    // arrange: order matters
    if (b < 0 || a < 0 || b > a) return 0;
    return f[a] * fi[a - b] % MOD;
  }
  static int C2(int a) { return a < 2 ? 0 : a % 2 ? a * ((a - 1) / 2) : (a / 2) * (a - 1); }

  // ---- named families ----
  // Catalan: balanced bracket sequences of length 2k, binary trees with k
  // internal nodes, triangulations of a (k+2)-gon, monotone lattice paths that
  // stay weakly below the diagonal, stack-sortable permutations. NEEDS n >= 2k.
  int catalan(int k) const { return ((C(2 * k, k) - C(2 * k, k + 1)) % MOD + MOD) % MOD; }

  // Permutations of k items with NO fixed point. der[k] = round(k!/e).
  int derange(int k) const { return der[k]; }

  // Ways to split k identical items into `boxes` ordered NON-NEGATIVE parts.
  int stars(int k, int boxes) const { return C(k + boxes - 1, boxes - 1); }
  // ... into POSITIVE parts (each box gets at least one).
  int starsPos(int k, int boxes) const { return C(k - 1, boxes - 1); }
  // ... with an upper bound of cap on every box: inclusion-exclusion over how
  // many boxes overflow.   sum_j (-1)^j C(boxes,j) C(k - j*(cap+1) + boxes-1, boxes-1)
  int starsCapped(int k, int boxes, int cap) const {
    int res = 0;
    rep(j, 0, boxes + 1) {
      int rem = k - j * (cap + 1);
      if (rem < 0) break;
      int t = C(boxes, j) * C(rem + boxes - 1, boxes - 1) % MOD;
      res = (j % 2 ? res - t : res + t) % MOD;
    }
    return (res + MOD) % MOD;
  }

  // Ways to arrange a multiset with the given counts: (sum c)! / prod c_i!
  int multinomial(const vi &c) const {
    int s = 0, r = 1;
    for (int x : c) { s += x; r = r * fi[x] % MOD; }
    return r * f[s] % MOD;
  }

  // # of surjections from an a-set onto a b-set (all boxes non-empty, labelled)
  int surjections(int a, int b) const {
    int res = 0;
    rep(j, 0, b + 1) {
      int t = C(b, j) * powmod(b - j, a, MOD) % MOD;
      res = (j % 2 ? res - t : res + t) % MOD;
    }
    return (res + MOD) % MOD;
  }
};

// ---- Lucas: C(n, k) mod a SMALL prime p ------------------------------------
// C(n,k) = prod C(n_i, k_i) mod p over the base-p digits of n and k.
// Build cost O(p), so p up to about 1e7. n and k may be up to 1e18.
struct Lucas {
  int p;
  vi f, fi;
  Lucas(int p) : p(p), f(p), fi(p) {
    f[0] = 1;
    rep(i, 1, p) f[i] = f[i - 1] * i % p;
    fi[p - 1] = powmod(f[p - 1], p - 2, p);
    for (int i = p - 1; i > 0; i--) fi[i - 1] = fi[i] * i % p;
  }
  int small(int a, int b) const {
    if (b < 0 || b > a) return 0;
    return f[a] * fi[b] % p * fi[a - b] % p;
  }
  int C(int a, int b) const {
    int r = 1;
    while ((a || b) && r) { r = r * small(a % p, b % p) % p; a /= p; b /= p; }
    return r;
  }
};

// C(n, k) mod prime p when n is ASTRONOMICAL but k is small (k < p, k <= ~1e6).
// Multiplies out the k-term falling factorial directly.
int Cbig(int n, int k, int p) {
  if (k < 0 || n < 0 || k > n) return 0;
  int num = 1, den = 1;
  rep(i, 0, k) { num = num * ((n - i) % p) % p; den = den * ((i + 1) % p) % p; }
  return num * powmod(den, p - 2, p) % p;
}

// Pascal's triangle mod ANY modulus (works when m is composite and inverse
// factorials are unavailable). O(n^2) -- so n up to a few thousand.
vvi pascal(int n, int m) {
  vvi c(n + 1, vi(n + 1, 0));
  rep(i, 0, n + 1) {
    c[i][0] = 1 % m;
    rep(j, 1, i + 1) c[i][j] = (c[i - 1][j - 1] + c[i - 1][j]) % m;
  }
  return c;
}

// ---- Stirling / Bell -------------------------------------------------------
// S2[n][k] = ways to partition an n-set into k NON-EMPTY UNLABELLED blocks.
//   S2[n][k] = k*S2[n-1][k] + S2[n-1][k-1]
// S1[n][k] = permutations of n with exactly k cycles (unsigned Stirling 1st).
//   S1[n][k] = (n-1)*S1[n-1][k] + S1[n-1][k-1]
// Bell[n] = sum_k S2[n][k] = total partitions of an n-set.
vvi stirling2(int n) {
  vvi s(n + 1, vi(n + 1, 0));
  s[0][0] = 1;
  rep(i, 1, n + 1) rep(j, 1, i + 1)
    s[i][j] = (j * s[i - 1][j] + s[i - 1][j - 1]) % MOD;
  return s;
}
vvi stirling1(int n) {
  vvi s(n + 1, vi(n + 1, 0));
  s[0][0] = 1;
  rep(i, 1, n + 1) rep(j, 1, i + 1)
    s[i][j] = ((i - 1) * s[i - 1][j] + s[i - 1][j - 1]) % MOD;
  return s;
}
// Closed form, O(k log n), for a single S2(n, k) with huge n:
//   S2(n,k) = (1/k!) * sum_{j=0}^{k} (-1)^j C(k,j) (k-j)^n
int stirling2Single(int n, int k, const Comb &c) {
  int res = 0;
  rep(j, 0, k + 1) {
    int t = c.C(k, j) * powmod(k - j, n, MOD) % MOD;
    res = (j % 2 ? res - t : res + t) % MOD;
  }
  return (res + MOD) % MOD * c.fi[k] % MOD;
}

// p(n) = integer partitions of n (unordered summands), O(n sqrt n) via Euler's
// pentagonal number theorem:  p(n) = sum_k (-1)^(k+1) [p(n - k(3k-1)/2) + p(n - k(3k+1)/2)]
vi partitions(int n) {
  vi p(n + 1, 0);
  p[0] = 1;
  rep(i, 1, n + 1) {
    for (int k = 1;; k++) {
      int g1 = k * (3 * k - 1) / 2, g2 = k * (3 * k + 1) / 2;
      if (g1 > i && g2 > i) break;
      int sgn = (k % 2) ? 1 : -1;
      if (g1 <= i) p[i] = (p[i] + sgn * p[i - g1]) % MOD;
      if (g2 <= i) p[i] = (p[i] + sgn * p[i - g2]) % MOD;
    }
    p[i] = (p[i] % MOD + MOD) % MOD;
  }
  return p;
}

// ---- counting with symmetry ------------------------------------------------
//  BURNSIDE / POLYA: the number of distinct objects under a group G of
//  symmetries is (1/|G|) * sum over g in G of (# objects fixed by g).
//  Necklaces with k colours and n beads under rotation:
//      (1/n) * sum_{d | n} phi(n/d) * k^d
//  ... under rotation AND reflection (bracelets), add the reflections:
//      n odd : (1/2n)[ sum_{d|n} phi(n/d) k^d + n * k^((n+1)/2) ]
//      n even: (1/2n)[ sum_{d|n} phi(n/d) k^d + (n/2)(k^(n/2) + k^(n/2+1)) ]
//
//  INCLUSION-EXCLUSION: |A1 u ... u An| = sum |Ai| - sum |Ai n Aj| + ...
//  With n <= 20 properties, iterate over subsets and weight by (-1)^popcount.

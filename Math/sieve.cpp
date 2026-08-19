// ============================================================================
//  SIEVE  -- primes, smallest-prime-factor, phi, mobius, divisors up to N
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    One linear sieve gives you the smallest prime factor of every number up to
//    N, and spf is the key to everything else: factorising any x <= N in
//    O(log x), listing its divisors, and (optionally, built in the same pass)
//    Euler's phi and the Mobius function.
//    For a SINGLE number bigger than N use factor.cpp (Miller-Rabin + rho).
//
//  API
//    Sieve S(N);                    // spf + primes only
//    Sieve S(N, true, true);        // also build phi[] and mu[]
//    S.isPrime(x)                   // O(1)
//    S.primes                       // vector of all primes <= N
//    S.factor(x)   -> {{p, e}, ...}     S.factorList(x) -> flat with multiplicity
//    S.divisors(x) -> sorted vector      S.numDivisors(x)   S.sumDivisors(x)
//    S.phi[x]  S.mu[x]              // only if you asked for them
//    S.phiOf(x)  S.muOf(x)          // computed from spf, no table needed
//
//    segmentedSieve(L, R)           // primality of every x in [L, R], R-L <= 1e7
//    primesInRange(L, R)            // the actual list
//    countPrimes(n)                 // pi(n) for n up to ~1e11
//    DivTable T(N);  for (int d : T[i])  // ALL divisors of EVERY i <= N
//    divisorTable(N)                     // same, as vector<vector<int32_t>>
//    divisorCountTable(N) / divisorSumTable(N)   // d(i), sigma(i) for all i
//    divisorZeta / divisorMobius / multipleZeta / multipleMobius
//    divisorBlockSum(n)                  // sum_{i<=n} floor(n/i), O(sqrt n)
//    CoprimeCounter C(N);  C(m)          // pairs in [1,m]^2 with gcd 1
//                          C.gcdExactly(m, k)
//    FixedCoprime F(primesOf(x));        // x fixed; reuse across queries
//      F.upTo(N)  F.inRange(A, B)        // # y coprime to x
//      F.sumUpTo(N)  F.sumInRange(A, B)  // their sum
//      F.nth(k)                          // k-th y coprime to x
//    countCoprimeTo(N, primesOf(x))      // one-shot wrappers of the above
//    countCoprimeToRange(A, B, primesOf(x))
//    countGcdWith(N, g, primesOf(x/g))   // y in [1,N] with gcd(y, x) == g
//    countGcdWithRange(A, B, g, primesOf(x/g))
//
//  COMPLEXITY  build O(N) time. Memory is 4 bytes/number for spf, plus 4 for
//              phi and 1 for mu. N = 1e7 costs 40MB / 80MB / 90MB -- watch the
//              memory limit, and prefer a plain bitset sieve if you only need
//              "is it prime".
//
//  PITFALLS
//    * The big arrays are DELIBERATELY int32_t/int8_t, not `int`. With
//      `#define int long long` in your template a vector<int> of size 1e7 is
//      80MB and will MLE. Don't "tidy" those types away.
//    * phi and mu are only filled when you pass the flags; otherwise they are
//      empty vectors and indexing them is UB.
//    * factor(x) needs x <= N and x >= 1. factor(1) is empty (correct: 1 has
//      no prime factors, and its only divisor is itself).
//    * segmentedSieve indexes by (x - L); L may be 1 (handled) but not 0.
//
//  FACTS
//    phi(n) = n * prod_{p | n} (1 - 1/p),  sum_{d | n} phi(d) = n
//    mu(n) = 0 if n has a squared factor, else (-1)^(#distinct primes)
//    Mobius inversion: g(n) = sum_{d|n} f(d)  <=>  f(n) = sum_{d|n} mu(d) g(n/d)
//    d(n) = prod (e_i + 1),  sigma(n) = prod (p^(e_i+1) - 1)/(p - 1)
//    n <= 1e18 has at most 103680 divisors; n <= 1e9 at most 1344.
//    pi(n) ~ n / ln n:  pi(1e6) = 78498, pi(1e7) = 664579, pi(1e9) = 50847534.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define all(x) x.begin(), x.end()
typedef vector<int> vi;
typedef pair<int, int> pii;

struct Sieve {
  int n;
  vector<int32_t> spf;      // spf[x] = smallest prime factor of x (0 for 0,1)
  vector<int32_t> primes;   // every prime <= n, ascending
  vector<int32_t> phi;      // only if wantPhi
  vector<int8_t> mu;        // only if wantMu

  Sieve(int n, bool wantPhi = false, bool wantMu = false) : n(n), spf(n + 1, 0) {
    if (wantPhi) { phi.assign(n + 1, 0); if (n >= 1) phi[1] = 1; }
    if (wantMu)  { mu.assign(n + 1, 0);  if (n >= 1) mu[1] = 1;  }
    for (int i = 2; i <= n; i++) {
      if (!spf[i]) {                            // i is prime
        spf[i] = i; primes.push_back(i);
        if (wantPhi) phi[i] = i - 1;
        if (wantMu)  mu[i] = -1;
      }
      // Every composite is struck exactly once, by its smallest prime factor.
      for (int32_t p : primes) {
        if (p > spf[i] || i * p > n) break;
        spf[i * p] = p;
        if (p == spf[i]) {                      // p^2 divides i*p
          if (wantPhi) phi[i * p] = phi[i] * p;
          if (wantMu)  mu[i * p] = 0;
        } else {
          if (wantPhi) phi[i * p] = phi[i] * (p - 1);
          if (wantMu)  mu[i * p] = -mu[i];
        }
      }
    }
  }

  bool isPrime(int x) const { return x >= 2 && spf[x] == x; }

  // {prime, exponent} pairs, ascending. O(log x).
  vector<pii> factor(int x) const {
    vector<pii> r;
    while (x > 1) {
      int p = spf[x], e = 0;
      while (x % p == 0) x /= p, e++;
      r.push_back({p, e});
    }
    return r;
  }
  // flat list with multiplicity: 12 -> {2, 2, 3}
  vi factorList(int x) const {
    vi r;
    while (x > 1) { r.push_back(spf[x]); x /= spf[x]; }
    return r;
  }
  // every divisor, sorted. O(d(x) log d(x)).
  vi divisors(int x) const {
    vi d = {1};
    for (auto [p, e] : factor(x)) {
      int sz = d.size(), pk = 1;
      rep(k, 0, e) { pk *= p; rep(i, 0, sz) d.push_back(d[i] * pk); }
    }
    sort(all(d));
    return d;
  }
  int numDivisors(int x) const {
    int r = 1;
    for (auto [p, e] : factor(x)) r *= e + 1;
    return r;
  }
  int sumDivisors(int x) const {
    int r = 1;
    for (auto [p, e] : factor(x)) {
      int t = 1, pk = 1;
      rep(k, 0, e) { pk *= p; t += pk; }
      r *= t;
    }
    return r;
  }
  int phiOf(int x) const {
    int r = x;
    for (auto [p, e] : factor(x)) { (void)e; r -= r / p; }
    return r;
  }
  int muOf(int x) const {
    int r = 1;
    for (auto [p, e] : factor(x)) { (void)p; if (e > 1) return 0; r = -r; }
    return r;
  }
};

// ---- range sieving ---------------------------------------------------------
// Primality of every integer in [L, R], indexed by (x - L). Needs only
// O(sqrt R) memory for the base primes plus O(R - L) for the answer, so it
// handles R up to ~1e12 as long as the WIDTH R-L stays around 1e7.
vector<char> segmentedSieve(int L, int R) {
  int lim = sqrtl((long double)R) + 1;
  vector<char> mark(lim + 1, false);
  vi base;
  for (int i = 2; i <= lim; i++)
    if (!mark[i]) {
      base.push_back(i);
      for (int j = i * i; j <= lim; j += i) mark[j] = true;
    }
  vector<char> isPrime(R - L + 1, true);
  for (int p : base)
    for (int j = max(p * p, (L + p - 1) / p * p); j <= R; j += p) isPrime[j - L] = false;
  if (L <= 1) for (int x = L; x <= min(R, 1LL); x++) isPrime[x - L] = false;
  return isPrime;
}
vi primesInRange(int L, int R) {
  auto ip = segmentedSieve(L, R);
  vi res;
  rep(i, 0, R - L + 1) if (ip[i]) res.push_back(L + i);
  return res;
}

// pi(n) = how many primes are <= n, by sieving n in blocks of S. O(n log log n)
// time but only O(sqrt n + S) memory, so n ~ 1e11 is reachable (slowly).
int countPrimes(int n) {
  const int S = 10000;
  vi base;
  int nsqrt = sqrtl((long double)n);
  vector<char> ip(nsqrt + 2, true);
  for (int i = 2; i <= nsqrt; i++)
    if (ip[i]) {
      base.push_back(i);
      for (int j = i * i; j <= nsqrt; j += i) ip[j] = false;
    }
  int result = 0;
  vector<char> block(S);
  for (int k = 0; k * S <= n; k++) {
    fill(all(block), true);
    int start = k * S;
    for (int p : base) {
      int startIdx = (start + p - 1) / p;
      int j = max(startIdx, p) * p - start;
      for (; j < S; j += p) block[j] = false;
    }
    if (k == 0) block[0] = block[1] = false;
    for (int i = 0; i < S && start + i <= n; i++) result += block[i];
  }
  return result;
}

// ---- ALL divisors of EVERY i <= n ------------------------------------------
// When you need the actual divisor LISTS (not just counts) for every number up
// to n, don't call divisors(x) in a loop -- sieve them all at once. Each d is
// pushed onto its own multiples, so the total work is the harmonic series:
//   sum_{d<=n} n/d  ~  n ln n     (n = 1e6 -> 1.4e7 entries, not 1e12)
//
// This flat/CSR layout is the one to use: the divisors of i live contiguously
// in dv[st[i] .. st[i+1]). One pass counts the sizes, a prefix sum turns those
// into offsets, and a second pass fills. Because d ascends, each row comes out
// already SORTED -- no sorting anywhere.
//
//   DivTable T(1000000);
//   for (int d : T[12]) ...        // 1 2 3 4 6 12
//   T[12].size()                   // 6
//
// Measured at n = 1e6:  61 ms, 60 MB.
// The obvious `vector<vector<int>> d(n+1); ... d[j].push_back(i);` one-liner is
// 491 ms and 92 MB -- 1e6 vectors each doubling their capacity as they grow.
// Pre-reserving exact sizes (divisorTable below) gets that to 95 ms / 80 MB.
struct DivTable {
  vector<int32_t> st;    // st[i] = where i's divisors start in dv
  vector<int32_t> dv;    // all divisor lists, concatenated

  DivTable(int n) : st(n + 2, 0) {
    for (int i = 1; i <= n; i++)
      for (int j = i; j <= n; j += i) st[j + 1]++;          // count
    for (int i = 1; i <= n + 1; i++) st[i] += st[i - 1];    // -> offsets
    dv.resize(st[n + 1]);
    vector<int32_t> pos(st.begin(), st.begin() + n + 1);
    for (int i = 1; i <= n; i++)
      for (int j = i; j <= n; j += i) dv[pos[j]++] = i;     // fill, ascending
  }

  // a lightweight view so `for (int d : T[i])` and `T[i].size()` both work
  struct Row {
    const int32_t *b, *e;
    const int32_t *begin() const { return b; }
    const int32_t *end() const { return e; }
    int size() const { return e - b; }
    int32_t operator[](int k) const { return b[k]; }
  };
  Row operator[](int i) const { return {dv.data() + st[i], dv.data() + st[i + 1]}; }
  int count(int i) const { return st[i + 1] - st[i]; }
};

// Same data as vector-of-vectors, if you need to MUTATE the rows. Sizes are
// counted first so no vector ever reallocates -- that alone is a 5x speedup
// over the naive push_back version. Rows come out sorted.
vector<vector<int32_t>> divisorTable(int n) {
  vector<int32_t> cnt(n + 1, 0);
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j += i) cnt[j]++;
  vector<vector<int32_t>> d(n + 1);
  for (int i = 1; i <= n; i++) d[i].reserve(cnt[i]);
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j += i) d[j].push_back(i);
  return d;
}
// MEMORY is the binding constraint here, not time: the entry count grows like
// n ln n, so n = 1e6 costs ~60MB and n = 1e7 costs ~750MB (hopeless). If you
// only need counts or sums use the two tables below (4 bytes per number); if
// you need the divisors of just a FEW numbers, use Sieve::divisors(x) or the
// sqrt-trial-division divisorsSqrt() in factor.cpp instead.

// ---- divisor tables --------------------------------------------------------
// d[i] = number of divisors of i, sigma[i] = sum of divisors, for ALL i <= n.
// The harmonic-series loop is O(n log n) -- much simpler than factorising each
// number, and fast enough for n = 1e7.
vector<int32_t> divisorCountTable(int n) {
  vector<int32_t> d(n + 1, 0);
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j += i) d[j]++;
  return d;
}
vi divisorSumTable(int n) {
  vi s(n + 1, 0);
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j += i) s[j] += i;
  return s;
}

// ---- divisor / multiple transforms ----------------------------------------
// These four turn "sum over divisors" and "sum over multiples" into O(n log log
// n) passes -- the number-theoretic analogue of a SOS DP, and the standard way
// to answer "how many pairs have gcd exactly k".
//
//   divisorZeta:    g[i] = sum_{d | i} f[d]
//   divisorMobius:  the exact inverse of divisorZeta
//   multipleZeta:   g[i] = sum_{i | m} f[m]
//   multipleMobius: the exact inverse of multipleZeta
//
// Recipe -- count pairs with gcd exactly k:
//   c[d] = # of array elements divisible by d;  exact[d] = C(c[d], 2);
//   multipleMobius(exact, primes);  // now exact[d] counts pairs with gcd == d
void divisorZeta(vi &f, const vector<int32_t> &primes) {
  int n = f.size() - 1;
  for (int p : primes) for (int i = 1; i * p <= n; i++) f[i * p] += f[i];
}
void divisorMobius(vi &f, const vector<int32_t> &primes) {
  int n = f.size() - 1;
  for (int p : primes) for (int i = n / p; i >= 1; i--) f[i * p] -= f[i];
}
void multipleZeta(vi &f, const vector<int32_t> &primes) {
  int n = f.size() - 1;
  for (int p : primes) for (int i = n / p; i >= 1; i--) f[i] += f[i * p];
}
void multipleMobius(vi &f, const vector<int32_t> &primes) {
  int n = f.size() - 1;
  for (int p : primes) for (int i = 1; i * p <= n; i++) f[i] -= f[i * p];
}

// ---- divisor-block enumeration --------------------------------------------
// floor(n/i) only takes O(sqrt n) distinct values. This loop walks the maximal
// blocks where it is constant, which is how you evaluate sums like
//   sum_{i=1}^{n} floor(n/i)   or   sum_{i=1}^{n} mu(i) * floor(n/i)
// in O(sqrt n) instead of O(n).
//
// The weighted version needs one extra ingredient: to collapse a whole block to
// a single term you must be able to sum the weight over [l, r] in O(1), i.e.
// you need PREFIX SUMS of the weight. Unweighted (weight 1) that is just
// (r - l + 1); for mu it is the Mertens function -- see CoprimeCounter below.
int divisorBlockSum(int n) {                 // sum_{i=1}^{n} floor(n/i)
  int res = 0;
  for (int l = 1, r; l <= n; l = r + 1) {
    int q = n / l;
    r = n / q;                               // largest i with n/i == q
    res += q * (r - l + 1);
  }
  return res;
}

// ---- counting by gcd -------------------------------------------------------
// # of pairs (i, j) in [1,n]^2 with gcd(i,j) == 1 is
//   sum_{d=1}^{n} mu(d) * floor(n/d)^2
// because sum_{d | g} mu(d) = [g == 1] turns "gcd is 1" into inclusion-exclusion
// over the common divisor d, and floor(n/d)^2 counts the pairs both of whose
// entries are multiples of d.
//
// Evaluating that sum term by term is O(n). Block-walking floor(n/d) and using
// Mertens M(x) = sum_{k<=x} mu(k) collapses each block to one term:
//   block [l, r] where floor(n/d) == q contributes (M(r) - M(l-1)) * q^2
// so each QUERY is O(sqrt n). Building M is a linear sieve, O(N).
//
// Use it when you answer many queries, or need gcd == k for several k. For one
// query at one n the plain O(n) loop over mu[] is simpler and just as fast.
//
//   CoprimeCounter C(1e7);
//   C(n);            // pairs in [1,n]^2 with gcd exactly 1
//   C.gcdExactly(n, k);
//
// MEMORY: peak is during construction -- spf (4 bytes/number) + mu (1) + M (4),
// so ~9N bytes. N = 1e7 is 90MB; the Sieve is freed on the way out, leaving 4N.
// M is int32_t on purpose: |M(x)| stays tiny (well under 1e5 for x <= 1e16), so
// long long would just double the footprint for nothing.
struct CoprimeCounter {
  int n;
  vector<int32_t> M;                         // M[i] = sum_{k=1}^{i} mu(k)

  CoprimeCounter(int n) : n(n), M(n + 1, 0) {
    Sieve S(n, false, true);
    for (int i = 1; i <= n; i++) M[i] = M[i - 1] + S.mu[i];
  }

  // pairs (i, j) in [1,m]^2 with gcd(i,j) == 1. Requires m <= n.
  int operator()(int m) const {
    int res = 0;
    for (int l = 1, r; l <= m; l = r + 1) {
      int q = m / l;
      r = m / q;
      res += (int)(M[r] - M[l - 1]) * q * q;
    }
    return res;
  }

  // pairs (i, j) in [1,m]^2 with gcd(i,j) == k exactly. Both entries must be
  // multiples of k, and dividing through by k makes them coprime.
  int gcdExactly(int m, int k) const {
    return k > m ? 0 : (*this)(m / k);
  }
};

// ---- ONE side fixed --------------------------------------------------------
// CoprimeCounter lets both entries range over [1,n]. If instead x is FIXED and
// you want   # of y in [1,N] with gcd(y, x) == 1,   the inclusion-exclusion runs
// over the distinct primes of x rather than over every d <= N:
//   sum_{d | rad(x)} mu(d) * floor(N / d)
// That is O(2^w) with w = #distinct primes -- at most 15 for x <= 1e18, 9 for
// x <= 1e9. NOTHING is sieved; you only need x's prime factors, from
// Sieve::factor (needs x <= N) or factor.cpp's Pollard rho (x up to 1e18).
//
// Sanity checks: x == 1 gives N; N == x gives exactly phi(x).
//
// Everything here is a PREFIX count, so a half-open range [A,B] is just
// upTo(B) - upTo(A-1). Same trick gives sums, and inverting upTo by binary
// search gives the k-th coprime.
//
// Build the struct once per fixed x and reuse it -- the 2^w divisor list is the
// only real work. The free functions below are one-shot wrappers.
struct FixedCoprime {
  vi sd;                    // squarefree divisors of rad(x), signed by their mu

  // primesOfX may contain repeats: Sieve::factorList and factor.cpp's factor()
  // both return primes WITH multiplicity, and {2,2,3} would otherwise
  // double-subtract the multiples of 2. Dedupe so either can be passed straight
  // in. Entries <= 1 are dropped -- a stray 1 would put +1 and -1 in sd and
  // silently zero every answer.
  FixedCoprime(vi primesOfX) {
    sd.push_back(1);                         // the empty product, mu(1) = +1
    sort(all(primesOfX));
    primesOfX.erase(unique(all(primesOfX)), primesOfX.end());
    for (int p : primesOfX) {
      if (p <= 1) continue;
      int sz = sd.size();
      rep(i, 0, sz) sd.push_back(-sd[i] * p);
    }
  }

  // # of y in [1,N] with gcd(y, x) == 1.
  int upTo(int N) const {
    if (N <= 0) return 0;                    // truncation makes N<0 nonsense
    int res = 0;
    for (int d : sd) res += d > 0 ? N / d : -(N / -d);
    return res;
  }
  // # of y in [A,B] with gcd(y, x) == 1. A may be <= 0; empty range gives 0.
  int inRange(int A, int B) const {
    return A > B ? 0 : upTo(B) - upTo(A - 1);
  }

  // SUM of the y in [1,N] with gcd(y, x) == 1, via
  //   sum_{d | rad(x)} mu(d) * d * T(floor(N/d)),   T(k) = k(k+1)/2
  // since the multiples of d contribute d * (1 + 2 + ... + floor(N/d)).
  // OVERFLOW: the answer is ~0.3 * N^2, so it only fits in a long long for
  // N up to ~4e9. The internal __int128 protects the intermediate T(), not the
  // return value -- past that you want this reduced mod p instead.
  int sumUpTo(int N) const {
    if (N <= 0) return 0;
    __int128 res = 0;
    for (int d : sd) {
      int a = d > 0 ? d : -d, q = N / a;
      __int128 t = (__int128)q * (q + 1) / 2 * a;
      res += d > 0 ? t : -t;
    }
    return (int)res;
  }
  int sumInRange(int A, int B) const {
    return A > B ? 0 : sumUpTo(B) - sumUpTo(A - 1);
  }

  // The k-th smallest positive y with gcd(y, x) == 1, 1-indexed. upTo is
  // non-decreasing, so binary search inverts it. 1 is always coprime, so the
  // count grows without bound and the doubling terminates.
  int nth(int k) const {
    if (k <= 0) return 0;
    int lo = 1, hi = 1;
    while (upTo(hi) < k) hi *= 2;
    while (lo < hi) {                        // smallest y with upTo(y) >= k
      int mid = lo + (hi - lo) / 2;
      if (upTo(mid) >= k) hi = mid; else lo = mid + 1;
    }
    return lo;
  }
};

// One-shot wrappers. If you query the same x more than once, build a
// FixedCoprime instead of paying for the divisor list every call.
int countCoprimeTo(int N, const vi &primesOfX) {
  return FixedCoprime(primesOfX).upTo(N);
}
int countCoprimeToRange(int A, int B, const vi &primesOfX) {
  return FixedCoprime(primesOfX).inRange(A, B);
}

// # of y in [1,N] with gcd(y, x) == g exactly, for fixed x. Zero unless g | x:
// y must be a multiple of g, and dividing through leaves y/g coprime to x/g.
// Pass the primes of x/g (repeats fine), NOT of x.
int countGcdWith(int N, int g, const vi &primesOfXoverG) {
  return g > N ? 0 : countCoprimeTo(N / g, primesOfXoverG);
}
// Same over [A,B]: y = g*z with z in [ceil(A/g), floor(B/g)] and gcd(z, x/g) 1.
int countGcdWithRange(int A, int B, int g, const vi &primesOfXoverG) {
  if (A > B) return 0;
  int lo = A <= 0 ? 1 : (A + g - 1) / g;     // ceil, guarding A <= 0
  return FixedCoprime(primesOfXoverG).inRange(lo, B / g);
}

// ============================================================================
//  FACTORISATION OF ONE BIG NUMBER  -- Miller-Rabin + Pollard's rho
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Primality and full factorisation for a SINGLE n up to 1e18, where a sieve
//    is hopeless. Miller-Rabin with the standard 12 bases is DETERMINISTIC for
//    every 64-bit n (not probabilistic -- there are no false answers).
//    If you need to factor many numbers all below ~1e7, build a Sieve instead
//    (sieve.cpp): its spf table is far faster per query.
//
//  API
//    isPrimeMR(n)            deterministic for all n < 2^64
//    factor(n)   -> sorted vector<int> with multiplicity, {} for n == 1
//    factorMap(n)-> {{p, e}, ...}
//    divisorsSqrt(n) -> EVERY divisor (not just primes), O(sqrt n) -- the default
//    divisors(n) -> every divisor via factorisation, O(n^1/4), for huge n
//    phi(n)  numDivisors(n)  sumDivisors(n)
//    fermat(n)               finds a factor fast IF n = a*b with a, b close
//
//  COMPLEXITY  isPrimeMR O(12 log^3 n) bit ops -- microseconds.
//              factor O(n^(1/4)) expected; ~1e18 factors in a few ms.
//
//  PITFALLS
//    * factor(n) is RECURSIVE and randomised in effect; it is expected-time,
//      not worst-case. Fine for 1e5 calls on 1e18 inputs.
//    * n must be > 0. factor(1) = {} by definition.
//    * pollard() returns SOME non-trivial factor, not necessarily prime --
//      that's why factor() recurses on both halves.
//    * Everything goes through __int128 mulmod, so it is correct but ~3x
//      slower than plain multiplication. If n < 1e9, trial division to sqrt(n)
//      is simpler and quicker.
//    * fermat() loops forever-ish on a prime -- only call it when you know n is
//      composite, or bound the iterations.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define all(x) x.begin(), x.end()
typedef vector<int> vi;
typedef pair<int, int> pii;
typedef unsigned long long ull;

ull mulmodu(ull a, ull b, ull m) { return (__uint128_t)a * b % m; }
ull powmodu(ull a, ull e, ull m) {
  ull r = 1 % m; a %= m;
  while (e) { if (e & 1) r = mulmodu(r, a, m); a = mulmodu(a, a, m); e >>= 1; }
  return r;
}

// ---- Miller-Rabin ----------------------------------------------------------
// Write n-1 = d * 2^s. n is prime iff for every base a either a^d = 1 or
// a^(d*2^r) = -1 for some r < s. These 12 bases (all primes up to 37) are a
// PROVEN witness set for every n < 2^64, so no randomness is needed.
bool isPrimeMR(ull n) {
  if (n < 2) return false;
  for (ull p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37})
    if (n % p == 0) return n == p;
  ull d = n - 1;
  int s = 0;
  while (!(d & 1)) d >>= 1, s++;
  for (ull a : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
    ull x = powmodu(a, d, n);
    if (x == 1 || x == n - 1) continue;
    bool composite = true;
    rep(r, 1, s) { x = mulmodu(x, x, n); if (x == n - 1) { composite = false; break; } }
    if (composite) return false;
  }
  return true;
}

// ---- Pollard's rho (Brent's variant) --------------------------------------
// Iterate x -> x^2 + 1 mod n and look for a collision mod a hidden factor.
// The products are batched (prd) so we only pay one gcd every 40 steps.
// Returns a non-trivial divisor of a COMPOSITE n.
ull pollard(ull n) {
  auto f = [n](ull x) { return mulmodu(x, x, n) + 1; };
  ull x = 0, y = 0, t = 30, prd = 2, i = 1, q;
  while (t++ % 40 || gcd(prd, n) == 1) {
    if (x == y) x = ++i, y = f(x);
    if ((q = mulmodu(prd, x > y ? x - y : y - x, n))) prd = q;
    x = f(x), y = f(f(y));
  }
  return gcd(prd, n);
}

// Full prime factorisation, sorted, with multiplicity.  360 -> {2,2,2,3,3,5}
vi factor(int n) {
  if (n <= 1) return {};
  if (isPrimeMR(n)) return {n};
  // strip the small primes first -- rho is slow on tiny factors
  for (int p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37})
    if (n % p == 0) { vi r = factor(n / p); r.push_back(p); sort(all(r)); return r; }
  ull d = pollard(n);
  vi l = factor(d), r = factor(n / d);
  l.insert(l.end(), all(r));
  sort(all(l));
  return l;
}

// {prime, exponent} pairs, ascending
vector<pii> factorMap(int n) {
  vi f = factor(n);
  vector<pii> r;
  for (int p : f) {
    if (!r.empty() && r.back().first == p) r.back().second++;
    else r.push_back({p, 1});
  }
  return r;
}

// ---- divisors by trial division -------------------------------------------
// THE DEFAULT for a single n: divisors come in pairs (i, n/i) straddling
// sqrt(n), so one scan to sqrt(n) collects both halves. No factorisation, no
// dependencies, four lines.
//
// WHEN TO USE WHICH
//   n <= ~1e14, one-off        -> divisorsSqrt   (25 ms at 1e14, 9 ms at 1e12)
//   n bigger than that         -> divisors       (Pollard rho, O(n^1/4))
//   many x, all <= 1e7         -> Sieve::divisors in sieve.cpp
//   every i <= N               -> DivTable in sieve.cpp
//
// TWO TRAPS, both handled below:
//   * `i * i <= n` OVERFLOWS once n approaches 9e18. Written as i <= n/i it is
//     the same test with no product at all.
//   * a PERFECT SQUARE would otherwise emit sqrt(n) twice -- hence the
//     `i != n/i` guard. divisorsSqrt(36) must give 1 2 3 4 6 9 12 18 36.
vi divisorsSqrt(int n) {
  vi d;
  for (int i = 1; i <= n / i; i++)
    if (n % i == 0) {
      d.push_back(i);
      if (i != n / i) d.push_back(n / i);
    }
  sort(all(d));
  return d;
}
// Same thing without the final sort: the small half is already ascending and
// the large half is descending, so walking it backwards concatenates in order.
// Saves the O(sqrt(n) log) sort -- rarely matters, but it's free.
vi divisorsSqrtNoSort(int n) {
  vi lo, hi;
  for (int i = 1; i <= n / i; i++)
    if (n % i == 0) {
      lo.push_back(i);
      if (i != n / i) hi.push_back(n / i);
    }
  lo.insert(lo.end(), hi.rbegin(), hi.rend());
  return lo;
}

// every divisor of n, sorted, VIA FACTORISATION -- O(n^1/4), for n too big for
// divisorsSqrt. d(n) <= 103680 for n <= 1e18.
vi divisors(int n) {
  vi d = {1};
  for (auto [p, e] : factorMap(n)) {
    int sz = d.size(), pk = 1;
    rep(k, 0, e) { pk *= p; rep(i, 0, sz) d.push_back(d[i] * pk); }
  }
  sort(all(d));
  return d;
}

int phi(int n) {                       // Euler totient of a big n
  int r = n;
  for (auto [p, e] : factorMap(n)) { (void)e; r -= r / p; }
  return r;
}
int numDivisors(int n) {
  int r = 1;
  for (auto [p, e] : factorMap(n)) { (void)p; r *= e + 1; }
  return r;
}
int sumDivisors(int n) {
  int r = 1;
  for (auto [p, e] : factorMap(n)) {
    int t = 1, pk = 1;
    rep(k, 0, e) { pk *= p; t += pk; }
    r *= t;
  }
  return r;
}

// ---- Fermat's method -------------------------------------------------------
// n = a^2 - b^2 = (a-b)(a+b). Walks a upward from ceil(sqrt(n)) looking for a
// perfect-square b^2 = a^2 - n, so it is FAST when n has two nearly equal
// factors (e.g. a semiprime with balanced halves) and hopeless otherwise.
// Returns the smaller factor. n must be ODD and COMPOSITE.
int fermat(int n) {
  int a = ceill(sqrtl((long double)n));
  while (a * a < n) a++;
  while (true) {
    int b2 = a * a - n;
    int b = sqrtl((long double)b2);
    while (b > 0 && b * b > b2) b--;
    while ((b + 1) * (b + 1) <= b2) b++;
    if (b * b == b2) return a - b;
    a++;
  }
}

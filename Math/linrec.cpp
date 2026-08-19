// ============================================================================
//  LINEAR RECURRENCES  -- Berlekamp-Massey + Kitamasa
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Two halves of one very strong trick:
//      1. berlekampMassey(S) takes the first terms of a sequence and RECOVERS
//         the shortest linear recurrence they satisfy, mod a prime.
//      2. linearRec(S, tr, k) evaluates the k-th term of that recurrence for k
//         up to 1e18 in O(d^2 log k).
//    Together: brute-force the first ~2d terms of ANY sequence you suspect is
//    linearly recurrent (counting DPs, matrix powers, walk counts), feed them
//    in, and get the 1e18-th term without ever deriving the recurrence by hand.
//    That includes sequences whose recurrence you have no idea how to prove.
//
//  API
//    vl tr = berlekampMassey(S);     // S = first terms, mod MOD
//    ll x  = linearRec(S, tr, k);    // k-th term (0-indexed)
//
//  COMPLEXITY  BM O(n^2) on n given terms; linearRec O(d^2 log k) for a
//              degree-d recurrence. Beats the O(d^3 log k) matrix power in
//              matrix.cpp as soon as d is more than about 10.
//
//  PITFALLS
//    * MOD MUST BE PRIME (there is a modular inverse inside BM).
//    * Give BM at least 2*d terms, where d is the true recurrence order --
//      ideally more. Too few terms and it happily returns a SHORTER wrong
//      recurrence that fits your samples. When in doubt, generate 2d + 10 and
//      verify the extra terms.
//    * The convention is  S[i] = sum_j tr[j] * S[i-1-j].  tr[0] multiplies the
//      MOST RECENT term.
//    * linearRec needs sz(S) >= sz(tr); it reads S[0 .. d-1].
//    * If the sequence is not linearly recurrent at all, BM returns garbage
//      silently. Always sanity-check against a brute force.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define all(x) x.begin(), x.end()
#define sz(x) (int)(x).size()
typedef long long ll;
typedef vector<ll> vl;

const ll LMOD = 1000000007;      // MUST be prime

ll powmodL(ll a, ll e, ll m = LMOD) {
  ll r = 1 % m; a %= m; if (a < 0) a += m;
  while (e > 0) { if (e & 1) r = r * a % m; a = a * a % m; e >>= 1; }
  return r;
}

// Shortest linear recurrence for S, as coefficients tr with
//     S[i] = tr[0]*S[i-1] + tr[1]*S[i-2] + ... + tr[d-1]*S[i-d].
vl berlekampMassey(vl S) {
  int n = sz(S), L = 0, m = 0;
  vl C(n), B(n), T;
  C[0] = B[0] = 1;
  ll b = 1;
  rep(i, 0, n) {
    ++m;
    ll d = S[i] % LMOD;
    rep(j, 1, L + 1) d = (d + C[j] * S[i - j]) % LMOD;
    if (!d) continue;                       // current guess still predicts S[i]
    T = C;
    ll coef = d * powmodL(b, LMOD - 2) % LMOD;
    rep(j, m, n) C[j] = (C[j] - coef * B[j - m]) % LMOD;
    if (2 * L > i) continue;
    L = i + 1 - L; B = T; b = d; m = 0;     // the recurrence had to grow
  }
  C.resize(L + 1);
  C.erase(C.begin());
  for (ll &x : C) x = (LMOD - x) % LMOD;
  return C;
}

// k-th term (0-indexed) of the recurrence `tr` with initial values S.
// Kitamasa: compute x^k mod the characteristic polynomial, then dot it with S.
ll linearRec(const vl &S, const vl &tr, ll k) {
  int n = sz(tr);
  auto combine = [&](vl a, vl b) {
    vl res(n * 2 + 1);
    rep(i, 0, n + 1) rep(j, 0, n + 1) res[i + j] = (res[i + j] + a[i] * b[j]) % LMOD;
    for (int i = 2 * n; i > n; --i)          // reduce x^i using the recurrence
      rep(j, 0, n) res[i - 1 - j] = (res[i - 1 - j] + res[i] * tr[j]) % LMOD;
    res.resize(n + 1);
    return res;
  };
  vl pol(n + 1), e(pol);
  pol[0] = e[1] = 1;
  for (++k; k; k /= 2) {
    if (k % 2) pol = combine(pol, e);
    e = combine(e, e);
  }
  ll res = 0;
  rep(i, 0, n) res = (res + pol[i + 1] * S[i]) % LMOD;
  return res;
}

// ---- worked example --------------------------------------------------------
//   // I have a DP but n is 1e18. Brute-force the small cases first:
//   vl S;
//   for (int i = 0; i < 60; i++) S.push_back(bruteForce(i) % LMOD);
//   vl tr = berlekampMassey(S);
//   cout << tr.size() << "\n";              // <- eyeball the order; if it is
//                                           //    close to 30 you need more terms
//   cout << linearRec(S, tr, 1e18) << "\n";
//
// Fibonacci sanity check: S = {0,1,1,2,3,5,8,13} gives tr = {1,1} and
// linearRec(S, tr, 10) = 55.

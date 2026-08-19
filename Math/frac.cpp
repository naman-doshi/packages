// ============================================================================
//  FRAC  -- exact rational arithmetic (no floating point anywhere)
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Use it whenever comparing two fractions with doubles would be a precision
//    coin-flip: sorting by slope or by ratio, comparing a/b vs c/d, geometry
//    intersections on a grid, probability answers that must be exact, or any
//    "output p/q in lowest terms" problem. Always normalised: q > 0 and
//    gcd(|p|, q) = 1, so == and < mean exactly what you expect.
//
//  API
//    Frac f(3, 6);        // -> 1/2
//    f + g   f - g   f * g   f / g   -f
//    f < g   f == g       // exact, via 128-bit cross multiplication
//    f.val()              // double, only for printing/debugging
//    cout << f            // prints "p/q" (or just "p" when q == 1)
//    bestApprox(x, N)     // closest fraction to x with denominator <= N
//    contFrac(p, q)       // continued-fraction expansion of p/q
//
//  COMPLEXITY  every operation is one gcd, O(log). bestApprox is O(log N).
//
//  PITFALLS
//    * OVERFLOW is the real risk, not precision. p and q grow multiplicatively:
//      adding n fractions with distinct denominators can hit a denominator of
//      lcm(all of them). Normalising after every op helps a lot but is not a
//      guarantee -- if you're summing many terms, think about whether the final
//      denominator can fit in 64 bits, or use __int128 members.
//    * Comparison uses __int128 so it is safe even when p*d would overflow, but
//      arithmetic (+, *) is NOT protected -- it reduces first, then multiplies.
//    * Frac(1, 0) is not checked. Don't divide by zero.
//    * val() is for output only. Never sort by val().
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long

#define rep(i, a, b) for(int i = a; i < (b); ++i)
typedef vector<int> vi;
typedef pair<int, int> pii;

struct Frac {
  int p, q;                                // always q > 0, gcd(|p|, q) = 1

  Frac(int p = 0, int q = 1) : p(p), q(q) { norm(); }
  void norm() {
    if (q < 0) { p = -p; q = -q; }
    int g = gcd(p < 0 ? -p : p, q);
    if (g > 1) { p /= g; q /= g; }
    if (p == 0) q = 1;
  }

  Frac operator+(const Frac &o) const { return Frac(p * o.q + o.p * q, q * o.q); }
  Frac operator-(const Frac &o) const { return Frac(p * o.q - o.p * q, q * o.q); }
  Frac operator*(const Frac &o) const {
    // cross-reduce BEFORE multiplying to keep the numbers as small as possible
    int g1 = gcd(p < 0 ? -p : p, o.q), g2 = gcd(o.p < 0 ? -o.p : o.p, q);
    if (!g1) g1 = 1;
    if (!g2) g2 = 1;
    return Frac((p / g1) * (o.p / g2), (q / g2) * (o.q / g1));
  }
  Frac operator/(const Frac &o) const { return *this * Frac(o.q, o.p); }
  Frac operator-() const { return Frac(-p, q); }
  Frac &operator+=(const Frac &o) { return *this = *this + o; }
  Frac &operator-=(const Frac &o) { return *this = *this - o; }
  Frac &operator*=(const Frac &o) { return *this = *this * o; }
  Frac &operator/=(const Frac &o) { return *this = *this / o; }

  // exact comparison: p/q < o.p/o.q  <=>  p*o.q < o.p*q  (q, o.q > 0)
  bool operator<(const Frac &o) const { return (__int128)p * o.q < (__int128)o.p * q; }
  bool operator>(const Frac &o) const { return o < *this; }
  bool operator<=(const Frac &o) const { return !(o < *this); }
  bool operator>=(const Frac &o) const { return !(*this < o); }
  bool operator==(const Frac &o) const { return p == o.p && q == o.q; }
  bool operator!=(const Frac &o) const { return !(*this == o); }

  int floorVal() const { return p >= 0 ? p / q : -((-p + q - 1) / q); }
  int ceilVal()  const { return p >= 0 ? (p + q - 1) / q : -((-p) / q); }
  Frac abs() const { return Frac(p < 0 ? -p : p, q); }
  double val() const { return (double)p / q; }

  friend ostream &operator<<(ostream &os, const Frac &f) {
    return f.q == 1 ? os << f.p : os << f.p << "/" << f.q;
  }
};

// ---- continued fractions ---------------------------------------------------
// p/q = a0 + 1/(a1 + 1/(a2 + ...)). The convergents (partial evaluations) are
// the best rational approximations of p/q with small denominators.
vi contFrac(int p, int q) {
  vi a;
  while (q) { a.push_back(p / q); int t = p % q; p = q; q = t; }
  return a;
}

// Closest fraction to x with denominator <= N.
// Walks the continued fraction of x, keeping the convergents h/k, and stops at
// the first term that would push the denominator past N. The answer is either
// that last convergent or the best SEMICONVERGENT below it (the largest a' with
// a'*k1 + k0 <= N) -- checking both is what makes this exactly optimal rather
// than merely close. O(log N): each step is one continued-fraction term.
//   bestApprox(M_PI, 100)  = 311/99      (closer than 22/7)
//   bestApprox(M_PI, 1000) = 355/113
Frac bestApprox(double x, int N) {
  int sgn = x < 0 ? -1 : 1;
  x = fabs(x);
  int p0 = 0, q0 = 1, p1 = 1, q1 = 0;      // h(-2)/k(-2) = 0/1, h(-1)/k(-1) = 1/0
  double v = x;
  Frac best(0, 1);
  rep(iter, 0, 64) {
    int a = (int)floor(v);
    int p2 = a * p1 + p0, q2 = a * q1 + q0;
    if (q2 > N) {                          // this term overshoots the bound
      int amax = (N - q0) / q1;            // q1 > 0 once we're past iter 0
      int ps = amax * p1 + p0, qs = amax * q1 + q0;
      double e1 = fabs(x - (double)p1 / q1);
      double es = fabs(x - (double)ps / qs);
      best = es < e1 ? Frac(ps, qs) : Frac(p1, q1);
      return Frac(sgn * best.p, best.q);
    }
    p0 = p1; q0 = q1; p1 = p2; q1 = q2;
    double f = v - a;
    if (f < 1e-15) break;                  // x was exactly rational
    v = 1.0 / f;
  }
  return Frac(sgn * p1, q1);
}

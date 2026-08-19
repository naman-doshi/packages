// ============================================================================
//  INTEGER BASICS  -- safe division, exact roots, powers, series, floor sums
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    The small helpers that sit inside every other math routine: division that
//    rounds the way you actually want on negatives, EXACT integer sqrt/cbrt
//    (never trust sqrt() on a long long), 128-bit multiply-mod, fast-doubling
//    Fibonacci, closed-form series, and floorSum -- an O(log) evaluator for
//    sum floor((a*i+b)/m) that collapses most "count lattice points under a
//    line" problems into a single call.
//
//  WHAT'S HERE
//    fdiv / cdiv       floor / ceil division, correct for negative operands
//    pmod              a % m always in [0, m)
//    (gcd / lcm)       use std::gcd / std::lcm -- see the note below
//    isqrt / icbrt     exact integer roots; isSquare, ilog
//    mulmod / powmod   overflow-free a*b%m and a^e%m
//    binpow            plain a^e (no modulus)
//    fib               {F(n), F(n+1)} mod m in O(log n)
//    sumN / sumSq / sumCube / geoSum      closed forms, mod-safe variants
//    floorSum          sum_{i=0}^{n-1} floor((a*i+b)/m), O(log)
//    countMultiples    how many multiples of k lie in [l, r]
//    digitSum / numDigits / digitalRoot / toBase / fromBase
//    josephus          survivor of the k-th elimination circle (O(n) & O(k log n))
//
//  PITFALLS
//    * DON'T define your own gcd/lcm: std::gcd/std::lcm (C++17) already exist,
//      and a global gcd() becomes AMBIGUOUS with std::gcd under
//      `using namespace std;` -- a compile error at every call site.
//    * a/b in C++ truncates toward ZERO, so -7/2 == -3, not -4. Any time an
//      index or a count can go negative you want fdiv/cdiv.
//    * (int)sqrt(n) is off by one near 9e18 -- always isqrt().
//    * mulmod needs m < 2^63. powmod needs e >= 0.
//    * fib/geoSum add two residues, so they want m < 4.6e18.
//    * sumSq/sumCube overflow past n ~ 3e6 / 5e4 -- use the mod versions.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define all(x) x.begin(), x.end()
typedef vector<int> vi;
typedef pair<int, int> pii;

const int INF = 9e18;

// ---- division / rounding ---------------------------------------------------
// fdiv(-7,2) = -4, cdiv(-7,2) = -3.  Both undefined for b == 0.
int fdiv(int a, int b) { return a / b - ((a % b != 0) && ((a < 0) != (b < 0))); }
int cdiv(int a, int b) { return a / b + ((a % b != 0) && ((a < 0) == (b < 0))); }
int pmod(int a, int m) { return ((a % m) + m) % m; }   // result always in [0, m)

// round x down / up to a multiple of k (k > 0)
int floorTo(int x, int k) { return fdiv(x, k) * k; }
int ceilTo (int x, int k) { return cdiv(x, k) * k; }

// ---- gcd / lcm -------------------------------------------------------------
// Use std::gcd / std::lcm from <numeric> (C++17). Both take the absolute value
// of their arguments, and std::lcm already divides before multiplying, so it
// does NOT overflow the way a*b/g does. GCC also has __gcd(a,b), but Apple's
// libc++ static-asserts on signed types there -- prefer std::gcd everywhere.
//   gcd(a, b)   lcm(a, b)
// gcd of a whole array:  int g = 0; for (int x : v) g = gcd(g, x);
// (gcd(0, x) = x, which makes 0 the correct identity to start from.)

// ---- exact integer roots ---------------------------------------------------
// floor(sqrt(n)) with no floating-point error. The comparisons are written as
// divisions so nothing overflows even at n = 9e18.
int isqrt(int n) {
  if (n < 0) return -1;
  int x = sqrtl((long double)n);
  while (x > 0 && x > n / x) x--;
  while ((x + 1) <= n / (x + 1)) x++;
  return x;
}
bool isSquare(int n) { if (n < 0) return false; int r = isqrt(n); return r * r == n; }

// floor(cbrt(n)) for 0 <= n <= 9e18
int icbrt(int n) {
  int x = cbrtl((long double)n);
  while (x > 0 && x * x * x > n) x--;
  while ((x + 1) * (x + 1) * (x + 1) <= n) x++;
  return x;
}

// floor(log_b n) for n >= 1, b >= 2.  (For b == 2 use 63 - __builtin_clzll(n).)
int ilog(int n, int b) { int r = 0; while (n >= b) n /= b, r++; return r; }

// ---- powers ----------------------------------------------------------------
// a^e with NO modulus -- overflows silently, only for tiny exponents.
int binpow(int a, int e) {
  int r = 1;
  while (e > 0) { if (e & 1) r *= a; a *= a; e >>= 1; }
  return r;
}

// a*b % m and a^e % m for any m < 2^63. __int128 keeps the product exact.
// If m < 2^31 you can drop the __int128 and just write a * b % m.
int mulmod(int a, int b, int m) { return (__int128)a * b % m; }
int powmod(int a, int e, int m) {
  int r = 1 % m; a = pmod(a, m);
  while (e > 0) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
  return r;
}

// true iff a*b overflows long long (GCC builtin, res holds the product if not)
bool mulSafe(int a, int b, int &res) { return !__builtin_mul_overflow(a, b, &res); }

// ---- Fibonacci -------------------------------------------------------------
// Fast doubling: F(2k) = F(k)*(2*F(k+1) - F(k)),  F(2k+1) = F(k)^2 + F(k+1)^2.
// Returns {F(n), F(n+1)} mod m in O(log n).  F(93) already overflows int64, so
// there is deliberately no modulus-free version -- pass m = INF if you must.
//   Identities worth remembering:
//     gcd(F(a), F(b)) = F(gcd(a,b));  F(n) | F(kn);  sum_{i<=n} F(i) = F(n+2)-1
pii fib(int n, int m) {
  if (n == 0) return {0, 1 % m};
  auto [a, b] = fib(n >> 1, m);
  int t = pmod(b - a, m);                     // (2b - a) mod m, split to avoid
  t = (t + b) % m;                            // overflowing on 2*b
  int c = mulmod(a, t, m);                    // F(2k)
  int d = (mulmod(a, a, m) + mulmod(b, b, m)) % m;   // F(2k+1)
  return (n & 1) ? pii{d, (c + d) % m} : pii{c, d};
}

// ---- series ----------------------------------------------------------------
int sumN(int n) { return (n % 2 == 0) ? (n / 2) * (n + 1) : n * ((n + 1) / 2); } // 1+..+n
int sumSq(int n)   { return n * (n + 1) * (2 * n + 1) / 6; }   // overflows past n ~ 3e6
int sumCube(int n) { int s = sumN(n); return s * s; }          // (1+..+n)^2

// 1 + a + a^2 + ... + a^(n-1) mod m, in O(log n). Works for any m (no inverse
// needed), which is exactly when you can't just use (a^n - 1)/(a - 1).
int geoSum(int a, int n, int m) {
  if (n == 0) return 0;
  if (n & 1) return (mulmod(a, geoSum(a, n - 1, m), m) + 1) % m;
  int h = geoSum(a, n / 2, m);
  return mulmod(h, (1 + powmod(a, n / 2, m)) % m, m);
}

// ---- floor sums ------------------------------------------------------------
// sum_{i=0}^{n-1} floor((a*i + b) / m) in O(log(a+m)) -- a Euclidean-style
// recursion. n >= 0, m >= 1; a and b may be negative.
// Use it for: lattice points under the line y = (a x + b)/m, counting i with
// (a*i+b) mod m < k (= floorSum(...,+k) - floorSum(...)), sums of floor(n/i)
// style expressions, and "how many multiples of m are hit by an AP".
int floorSum(int n, int m, int a, int b) {
  int ans = 0;
  if (a < 0) { int a2 = pmod(a, m); ans -= n * (n - 1) / 2 * ((a2 - a) / m); a = a2; }
  if (b < 0) { int b2 = pmod(b, m); ans -= n * ((b2 - b) / m);               b = b2; }
  while (true) {
    if (a >= m) { ans += n * (n - 1) / 2 * (a / m); a %= m; }
    if (b >= m) { ans += n * (b / m);               b %= m; }
    int ymax = a * n + b;
    if (ymax < m) break;
    n = ymax / m; b = ymax % m; swap(m, a);
  }
  return ans;
}

// how many multiples of k lie in [l, r]  (k > 0, l <= r, negatives fine)
int countMultiples(int l, int r, int k) { return fdiv(r, k) - fdiv(l - 1, k); }

// ---- digits ----------------------------------------------------------------
int digitSum(int n)  { n = llabs(n); int s = 0; while (n) s += n % 10, n /= 10; return s; }
int numDigits(int n) { n = llabs(n); int c = 0; do c++, n /= 10; while (n); return c; }
// repeated digit sum until one digit; equals 1 + (n-1) % 9 for n > 0
int digitalRoot(int n) { return n == 0 ? 0 : 1 + (n - 1) % 9; }

string toBase(int n, int b) {                  // b in [2, 36]
  if (n == 0) return "0";
  bool neg = n < 0; n = llabs(n);
  string s;
  while (n) { int d = n % b; s += (d < 10 ? char('0' + d) : char('A' + d - 10)); n /= b; }
  if (neg) s += '-';
  reverse(all(s));
  return s;
}
int fromBase(const string &s, int b) {
  int n = 0, i = 0, sgn = 1;
  if (s[0] == '-') sgn = -1, i = 1;
  for (; i < (int)s.size(); i++) {
    int d = isdigit(s[i]) ? s[i] - '0' : toupper(s[i]) - 'A' + 10;
    n = n * b + d;
  }
  return sgn * n;
}

// ---- Josephus --------------------------------------------------------------
// n people in a circle, every k-th is removed. Returns the 0-indexed survivor.
int josephus(int n, int k) {                   // O(n)
  int r = 0;
  rep(i, 2, n + 1) r = (r + k) % i;
  return r;
}
int josephusFast(int n, int k) {               // O(k log n) -- for huge n
  if (n == 1) return 0;
  if (k == 1) return n - 1;
  if (k > n) return (josephusFast(n - 1, k) + k) % n;
  int res = josephusFast(n - n / k, k) - n % k;
  if (res < 0) res += n; else res += res / (k - 1);
  return res;
}

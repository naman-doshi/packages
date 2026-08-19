// ============================================================================
//  BIGINT  -- arbitrary-precision integers, base 1e9
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    For the rare problem where the answer genuinely doesn't fit in 64 bits and
//    isn't asked for modulo anything: exact factorials, huge Fibonacci /
//    Catalan numbers, big-number I/O.
//
//    CHECK THESE FIRST -- they solve most "it overflows" moments and are far
//    less code:
//      * __int128         (up to ~1.7e38, prints via a helper, no I/O support)
//      * long double      (64-bit mantissa; fine when you only need magnitude)
//      * work mod a prime (almost every counting problem wants this anyway)
//      * work in logs     (comparing products/powers without computing them)
//
//  API
//    BigInt a("123456789012345678901234567890"), b = 12345LL;
//    a + b   a - b   a * b   a / b   a % b   -a
//    a < b   a == b        a.pow(k)      a.isZero()
//    a / 1000LL   a % 1000LL           // fast small-integer forms
//    cout << a;   cin >> a;   a.toString()
//
//  COMPLEXITY  add/sub O(n).  multiply O(n*m) schoolbook -- 10^4 digits is
//              fine, 10^6 is not (that needs FFT, see fft.cpp).
//              BigInt / BigInt is O(n*m*30) via per-digit binary search: fine
//              for a few thousand digits, deliberately simple over fast.
//              / and % by a machine integer are O(n).
//
//  PITFALLS
//    * Division truncates TOWARD ZERO and the remainder takes the sign of the
//      dividend, matching C++ -- NOT floor division.
//    * The internal limb type is int32 in base 1e9; a limb product is taken in
//      long long. Don't change BASE without rechecking every carry.
//    * Zero is stored as an EMPTY limb vector with sign = +1, so `a.a.back()`
//      is invalid on zero -- use isZero().
//    * No square root, no gcd. Add them only if a problem actually needs them.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;

typedef long long ll;

struct BigInt {
  static const int BASE = 1000000000, DIG = 9;
  vector<int> a;                     // little-endian base-1e9 limbs
  int sign = 1;                      // +1 or -1; zero is (empty, +1)

  BigInt() {}
  BigInt(ll v) { *this = v; }
  BigInt(const string &s) { read(s); }

  BigInt &operator=(ll v) {
    sign = 1;
    if (v < 0) { sign = -1; v = -v; }
    a.clear();
    for (; v > 0; v /= BASE) a.push_back(v % BASE);
    return *this;
  }
  void trim() {
    while (!a.empty() && !a.back()) a.pop_back();
    if (a.empty()) sign = 1;
  }
  bool isZero() const { return a.empty(); }

  void read(const string &s) {
    sign = 1; a.clear();
    int pos = 0;
    while (pos < (int)s.size() && (s[pos] == '-' || s[pos] == '+')) {
      if (s[pos] == '-') sign = -sign;
      pos++;
    }
    for (int i = (int)s.size() - 1; i >= pos; i -= DIG) {
      int x = 0;
      for (int j = max(pos, i - DIG + 1); j <= i; j++) x = x * 10 + (s[j] - '0');
      a.push_back(x);
    }
    trim();
  }
  string toString() const {
    string s = (sign < 0 && !isZero()) ? "-" : "";
    s += a.empty() ? "0" : to_string(a.back());
    for (int i = (int)a.size() - 2; i >= 0; i--) {
      string t = to_string(a[i]);
      s += string(DIG - t.size(), '0') + t;
    }
    return s;
  }

  // ---- comparison ----
  static int cmpAbs(const BigInt &x, const BigInt &y) {
    if (x.a.size() != y.a.size()) return x.a.size() < y.a.size() ? -1 : 1;
    for (int i = (int)x.a.size() - 1; i >= 0; i--)
      if (x.a[i] != y.a[i]) return x.a[i] < y.a[i] ? -1 : 1;
    return 0;
  }
  bool operator==(const BigInt &o) const { return sign == o.sign && a == o.a; }
  bool operator!=(const BigInt &o) const { return !(*this == o); }
  bool operator<(const BigInt &o) const {
    if (sign != o.sign) return sign < o.sign;
    int c = cmpAbs(*this, o);
    return sign > 0 ? c < 0 : c > 0;
  }
  bool operator>(const BigInt &o) const { return o < *this; }
  bool operator<=(const BigInt &o) const { return !(o < *this); }
  bool operator>=(const BigInt &o) const { return !(*this < o); }

  // ---- addition / subtraction ----
  static BigInt addAbs(const BigInt &x, const BigInt &y) {
    BigInt r; r.a = x.a;
    int carry = 0;
    for (size_t i = 0; i < y.a.size() || carry; i++) {
      if (i == r.a.size()) r.a.push_back(0);
      ll cur = (ll)r.a[i] + carry + (i < y.a.size() ? y.a[i] : 0);
      r.a[i] = cur % BASE;
      carry = cur / BASE;
    }
    return r;
  }
  static BigInt subAbs(const BigInt &x, const BigInt &y) {   // needs |x| >= |y|
    BigInt r; r.a = x.a;
    int borrow = 0;
    for (size_t i = 0; i < y.a.size() || borrow; i++) {
      ll cur = (ll)r.a[i] - borrow - (i < y.a.size() ? y.a[i] : 0);
      if (cur < 0) { cur += BASE; borrow = 1; } else borrow = 0;
      r.a[i] = cur;
    }
    r.trim();
    return r;
  }
  BigInt operator-() const {
    BigInt r = *this;
    if (!r.isZero()) r.sign = -r.sign;
    return r;
  }
  BigInt operator+(const BigInt &o) const {
    BigInt r;
    if (sign == o.sign) { r = addAbs(*this, o); r.sign = sign; }
    else if (cmpAbs(*this, o) >= 0) { r = subAbs(*this, o); r.sign = sign; }
    else { r = subAbs(o, *this); r.sign = o.sign; }
    r.trim();
    return r;
  }
  BigInt operator-(const BigInt &o) const { return *this + (-o); }

  // ---- multiplication ----
  BigInt operator*(const BigInt &o) const {
    if (isZero() || o.isZero()) return BigInt();
    BigInt r;
    r.a.assign(a.size() + o.a.size(), 0);
    for (size_t i = 0; i < a.size(); i++) {
      ll carry = 0;
      for (size_t j = 0; j < o.a.size() || carry; j++) {
        ll cur = r.a[i + j] + carry + (j < o.a.size() ? (ll)a[i] * o.a[j] : 0);
        r.a[i + j] = cur % BASE;
        carry = cur / BASE;
      }
    }
    r.sign = sign * o.sign;
    r.trim();
    return r;
  }

  // ---- division by a machine integer (|v| < BASE) ----
  BigInt operator/(ll v) const {
    BigInt r = *this;
    if (v < 0) { r.sign = -r.sign; v = -v; }
    ll rem = 0;
    for (int i = (int)a.size() - 1; i >= 0; i--) {
      ll cur = a[i] + rem * BASE;
      r.a[i] = cur / v;
      rem = cur % v;
    }
    r.trim();
    return r;
  }
  ll operator%(ll v) const {
    if (v < 0) v = -v;
    ll rem = 0;
    for (int i = (int)a.size() - 1; i >= 0; i--) rem = (a[i] + rem * BASE) % v;
    return rem * sign;
  }

  // ---- full division: {quotient, remainder}, truncating toward zero ----
  // Long division, one base-1e9 digit at a time, each found by binary search.
  static pair<BigInt, BigInt> divmod(const BigInt &A, const BigInt &B) {
    BigInt babs = B; babs.sign = 1;
    BigInt q, r;
    q.a.assign(A.a.size(), 0);
    for (int i = (int)A.a.size() - 1; i >= 0; i--) {
      r.a.insert(r.a.begin(), A.a[i]);       // r = r * BASE + digit
      r.trim();
      int lo = 0, hi = BASE - 1, d = 0;
      while (lo <= hi) {                     // largest d with babs*d <= r
        int mid = lo + (hi - lo) / 2;
        if (cmpAbs(babs * BigInt((ll)mid), r) <= 0) { d = mid; lo = mid + 1; }
        else hi = mid - 1;
      }
      q.a[i] = d;
      r = subAbs(r, babs * BigInt((ll)d));
    }
    q.sign = A.sign * B.sign; q.trim();
    r.sign = A.sign; r.trim();
    return {q, r};
  }
  BigInt operator/(const BigInt &o) const { return divmod(*this, o).first; }
  BigInt operator%(const BigInt &o) const { return divmod(*this, o).second; }

  BigInt pow(ll e) const {
    BigInt r(1), b = *this;
    while (e > 0) { if (e & 1) r = r * b; b = b * b; e >>= 1; }
    return r;
  }

  friend ostream &operator<<(ostream &os, const BigInt &x) { return os << x.toString(); }
  friend istream &operator>>(istream &is, BigInt &x) { string s; is >> s; x.read(s); return is; }
};

// ============================================================================
//  MINT  -- a number that carries its modulus around with it
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Drop-in replacement for `int` in any counting DP or combinatorics problem
//    with a "print the answer mod 1e9+7" tail. Every operator reduces for you,
//    so the entire class of "forgot one % MOD" and "went negative after a
//    subtraction" bugs disappears. Division works too (it multiplies by the
//    modular inverse), which is what makes nCr read like actual maths.
//
//  API
//    using mint = Mint<1000000007>;    // or Mint<998244353> for NTT problems
//    mint a = 5, b = 7;
//    a + b   a - b   a * b   a / b     // and the compound forms
//    a.pow(k)   a.inv()   -a           // pow takes negative exponents too
//    cin >> a;  cout << a;             // reads/prints the raw value
//    (long long)a                      // explicit cast back out
//    vector<mint> dp(n);               // works in containers, default is 0
//
//  COMPLEXITY  +,-,* are O(1);  / and inv() are O(log MOD) (one powmod).
//
//  PITFALLS
//    * MOD MUST BE PRIME for / and inv() -- they use Fermat. For a composite
//      modulus use modinv() from modular.cpp (extended Euclid) instead.
//    * MOD must be < 3.03e9 so that v * o.v fits in a signed 64-bit product.
//    * NEVER divide by zero-mod-p: mint(MOD) / mint(MOD) silently gives 0, not
//      an error. Check the denominator when it can vanish (e.g. C(n,k)/k).
//    * Comparisons compare RESIDUES, so a < b is only meaningful for sorting
//      or use in a std::set -- it is not the ordering of the original numbers.
//    * `mint x = 1e18;` works (the constructor reduces), but `mint x = a * b;`
//      where a, b are raw long longs overflows BEFORE the constructor runs.
//      Convert first: mint(a) * mint(b).
//    * Two different Mint<> instantiations are different types and won't mix.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;

template <long long MOD> struct Mint {
  long long v;

  Mint() : v(0) {}
  Mint(long long x) : v(x % MOD) { if (v < 0) v += MOD; }

  Mint &operator+=(const Mint &o) { v += o.v; if (v >= MOD) v -= MOD; return *this; }
  Mint &operator-=(const Mint &o) { v -= o.v; if (v < 0) v += MOD; return *this; }
  Mint &operator*=(const Mint &o) { v = v * o.v % MOD; return *this; }
  Mint &operator/=(const Mint &o) { return *this *= o.inv(); }

  friend Mint operator+(Mint a, const Mint &b) { return a += b; }
  friend Mint operator-(Mint a, const Mint &b) { return a -= b; }
  friend Mint operator*(Mint a, const Mint &b) { return a *= b; }
  friend Mint operator/(Mint a, const Mint &b) { return a /= b; }
  Mint operator-() const { return Mint(0) - *this; }

  // a^e, e may be negative (then it inverts first)
  Mint pow(long long e) const {
    if (e < 0) return inv().pow(-e);
    Mint r = 1, b = *this;
    while (e > 0) { if (e & 1) r *= b; b *= b; e >>= 1; }
    return r;
  }
  Mint inv() const { return pow(MOD - 2); }   // MOD MUST BE PRIME

  bool operator==(const Mint &o) const { return v == o.v; }
  bool operator!=(const Mint &o) const { return v != o.v; }
  bool operator<(const Mint &o) const { return v < o.v; }   // for set/map only
  explicit operator long long() const { return v; }

  friend ostream &operator<<(ostream &os, const Mint &m) { return os << m.v; }
  friend istream &operator>>(istream &is, Mint &m) { long long x; is >> x; m = Mint(x); return is; }
};

using mint = Mint<1000000007>;
using mint2 = Mint<998244353>;      // the NTT-friendly prime, = 119 * 2^23 + 1

// ---- runtime modulus -------------------------------------------------------
// When the modulus is read from input you can't use a template parameter.
// Set MDYN once at the top of main, then use dmint exactly like mint.
long long MDYN = 1000000007;
struct DMint {
  long long v;
  DMint() : v(0) {}
  DMint(long long x) : v(x % MDYN) { if (v < 0) v += MDYN; }
  DMint &operator+=(const DMint &o) { v += o.v; if (v >= MDYN) v -= MDYN; return *this; }
  DMint &operator-=(const DMint &o) { v -= o.v; if (v < 0) v += MDYN; return *this; }
  DMint &operator*=(const DMint &o) { v = (__int128)v * o.v % MDYN; return *this; }
  DMint &operator/=(const DMint &o) { return *this *= o.inv(); }
  friend DMint operator+(DMint a, const DMint &b) { return a += b; }
  friend DMint operator-(DMint a, const DMint &b) { return a -= b; }
  friend DMint operator*(DMint a, const DMint &b) { return a *= b; }
  friend DMint operator/(DMint a, const DMint &b) { return a /= b; }
  DMint operator-() const { return DMint(0) - *this; }
  DMint pow(long long e) const {
    if (e < 0) return inv().pow(-e);
    DMint r = 1, b = *this;
    while (e > 0) { if (e & 1) r *= b; b *= b; e >>= 1; }
    return r;
  }
  DMint inv() const { return pow(MDYN - 2); }   // MDYN MUST BE PRIME
  bool operator==(const DMint &o) const { return v == o.v; }
  friend ostream &operator<<(ostream &os, const DMint &m) { return os << m.v; }
};

// ---- recipes ---------------------------------------------------------------
//  factorials + nCr (see combinatorics.cpp for the full struct):
//    vector<mint> f(N+1, 1), fi(N+1);
//    for (int i = 1; i <= N; i++) f[i] = f[i-1] * i;
//    fi[N] = f[N].inv();
//    for (int i = N; i > 0; i--) fi[i-1] = fi[i] * i;     // one inverse total
//    auto C = [&](int n, int k) { return k < 0 || k > n ? mint(0) : f[n]*fi[k]*fi[n-k]; };
//
//  matrix of mint: Matrix in matrix.cpp works verbatim if you swap `int` for
//  `mint` and delete every `% MOD` -- the operators already do it.
//
//  geometric series 1 + r + ... + r^(n-1):  r == 1 ? mint(n) : (r.pow(n)-1)/(r-1)

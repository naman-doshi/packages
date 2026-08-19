// ============================================================================
//  CONVOLUTIONS  -- FFT, NTT, arbitrary-mod FFT, and the bitwise transforms
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    c[k] = sum_{i+j=k} a[i]*b[j] in O(n log n) instead of O(n^2). Reach for it
//    whenever an answer is "for every k, count pairs/ways that sum to k":
//      - polynomial multiplication
//      - counting pairs (i, j) with a[i] + b[j] == k  (indicator arrays)
//      - distributions of a sum of independent random variables
//      - string matching with wildcards, and shifted correlations (reverse one
//        array first: correlation is convolution with b reversed)
//      - big-integer multiplication (see bigint.cpp)
//    The bitwise transforms at the bottom do the same for i|j, i&j and i^j.
//
//  WHICH ONE
//    conv(a, b)          doubles. Fast, but PRECISION-limited (see pitfalls).
//    convNTT(a, b)       exact, mod 998244353. Use this whenever the problem
//                        already asks for that modulus -- it is the safest.
//    convMod<M>(a, b)    exact for ANY modulus M < ~1e9, 3 FFTs, ~2x slower.
//    xorConv / orConv / andConv    bitwise, arrays must have length 2^k.
//
//  COMPLEXITY  O(n log n) with n rounded up to a power of two.
//              n = 1e6 is a few hundred ms; 1e5 is instant.
//
//  PITFALLS
//    * DOUBLE FFT PRECISION is the classic wrong-answer source. The rule of
//      thumb: it is safe while n * max|a| * max|b| stays under about 1e14 (an
//      error of ~1e-9 * that product must round to the right integer). Values
//      near 1e9 with n = 1e5 are ALREADY unsafe -- use convNTT or convMod.
//    * Always round with llround / (ll)(x + 0.5), never a bare cast.
//    * The result has size a.size() + b.size() - 1. Empty input -> empty output.
//    * fft()/ntt() work in place on a POWER-OF-TWO length; the conv wrappers
//      handle the padding, the raw transforms do not.
//    * NTT's modulus 998244353 = 119 * 2^23 + 1 supports lengths up to 2^23.
//      Other usable primes: 5<<25|1, 7<<26|1, 479<<21|1, 483<<21|1 (root 3).
//    * The bitwise transforms need EXACTLY length 2^k, and xorConv's inverse
//      divides by n -- do that in the ring you're working in (for a modulus,
//      multiply by inv(n) instead).
// ============================================================================
#include <bits/stdc++.h>
using namespace std;

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define all(x) x.begin(), x.end()
#define sz(x) (int)(x).size()
typedef long long ll;
typedef vector<ll> vl;
typedef vector<double> vd;
typedef vector<int> vi;
typedef complex<double> C;

// ---- complex FFT -----------------------------------------------------------
// In-place, a.size() must be a power of two. Roots are cached across calls in
// long double for accuracy, then stored as doubles for speed.
void fft(vector<C> &a) {
  int n = sz(a), L = 31 - __builtin_clz(n);
  static vector<complex<long double>> R(2, 1);
  static vector<C> rt(2, 1);
  for (static int k = 2; k < n; k *= 2) {
    R.resize(n); rt.resize(n);
    auto x = polar(1.0L, acos(-1.0L) / k);
    rep(i, k, 2 * k) rt[i] = R[i] = i & 1 ? R[i / 2] * x : R[i / 2];
  }
  vi rev(n);
  rep(i, 0, n) rev[i] = (rev[i / 2] | (i & 1) << L) / 2;
  rep(i, 0, n) if (i < rev[i]) swap(a[i], a[rev[i]]);
  for (int k = 1; k < n; k *= 2)
    for (int i = 0; i < n; i += 2 * k) rep(j, 0, k) {
      C z = rt[j + k] * a[i + j + k];
      a[i + j + k] = a[i + j] - z;
      a[i + j] += z;
    }
}

// Convolution of two real sequences. Packs a into the real part and b into the
// imaginary part, so it costs only TWO transforms instead of three.
vd conv(const vd &a, const vd &b) {
  if (a.empty() || b.empty()) return {};
  vd res(sz(a) + sz(b) - 1);
  int L = 32 - __builtin_clz(sz(res)), n = 1 << L;
  vector<C> in(n), out(n);
  copy(all(a), begin(in));
  rep(i, 0, sz(b)) in[i].imag(b[i]);
  fft(in);
  for (C &x : in) x *= x;
  rep(i, 0, n) out[i] = in[-i & (n - 1)] - conj(in[i]);
  fft(out);
  rep(i, 0, sz(res)) res[i] = imag(out[i]) / (4 * n);
  return res;
}
// Integer wrapper: rounds properly. Check the precision rule above first.
vl convInt(const vl &a, const vl &b) {
  vd x(all(a)), y(all(b));
  vd z = conv(x, y);
  vl res(z.size());
  rep(i, 0, sz(z)) res[i] = llround(z[i]);
  return res;
}

// ---- NTT (exact, mod 998244353) -------------------------------------------
const ll NMOD = 998244353, NROOT = 62;   // 998244353 = (119 << 23) + 1
ll powmodN(ll a, ll e, ll m) {
  ll r = 1 % m; a %= m;
  while (e > 0) { if (e & 1) r = r * a % m; a = a * a % m; e >>= 1; }
  return r;
}
void ntt(vl &a) {
  int n = sz(a), L = 31 - __builtin_clz(n);
  static vl rt(2, 1);
  for (static int k = 2, s = 2; k < n; k *= 2, s++) {
    rt.resize(n);
    vl z = {1, powmodN(NROOT, NMOD >> s, NMOD)};
    rep(i, k, 2 * k) rt[i] = rt[i / 2] * z[i & 1] % NMOD;
  }
  vi rev(n);
  rep(i, 0, n) rev[i] = (rev[i / 2] | (i & 1) << L) / 2;
  rep(i, 0, n) if (i < rev[i]) swap(a[i], a[rev[i]]);
  for (int k = 1; k < n; k *= 2)
    for (int i = 0; i < n; i += 2 * k) rep(j, 0, k) {
      ll z = rt[j + k] * a[i + j + k] % NMOD, &ai = a[i + j];
      a[i + j + k] = ai - z + (z > ai ? NMOD : 0);
      ai += (ai + z >= NMOD ? z - NMOD : z);
    }
}
vl convNTT(const vl &a, const vl &b) {
  if (a.empty() || b.empty()) return {};
  int s = sz(a) + sz(b) - 1, B = 32 - __builtin_clz(s), n = 1 << B;
  ll inv = powmodN(n, NMOD - 2, NMOD);
  vl L(a), R(b), out(n);
  L.resize(n); R.resize(n);
  ntt(L); ntt(R);
  rep(i, 0, n) out[-i & (n - 1)] = L[i] * R[i] % NMOD * inv % NMOD;
  ntt(out);
  return {out.begin(), out.begin() + s};
}

// ---- arbitrary-modulus convolution ----------------------------------------
// Splits every coefficient as a = hi*cut + lo with cut ~ sqrt(M), so each half
// is small enough for exact double arithmetic, then recombines. Correct for any
// M below ~1e9 and n up to ~1e6.
template <int M> vl convMod(const vl &a, const vl &b) {
  if (a.empty() || b.empty()) return {};
  vl res(sz(a) + sz(b) - 1);
  int B = 32 - __builtin_clz(sz(res)), n = 1 << B, cut = int(sqrt(M));
  vector<C> L(n), R(n), outs(n), outl(n);
  rep(i, 0, sz(a)) L[i] = C((int)(a[i] / cut), (int)(a[i] % cut));
  rep(i, 0, sz(b)) R[i] = C((int)(b[i] / cut), (int)(b[i] % cut));
  fft(L); fft(R);
  rep(i, 0, n) {
    int j = -i & (n - 1);
    outl[j] = (L[i] + conj(L[j])) * R[i] / (2.0 * n);
    outs[j] = (L[i] - conj(L[j])) * R[i] / (2.0 * n) * C(0, -1);  // "/ i"
  }
  fft(outl); fft(outs);
  rep(i, 0, sz(res)) {
    ll av = llround(real(outl[i])), cv = llround(imag(outs[i]));
    ll bv = llround(imag(outl[i])) + llround(real(outs[i]));
    res[i] = ((av % M * cut + bv) % M * cut + cv) % M;
  }
  return res;
}

// ---- bitwise convolutions --------------------------------------------------
// c[k] = sum_{i OP j == k} a[i] * b[j], for OP in {xor, or, and}.
// All arrays must have length exactly 2^B (index space = subsets of B bits).
//
// xor uses the Walsh-Hadamard transform; or/and use the subset/superset zeta
// transform (a "SOS DP") and its Mobius inverse. Same shape every time:
// transform both, multiply pointwise, transform back.
void wht(vl &a, bool inverse) {                  // XOR
  int n = sz(a);
  for (int len = 1; 2 * len <= n; len <<= 1)
    for (int i = 0; i < n; i += 2 * len)
      rep(j, i, i + len) {
        ll u = a[j], v = a[j + len];
        a[j] = u + v; a[j + len] = u - v;
      }
  if (inverse) for (ll &x : a) x /= n;           // mod: multiply by inv(n)
}
void subsetZeta(vl &a, bool inverse) {           // OR: a[S] = sum_{T subset S}
  int n = sz(a);
  for (int len = 1; 2 * len <= n; len <<= 1)
    for (int i = 0; i < n; i += 2 * len)
      rep(j, i, i + len) a[j + len] += inverse ? -a[j] : a[j];
}
void supersetZeta(vl &a, bool inverse) {         // AND: a[S] = sum_{T superset S}
  int n = sz(a);
  for (int len = 1; 2 * len <= n; len <<= 1)
    for (int i = 0; i < n; i += 2 * len)
      rep(j, i, i + len) a[j] += inverse ? -a[j + len] : a[j + len];
}
vl bitConv(vl a, vl b, int type) {   // 0 = xor, 1 = or, 2 = and
  int n = sz(a);
  auto fwd = [&](vl &v, bool inv) {
    if (type == 0) wht(v, inv);
    else if (type == 1) subsetZeta(v, inv);
    else supersetZeta(v, inv);
  };
  fwd(a, false); fwd(b, false);
  rep(i, 0, n) a[i] *= b[i];
  fwd(a, true);
  return a;
}

// ---- recipes ---------------------------------------------------------------
//  polynomial power:      repeated squaring with conv, or exp(k * log(P)).
//  a[i] + b[j] == k:      indicator arrays, then convolve.
//  a[i] - b[j] == k:      reverse b, convolve, and shift the index by sz(b)-1.
//  count distinct sums of a multiset: convolve the indicator with itself.
//  string match with '?':  treat '?' as 0 and match sum (s[i]-t[j])^2 s[i] t[j]
//                          -- expand into three convolutions.
//  large power of a single polynomial mod x^n: use convNTT in a loop; there is
//  no need for a full polynomial-inverse library at contest scale.

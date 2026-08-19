// ============================================================================
//  MODULAR ARITHMETIC  -- inverses, CRT, Diophantine, discrete log & sqrt
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Everything built on the extended Euclidean algorithm plus the two "hard"
//    modular questions: given a^x = b (mod m) find x (discrete log), and given
//    x^2 = a (mod p) find x (modular square root).
//
//  WHAT'S HERE
//    extgcd(a,b,x,y)      a*x + b*y = gcd(a,b)
//    modinv(a, m)         inverse for ANY m with gcd(a,m)=1 (-1 if none)
//    modinvPrime(a, p)    inverse when p is prime (Fermat), one powmod
//    invTable(n, p)       inverses of 1..n in O(n)
//    modLinear(a,b,m)     ALL x with a*x = b (mod m)
//    diophantine(...)     one solution of a*x + b*y = c, plus how to shift it
//    crt2 / crt           merge congruences x = r_i (mod m_i), non-coprime OK
//    discreteLog(a,b,m)   smallest x >= 0 with a^x = b (mod m), any m
//    sqrtMod(a, p)        x with x^2 = a (mod p), p odd prime
//    primitiveRoot(p)     generator of Z_p*
//    discreteRoot(k,a,p)  all x with x^k = a (mod p)
//    eulerPhi(n)          phi of a single n by trial division
//
//  KEY FACTS
//    * Euler:  a^phi(m) = 1 (mod m) when gcd(a,m)=1, so a^-1 = a^(phi(m)-1).
//      Fermat is the case m = p prime, phi = p-1, so a^-1 = a^(p-2).
//    * Exponents reduce mod phi(m), NOT mod m: a^e = a^(e mod phi(m)) (mod m)
//      -- only valid when gcd(a,m)=1. Otherwise use the lifting-the-exponent
//      form: for e >= log2(m), a^e = a^(phi(m) + e mod phi(m)) (mod m).
//    * a/b (mod m) is a * modinv(b, m); it only exists when gcd(b,m)=1.
//
//  PITFALLS
//    * modinv returns -1 when the inverse doesn't exist -- check it.
//    * Everything assumes m >= 1 and normalises inputs, but exponents must be
//      non-negative.
//    * discreteLog is O(sqrt(m)) with a hash map; fine to m ~ 1e12.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define all(x) x.begin(), x.end()
typedef vector<int> vi;
typedef pair<int, int> pii;

int pmod(int a, int m) { return ((a % m) + m) % m; }
int mulmod(int a, int b, int m) { return (__int128)a * b % m; }
int powmod(int a, int e, int m) {
  int r = 1 % m; a = pmod(a, m);
  while (e > 0) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
  return r;
}

// ---- extended Euclid -------------------------------------------------------
// Returns g = gcd(a,b) and fills x, y with a*x + b*y = g.
// |x| <= b/(2g), |y| <= a/(2g), so nothing overflows for a, b < 9e18.
int extgcd(int a, int b, int &x, int &y) {
  if (b == 0) { x = 1; y = 0; return a; }
  int x1, y1, g = extgcd(b, a % b, x1, y1);
  x = y1; y = x1 - (a / b) * y1;
  return g;
}

// ---- inverses --------------------------------------------------------------
// Inverse of a mod m for ANY modulus, provided gcd(a, m) = 1. Returns -1 if not.
int modinv(int a, int m) {
  int x, y, g = extgcd(pmod(a, m), m, x, y);
  return g == 1 ? pmod(x, m) : -1;
}
// Inverse mod a PRIME p (Fermat). Slightly shorter, same speed class.
int modinvPrime(int a, int p) { return powmod(a, p - 2, p); }

// inv[i] = i^-1 mod p for i in [1, n], built in O(n). p must be prime and > n.
// Derivation: p = q*i + r  =>  0 = q*i + r  =>  i^-1 = -q * r^-1.
vi invTable(int n, int p) {
  vi inv(n + 1);
  inv[1] = 1;
  rep(i, 2, n + 1) inv[i] = (p - (p / i) * inv[p % i] % p) % p;
  return inv;
}

// ---- linear congruences ----------------------------------------------------
// All solutions of a*x = b (mod m), returned as residues in [0, m).
// There are g = gcd(a, m) of them (spaced m/g apart) when g | b, else none.
vi modLinear(int a, int b, int m) {
  a = pmod(a, m); b = pmod(b, m);
  int x, y, g = extgcd(a, m, x, y);
  if (b % g) return {};
  vi res;
  int x0 = mulmod(pmod(x, m), b / g, m);
  rep(i, 0, g) res.push_back((x0 + i * (m / g)) % m);
  sort(all(res));
  return res;
}

// One integer solution of a*x + b*y = c (a, b not both 0). Returns false if
// gcd(a,b) does not divide c. The full family is
//     x + k*(b/g),  y - k*(a/g)   for every integer k,
// which is how you search for the solution with x >= 0, or minimal |x| + |y|.
bool diophantine(int a, int b, int c, int &x, int &y, int &g) {
  g = extgcd(llabs(a), llabs(b), x, y);
  if (c % g) return false;
  x *= c / g; y *= c / g;
  if (a < 0) x = -x;
  if (b < 0) y = -y;
  return true;
}

// ---- Chinese Remainder Theorem --------------------------------------------
// Merge x = r1 (mod m1) and x = r2 (mod m2). m1, m2 NEED NOT be coprime.
// Returns {remainder, lcm(m1,m2)}, or {0, -1} if the two are inconsistent.
pii crt2(int r1, int m1, int r2, int m2) {
  int x, y, g = extgcd(m1, m2, x, y);
  if ((r2 - r1) % g) return {0, -1};
  int l = m1 / g * m2;                       // lcm, assumed to fit in int64
  int t = mulmod(pmod((r2 - r1) / g, m2 / g), pmod(x, m2 / g), m2 / g);
  return {pmod(r1 + (__int128)t % l * m1 % l, l), l};
}
// Fold a whole list. Same contract: {0, -1} means no solution exists.
pii crt(const vi &r, const vi &m) {
  pii cur = {0, 1};
  rep(i, 0, (int)r.size()) {
    cur = crt2(cur.first, cur.second, pmod(r[i], m[i]), m[i]);
    if (cur.second == -1) return cur;
  }
  return cur;
}

// ---- Euler's totient -------------------------------------------------------
// phi(n) = # of x in [1, n] with gcd(x, n) = 1, by trial division O(sqrt n).
// For every n up to N at once use the linear sieve in sieve.cpp; for n up to
// 1e18 factor with Pollard rho (factor.cpp) and use phi = n * prod (1 - 1/p).
int eulerPhi(int n) {
  int res = n;
  for (int p = 2; p * p <= n; p++)
    if (n % p == 0) { while (n % p == 0) n /= p; res -= res / p; }
  if (n > 1) res -= res / n;
  return res;
}

// ---- discrete logarithm ----------------------------------------------------
// Baby-step giant-step: smallest x >= 0 with a^x = b (mod m), or -1.
// Works for ANY m (the leading gcd-stripping loop handles gcd(a,m) > 1).
// O(sqrt(m)) time and memory.
int discreteLog(int a, int b, int m) {
  a = pmod(a, m); b = pmod(b, m);
  int k = 1, add = 0, g;
  while ((g = gcd(a, m)) > 1) {            // strip common factors
    if (b == k) return add;
    if (b % g) return -1;
    b /= g; m /= g; add++;
    k = mulmod(k, a / g, m);
  }
  int n = sqrtl((long double)m) + 1;
  unordered_map<int, int> vals;
  int cur = b;
  rep(q, 0, n + 1) { vals[cur] = q; cur = mulmod(cur, a, m); }   // baby steps
  int an = powmod(a, n, m);
  cur = k;
  rep(p, 1, n + 1) {                                            // giant steps
    cur = mulmod(cur, an, m);
    if (vals.count(cur)) {
      int ans = n * p - vals[cur] + add;
      if (ans >= 0) return ans;
    }
  }
  return -1;
}

// ---- modular square root ---------------------------------------------------
// Legendre symbol: 1 if a is a non-zero QR mod p, p-1 (= -1) if not, 0 if p | a.
int legendre(int a, int p) { return powmod(a, (p - 1) / 2, p); }

// Tonelli-Shanks: some x with x^2 = a (mod p) for an ODD PRIME p, or -1 if a is
// not a quadratic residue. The other root is p - x.
int sqrtMod(int a, int p) {
  a = pmod(a, p);
  if (a == 0) return 0;
  if (p == 2) return a;
  if (legendre(a, p) != 1) return -1;
  if (p % 4 == 3) return powmod(a, (p + 1) / 4, p);       // the common fast case
  int s = p - 1, e = 0;
  while (s % 2 == 0) s /= 2, e++;
  int n = 2;
  while (legendre(n, p) != p - 1) n++;                    // any non-residue
  int x = powmod(a, (s + 1) / 2, p), b = powmod(a, s, p), g = powmod(n, s, p);
  for (int r = e;;) {
    int t = b, m = 0;
    for (; m < r && t != 1; m++) t = mulmod(t, t, p);
    if (m == 0) return x;
    int gs = powmod(g, 1LL << (r - m - 1), p);
    g = mulmod(gs, gs, p);
    x = mulmod(x, gs, p);
    b = mulmod(b, g, p);
    r = m;
  }
}

// ---- primitive roots and k-th roots ---------------------------------------
// Smallest generator g of Z_p* (p prime): g^0..g^(p-2) hits every non-zero
// residue. Exists for p, p^k, 2p^k and 4. O(ans * log p * log p).
int primitiveRoot(int p) {
  if (p == 2) return 1;
  vi fac;
  int phi = p - 1, n = phi;
  for (int i = 2; i * i <= n; i++)
    if (n % i == 0) { fac.push_back(i); while (n % i == 0) n /= i; }
  if (n > 1) fac.push_back(n);
  for (int g = 2; g <= p; g++) {
    bool ok = true;
    for (int f : fac) if (powmod(g, phi / f, p) == 1) { ok = false; break; }
    if (ok) return g;
  }
  return -1;
}

// All x with x^k = a (mod p), p prime. Reduces to a discrete log in the
// exponent: write x = g^y and solve k*y = log_g(a) (mod p-1).
vi discreteRoot(int k, int a, int p) {
  a = pmod(a, p);
  if (a == 0) return {0};
  int g = primitiveRoot(p);
  int l = discreteLog(powmod(g, k, p), a, p);            // (g^k)^l = a
  if (l == -1) return {};
  vi res;
  int x0 = powmod(g, l, p), step = powmod(g, (p - 1) / gcd(k, p - 1), p);
  int cnt = gcd(k, p - 1);
  int cur = x0;
  rep(i, 0, cnt) { res.push_back(cur); cur = mulmod(cur, step, p); }
  sort(all(res));
  res.erase(unique(all(res)), res.end());
  return res;
}

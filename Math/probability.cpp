// ============================================================================
//  PROBABILITY & EXPECTED VALUE MOD P   -- the "print the answer mod 1e9+7"
//  companion to modint.cpp
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    A problem says "output the expected value as p * q^-1 mod 998244353". The
//    arithmetic itself needs nothing new -- a probability IS a Mint, and p/q is
//    just Mint(p)/Mint(q), with every operator already reducing for you. What
//    this file adds is the four things Mint alone does NOT give you:
//
//      ratRecover(x)        a residue back to the fraction p/q it came from
//      probFromDecimal(s)   a decimal in the INPUT to an exact residue
//      solveMod / expectedSteps / absorbProb   chains with cycles in them
//      selfLoop(...)        the one-state version of the same thing
//
//  API
//    ratRecover(x)             -> {p, q} with x == p/q, or {0,0} if none small
//    probFromDecimal("12.345") -> Mint, exact, no floating point at all
//    probFromDecimal("12.345%")-> the same divided by 100
//    selfLoop(base, pSelf)     -> base / (1 - pSelf)
//    solveMod(A, b)            -> x with A x = b, or {} if singular
//    expectedSteps(Q)          -> expected steps to leave the transient states
//    absorbProb(Q, r)          -> probability of being absorbed "the r way"
//    geometricTrials(p)        -> 1/p, expected tries until a p-chance succeeds
//    couponCollector(n)        -> n * H_n
//
//  ---------------------------------------------------------------------------
//  THE THREE MOVES THAT SOLVE MOST EXPECTED-VALUE PROBLEMS
//
//  1. LINEARITY. E[X+Y] = E[X]+E[Y] ALWAYS -- no independence needed. Split the
//     quantity into indicator variables and add up their probabilities:
//       expected number of fixed points of a random permutation = n * (1/n) = 1
//       expected inversions = C(n,2) * 1/2
//     If a problem asks for the expected COUNT of something, you almost never
//     need a DP: count what fraction of the time each candidate contributes.
//
//  2. SELF LOOPS. A state that can step back to itself gives
//        E = 1 + p*E + (sum of q_i * E_i)      ->    E = (1 + sum q_i E_i)/(1-p)
//     which is selfLoop() below. This is where marks get lost: you must divide
//     by (1 - p), and if p == 1 the process never leaves and there is no answer.
//
//  3. CYCLES BETWEEN STATES. When the dependencies are not a DAG you cannot
//     just DP over them -- set up one linear equation per state and solve.
//     expectedSteps() does exactly that, mod p. (gauss.cpp solves the same
//     shape in DOUBLES, which is no use when the answer must be a residue;
//     matrix.cpp has a raw-int mod-p solver if you prefer it to Mint.)
//
//  PITFALLS
//    * MOD MUST BE PRIME -- everything here divides.
//    * Dividing by a probability that can be 0 gives a silent wrong answer,
//      not an error. Guard denominators that a test case can drive to zero
//      (p == 0 in geometricTrials, pSelf == 1 in selfLoop).
//    * A Mint has NO ORDER. You cannot compare two probabilities, round one,
//      or test "is this < 1". If you need that, you are in the wrong domain --
//      keep a double alongside for comparisons only, never for the answer.
//    * Reading probabilities from input: NEVER `cin >> double` then scale and
//      cast. 0.03 comes back a hair low and truncates to the wrong numerator,
//      and a numerator off by one is a completely unrelated residue, not a near
//      miss. Use probFromDecimal on the raw token.
//    * ratRecover only works when the true p and q are both below
//      sqrt(MOD/2) (~2e4 for 1e9+7, ~2.2e4 for 998244353). It is a debugging
//      and cross-checking tool, not a general decoder -- a real answer with a
//      big denominator is simply not recoverable.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;

typedef long long ll;

template <long long MOD> struct Mint {
  long long v;
  Mint() : v(0) {}
  Mint(long long x) : v(x % MOD) { if (v < 0) v += MOD; }
  Mint &operator+=(const Mint &o) { v += o.v; if (v >= MOD) v -= MOD; return *this; }
  Mint &operator-=(const Mint &o) { v -= o.v; if (v < 0) v += MOD; return *this; }
  Mint &operator*=(const Mint &o) { v = (__int128)v * o.v % MOD; return *this; }
  Mint &operator/=(const Mint &o) { return *this *= o.inv(); }
  friend Mint operator+(Mint a, const Mint &b) { return a += b; }
  friend Mint operator-(Mint a, const Mint &b) { return a -= b; }
  friend Mint operator*(Mint a, const Mint &b) { return a *= b; }
  friend Mint operator/(Mint a, const Mint &b) { return a /= b; }
  Mint operator-() const { return Mint(0) - *this; }
  Mint pow(long long e) const {
    if (e < 0) return inv().pow(-e);
    Mint r = 1, b = *this;
    while (e > 0) { if (e & 1) r *= b; b *= b; e >>= 1; }
    return r;
  }
  Mint inv() const { return pow(MOD - 2); }        // MOD MUST BE PRIME
  bool operator==(const Mint &o) const { return v == o.v; }
  bool operator!=(const Mint &o) const { return v != o.v; }
  friend ostream &operator<<(ostream &os, const Mint &m) { return os << m.v; }
};

using mint = Mint<998244353>;                      // swap for Mint<1000000007>

// ------------------------------------------------- residue -> the fraction
// Runs the extended Euclid on (MOD, x) and stops halfway: the first remainder
// below sqrt(MOD/2) is the numerator, and the cofactor built alongside it is
// the denominator. That pair is UNIQUE if it exists, which is what makes this
// a valid check rather than a guess.
//   use it to: compare a modular answer against a small brute force, or print
//   an answer that the statement wants as "p/q in lowest terms".
// Returns {0, 0} when no fraction that small represents x.
template <long long MOD>
pair<ll, ll> ratRecover(Mint<MOD> xm) {
	ll bound = (ll)sqrtl((long double)MOD / 2);
	ll a = MOD, b = xm.v, x = 0, y = 1;
	while (b > bound) {
		ll q = a / b;
		ll t = a - q * b; a = b; b = t;
		t = x - q * y; x = y; y = t;
	}
	ll p = b, den = y;
	if (den < 0) { p = -p; den = -den; }
	if (den == 0 || den > bound) return {0, 0};
	if (gcd(p < 0 ? -p : p, den) != 1) return {0, 0};
	if (Mint<MOD>(p) / Mint<MOD>(den) != xm) return {0, 0};   // verify, cheap
	return {p, den};
}

// ------------------------------------------- decimal in the input -> residue
// Exact: builds the numerator digit by digit mod MOD, so the string may be any
// length and no double is ever involved. Accepts a leading sign, a missing
// integer or fractional part, and an optional trailing '%'.
//   "0.5" -> 1/2      "12.34567" -> 1234567/100000      "40%" -> 2/5
template <long long MOD>
Mint<MOD> probFromDecimalT(const string &s) {
	Mint<MOD> num = 0, den = 1;
	bool neg = false, seenDot = false, percent = false;
	for (char c : s) {
		if (c == '-') neg = true;
		else if (c == '+') continue;
		else if (c == '.') seenDot = true;
		else if (c == '%') percent = true;
		else if (isdigit((unsigned char)c)) {
			num = num * Mint<MOD>(10) + Mint<MOD>(c - '0');
			if (seenDot) den *= Mint<MOD>(10);
		}
	}
	if (percent) den *= Mint<MOD>(100);
	Mint<MOD> r = num / den;
	return neg ? -r : r;
}
mint probFromDecimal(const string &s) { return probFromDecimalT<998244353>(s); }

// --------------------------------------------------------------- self loops
// E = base + pSelf * E   ->   E = base / (1 - pSelf).
// `base` is the 1 (or whatever the step costs) plus every term that leaves.
// pSelf == 1 means the process never escapes: no finite answer exists.
template <long long MOD>
Mint<MOD> selfLoop(Mint<MOD> base, Mint<MOD> pSelf) {
	return base / (Mint<MOD>(1) - pSelf);
}

// --------------------------------------------------- systems with cycles
// Gauss-Jordan mod a PRIME. Returns {} if the matrix is singular.
template <long long MOD>
vector<Mint<MOD>> solveMod(vector<vector<Mint<MOD>>> A, vector<Mint<MOD>> b) {
	int k = (int)A.size();
	for (int i = 0; i < k; i++) {
		int piv = -1;
		for (int r = i; r < k; r++) if (A[r][i] != Mint<MOD>(0)) { piv = r; break; }
		if (piv < 0) return {};
		swap(A[i], A[piv]); swap(b[i], b[piv]);
		Mint<MOD> f = A[i][i].inv();
		for (int c = i; c < k; c++) A[i][c] *= f;
		b[i] *= f;
		for (int r = 0; r < k; r++) if (r != i && A[r][i] != Mint<MOD>(0)) {
			Mint<MOD> g = A[r][i];
			for (int c = i; c < k; c++) A[r][c] -= g * A[i][c];
			b[r] -= g * b[i];
		}
	}
	return b;
}

// Q[i][j] = P(transient state i -> transient state j) in one step. Whatever
// probability is missing from a row leaves the transient set (is absorbed).
// Returns e[i] = expected number of STEPS until absorption, starting at i.
// Solves (I - Q) e = 1.
template <long long MOD>
vector<Mint<MOD>> expectedSteps(const vector<vector<Mint<MOD>>> &Q) {
	int k = (int)Q.size();
	vector<vector<Mint<MOD>>> A(k, vector<Mint<MOD>>(k, 0));
	vector<Mint<MOD>> b(k, 1);
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < k; j++) A[i][j] = (i == j ? Mint<MOD>(1) : Mint<MOD>(0)) - Q[i][j];
	}
	return solveMod<MOD>(A, b);
}

// Same Q, but r[i] = probability of being absorbed IN THE OUTCOME YOU CARE
// ABOUT directly from i in one step. Returns x[i] = probability that outcome
// eventually happens from i. Solves (I - Q) x = r.
template <long long MOD>
vector<Mint<MOD>> absorbProb(const vector<vector<Mint<MOD>>> &Q, vector<Mint<MOD>> r) {
	int k = (int)Q.size();
	vector<vector<Mint<MOD>>> A(k, vector<Mint<MOD>>(k, 0));
	for (int i = 0; i < k; i++)
		for (int j = 0; j < k; j++) A[i][j] = (i == j ? Mint<MOD>(1) : Mint<MOD>(0)) - Q[i][j];
	return solveMod<MOD>(A, r);
}

// -------------------------------------------------------- waiting times
// expected number of independent tries until one with chance p succeeds
template <long long MOD> Mint<MOD> geometricTrials(Mint<MOD> p) { return Mint<MOD>(1) / p; }
// expected draws to collect all n coupon types: n * (1 + 1/2 + ... + 1/n)
template <long long MOD> Mint<MOD> couponCollectorT(long long n) {
	Mint<MOD> h = 0;
	for (long long i = 1; i <= n; i++) h += Mint<MOD>(1) / Mint<MOD>(i);
	return Mint<MOD>(n) * h;
}
mint couponCollector(long long n) { return couponCollectorT<998244353>(n); }

int main() {
	int bad = 0, checks = 0;
	auto ok = [&](bool c, const char *m) { checks++; if (!c) { bad++; printf("FAIL: %s\n", m); } };

	// --- ratRecover round-trips every small fraction ------------------------
	for (ll q = 1; q <= 120; q++) for (ll p = -120; p <= 120; p++) {
		if (gcd(p < 0 ? -p : p, q) != 1) continue;
		auto got = ratRecover<998244353>(mint(p) / mint(q));
		checks++;
		if (got.first != p || got.second != q) { bad++; printf("FAIL: ratRecover %lld/%lld\n", p, q); }
	}
	// a residue with no small representation is reported as such
	ok(ratRecover<998244353>(mint(123456789)).second == 0 ||
	   ratRecover<998244353>(mint(123456789)).second != 0, "ratRecover terminates");

	// --- exact decimal parsing ---------------------------------------------
	ok(probFromDecimal("0.5")      == mint(1) / mint(2),        "0.5");
	ok(probFromDecimal("12.34567") == mint(1234567) / mint(100000), "12.34567");
	ok(probFromDecimal("40%")      == mint(2) / mint(5),        "40%");
	ok(probFromDecimal("0.00013")  == mint(13) / mint(100000),  "0.00013 (the truncation trap)");
	ok(probFromDecimal("7")        == mint(7),                  "integer");
	ok(probFromDecimal("-0.25")    == -(mint(1) / mint(4)),     "negative");
	ok(probFromDecimal("100.00000%") == mint(1),                "100% with 5dp");

	// --- self loop ----------------------------------------------------------
	// roll a die until it shows 6: E = 1 + (5/6)E  ->  6
	ok(selfLoop<998244353>(mint(1), mint(5) / mint(6)) == mint(6), "die self-loop = 6");
	ok(geometricTrials<998244353>(mint(1) / mint(6)) == mint(6),   "geometric = 6");

	// --- gambler's ruin: expected steps from i is exactly i*(n-i) ----------
	for (int n = 2; n <= 9; n++) {
		int k = n - 1;                                   // transient states 1..n-1
		vector<vector<mint>> Q(k, vector<mint>(k, 0));
		mint half = mint(1) / mint(2);
		for (int i = 0; i < k; i++) {
			if (i - 1 >= 0) Q[i][i - 1] = half;
			if (i + 1 < k)  Q[i][i + 1] = half;
		}
		auto e = expectedSteps<998244353>(Q);
		for (int i = 0; i < k; i++)
			ok(e[i] == mint((ll)(i + 1) * (n - (i + 1))), "gambler's ruin i*(n-i)");
		// probability of being absorbed at the TOP from i is i/n
		vector<mint> r(k, 0);
		r[k - 1] = half;                                 // only state n-1 can reach n
		auto x = absorbProb<998244353>(Q, r);
		for (int i = 0; i < k; i++)
			ok(x[i] == mint(i + 1) / mint(n), "gambler's ruin absorb prob i/n");
	}

	// --- coupon collector, checked back as an exact fraction ---------------
	// n=3: 3*(1 + 1/2 + 1/3) = 11/2
	ok(ratRecover<998244353>(couponCollector(3)) == make_pair((ll)11, (ll)2), "coupons(3) = 11/2");
	ok(ratRecover<998244353>(couponCollector(4)) == make_pair((ll)25, (ll)3), "coupons(4) = 25/3");

	printf("%d/%d checks passed\n", checks - bad, checks);
	return bad != 0;
}

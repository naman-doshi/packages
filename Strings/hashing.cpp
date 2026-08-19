// ============================================================================
//  POLYNOMIAL STRING HASHING  -- O(1) substring equality after O(n) prep
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    hash(s[l..r]) in O(1), so you can compare any two substrings for equality
//    in constant time. That gives you, cheaply:
//      * "are these two substrings equal?"          -> compare hashes
//      * longest common prefix of two suffixes      -> binary search on length
//      * find all occurrences of a pattern          -> compare every window
//      * count distinct substrings of a length      -> hash them into a set
//      * palindrome check                           -> hash vs reversed hash
//      * comparing substrings lexicographically     -> LCP, then one char
//
//  API
//    StrHash h(s);
//    h.get(l, r)            // hash of s[l..r] INCLUSIVE, as a pair
//    h.get(l, r) == h.get(a, b)     // substring equality
//    h.pref(len)            // hash of the first `len` characters
//    StrHash::combine(a, la, b, lb) // hash of the concatenation
//
//  WHY TWO MODULI
//    A single 1e9 modulus loses to the birthday paradox: comparing q pairs gives
//    a collision chance around q^2 / 1e9, which at q = 1e5 is already ~1%, and
//    on Codeforces someone will hand you an ANTI-HASH TEST built against the
//    common bases. Two independent moduli square the space to ~1e18, and the
//    base is randomised per run so a fixed adversarial test can't target it.
//
//  COMPLEXITY  O(n) build, O(1) per query, O(n) memory (4 arrays).
//
//  PITFALLS
//    * get(l, r) is INCLUSIVE on both ends. get(l, l-1) (empty) returns {0,0},
//      which is fine, but l > r+1 is nonsense -- don't.
//    * Add 1 to each character (`s[i] - 'a' + 1`), never 0: a leading 'a'
//      mapped to 0 makes "a", "aa", "aaa" all hash the same.
//    * The base is random per RUN, so hashes are NOT comparable across runs and
//      must not be stored/printed. Two StrHash objects in the same run DO share
//      the base (it's static), so cross-string comparison works.
//    * Hashing proves equality only PROBABILISTICALLY. If the problem is small
//      enough for Z-function / KMP / suffix automaton, prefer those -- they are
//      exact.
//    * Everything is mod ~1e9 and multiplied in long long. Don't switch the
//      moduli to something above 3e9 or the products overflow.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef pair<ll, ll> pll;

struct StrHash {
  static const ll M1 = 1000000007, M2 = 1000000009;
  static ll B1, B2;                 // randomised once per run
  static vector<ll> p1, p2;         // p1[i] = B1^i, grown on demand

  vector<ll> h1, h2;                // prefix hashes, h[i] covers s[0..i-1]

  static void grow(size_t n) {
    while (p1.size() <= n) {
      p1.push_back(p1.back() * B1 % M1);
      p2.push_back(p2.back() * B2 % M2);
    }
  }

  StrHash() {}
  StrHash(const string &s) { build(s); }

  void build(const string &s) {
    int n = s.size();
    grow(n);
    h1.assign(n + 1, 0); h2.assign(n + 1, 0);
    for (int i = 0; i < n; i++) {
      // +1 so that a character mapping to 0 can't make "a" == "aa"
      ll c = (unsigned char)s[i] + 1;
      h1[i + 1] = (h1[i] * B1 + c) % M1;
      h2[i + 1] = (h2[i] * B2 + c) % M2;
    }
  }

  // hash of s[l..r], INCLUSIVE. Empty range (l > r) hashes to {0, 0}.
  pll get(int l, int r) const {
    if (l > r) return {0, 0};
    int len = r - l + 1;
    ll a = ((h1[r + 1] - h1[l] * p1[len]) % M1 + M1) % M1;
    ll b = ((h2[r + 1] - h2[l] * p2[len]) % M2 + M2) % M2;
    return {a, b};
  }
  pll pref(int len) const { return len <= 0 ? pll{0, 0} : get(0, len - 1); }
  size_t size() const { return h1.empty() ? 0 : h1.size() - 1; }

  // hash of (a concatenated with b), given each piece's hash and LENGTH
  static pll combine(pll a, int la, pll b, int lb) {
    (void)la;
    grow(lb);
    return {(a.first * p1[lb] + b.first) % M1, (a.second * p2[lb] + b.second) % M2};
  }
};
ll StrHash::B1 = 0, StrHash::B2 = 0;
vector<ll> StrHash::p1 = {1}, StrHash::p2 = {1};

// Call ONCE at the top of main, before building any StrHash.
// (Without it B1 = B2 = 0 and every hash collapses to the last character.)
void initHash() {
  mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
  StrHash::B1 = 131 + rng() % 1000000;
  StrHash::B2 = 137 + rng() % 1000000;
}

// ---- recipes ---------------------------------------------------------------
//  longest common prefix of s[i..] and s[j..]  (binary search on the length):
//    int lo = 0, hi = min(n-i, n-j);
//    while (lo < hi) { int m = (lo+hi+1)/2;
//      if (h.get(i,i+m-1) == h.get(j,j+m-1)) lo = m; else hi = m-1; }
//    // lo = LCP length
//
//  find every occurrence of pattern p in text t:
//    StrHash ht(t), hp(p);
//    pll target = hp.get(0, p.size()-1);
//    for (int i = 0; i + p.size() <= t.size(); i++)
//      if (ht.get(i, i+p.size()-1) == target) found at i;
//
//  is s[l..r] a palindrome:
//    build a second StrHash on the reversed string; the reverse of s[l..r] sits
//    at [n-1-r, n-1-l] there, so compare those two hashes.
//
//  count distinct substrings of length k:
//    set<pll> seen; for (i) seen.insert(h.get(i, i+k-1)); -> seen.size()

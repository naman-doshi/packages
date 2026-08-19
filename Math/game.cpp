// ============================================================================
//  COMBINATORIAL GAME THEORY  -- nim, Grundy numbers, Sprague-Grundy
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Two players, alternating moves, perfect information, no draws, and the
//    player unable to move LOSES (normal play). Then:
//
//      SPRAGUE-GRUNDY: every such game is equivalent to a nim pile of size
//      g(state), where  g(x) = mex{ g(y) : y reachable from x }.
//      A position is LOSING for the player to move iff g == 0.
//      A game that splits into INDEPENDENT sub-games has grundy value equal to
//      the XOR of the sub-games' values.
//
//    So the whole method is: define the state, compute grundy values for small
//    states by brute force, XOR the components, and look for the pattern.
//
//  API
//    mex(v)                       smallest non-negative integer not in v
//    Grundy G(n, movesFn);        G.g[x] for every state 0..n
//    subtractionGame(n, S)        grundy of a pile of size x, moves = remove s in S
//    nimWin(piles)                true if the player to move wins
//    nimMove(piles)               a winning move as {pile index, new size}
//    misereNimWin(piles)          normal nim, but taking the LAST object loses
//
//  PITFALLS
//    * Grundy values only XOR for games played SIMULTANEOUSLY and
//      INDEPENDENTLY (a move touches exactly one component). A game where a
//      move affects two piles at once is not a sum, and the theorem says
//      nothing.
//    * The theory is for NORMAL play (can't move => you lose). Misere games
//      need separate analysis; only for nim is the misere rule this simple.
//    * Loops (a state reachable from itself) break Grundy entirely -- that's a
//      "loopy game", solved with retrograde analysis instead.
//    * Brute-force grundy tables are for FINDING the pattern. Print the first
//      50 values, spot the periodicity, then prove/assume it for the real n.
//
//  PATTERNS WORTH KNOWING
//    Nim:               g(pile of n) = n;  win iff XOR of piles != 0.
//    Misere nim:        if every pile is 1, win iff the count is EVEN;
//                       otherwise the normal rule (XOR != 0) applies.
//    Staircase nim:     tokens on steps 1..n, move any number down one step.
//                       Only the ODD-indexed steps matter: XOR those counts.
//    Subtraction {1..k}: g(n) = n mod (k+1).
//    Take-away halving / "divide" games: almost always eventually periodic.
//    Wythoff: two piles, remove any amount from one OR the same from both.
//             Losing positions are ( floor(k*phi), floor(k*phi^2) ), phi the
//             golden ratio -- equivalently the pairs whose difference is k.
//    Turning games / Mock Turtles etc: grundy of a coin at position i, XORed.
//    Green Hackenbush: a stalk of length n is a nim pile of n; fuse cycles.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define all(x) x.begin(), x.end()
typedef vector<int> vi;
typedef pair<int, int> pii;

// smallest non-negative integer NOT present in v
int mex(const vi &v) {
  vector<char> seen(v.size() + 1, false);
  for (int x : v) if (x >= 0 && x <= (int)v.size()) seen[x] = true;
  rep(i, 0, (int)v.size() + 1) if (!seen[i]) return i;
  return v.size();
}

// Grundy values for states 0..n, where moves(x) lists the states reachable
// from x. Requires an acyclic game whose states are integers, each move going
// to a SMALLER state (that's what makes the bottom-up loop valid).
struct Grundy {
  vi g;
  Grundy(int n, const function<vi(int)> &moves) : g(n + 1, 0) {
    rep(x, 0, n + 1) {
      vi vals;
      for (int y : moves(x)) vals.push_back(g[y]);
      g[x] = mex(vals);
    }
  }
};

// One pile of x objects, a legal move removes some s in S. O(n |S|).
vi subtractionGame(int n, const vi &S) {
  vi g(n + 1, 0);
  rep(x, 1, n + 1) {
    vi vals;
    for (int s : S) if (s <= x) vals.push_back(g[x - s]);
    g[x] = mex(vals);
  }
  return g;
}

// ---- nim -------------------------------------------------------------------
bool nimWin(const vi &piles) {
  int x = 0;
  for (int p : piles) x ^= p;
  return x != 0;
}
// A winning move: shrink some pile p to (p ^ total). Returns {-1,-1} if lost.
// Such a pile always exists when total != 0: take the one whose highest set bit
// matches the highest set bit of total.
pii nimMove(const vi &piles) {
  int x = 0;
  for (int p : piles) x ^= p;
  if (!x) return {-1, -1};
  rep(i, 0, (int)piles.size()) if ((piles[i] ^ x) < piles[i]) return {i, piles[i] ^ x};
  return {-1, -1};
}
// Misere nim: taking the LAST object loses.
bool misereNimWin(const vi &piles) {
  int x = 0, big = 0;
  for (int p : piles) { x ^= p; if (p > 1) big++; }
  if (big) return x != 0;                    // some pile > 1: normal rule
  int ones = 0;
  for (int p : piles) ones += (p == 1);
  return ones % 2 == 0;                      // all piles are 1: parity flips
}

// ---- worked example --------------------------------------------------------
//   // "a pile of n stones; a move removes 1, 3 or 4. Who wins?"
//   vi g = subtractionGame(1000, {1, 3, 4});
//   // print g[0..20] -> 0 1 0 1 2 3 2 0 1 0 1 2 3 2 ... period 7
//   // For several piles, XOR their grundy values and test != 0.
//
//   // A game on a DAG with arbitrary states: memoise instead.
//   //   map<State, int> memo;
//   //   int grundy(State s) { if (memo.count(s)) return memo[s];
//   //     vi v; for (auto t : moves(s)) v.push_back(grundy(t));
//   //     return memo[s] = mex(v); }

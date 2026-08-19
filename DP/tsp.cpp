// ============================================================================
//  BITMASK DP over SUBSETS  /  TRAVELLING SALESMAN  -- O(2^n * n^2)
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    dp[mask][j] = cheapest way to have visited exactly the set `mask`, standing
//    on node j. This is THE shape for "visit every one of n things once, order
//    matters, n is small (<= 20)". TSP is the named case, but the same table
//    solves: shortest Hamiltonian path, minimum-cost assignment when n is tiny,
//    "cheapest order to process jobs with pairwise switch costs", and counting
//    Hamiltonian paths (swap min for +).
//
//    The loop order needs no thought: iterating mask upward is automatically a
//    valid topological order, because adding a bit always makes mask LARGER, so
//    every state it depends on was computed earlier.
//
//  API
//    TSP t(cost);            // free start: path may begin anywhere
//    TSP t(cost, 0);         // start pinned to node 0 (needed for cycles)
//    t.bestPath()            // cheapest Hamiltonian PATH cost, INF if none
//    t.bestCycle()           // cheapest Hamiltonian CYCLE through the start
//    t.pathOrder()           // the visiting order achieving bestPath()
//    t.cycleOrder()          // ... achieving bestCycle(), start node first
//    t.dp[mask][j]           // the raw table, if you want to read it yourself
//
//  COMPLEXITY  O(2^n * n^2) time, O(2^n * n) memory.
//    n = 15 -> 7e6, instant.      n = 18 -> 8.5e7, fine.
//    n = 20 -> 4e8 ops AND 168MB with 8-byte cells -- too slow AND too big.
//    Past ~20 you need a different idea entirely (Held-Karp is optimal for
//    exact TSP; consider DP over profiles, or an approximation).
//
//  PITFALLS
//    * cost[i][j] == INF means "no edge"; the code checks for it, so don't
//      let INF + INF silently wrap. Keep INF at 9e18 / 4 or use a smaller
//      sentinel if you ever add two together.
//    * A Hamiltonian CYCLE needs the start pinned (a cycle has no natural
//      start, so fixing one loses nothing and halves the work). bestCycle()
//      asserts on a free-start instance.
//    * The graph may be DIRECTED -- cost need not be symmetric. Nothing here
//      assumes it is.
//    * INF is returned when no Hamiltonian path/cycle exists. Check for it.
//    * n = 0 is not handled (1 << 0 = 1 mask, no nodes). Guard upstream.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long
#define ll long long

const int INF = 2e18;

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define trav(a, x) for(auto& a : x)
#define all(x) x.begin(), x.end()
typedef vector<vector<int>> vvi;
typedef vector<pair<int, int>> vpi;
typedef vector<int> vi;
typedef pair<int, int> pii;

struct TSP {
  int n, start;
  vvi cost;        // cost[i][j] = weight of edge i -> j, INF if absent
  vvi dp;          // dp[mask][j] = min cost visiting exactly `mask`, ending at j
  vvi par;         // par[mask][j] = the node we came from, for reconstruction

  // start = -1 lets the path begin at any node; pin it for cycles.
  TSP(const vvi &c, int start = -1)
      : n(c.size()), start(start), cost(c),
        dp(1LL << n, vi(n, INF)), par(1LL << n, vi(n, -1)) {
    if (start < 0) rep(i, 0, n) dp[1LL << i][i] = 0;
    else dp[1LL << start][start] = 0;

    rep(mask, 0, 1LL << n) rep(j, 0, n) {
      if (dp[mask][j] == INF || !((mask >> j) & 1)) continue;
      rep(k, 0, n) {
        if ((mask >> k) & 1) continue;            // already visited
        if (cost[j][k] >= INF) continue;          // no edge
        int nm = mask | (1LL << k), nv = dp[mask][j] + cost[j][k];
        if (nv < dp[nm][k]) { dp[nm][k] = nv; par[nm][k] = j; }
      }
    }
  }

  int full() const { return (1LL << n) - 1; }

  // Cheapest Hamiltonian path. endNode (if given) receives the final node.
  int bestPath(int *endNode = nullptr) const {
    int best = INF, at = -1;
    rep(j, 0, n) if (dp[full()][j] < best) best = dp[full()][j], at = j;
    if (endNode) *endNode = at;
    return best;
  }
  // Cheapest Hamiltonian cycle returning to `start`. Requires a pinned start.
  int bestCycle() const {
    assert(start >= 0);
    int best = INF;
    rep(j, 0, n) {
      if (dp[full()][j] >= INF || cost[j][start] >= INF) continue;
      best = min(best, dp[full()][j] + cost[j][start]);
    }
    return best;
  }

  // Walk `par` backwards from (full mask, end) to recover the order.
  vi orderEndingAt(int end) const {
    vi order;
    int mask = full(), at = end;
    while (at != -1) { order.push_back(at); int p = par[mask][at]; mask ^= 1LL << at; at = p; }
    reverse(all(order));
    return order;
  }
  vi pathOrder() const {
    int end; if (bestPath(&end) >= INF) return {};
    return orderEndingAt(end);
  }
  // Order for the best cycle; the start node is first and is NOT repeated.
  vi cycleOrder() const {
    assert(start >= 0);
    int best = INF, at = -1;
    rep(j, 0, n) {
      if (dp[full()][j] >= INF || cost[j][start] >= INF) continue;
      if (dp[full()][j] + cost[j][start] < best) best = dp[full()][j] + cost[j][start], at = j;
    }
    return at < 0 ? vi{} : orderEndingAt(at);
  }
};

// ---- the same skeleton, without the struct ---------------------------------
// When you need a different aggregate (count paths, maximise, track a second
// dimension) it's usually faster to inline the three loops than to bend TSP:
//
//   vvi dp(1 << n, vi(n, INF));
//   rep(i, 0, n) dp[1 << i][i] = 0;
//   rep(mask, 0, 1 << n) rep(j, 0, n) {
//     if (dp[mask][j] == INF || !(mask >> j & 1)) continue;
//     rep(k, 0, n) {
//       if (mask >> k & 1) continue;
//       int nm = mask | 1 << k;
//       dp[nm][k] = min(dp[nm][k], dp[mask][j] + cost[j][k]);
//     }
//   }
//
// COUNTING instead of minimising: start dp[1<<i][i] = 1 and use += .
// OTHER SUBSET-DP SHAPES worth recognising:
//   * dp[mask] alone (no "current node") when order within the set is free --
//     e.g. partition n items into groups: iterate submasks, O(3^n).
//     for (int sub = mask; sub; sub = (sub - 1) & mask) { ... }
//   * SOS / superset sums in O(2^n * n) -- see Bitwise/bitwise.md.

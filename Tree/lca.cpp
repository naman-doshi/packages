#include <bits/stdc++.h>
using namespace std;

const int INF = 9000372036854775800;

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define trav(a, x) for(auto& a : x)
#define all(x) x.begin(), x.end()
#define print(arr) for (auto i : arr) cout << i << " "; cout << endl;
typedef vector<vector<int>> vvi;
typedef vector<pair<int, int>> vpi;
typedef vector<int> vi;
typedef pair<int, int> pii;

// Binary lifting LCA + Euler tour. O(n log n) build, O(log n) per query.
// Static structure -- if edge weights change between queries use HLD, if the
// tree shape changes use LCT.
//
// Usage:  LCA t(n); t.addEdge(u, v, w); ...; t.build();
//   lca(a,b)  dist(a,b)  wdist(a,b)  pathAgg(a,b)  ancestor(v,k)  jump(a,b,k)
//   isAncestor(a,b)  subtreeSize(v)  onPath(u,v,w)  meetPoint(a,b,c)
//   rootedLca(a,b,r)  climbWhile(v,pred)
//
// tin/tout give subtree v the contiguous range [tin[v], tout[v]), so a BIT over
// that array answers subtree queries. Path update "add x to every edge on u-v"
// offline: d[u]+=x, d[v]+=x, d[lca]-=2x, then each edge's value is the subtree
// sum below it.
struct LCA {
  // ---- customise the aggregate here ----------------------------------------
  // op must be ASSOCIATIVE and COMMUTATIVE: the two halves of a path are folded
  // in from both endpoints, so order is not preserved. max/min/sum/gcd/xor/count
  // are all fine; "first edge", "last edge", matrix product are NOT.
  // ID must satisfy op(ID, x) == x for every weight x you actually insert.
  // For path SUM prefer wdist() -- it is O(1) and leaves op free for max/min/gcd.
  static int op(int a, int b) { return max(a, b); }
  static const int ID = 0;     // max: 0 (or -INF if weights can be negative)
                               // sum: 0    min: INF    gcd: 0    xor: 0
  // -------------------------------------------------------------------------

  int n, LOG;
  vector<vpi> adj;   // adj[u] = {v, w}
  vvi up;            // up[j][v]  = 2^j-th ancestor of v, -1 if none
  vvi agg;           // agg[j][v] = op over the edges from v up to up[j][v]
  vi depth;          // edges from root
  vi wdep;           // summed edge weights from root
  vi tin, tout;      // euler tour: subtree of v is [tin[v], tout[v])

  LCA(int n) : n(n), LOG(1) {
    while ((1 << LOG) < n) LOG++;
    LOG++;                          // so 2^(LOG-1) >= n-1 = max possible depth
    adj.assign(n, {});
    depth.assign(n, 0);
    wdep.assign(n, 0);
    tin.assign(n, 0);
    tout.assign(n, 0);
    up.assign(LOG, vi(n, -1));
    agg.assign(LOG, vi(n, ID));
  }

  void addEdge(int u, int v, int w = 0) {
    adj[u].push_back({v, w});
    adj[v].push_back({u, w});
  }

  // Roots the tree at `root` (iterative, so no stack overflow), then builds the
  // jump tables. Extra components get rooted at their smallest node, so this is
  // safe on forests too.
  void build(int root = 0) {
    vector<char> seen(n, 0);
    vi it(n, 0);                    // how far into adj[u] the dfs has walked
    int timer = 0;
    rep(s, 0, n) {
      int start = (s == 0 ? root : s);
      if (seen[start]) continue;
      seen[start] = 1;
      depth[start] = wdep[start] = 0;
      up[0][start] = -1;
      agg[0][start] = ID;
      tin[start] = timer++;
      vi st = {start};
      while (!st.empty()) {
        int u = st.back();
        if (it[u] == (int)adj[u].size()) {
          tout[u] = timer;
          st.pop_back();
          continue;
        }
        pii e = adj[u][it[u]++];
        int v = e.first;
        if (seen[v]) continue;
        seen[v] = 1;
        up[0][v] = u;
        agg[0][v] = e.second;       // the single edge v -- parent(v)
        depth[v] = depth[u] + 1;
        wdep[v] = wdep[u] + e.second;
        tin[v] = timer++;
        st.push_back(v);
      }
    }
    rep(j, 1, LOG) rep(v, 0, n) {
      int mid = up[j-1][v];
      if (mid == -1) {
        up[j][v] = -1;
        agg[j][v] = ID;
      } else {
        up[j][v] = up[j-1][mid];
        agg[j][v] = op(agg[j-1][v], agg[j-1][mid]);
      }
    }
  }

  // k-th ancestor of v, or -1 if it doesn't exist
  int ancestor(int v, int k) {
    rep(j, 0, LOG) {
      if (v == -1) break;
      if (k >> j & 1) v = up[j][v];
    }
    return v;
  }

  int lca(int a, int b) {
    if (depth[a] < depth[b]) swap(a, b);
    a = ancestor(a, depth[a] - depth[b]);
    if (a == b) return a;
    for (int j = LOG - 1; j >= 0; j--) {
      if (up[j][a] != up[j][b]) a = up[j][a], b = up[j][b];
    }
    return up[0][a];
  }

  // op folded over every edge on the path a -- b. Returns ID if a == b.
  int pathAgg(int a, int b) {
    if (depth[a] < depth[b]) swap(a, b);
    int acc = ID, k = depth[a] - depth[b];
    rep(j, 0, LOG) if (k >> j & 1) {          // lift the deeper node to b's level
      acc = op(acc, agg[j][a]);
      a = up[j][a];
    }
    if (a == b) return acc;
    for (int j = LOG - 1; j >= 0; j--) {      // climb in lockstep, stopping just
      if (up[j][a] != up[j][b]) {             // below the lca
        acc = op(acc, op(agg[j][a], agg[j][b]));
        a = up[j][a], b = up[j][b];
      }
    }
    return op(acc, op(agg[0][a], agg[0][b])); // the two edges into the lca itself
  }

  // Highest ancestor u of v with pred(op over the edges v..u) still true.
  // pred MUST be monotone -- once it fails as the path grows it stays failed.
  // Returns v itself if even the first edge up fails. pred(ID) is assumed true.
  // e.g. farthest ancestor reachable using only edges <= X (with op = max):
  //   t.climbWhile(v, [&](int m) { return m <= X; })
  template <class P>
  int climbWhile(int v, P pred) {
    int acc = ID;
    for (int j = LOG - 1; j >= 0; j--) {
      if (up[j][v] == -1) continue;
      int cand = op(acc, agg[j][v]);
      if (pred(cand)) acc = cand, v = up[j][v];
    }
    return v;
  }

  int dist(int a, int b) {          // number of edges between a and b
    return depth[a] + depth[b] - 2 * depth[lca(a, b)];
  }

  int wdist(int a, int b) {         // summed edge weights between a and b
    return wdep[a] + wdep[b] - 2 * wdep[lca(a, b)];
  }

  bool isAncestor(int a, int b) {   // is a an ancestor of b (or a == b)?  O(1)
    return tin[a] <= tin[b] && tout[b] <= tout[a];
  }

  int subtreeSize(int v) { return tout[v] - tin[v]; }

  bool onPath(int u, int v, int w) {          // does w lie on the path u -- v?
    return dist(u, w) + dist(w, v) == dist(u, v);
  }

  // The one node where the paths between a, b, c all meet. Of the three pairwise
  // lcas two coincide and the third is this node, so lca(a,b)^lca(b,c)^lca(a,c)
  // also works; taking the deepest is clearer.
  int meetPoint(int a, int b, int c) {
    int x = lca(a, b), y = lca(b, c), z = lca(a, c);
    if (depth[y] > depth[x]) x = y;
    if (depth[z] > depth[x]) x = z;
    return x;
  }

  // lca of a and b if the tree were rerooted at r -- the same node as the meet
  // point of the triple.
  int rootedLca(int a, int b, int r) { return meetPoint(a, b, r); }

  // the node k steps along the path a -> b (k = 0 gives a)
  int jump(int a, int b, int k) {
    int l = lca(a, b), up_len = depth[a] - depth[l];
    if (k <= up_len) return ancestor(a, k);
    return ancestor(b, dist(a, b) - k);
  }
};


signed main() {
  //          0
  //     5 /  |2 \ 7
  //      3   4   1
  //   1 / \4      \ 3
  //    2   6       5
  //         \ 9
  //          7
  LCA t(8);
  t.addEdge(0, 3, 5), t.addEdge(0, 4, 2), t.addEdge(0, 1, 7);
  t.addEdge(3, 2, 1), t.addEdge(3, 6, 4);
  t.addEdge(1, 5, 3);
  t.addEdge(6, 7, 9);
  t.build(0);

  cout << t.lca(2, 7) << endl;              // 3
  cout << t.dist(5, 7) << endl;             // 5   path 5-1-0-3-6-7
  cout << t.wdist(5, 7) << endl;            // 28  = 3+7+5+4+9
  cout << t.jump(5, 7, 2) << endl;          // 0
  cout << t.pathAgg(2, 7) << endl;          // 9   edges 1, 4, 9
  cout << t.pathAgg(2, 4) << endl;          // 5   edges 1, 5, 2
  cout << t.pathAgg(3, 3) << endl;          // 0   ID
  cout << t.isAncestor(3, 7) << endl;       // 1
  cout << t.isAncestor(1, 7) << endl;       // 0
  cout << t.subtreeSize(3) << endl;         // 4   {3, 2, 6, 7}
  cout << t.onPath(5, 7, 0) << endl;        // 1
  cout << t.onPath(5, 7, 4) << endl;        // 0
  cout << t.meetPoint(2, 7, 5) << endl;     // 3
  cout << t.rootedLca(2, 4, 7) << endl;     // 3   rerooted at 7
  // farthest ancestor of 7 reachable over edges <= X
  cout << t.climbWhile(7, [](int m) { return m <= 9; }) << endl;  // 0
  cout << t.climbWhile(7, [](int m) { return m <= 8; }) << endl;  // 7
  cout << t.climbWhile(2, [](int m) { return m <= 4; }) << endl;  // 3
}

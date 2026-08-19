// ============================================================================
//  DAGs  --  topological sort, cycle detection, and DP on the resulting order
// ----------------------------------------------------------------------------
//  A topological order lists the vertices so that every edge u -> v has u
//  before v. It exists iff the directed graph is acyclic.
//
//  WHY YOU CARE: a DAG gives you a valid evaluation order for a recurrence, so
//  "DP on a DAG" and "longest path in a DAG" are the same thing. Any DP is a
//  DAG of subproblems; if your dependency graph has a cycle, the DP is ill-
//  defined and you need something else (Dijkstra, SCC condensation, ...).
//
//  CONTENTS
//    1.  topoDFS        reverse postorder, O(V+E)
//    2.  kahn           indegree/queue version -- also detects cycles, and
//                       gives the lexicographically smallest order with a heap
//    3.  hasCycle       directed cycle detection (3-colour DFS)
//    4.  findCycle      recover an actual directed cycle
//    5.  longestPath    longest path in a DAG (the lecture's Water Falls)
//    6.  countPaths     number of paths s -> t in a DAG
//    7.  NOTES
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

// ============================================================================
//  1. TOPOLOGICAL SORT BY DFS -- reverse postorder
// ----------------------------------------------------------------------------
//  A vertex is appended once all of its descendants are finished, so reversing
//  the finish order puts every vertex before everything it can reach.
//  Does NOT check for cycles: on a cyclic graph it silently returns garbage.
//  Pair it with hasCycle(), or use Kahn's below which detects cycles for free.
//
//  Recursion depth is O(V). At V = 1e5+ raise the stack or use Kahn's.
// ============================================================================
vector<int> topoDFS(const vector<vector<int>>& adj) {
	int n = (int)adj.size();
	vector<char> seen(n, 0);
	vector<int> order;
	// iterative DFS so we cannot blow the stack
	vector<pair<int, int>> st;                  // (node, next child index)
	st.reserve(n);                              // so push_back never reallocates
	for (int s = 0; s < n; s++) {
		if (seen[s]) continue;
		st.clear();
		st.push_back({s, 0});
		seen[s] = 1;
		while (!st.empty()) {
			auto& [u, i] = st.back();
			if (i < (int)adj[u].size()) {
				int v = adj[u][i++];
				if (!seen[v]) { seen[v] = 1; st.push_back({v, 0}); }
			} else {
				order.push_back(u);                // postorder
				st.pop_back();
			}
		}
	}
	reverse(order.begin(), order.end());
	return order;
}

// ============================================================================
//  2. KAHN'S ALGORITHM -- repeatedly remove a vertex of indegree 0
// ----------------------------------------------------------------------------
//  Returns an empty vector iff the graph has a cycle (fewer than n vertices
//  ever reach indegree 0). This is the version to write under exam pressure:
//  no recursion, and cycle detection comes free.
//
//  Swap the queue for a priority_queue<int, vector<int>, greater<int>> to get
//  the LEXICOGRAPHICALLY SMALLEST topological order (a common extra
//  requirement). Note that is not the same as "smallest at each greedy step"
//  for the reversed problem -- for lexicographically largest, run on the
//  reversed graph with a max-heap and reverse the result.
//
//  If at every step exactly one vertex has indegree 0, the order is UNIQUE
//  (equivalently, the DAG has a Hamiltonian path).
// ============================================================================
vector<int> kahn(const vector<vector<int>>& adj, bool lexSmallest = false) {
	int n = (int)adj.size();
	vector<int> indeg(n, 0);
	for (int u = 0; u < n; u++) for (int v : adj[u]) indeg[v]++;

	vector<int> order;
	if (lexSmallest) {
		priority_queue<int, vector<int>, greater<int>> pq;
		for (int i = 0; i < n; i++) if (!indeg[i]) pq.push(i);
		while (!pq.empty()) {
			int u = pq.top();
			pq.pop();
			order.push_back(u);
			for (int v : adj[u]) if (--indeg[v] == 0) pq.push(v);
		}
	} else {
		queue<int> q;
		for (int i = 0; i < n; i++) if (!indeg[i]) q.push(i);
		while (!q.empty()) {
			int u = q.front();
			q.pop();
			order.push_back(u);
			for (int v : adj[u]) if (--indeg[v] == 0) q.push(v);
		}
	}
	if ((int)order.size() != n) return {};     // cycle
	return order;
}

// ============================================================================
//  3. DIRECTED CYCLE DETECTION -- the 3-colour DFS
// ----------------------------------------------------------------------------
//  UNSEEN -> ACTIVE (on the recursion stack) -> COMPLETE (returned from).
//  An edge into an ACTIVE vertex is a back edge, hence a cycle.
//
//  The classic mistake is testing "have I seen v before" instead of "is v still
//  on the stack". On the diamond 1->2->3<-4<-1 you see 3 twice with no cycle.
//  (For UNDIRECTED graphs it really is that simple: any edge to a seen vertex
//  other than the one you came in on closes a cycle.)
// ============================================================================
bool hasCycle(const vector<vector<int>>& adj) {
	int n = (int)adj.size();
	const int UNSEEN = 0, ACTIVE = 1, COMPLETE = 2;
	vector<int> status(n, UNSEEN);
	// iterative to avoid stack overflow
	vector<pair<int, int>> st;
	st.reserve(n);                              // so push_back never reallocates
	for (int s = 0; s < n; s++) {
		if (status[s] != UNSEEN) continue;
		st.clear();
		st.push_back({s, 0});
		status[s] = ACTIVE;
		while (!st.empty()) {
			auto& [u, i] = st.back();
			if (i < (int)adj[u].size()) {
				int v = adj[u][i++];
				if (status[v] == ACTIVE) return true;
				if (status[v] == UNSEEN) { status[v] = ACTIVE; st.push_back({v, 0}); }
			} else {
				status[u] = COMPLETE;
				st.pop_back();
			}
		}
	}
	return false;
}

// Recover an actual directed cycle (node list, in order). Empty if acyclic.
vector<int> findCycle(const vector<vector<int>>& adj) {
	int n = (int)adj.size();
	vector<int> status(n, 0), par(n, -1);
	int cs = -1, ce = -1;
	function<bool(int)> dfs = [&](int u) -> bool {
		status[u] = 1;
		for (int v : adj[u]) {
			if (status[v] == 1) { cs = v; ce = u; return true; }
			if (status[v] == 0) { par[v] = u; if (dfs(v)) return true; }
		}
		status[u] = 2;
		return false;
	};
	for (int i = 0; i < n; i++)
		if (status[i] == 0 && dfs(i)) break;
	if (cs == -1) return {};
	vector<int> cyc;
	for (int v = ce; v != cs; v = par[v]) cyc.push_back(v);
	cyc.push_back(cs);
	reverse(cyc.begin(), cyc.end());
	return cyc;
}

// ============================================================================
//  5. LONGEST PATH IN A DAG  (== the lecture's Water Falls problem)
// ----------------------------------------------------------------------------
//  NP-hard on general graphs, trivial on a DAG: process vertices in topological
//  order and relax forwards. This is the prototypical DP.
//
//  Returns the maximum path weight starting anywhere. best[u] = longest path
//  starting at u. For "longest path FROM s", initialise best[] to -INF except
//  best[s] = 0 and relax in topological order.
//
//  The same loop with min instead of max gives shortest paths in a DAG -- and
//  unlike Dijkstra it is correct with NEGATIVE weights, in O(V+E).
// ============================================================================
ll longestPathDAG(const vector<vector<pair<int, ll>>>& adj) {
	int n = (int)adj.size();
	vector<vector<int>> plain(n);
	for (int u = 0; u < n; u++) for (auto [v, w] : adj[u]) plain[u].push_back(v);
	vector<int> order = kahn(plain);
	if (order.empty() && n > 0) return -1;         // not a DAG

	vector<ll> best(n, 0);
	// relax in REVERSE topological order: everything u reaches is already done
	for (int i = n - 1; i >= 0; i--) {
		int u = order[i];
		for (auto [v, w] : adj[u]) best[u] = max(best[u], best[v] + w);
	}
	return *max_element(best.begin(), best.end());
}

// ============================================================================
//  6. COUNT PATHS s -> t IN A DAG
// ----------------------------------------------------------------------------
//  Also the template for "number of ways to ...". Counts can be astronomically
//  large -- take them mod p unless the problem bounds them.
// ============================================================================
const ll MOD = 1000000007;

ll countPathsDAG(const vector<vector<int>>& adj, int s, int t) {
	int n = (int)adj.size();
	vector<int> order = kahn(adj);
	vector<ll> ways(n, 0);
	ways[s] = 1;
	for (int u : order)
		if (ways[u])
			for (int v : adj[u]) ways[v] = (ways[v] + ways[u]) % MOD;
	return ways[t];
}

// ============================================================================
//  7. NOTES
// ----------------------------------------------------------------------------
//  RECOGNISING A DAG IN DISGUISE. Many problems hand you a general graph plus a
//  rule that secretly forbids cycles:
//     * "water flows to strictly lower peaks"  -> orient by height
//     * "you may only move right or down"      -> grid DAG
//     * "task i must finish before task j"      -> precedence DAG
//     * costs/times strictly increase along any move
//  Orient the edges by that quantity and you have a DAG; then it is just DP.
//  If two endpoints have EQUAL height/value, do not add an edge in either
//  direction -- that is exactly what would create the cycle.
//
//  NOT A DAG? CONDENSE IT. Contract each strongly connected component to a
//  single node (see scc.cpp) and the condensation is always a DAG. Then DP over
//  the condensation, with each super-node carrying an aggregate of its members
//  ("total value collectable in this SCC" -- once you enter you can take it
//  all). This is the standard "longest path in a general directed graph with
//  node values" solution.
//
//  TOPOLOGICAL ORDER IS NOT UNIQUE unless every step has exactly one indegree-0
//  vertex. If the problem asks for a specific one, say so with Kahn + heap.
//
//  DFS TOPO VS KAHN: reverse postorder is shorter to write and is what you want
//  when you already have a DFS; Kahn detects cycles, avoids recursion, and
//  supports tie-breaking. Under exam pressure, write Kahn.
//
//  RELATED: functional graphs (every node has exactly one outgoing edge) are
//  NOT DAGs -- they are "rho" shaped, a tree hanging off a cycle. Find the
//  cycle with Floyd/tortoise-hare or by walking with a visit stamp, then handle
//  the trees separately. See Graph/successor.cpp for binary lifting on those.
// ============================================================================

int main() {
	// 0 -> 1 -> 3,  0 -> 2 -> 3
	int n = 4;
	vector<vector<int>> adj(n);
	adj[0] = {1, 2};
	adj[1] = {3};
	adj[2] = {3};

	cout << "topo (kahn):";
	for (int v : kahn(adj)) cout << " " << v;
	cout << "\n";
	cout << "topo (dfs) :";
	for (int v : topoDFS(adj)) cout << " " << v;
	cout << "\n";
	cout << "has cycle: " << (hasCycle(adj) ? "yes" : "no") << "\n";
	cout << "paths 0 -> 3: " << countPathsDAG(adj, 0, 3) << "\n";

	vector<vector<pair<int, ll>>> wadj(n);
	wadj[0] = {{1, 5}, {2, 2}};
	wadj[1] = {{3, 1}};
	wadj[2] = {{3, 9}};
	cout << "longest path: " << longestPathDAG(wadj) << "\n";

	adj[3].push_back(0);                              // now cyclic
	cout << "after adding 3->0, has cycle: " << (hasCycle(adj) ? "yes" : "no") << "\n";
	cout << "cycle:";
	for (int v : findCycle(adj)) cout << " " << v;
	cout << "\n";
	return 0;
}

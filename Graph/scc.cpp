// ============================================================================
//  TARJAN'S STRONGLY CONNECTED COMPONENTS  (directed graph, O(V + E))
// ----------------------------------------------------------------------------
//  Finds SCCs in one DFS pass. After gen(): comp[v] = id of v's component (the
//  id is an arbitrary representative node of that component); comps lists one
//  representative per SCC in REVERSE topological order of the condensation DAG
//  -- exactly the order you want for DP on the condensation and for 2-SAT.
//
//  API  (nodes 0-based)
//    SCC s; s.init(N);
//    s.ae(u, v);            // add directed edge u -> v
//    s.gen();               // compute components
//    // then read s.comp[v], or iterate s.comps
//    // two nodes u,v are in the same SCC iff comp[u] == comp[v].
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long
#define ll long long

const int INF = 9 * 1e18;
 
#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define trav(a, x) for(auto& a : x)
#define all(x) x.begin(), x.end()
#define print(arr) for (auto i : arr) cout << i << " "; cout << endl;
typedef vector<vector<int>> vvi;
typedef vector<pair<int, int>> vpi;
typedef vector<int> vi;
typedef pair<int, int> pii;

struct SCC {
	int N, ti = 0;
	vector<vector<int>> adj;
	vector<int> disc, low, comp, st, comps;
	void init(int _N) {
		N = _N;
		adj.resize(N), disc.resize(N), low.resize(N), comp = vector<int>(N, -1);
	}
	void ae(int x, int y) { adj[x].push_back(y); }
	void dfs(int x) {
		disc[x] = low[x] = ++ti;
		st.push_back(x);  // disc[y] != 0 -> in stack
		for (int y : adj[x])
			if (comp[y] == -1) low[x] = min(low[x], disc[y] ?: (dfs(y), low[y]));
		if (low[x] == disc[x]) {  // make new SCC, pop off stack until you find x
			comps.push_back(x);
			for (int y = -1; y != x;) comp[y = st.back()] = x, st.pop_back();
		}
	}
	void gen() {
		for (int i = 0; i < N; i++)
			if (!disc[i]) dfs(i);
		reverse(begin(comps), end(comps));
	}
};
signed main() {
    SCC scc;
    int numNodes = 5;  // Number of nodes in the graph
    scc.init(numNodes);

    // Add edges to the graph
    scc.ae(0, 1);
    scc.ae(1, 2);
    scc.ae(2, 0);
    scc.ae(1, 3);
    scc.ae(3, 4);

    // Generate SCCs
    scc.gen();

    // Print the SCCs
    for (int i = 0; i < numNodes; i++) {
        cout << "Node " << i << " is in component " << scc.comp[i] << endl;
    }

    return 0;
}
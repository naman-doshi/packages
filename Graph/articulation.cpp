// ============================================================================
//  DFS TREE STRUCTURE  --  articulation points, bridges, biconnected and
//                          2-edge-connected components, bridge tree
// ----------------------------------------------------------------------------
//  All of it is ONE DFS, O(V + E), and all of it rests on the same two facts
//  about an UNDIRECTED graph:
//
//    (a) Every edge is either a TREE edge of the DFS tree or a BACK EDGE.
//        There are no "cross edges" -- if there were, the DFS would have
//        descended into that vertex instead.
//    (b) A back edge always points from a vertex to one of its ANCESTORS.
//
//  So: low[u] = the smallest discovery time reachable from u's subtree using
//  tree edges downwards plus at most one back edge. Then
//        edge (p, u) is a BRIDGE          iff  low[u] >  disc[p]
//        p is an ARTICULATION POINT       iff  low[u] >= disc[p] for some child
//                                             (and the root iff it has >= 2
//                                              children in the DFS tree)
//  The only difference between the two is > versus >=: an edge is a bridge when
//  the subtree cannot reach p AT ALL, while p is a cut vertex when the subtree
//  cannot reach ABOVE p.
//
//  CONTENTS
//    1.  CutPoints    articulation points + bridges together
//    2.  bridgeTree   contract 2-edge-connected components -> a TREE
//    3.  NOTES
//
//  (Bridges alone, with the parallel-edge handling spelled out, are also in
//   Graph/bridges.cpp -- this file is the superset.)
// ============================================================================
#include <bits/stdc++.h>
using namespace std;

// ============================================================================
//  1. ARTICULATION POINTS AND BRIDGES
// ----------------------------------------------------------------------------
//  PARALLEL EDGES: we skip the edge INDEX we arrived on, not the parent NODE.
//  If you skip the parent node instead, a doubled edge u==v is wrongly reported
//  as a bridge. Self-loops are never bridges and never create cut vertices.
//  Disconnected graphs are handled (we start a DFS from every unseen vertex).
//
//  API  (nodes 0-based)
//    CutPoints c; c.init(N);
//    c.ae(u, v);
//    c.gen();
//    c.isCut[v]      // true iff v is an articulation point
//    c.bridges       // vector of (u, v)
// ============================================================================
struct CutPoints {
	int N, timer = 0;
	vector<vector<pair<int, int>>> adj;    // (neighbour, edge id)
	vector<int> disc, low;
	vector<char> isCut;
	vector<pair<int, int>> bridges;
	int edgeCount = 0;

	void init(int _N) {
		N = _N;
		adj.assign(N, {});
		disc.assign(N, -1);
		low.assign(N, -1);
		isCut.assign(N, 0);
		bridges.clear();
		timer = edgeCount = 0;
	}
	void ae(int u, int v) {
		adj[u].push_back({v, edgeCount});
		adj[v].push_back({u, edgeCount});
		edgeCount++;
	}

	void dfs(int u, int inEdge) {
		disc[u] = low[u] = timer++;
		int children = 0;
		for (auto [v, id] : adj[u]) {
			if (id == inEdge) continue;                 // the edge we came in on
			if (disc[v] != -1) {                        // back edge
				low[u] = min(low[u], disc[v]);
			} else {                                     // tree edge
				children++;
				dfs(v, id);
				low[u] = min(low[u], low[v]);
				if (low[v] > disc[u]) bridges.push_back({u, v});
				if (low[v] >= disc[u] && inEdge != -1) isCut[u] = 1;
			}
		}
		if (inEdge == -1 && children > 1) isCut[u] = 1;  // root rule
	}

	void gen() {
		for (int i = 0; i < N; i++) if (disc[i] == -1) dfs(i, -1);
	}
};

// ============================================================================
//  2. BRIDGE TREE  (2-edge-connected component contraction)
// ----------------------------------------------------------------------------
//  Delete every bridge; each remaining connected piece is a 2-edge-connected
//  component (a maximal subgraph with no bridge -- between any two of its
//  vertices there are two edge-disjoint paths). Contract each to one node and
//  reconnect them with the bridges: the result is always a TREE (or forest).
//
//  This is the standard move for "any query about edges that must be crossed":
//  once it is a tree you have LCA, path queries, tree DP, diameter, ...
//    * "minimum number of edges whose removal disconnects u from v" on an
//      undirected graph = number of bridges on the u-v path in the bridge tree.
//    * "add the fewest edges to make the graph 2-edge-connected" =
//      ceil(number of leaves of the bridge tree / 2).
//
//  Returns comp[] (which 2ecc each vertex is in) and the tree's edge list.
// ============================================================================
struct BridgeTree {
	vector<int> comp;                 // comp[v] = id of v's 2-edge-connected comp
	int numComps = 0;
	vector<pair<int, int>> treeEdges;

	void build(CutPoints& c) {
		int n = c.N;
		set<pair<int, int>> isBridge;
		for (auto [u, v] : c.bridges) {
			isBridge.insert({u, v});
			isBridge.insert({v, u});
		}
		comp.assign(n, -1);
		numComps = 0;
		// flood fill without crossing bridges
		for (int s = 0; s < n; s++) {
			if (comp[s] != -1) continue;
			int id = numComps++;
			vector<int> st = {s};
			comp[s] = id;
			while (!st.empty()) {
				int u = st.back();
				st.pop_back();
				for (auto [v, eid] : c.adj[u]) {
					if (isBridge.count({u, v})) continue;
					if (comp[v] == -1) { comp[v] = id; st.push_back(v); }
				}
			}
		}
		treeEdges.clear();
		for (auto [u, v] : c.bridges) treeEdges.push_back({comp[u], comp[v]});
	}
};

// ============================================================================
//  3. NOTES
// ----------------------------------------------------------------------------
//  DEFINITIONS, so you answer the question actually asked:
//    * BRIDGE          an EDGE whose removal increases the component count.
//    * ARTICULATION    a VERTEX whose removal increases the component count.
//    * 2-EDGE-CONNECTED  no bridges. Contract -> bridge tree.
//    * 2-VERTEX-CONNECTED (biconnected) no articulation points. Contract ->
//      block-cut tree, whose nodes are the biconnected components plus the cut
//      vertices, and whose edges join a cut vertex to each block containing it.
//    A graph can have no bridges but still have a cut vertex (two triangles
//    sharing one vertex), so do not conflate the two.
//
//  BICOMPONENTS: push edges onto a stack as you traverse them; when a child v
//  of u satisfies low[v] >= disc[u], pop edges down to (u,v) -- that popped set
//  is one biconnected component. (Blocks partition the EDGES, not the vertices;
//  a cut vertex belongs to several blocks.)
//
//  DIRECTED GRAPHS: none of this applies. The DFS tree of a directed graph has
//  cross and forward edges too, and "strongly connected" is the right notion --
//  use scc.cpp (Tarjan) which is the same low-link idea adapted.
//
//  RECURSION DEPTH is O(V). At V >= ~1e5 with a path-shaped graph you will
//  stack overflow; either raise the stack or rewrite the DFS iteratively.
//
//  TYPICAL PHRASINGS THAT MEAN "BRIDGES"
//    * "roads whose closure would cut the country in two"
//    * "critical connections in a network"
//    * "edges lying on no cycle"      (a bridge is exactly such an edge)
//  AND THAT MEAN "ARTICULATION POINTS"
//    * "servers whose failure disconnects the network"
//    * "after removing city i, how many components remain"
// ============================================================================

int main() {
	// 0-1-2 triangle, bridge 2--3, then 3-4-5 triangle, bridge 5--6
	CutPoints c;
	c.init(7);
	c.ae(0, 1); c.ae(1, 2); c.ae(2, 0);
	c.ae(2, 3);
	c.ae(3, 4); c.ae(4, 5); c.ae(5, 3);
	c.ae(5, 6);
	c.gen();

	cout << "bridges:";
	for (auto [u, v] : c.bridges) cout << " (" << u << "," << v << ")";
	cout << "\n";
	cout << "articulation points:";
	for (int i = 0; i < c.N; i++) if (c.isCut[i]) cout << " " << i;
	cout << "\n";

	BridgeTree bt;
	bt.build(c);
	cout << "2-edge-connected components: " << bt.numComps << "\n";
	cout << "comp[] =";
	for (int i = 0; i < c.N; i++) cout << " " << bt.comp[i];
	cout << "\n";
	cout << "bridge tree edges:";
	for (auto [u, v] : bt.treeEdges) cout << " (" << u << "," << v << ")";
	cout << "\n";
	return 0;
}

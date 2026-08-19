#include <bits/stdc++.h>
using namespace std;
#define int long long

vector<vector<int>> graph;
vector<int> indegree;

// make sure nodes are 0 indexed
vector<int> topsort() {
  queue<int> q;
  for (int i = 0; i < indegree.size(); i++) {
    if (indegree[i] == 0) { q.push(i); }
  }
  vector<int> order;
  while (!q.empty()) {
    int curr = q.front();
    q.pop();
    order.push_back(curr);
    for (int next : graph[curr]) {
      indegree[next]--;
      if (indegree[next] == 0) { q.push(next); }
    }
  }
  if (order.size() != indegree.size()) {
    return vector<int>();
  }
  return order;
}

// find cycle starting at node n in directed graph
bool visited[(int)1e5 + 5], on_stack[(int)1e5 + 5];
vector<int> adj[(int)1e5 + 5];
vector<int> cycle;
int N, M;
bool dfs(int n) {
	visited[n] = on_stack[n] = true;
	for (int u : adj[n]) {
		if (on_stack[u]) {
			cycle.push_back(n);  // start cycle
			on_stack[n] = on_stack[u] = false;
			return true;
		} else if (!visited[u]) {
			if (dfs(u)) {  // continue cycle
				if (on_stack[n]) {
					cycle.push_back(n);
					on_stack[n] = false;
					return true;
				} else {  // found u again
					cycle.push_back(n);
					return false;
				}
			}
			if (!cycle.empty())  // finished with cycle
				return false;
		}
	}
	on_stack[n] = false;
	return false;
}
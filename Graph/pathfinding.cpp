#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <queue>
#include <deque>
#include <bitset>
#include <iterator>
#include <list>
#include <stack>
#include <map>
#include <set>
#include <functional>
#include <numeric>
#include <utility>
#include <limits>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cmath>

using namespace std;

const int INF = 2147483647;
const double PI =  3.1415926535897932384626433832795;

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define trav(a, x) for(auto& a : x)
#define all(x) x.begin(), x.end()
#define len(x) x.size()
// #define count(arr, elem) count(arr.begin(), arr.end(), elem)
#define in(arr, elem) find(arr.begin(), arr.end(), elem) != arr.end()
#define notin(arr, elem) find(arr.begin(), arr.end(), elem) == arr.end()
#define sum(arr) accumulate(arr.begin(), arr.end(), 0)
#define printVector(arr) for (auto elem : arr) cout << elem << " "; cout << endl;

typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;
typedef vector<string> vs;
typedef vector<bool> vb;
typedef vector<vector<int>> vvi;
typedef vector<vector<string>> vvs;
typedef vector<vector<char>> vvc;
typedef map<string, int> msi;
typedef map<string, string> mss;
typedef vector<pair<int, int>> vpi;

map<int, vi> graph;

vi shortestPath(int N, int origin, int dest) {
  vi dist(N + 1, INF), parent(N + 1);
  queue<int> q;
	dist[1] = 0;
	q.push(1);
	while (!q.empty()) {
		int x = q.front();
		q.pop();
		for (int t : graph[x])
			if (dist[t] == INF) {
				dist[t] = dist[x] + 1;
				parent[t] = x;
				q.push(t);
			}
	}
	if (dist[N] == INF) return vi();
	else {
		vi v{N};
		while (v.back() != 1) v.push_back(parent[v.back()]);
		reverse(begin(v), end(v));
		return v;
	}
}

int dijkstra(int start, int end, vector<vector<int>> weights) {
	vi dist(len(weights), INF);
	dist[start] = 0;
	priority_queue<pii, vector<pii>, greater<pii>> pq;
	pq.push({0, start});
	while (!pq.empty()) {
		pii curr = pq.top();
		pq.pop();
		int node = curr.second;
		int d = curr.first;
		if (d > dist[node]) continue;
		for (int i = 0; i < len(weights[node]); i++) {
			int neighbor = i;
			int weight = weights[node][i];
			if (dist[node] + weight < dist[neighbor]) {
				dist[neighbor] = dist[node] + weight;
				pq.push({dist[neighbor], neighbor});
			}
		}
	}
	return dist[end];
}

vector<vector<pair<int, int>>> adj;

// bellman
bool bellmanFord(int s, vector<int>& d) {
		int n = adj.size();
		d.assign(n, INF);
		d[s] = 0;
		for (int i = 0; i < n - 1; i++) {
				bool any = false;
				for (int u = 0; u < n; u++) {
						for (auto edge : adj[u]) {
								int v = edge.first;
								int len = edge.second;
								if (d[u] < INF && d[u] + len < d[v]) {
										d[v] = d[u] + len;
										any = true;
								}
						}
				}
				if (!any) break;
		}
		for (int u = 0; u < n; u++) {
				for (auto edge : adj[u]) {
						int v = edge.first;
						int len = edge.second;
						if (d[u] < INF && d[u] + len < d[v]) {
								return false;  // negative cycle
						}
				}
		}
		return true;
}

// spfa
bool spfa(int s, vector<int>& d) {
    int n = adj.size();
    d.assign(n, INF);
    vector<int> cnt(n, 0);
    vector<bool> inqueue(n, false);
    queue<int> q;

    d[s] = 0;
    q.push(s);
    inqueue[s] = true;
    while (!q.empty()) {
        int v = q.front();
        q.pop();
        inqueue[v] = false;

        for (auto edge : adj[v]) {
            int to = edge.first;
            int len = edge.second;

            if (d[v] + len < d[to]) {
                d[to] = d[v] + len;
                if (!inqueue[to]) {
                    q.push(to);
                    inqueue[to] = true;
                    cnt[to]++;
                    if (cnt[to] > n)
                        return false;  // negative cycle
                }
            }
        }
    }
    return true;
}


// floyd warshall
// good for all pairs shortest path instead of sssp
void floydWarshall(vector<vector<int>>& dist) {
	int n = dist.size();
	rep(i, 0, n) dist[i][i] = 0;
	for (int k = 0; k < n; k++) {
		for (int i = 0; i < n; i++) { 
			for (int j = 0; j < n; j++) {
				if (dist[i][k] < INF && dist[k][j] < INF) {
					dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
				}
			}
		}
	}
}


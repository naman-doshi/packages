#include <iostream>
#include <vector>

using namespace std;

const int maxn = 200010;

int n;
vector<int> adj[maxn];
int subtree_size[maxn];

struct DSU {
  vector<ll> comp;
  DSU (ll elements) : comp(elements, -1) {}
  bool Union(ll a, ll b) {
    a = repr(a), b = repr(b);
    if (a == b) return true;
    if (-comp[a] < -comp[b]) swap(a, b);
    comp[a] += comp[b];
    comp[b] = a;
    return false;
  }
  bool is_repr(ll x) {return comp[x] < 0; }
  ll repr(ll x) {
    if (comp[x] < 0) return x;
    ll par = comp[x];
    comp[x] = repr(par);
    return comp[x];
  }
  bool connected(ll a, ll b) {return repr(a) == repr(b); }
};

int get_subtree_size(int node, int parent = -1) {
	int &res = subtree_size[node];
	res = 1;
	for (int i : adj[node]) {
		if (i == parent) { continue; }
		res += get_subtree_size(i, node);
	}
	return res;
}

int get_centroid(int node, int parent = -1) {
	for (int i : adj[node]) {
		if (i == parent) { continue; }

		if (subtree_size[i] * 2 > n) { return get_centroid(i, node); }
	}
	return node;
}

vector<vector<int>> adj;
vector<bool> is_removed;
vector<int> subtree_size;

/** DFS to calculate the size of the subtree rooted at `node` */
int get_subtree_size(int node, int parent = -1) {
	subtree_size[node] = 1;
	for (int child : adj[node]) {
		if (child == parent || is_removed[child]) { continue; }
		subtree_size[node] += get_subtree_size(child, node);
	}
	return subtree_size[node];
}

/**
 * Returns a centroid (a tree may have two centroids) of the subtree
 * containing node `node` after node removals
 * @param node current node
 * @param tree_size size of current subtree after node removals
 * @param parent parent of u
 * @return first centroid found
 */
int get_centroid(int node, int tree_size, int parent = -1) {
	for (int child : adj[node]) {
		if (child == parent || is_removed[child]) { continue; }
		if (subtree_size[child] * 2 > tree_size) {
			return get_centroid(child, tree_size, node);
		}
	}
	return node;
}

/** Build up the centroid decomposition recursively */
void build_centroid_decomp(int node = 0) {
	int centroid = get_centroid(node, get_subtree_size(node));

	// do something

	is_removed[centroid] = true;

	for (int child : adj[centroid]) {
		if (is_removed[child]) { continue; }
		build_centroid_decomp(child);
	}
}

int main() {
	cin >> n;
	for (int i = 0; i < n - 1; i++) {
		int a, b;
		cin >> a >> b;
		a--;
		b--;
		adj[a].push_back(b);
		adj[b].push_back(a);
	}

	get_subtree_size(0);
	cout << get_centroid(0) + 1 << endl;
}
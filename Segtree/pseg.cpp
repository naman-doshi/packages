// ============================================================================
//  PERSISTENT SEGMENT TREE  -- keep every past version, query any of them
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    An ordinary segment tree, except each update COPIES only the O(log n) nodes
//    on the path it changes and returns a NEW root, sharing everything else with
//    the previous version. So all historical versions stay queryable in
//    O(log n), at only O(log n) extra memory per update.
//    Killer apps: k-th smallest in a subarray (build one version per prefix over
//    the value axis, then "subtract" version l-1 from version r); range "how many
//    values <= x"; querying past states; problems on a tree of versions.
//
//  THIS FILE = range-sum demo. roots[k] is the root of version k. main() shows:
//    (1) point-set on version k, (2) range-sum on version k, (3) clone a version.
//    Positions are 1-based over [1, n].
//
//  CORE API  (reuse these three; ignore main())
//    roots[0] = build();                         // version 0 from a[1..n]
//    roots[v2] = update(roots[v1], newVal, pos); // v2 = v1 but with a[pos]=newVal
//    ll s = query(roots[v], l, r);               // range sum on version v
//
//  COMPLEXITY  O(log n) per op; O(log n) new nodes per update; O(n + q log n) mem.
//
//  PITFALLS
//    * `new` is never freed -- fine for one CP run, but size arrays for the TOTAL
//      node count (~ 2n + q log n), not just n.
//    * `n` is a global read by default args -- set it before calling build().
//    * Positions here are 1-based; keep updates and queries consistent.
//    * roots[v] = update(roots[v], ...) OVERWRITES version v (in-place edit).
//      To BRANCH history, write the new root to a fresh index instead.
// ============================================================================
#include <bits/stdc++.h>
typedef long long ll;
using namespace std;

struct Node {
	ll val;
	Node *l, *r;

	Node(ll x) : val(x), l(nullptr), r(nullptr) {}
	Node(Node *ll, Node *rr) {
		l = ll, r = rr;
		val = 0;
		if (l) val += l->val;
		if (r) val += r->val;
	}
	Node(Node *cp) : val(cp->val), l(cp->l), r(cp->r) {}
};

int n, cnt = 1;
ll a[200001];
Node *roots[200001];

Node *build(int l = 1, int r = n) {
	if (l == r) return new Node(a[l]);
	int mid = (l + r) / 2;
	return new Node(build(l, mid), build(mid + 1, r));
}

Node *update(Node *node, int val, int pos, int l = 1, int r = n) {
	if (l == r) return new Node(val);
	int mid = (l + r) / 2;
	if (pos > mid)
		return new Node(node->l, update(node->r, val, pos, mid + 1, r));
	else return new Node(update(node->l, val, pos, l, mid), node->r);
}

ll query(Node *node, int a, int b, int l = 1, int r = n) {
	if (l > b || r < a) return 0;
	if (l >= a && r <= b) return node->val;
	int mid = (l + r) / 2;
	return query(node->l, a, b, l, mid) + query(node->r, a, b, mid + 1, r);
}

int main() {
	ios_base::sync_with_stdio(0);
	cin.tie(0);
	int q;
	cin >> n >> q;
	for (int i = 1; i <= n; i++) cin >> a[i];
	roots[cnt++] = build();

	while (q--) {
		int t;
		cin >> t;
		if (t == 1) {
			int k, i, x;
			cin >> k >> i >> x;
			roots[k] = update(roots[k], x, i);
		} else if (t == 2) {
			int k, l, r;
			cin >> k >> l >> r;
			cout << query(roots[k], l, r) << '\n';
		} else {
			int k;
			cin >> k;
			roots[cnt++] = new Node(roots[k]);
		}
	}
	return 0;
}
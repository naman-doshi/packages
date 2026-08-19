// ============================================================================
//  KD-TREE   -- REQUIRES point.cpp
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Static set of points, then many queries of the shape "what is near here":
//      * nearest neighbour of a query point            (the usual reason)
//      * k nearest neighbours
//      * how many points inside this axis-aligned box
//      * how many points within distance r of here
//    Build once, query many times. If you only need the single closest PAIR of
//    the set, closestPair.cpp is faster and simpler; if you only need box
//    COUNTS and the points are static, a merge-sort tree / offline BIT is
//    faster. The KD-tree wins when the queries are online and geometric.
//
//  API
//    KDTree T(pts);                       // pts: vector<P>, copied and reordered
//    auto [d2, q] = T.nearest(p);         // squared distance + the point
//    auto [d2, q] = T.nearest(p, true);   // ...ignoring a point equal to p
//    vector<pair<ll,P>> v = T.kNearest(p, k);      // sorted, closest first
//    ll c = T.countInRect(P(x0,y0), P(x1,y1));     // inclusive box
//    ll c = T.countInCircle(p, r);                 // r is a RADIUS, not r^2
//    vP v  = T.reportInRect(P(x0,y0), P(x1,y1));   // the points themselves
//
//  COMPLEXITY  build O(n log n). nearest ~O(log n) on random data, O(n) worst
//    case. countInRect O(sqrt n) for a balanced tree. kNearest ~O(k + log n).
//
//  PITFALLS
//    * All distances here are SQUARED and exact (long long). countInCircle
//      takes a plain radius and squares it internally -- pass r, not r*r.
//    * Coordinates up to 1e9 make a squared distance up to 8e18: it fits, but
//      do not add two of them together.
//    * nearest(p) with p in the tree returns p itself at distance 0. Pass
//      excludeSelf = true to skip points EQUAL to p (not "the same index" --
//      duplicates of p are skipped too).
//    * The tree is static. To delete points, rebuild, or mark and filter.
//    * Worst case really is O(n) per query (all points on a line, query far
//      away). If an adversary picks the points, randomly rotate them first.
// ============================================================================

struct KDTree {
	struct Node {
		ll x0, x1, y0, y1;                  // bounding box of this subtree
		int lo, hi;                         // pts[lo, hi)
		int l = -1, r = -1;
	};
	vP pts;
	vector<Node> nd;

	KDTree(vP p) : pts(std::move(p)) {
		if (!pts.empty()) build(0, (int)pts.size());
	}
	int build(int lo, int hi) {
		int id = (int)nd.size();
		nd.push_back(Node());
		Node cur;
		cur.lo = lo; cur.hi = hi;
		cur.x0 = cur.y0 = LLONG_MAX;
		cur.x1 = cur.y1 = LLONG_MIN;
		for (int i = lo; i < hi; i++) {
			cur.x0 = min(cur.x0, pts[i].x); cur.x1 = max(cur.x1, pts[i].x);
			cur.y0 = min(cur.y0, pts[i].y); cur.y1 = max(cur.y1, pts[i].y);
		}
		nd[id] = cur;
		if (hi - lo > 1) {
			bool byX = (cur.x1 - cur.x0) >= (cur.y1 - cur.y0);   // split the wider axis
			int mid = (lo + hi) / 2;
			nth_element(pts.begin() + lo, pts.begin() + mid, pts.begin() + hi,
			            [byX](const P &a, const P &b) { return byX ? a.x < b.x : a.y < b.y; });
			int L = build(lo, mid), R = build(mid, hi);
			nd[id].l = L; nd[id].r = R;
		}
		return id;
	}
	// squared distance from p to the node's bounding box (0 if inside)
	ll boxDist(int id, P p) const {
		ll dx = max({(ll)0, nd[id].x0 - p.x, p.x - nd[id].x1});
		ll dy = max({(ll)0, nd[id].y0 - p.y, p.y - nd[id].y1});
		return dx * dx + dy * dy;
	}

	// ------------------------------------------------------------- nearest
	void nearRec(int id, P p, bool excl, pair<ll, P> &best) const {
		if (nd[id].l < 0) {
			for (int i = nd[id].lo; i < nd[id].hi; i++) {
				if (excl && pts[i] == p) continue;
				ll d = dist2(pts[i], p);
				if (d < best.first) best = {d, pts[i]};
			}
			return;
		}
		int a = nd[id].l, b = nd[id].r;
		ll da = boxDist(a, p), db = boxDist(b, p);
		if (da > db) swap(a, b), swap(da, db);
		if (da < best.first) nearRec(a, p, excl, best);
		if (db < best.first) nearRec(b, p, excl, best);
	}
	pair<ll, P> nearest(P p, bool excludeSelf = false) const {
		pair<ll, P> best{LLONG_MAX, P()};
		if (!nd.empty()) nearRec(0, p, excludeSelf, best);
		return best;
	}

	// ------------------------------------------------------------ k nearest
	void kRec(int id, P p, int k, priority_queue<pair<ll, P>> &pq) const {
		if (nd[id].l < 0) {
			for (int i = nd[id].lo; i < nd[id].hi; i++) {
				pq.push({dist2(pts[i], p), pts[i]});
				if ((int)pq.size() > k) pq.pop();
			}
			return;
		}
		int a = nd[id].l, b = nd[id].r;
		ll da = boxDist(a, p), db = boxDist(b, p);
		if (da > db) swap(a, b), swap(da, db);
		auto worst = [&] { return (int)pq.size() < k ? LLONG_MAX : pq.top().first; };
		if (da <= worst()) kRec(a, p, k, pq);
		if (db <= worst()) kRec(b, p, k, pq);
	}
	vector<pair<ll, P>> kNearest(P p, int k) const {
		priority_queue<pair<ll, P>> pq;
		if (!nd.empty()) kRec(0, p, k, pq);
		vector<pair<ll, P>> res;
		while (!pq.empty()) res.push_back(pq.top()), pq.pop();
		reverse(res.begin(), res.end());
		return res;
	}

	// ----------------------------------------------------------- range query
	ll rectRec(int id, P lo, P hi) const {
		const Node &v = nd[id];
		if (v.x1 < lo.x || v.x0 > hi.x || v.y1 < lo.y || v.y0 > hi.y) return 0;
		if (v.x0 >= lo.x && v.x1 <= hi.x && v.y0 >= lo.y && v.y1 <= hi.y)
			return v.hi - v.lo;                          // fully inside: O(1)
		if (v.l < 0) {
			ll c = 0;
			for (int i = v.lo; i < v.hi; i++)
				c += (pts[i].x >= lo.x && pts[i].x <= hi.x && pts[i].y >= lo.y &&
				      pts[i].y <= hi.y);
			return c;
		}
		return rectRec(v.l, lo, hi) + rectRec(v.r, lo, hi);
	}
	ll countInRect(P lo, P hi) const { return nd.empty() ? 0 : rectRec(0, lo, hi); }

	void reportRec(int id, P lo, P hi, vP &out) const {
		const Node &v = nd[id];
		if (v.x1 < lo.x || v.x0 > hi.x || v.y1 < lo.y || v.y0 > hi.y) return;
		if (v.l < 0) {
			for (int i = v.lo; i < v.hi; i++)
				if (pts[i].x >= lo.x && pts[i].x <= hi.x && pts[i].y >= lo.y && pts[i].y <= hi.y)
					out.push_back(pts[i]);
			return;
		}
		reportRec(v.l, lo, hi, out);
		reportRec(v.r, lo, hi, out);
	}
	vP reportInRect(P lo, P hi) const {
		vP out;
		if (!nd.empty()) reportRec(0, lo, hi, out);
		return out;
	}

	ll circRec(int id, P c, ll r2) const {
		const Node &v = nd[id];
		if (boxDist(id, c) > r2) return 0;
		// the whole box is inside if its farthest corner is within r
		ll fx = max(llabs(v.x0 - c.x), llabs(v.x1 - c.x));
		ll fy = max(llabs(v.y0 - c.y), llabs(v.y1 - c.y));
		if (fx * fx + fy * fy <= r2) return v.hi - v.lo;
		if (v.l < 0) {
			ll cnt = 0;
			for (int i = v.lo; i < v.hi; i++) cnt += dist2(pts[i], c) <= r2;
			return cnt;
		}
		return circRec(v.l, c, r2) + circRec(v.r, c, r2);
	}
	ll countInCircle(P c, ll r) const { return nd.empty() ? 0 : circRec(0, c, r * r); }
};

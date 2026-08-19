// ============================================================================
//  SWEEP LINE   -- REQUIRES point.cpp, lines.cpp
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    The four sweeps that actually show up in contests:
//      rectUnionArea       area covered by >= 1 of n axis-aligned rectangles
//      rectUnionPerimeter  the outline length of that same union
//      skyline             the building outline, as (x, height) corners
//      anySegIntersect     do ANY two of n segments cross?  (Shamos-Hoey)
//      countHVIntersect    # crossing pairs of horizontal / vertical segments
//      maxOverlap          most intervals covering a single point
//
//  API
//    rects are {x1, y1, x2, y2} with x1 < x2 and y1 < y2 (half-open boxes)
//    ll a = rectUnionArea(rects);
//    ll p = rectUnionPerimeter(rects);
//    vector<pair<ll,ll>> sk = skyline(buildings);   // {left, right, height}
//    pair<int,int> hit = anySegIntersect(segs);     // {-1,-1} if none
//    ll c = countHVIntersect(segs);
//    int k = maxOverlap(intervals);
//
//  COMPLEXITY  all O(n log n).
//
//  PITFALLS
//    * rectUnionArea uses HALF-OPEN boxes [x1,x2) x [y1,y2). Two rectangles
//      that merely share an edge contribute no overlap -- which is what you
//      want for area. Degenerate rectangles (zero width or height) are skipped.
//    * The area can be huge: coordinates up to 1e9 give up to 4e18. It fits in
//      ll, but do not add anything else to it.
//    * rectUnionPerimeter counts the outline of the union, so internal shared
//      edges are not counted. Rectangles that touch only at a corner are fine.
//    * anySegIntersect returns SOME intersecting pair, not all of them, and not
//      necessarily the leftmost. Counting all crossings needs O(n log n + k)
//      Bentley-Ottmann, which is almost never the intended solution -- if the
//      problem wants a count, it is usually axis-parallel (countHVIntersect).
//    * anySegIntersect keeps a set ordered by the y at the current sweep x. Do
//      not insert into it from outside the sweep.
// ============================================================================

// segment tree over compressed y that tracks "total length covered at least once"
struct CoverSeg {
	int n;
	vector<ll> ys;
	vector<int> cnt;
	vector<ll> cov;
	CoverSeg(vector<ll> y) : ys(y) {
		n = (int)ys.size() - 1;
		if (n < 1) n = 1;
		cnt.assign(4 * n, 0);
		cov.assign(4 * n, 0);
	}
	void upd(int nd, int l, int r, int ql, int qr, int v) {
		if (qr <= l || r <= ql) return;
		if (ql <= l && r <= qr) cnt[nd] += v;
		else {
			int m = (l + r) / 2;
			upd(2 * nd, l, m, ql, qr, v);
			upd(2 * nd + 1, m, r, ql, qr, v);
		}
		if (cnt[nd] > 0) cov[nd] = ys[r] - ys[l];
		else if (r - l == 1) cov[nd] = 0;
		else cov[nd] = cov[2 * nd] + cov[2 * nd + 1];
	}
	void add(int l, int r, int v) { upd(1, 0, n, l, r, v); }
	ll covered() { return cov[1]; }
};

// total length of the y-projection swept, per unit x -- shared by area/perimeter
ll rectUnionArea(vector<array<ll, 4>> rs) {
	vector<ll> ys;
	for (auto &r : rs)
		if (r[0] < r[2] && r[1] < r[3]) { ys.push_back(r[1]); ys.push_back(r[3]); }
	sort(ys.begin(), ys.end());
	ys.erase(unique(ys.begin(), ys.end()), ys.end());
	if (ys.size() < 2) return 0;
	auto yi = [&](ll v) { return (int)(lower_bound(ys.begin(), ys.end(), v) - ys.begin()); };
	vector<array<ll, 4>> ev;                       // {x, delta, ylo, yhi}
	for (auto &r : rs) {
		if (r[0] >= r[2] || r[1] >= r[3]) continue;
		ev.push_back({r[0], 1, (ll)yi(r[1]), (ll)yi(r[3])});
		ev.push_back({r[2], -1, (ll)yi(r[1]), (ll)yi(r[3])});
	}
	if (ev.empty()) return 0;
	sort(ev.begin(), ev.end());
	CoverSeg T(ys);
	ll area = 0, prev = ev[0][0];
	for (auto &e : ev) {
		area += T.covered() * (e[0] - prev);
		prev = e[0];
		T.add((int)e[2], (int)e[3], (int)e[1]);
	}
	return area;
}

// sum over the sweep of |change in covered length| -- one half of the perimeter
ll _perimSide(vector<array<ll, 4>> rs) {
	vector<ll> ys;
	for (auto &r : rs)
		if (r[0] < r[2] && r[1] < r[3]) { ys.push_back(r[1]); ys.push_back(r[3]); }
	sort(ys.begin(), ys.end());
	ys.erase(unique(ys.begin(), ys.end()), ys.end());
	if (ys.size() < 2) return 0;
	auto yi = [&](ll v) { return (int)(lower_bound(ys.begin(), ys.end(), v) - ys.begin()); };
	vector<array<ll, 4>> ev;
	for (auto &r : rs) {
		if (r[0] >= r[2] || r[1] >= r[3]) continue;
		ev.push_back({r[0], 1, (ll)yi(r[1]), (ll)yi(r[3])});
		ev.push_back({r[2], -1, (ll)yi(r[1]), (ll)yi(r[3])});
	}
	if (ev.empty()) return 0;
	// ADDs before REMOVEs at the same x, and one event at a time. Both matter:
	// grouping a whole x together misses the case where the covered set changes
	// shape but not length (two boxes touching at a corner), and removing first
	// double-counts the shared edge of two boxes that abut.
	sort(ev.begin(), ev.end(), [](const array<ll, 4> &a, const array<ll, 4> &b) {
		return a[0] != b[0] ? a[0] < b[0] : a[1] > b[1];
	});
	CoverSeg T(ys);
	ll total = 0;
	for (auto &e : ev) {
		ll before = T.covered();
		T.add((int)e[2], (int)e[3], (int)e[1]);
		total += llabs(T.covered() - before);
	}
	return total;
}
ll rectUnionPerimeter(vector<array<ll, 4>> rs) {
	ll a = _perimSide(rs);
	for (auto &r : rs) { swap(r[0], r[1]); swap(r[2], r[3]); }
	return a + _perimSide(rs);
}

// building outline. Input {left, right, height}; output the corner points of
// the silhouette, i.e. (x, new height) at every place the height changes.
vector<pair<ll, ll>> skyline(const vector<array<ll, 3>> &bs) {
	vector<pair<ll, ll>> ev;                       // {x, +-height}
	for (auto &b : bs)
		if (b[0] < b[1] && b[2] > 0) {
			ev.push_back({b[0], -b[2]});           // negative = start
			ev.push_back({b[1], b[2]});
		}
	sort(ev.begin(), ev.end());
	multiset<ll> live{0};
	vector<pair<ll, ll>> res;
	for (size_t i = 0; i < ev.size();) {
		ll x = ev[i].first;
		while (i < ev.size() && ev[i].first == x) {
			if (ev[i].second < 0) live.insert(-ev[i].second);
			else live.erase(live.find(ev[i].second));
			i++;
		}
		ll h = *live.rbegin();
		if (res.empty() || res.back().second != h) res.push_back({x, h});
	}
	return res;
}

// ------------------------------------------------------- segment sweep
namespace shamos {
ld sweepX;
struct SLine {
	P a, b;
	int id;
	ld yAt(ld x) const {
		if (a.x == b.x) return (ld)min(a.y, b.y);
		return (ld)a.y + (ld)(b.y - a.y) * (x - (ld)a.x) / (ld)(b.x - a.x);
	}
	bool operator<(const SLine &o) const {
		ld p = yAt(sweepX), q = o.yAt(sweepX);
		if (fabs(p - q) > EPS) return p < q;
		return id < o.id;
	}
};
}  // namespace shamos

// {i, j} for SOME pair of intersecting segments, or {-1,-1} if none intersect.
pair<int, int> anySegIntersect(vector<pair<P, P>> segs) {
	using namespace shamos;
	int n = (int)segs.size();
	vector<array<ll, 3>> ev;                        // {x, type(0=open,1=close), id}
	for (int i = 0; i < n; i++) {
		if (segs[i].second < segs[i].first) swap(segs[i].first, segs[i].second);
		ev.push_back({segs[i].first.x, 0, (ll)i});
		ev.push_back({segs[i].second.x, 1, (ll)i});
	}
	sort(ev.begin(), ev.end(), [](const array<ll, 3> &a, const array<ll, 3> &b) {
		return a[0] != b[0] ? a[0] < b[0] : a[1] < b[1];   // opens before closes
	});
	set<SLine> S;
	vector<set<SLine>::iterator> where(n);
	auto cross2 = [&](int i, int j) {
		return i >= 0 && j >= 0 &&
		       segInter(segs[i].first, segs[i].second, segs[j].first, segs[j].second);
	};
	for (auto &e : ev) {
		int id = (int)e[2];
		sweepX = (ld)e[0];
		SLine cur{segs[id].first, segs[id].second, id};
		if (e[1] == 0) {
			auto it = S.insert(cur).first;
			where[id] = it;
			auto nxt = next(it);
			if (nxt != S.end() && cross2(id, nxt->id)) return {id, nxt->id};
			if (it != S.begin()) {
				auto prv = prev(it);
				if (cross2(id, prv->id)) return {id, prv->id};
			}
		} else {
			auto it = where[id];
			auto nxt = next(it);
			if (it != S.begin() && nxt != S.end()) {
				auto prv = prev(it);
				if (cross2(prv->id, nxt->id)) return {prv->id, nxt->id};
			}
			S.erase(it);
		}
	}
	return {-1, -1};
}

// number of crossing pairs among segments that are each horizontal or vertical
ll countHVIntersect(const vector<pair<P, P>> &segs) {
	vector<array<ll, 3>> hs;                        // {y, x1, x2}
	vector<array<ll, 3>> vs;                        // {x, y1, y2}
	vector<ll> xsAll;
	for (auto &s : segs) {
		P a = s.first, b = s.second;
		if (a.y == b.y) {
			hs.push_back({a.y, min(a.x, b.x), max(a.x, b.x)});
			xsAll.push_back(a.x); xsAll.push_back(b.x);
		} else if (a.x == b.x) {
			vs.push_back({a.x, min(a.y, b.y), max(a.y, b.y)});
			xsAll.push_back(a.x);
		}
	}
	sort(xsAll.begin(), xsAll.end());
	xsAll.erase(unique(xsAll.begin(), xsAll.end()), xsAll.end());
	int m = (int)xsAll.size();
	auto xi = [&](ll v) { return (int)(lower_bound(xsAll.begin(), xsAll.end(), v) - xsAll.begin()) + 1; };
	vector<ll> bit(m + 2, 0);
	auto add = [&](int i, ll v) { for (; i <= m; i += i & -i) bit[i] += v; };
	auto qry = [&](int i) { ll s = 0; for (; i > 0; i -= i & -i) s += bit[i]; return s; };
	// f(Y, x) = # horizontals with y <= Y whose x-span contains x, kept in a BIT
	// with range-add / point-query. Answer = sum over verticals of
	//   f(y2, x) - f(y1 - 1, x)
	// which counts exactly the horizontals crossing that vertical.
	sort(hs.begin(), hs.end());
	vector<array<ll, 3>> q;                          // {Y, sign, x}
	for (auto &v : vs) {
		q.push_back({v[1] - 1, -1, v[0]});
		q.push_back({v[2], 1, v[0]});
	}
	sort(q.begin(), q.end());
	ll ans = 0;
	size_t hp = 0;
	for (auto &e : q) {
		while (hp < hs.size() && hs[hp][0] <= e[0]) {
			add(xi(hs[hp][1]), 1);
			add(xi(hs[hp][2]) + 1, -1);
			hp++;
		}
		ans += e[1] * qry(xi(e[2]));
	}
	return ans;
}

// the largest number of intervals [l, r] covering one common point
int maxOverlap(vector<pair<ll, ll>> iv) {
	vector<pair<ll, int>> ev;
	for (auto &e : iv) { ev.push_back({e.first, 1}); ev.push_back({e.second + 1, -1}); }
	sort(ev.begin(), ev.end());
	int cur = 0, best = 0;
	for (auto &e : ev) best = max(best, cur += e.second);
	return best;
}

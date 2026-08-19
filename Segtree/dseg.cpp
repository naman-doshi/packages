// Dynamic (implicit) segtree over [0, SZ): point update, range max query.
// Index-pool based (no pointers) for a low constant factor and good cache use.
// Node 0 is the null sentinel (val = -INF); the root is node 1.
template <long long SZ>
struct dynseg {
  vi val{(int)-INF, (int)-INF};  // val[i] = max over node i's range
  vi lc{0, 0};         // left child index  (0 = none)
  vi rc{0, 0};         // right child index (0 = none)

  int nw() {  // allocate a fresh node, return its index
    val.push_back(-INF);
    lc.push_back(0);
    rc.push_back(0);
    return val.size() - 1;
  }

  void update(long long p, int v, int i = 1, long long l = 0, long long r = SZ - 1) {
    if (l == r) {
      val[i] = v;
      return;
    }
    long long m = (l + r) / 2;
    if (p <= m) {
      if (!lc[i]) lc[i] = nw();
      update(p, v, lc[i], l, m);
    } else {
      if (!rc[i]) rc[i] = nw();
      update(p, v, rc[i], m + 1, r);
    }
    val[i] = max(val[lc[i]], val[rc[i]]);
  }

  int query(long long lo, long long hi, int i = 1, long long l = 0, long long r = SZ - 1) {
    if (!i || lo > r || hi < l) return -INF;
    if (lo <= l && r <= hi) return val[i];
    long long m = (l + r) / 2;
    return max(query(lo, hi, lc[i], l, m), query(lo, hi, rc[i], m + 1, r));
  }
};

// usage: dynseg<1000000000> seg;  seg.update(x, y);  seg.query(lo, hi);

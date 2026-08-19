// ============================================================================
//  ITERATIVE SEGMENT TREE  -- point update, range query (O(log n))
// ----------------------------------------------------------------------------
//  The simplest segtree: flat array, leaves at [sz, 2*sz), no lazy / no range
//  updates. This copy stores range MAX (identity -INF). It's a MONOID: an
//  associative merge (max) plus its identity (-INF) -- change those two and you
//  change the operation. For range UPDATES too, use lazyseg.cpp instead.
//
//  API  (0-based; query range is INCLUSIVE)
//    segtree seg(n);          // n leaves, all -INF
//    seg.update(i, val);      // set position i to val
//    int m = seg.query(l, r); // max over [l, r]
//
//  CUSTOMISING  (full guide + examples in customSegtrees.cpp)
//    change BOTH the merge in pull()/query() AND the identity (seg init + the
//    out-of-range return):
//      sum : seg(2n, 0);   pull: seg[i]=seg[2i]+seg[2i+1];   query return 0
//      min : seg(2n, INF); pull: min(...);                   query return INF
//      gcd : seg(2n, 0);   pull: __gcd(...);                 query return 0
//  Needs `vi = vector<...>` and `INF` from the including file.
// ============================================================================
struct segtree{
 int sz;
 vi seg;
 segtree(int n) : sz(n), seg(2*n, -INF) {}
 void pull(int i){
  seg[i] = max(seg[2*i], seg[2*i+1]);
 }
 void update(int i, int val){
  seg[i+=sz]=val;
  while(i) i/=2, pull(i);
 }
 int query(int lo, int hi, int i=1, int l=0, int r=-1){
  if(r<0) r=sz-1;
  if(lo>r||hi<l) return -INF;
  if(lo<=l&&r<=hi) return seg[i];
  int m=(l+r)/2;
  return max(query(lo,hi,2*i,l,m), query(lo,hi,2*i+1,m+1,r));
 }
};
// usage: segtree seg(n);  seg.update(i, val);  seg.query(lo, hi);

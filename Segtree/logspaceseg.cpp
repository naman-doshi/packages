// Log-space (block-leaf) segment tree.
// Idea: stop the recursion when a node's range has size <= B (B ~ lg n), so the
// leaves are *blocks* of the original array instead of single cells. There are
// then O(n / lg n) leaves and O(n / lg n) nodes total -- a lg n factor less
// memory than a normal O(n)-node tree.
//
// Queries / updates stay O(lg n):
//  - internal nodes contribute their stored aggregate directly (as usual);
//  - a partial boundary leaf [l,r] is handled by scanning a[l..r], which costs
//    O(B) = O(lg n). Any range touches at most two such leaves.
//
// The win shows up when a problem needs many trees (e.g. one per bit): the total
// drops from O(n lg n) to O(n). For a single tree it rarely beats just shrinking
// the value type or using a Fenwick.
//
// This reference version: range-sum query + point update over backing array a[].
// (Relies on `vi` from the including file, like seg.cpp.)
struct logseg{
 int n, B;     // B = block size ~ lg n
 vi a;         // backing array, kept in sync with the tree
 vi sum;       // sum[i] = sum over node i's range

 logseg(const vi& init){
  a = init; n = a.size();
  B = max((int)1, (int)(63 - __builtin_clzll((unsigned long long)max((int)1, n))));  // ~ lg n
  int blocks = (n + B - 1) / B;
  sum.assign(4*blocks + 8, 0);               // O(n / lg n) nodes
  if(n) build(1, 0, n-1);
 }

 void build(int i, int l, int r){
  if(r - l < B){                             // leaf block: store its sum
   int s = 0; for(int k=l;k<=r;k++) s += a[k];
   sum[i] = s; return;
  }
  int m = (l + r) / 2;
  build(2*i, l, m); build(2*i+1, m+1, r);
  sum[i] = sum[2*i] + sum[2*i+1];
 }

 void update(int p, int val, int i=1, int l=0, int r=-1){
  if(r<0) r=n-1;
  if(r - l < B){                             // adjust this block by the delta
   sum[i] += val - a[p]; a[p] = val; return;
  }
  int m = (l + r) / 2;
  if(p <= m) update(p, val, 2*i, l, m);
  else       update(p, val, 2*i+1, m+1, r);
  sum[i] = sum[2*i] + sum[2*i+1];
 }

 int query(int lo, int hi, int i=1, int l=0, int r=-1){
  if(r<0) r=n-1;
  if(lo>r || hi<l) return 0;
  if(lo<=l && r<=hi) return sum[i];          // internal node fully covered
  if(r - l < B){                             // partial leaf: scan the block
   int s = 0; for(int k=max(l,lo);k<=min(r,hi);k++) s += a[k];
   return s;
  }
  int m = (l + r) / 2;
  return query(lo,hi,2*i,l,m) + query(lo,hi,2*i+1,m+1,r);
 }
};
// usage: logseg seg(vi{...});  seg.update(p, val);  seg.query(lo, hi);
//
// To adapt: swap `+`/`0` for your monoid (min/max/etc.). For range *updates*
// with lazy tags, a tag that fully covers a leaf block is stored as usual; a
// partial update of a block is applied by rewriting that block in O(B).

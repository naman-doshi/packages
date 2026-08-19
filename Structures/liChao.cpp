// ============================================================================
//  LI CHAO TREE  -- minimum of a set of lines, evaluated at integer x
// ----------------------------------------------------------------------------
//  (Despite the file name "lct", this is a Li Chao Tree, NOT a link-cut tree.)
//
//  WHAT / WHEN
//    Stores a growing set of lines y = k*x + b and answers "min y over all
//    inserted lines at a given x" in O(log C), C = domain size (maxn). It's the
//    Convex Hull Trick that works when slopes arrive in ARBITRARY order (plain
//    monotonic-stack CHT needs sorted slopes; Li Chao does not). Classic use:
//    speeding up a 1-D DP of the form  dp[i] = min_j (k_j * x_i + b_j).
//
//  ENCODING
//    A line is a complex<ll>: real = slope k, imag = intercept b.
//    f(line, x) = k*x + b.   add_line({k,b}) inserts;   get(x) returns the min.
//
//  API  (x must lie in [0, maxn))
//    add_line({k, b});        // insert line y = k*x + b
//    ll best = get(x);        // min over all inserted lines at x
//
//  MAXIMUM instead of minimum: insert {-k,-b} and negate get(x), OR flip the
//    `<` in add_line and the min() in get to max().
//
//  PITFALLS
//    * PHANTOM LINE: `line[]` is a global, so it starts zero-filled = the line
//      y=0 everywhere. If your true answer can be > 0 you'll wrongly get 0.
//      Fix: before use, fill line[] with a +INF line, e.g.
//          for (auto& L : line) L = point(0, INF);
//      (and reset it between test cases -- it is NOT cleared automatically).
//    * Domain is fixed integer x in [0, maxn). For negative / large / real x,
//      shift or coordinate-compress queries first; an x out of range reads junk.
//    * Overflow: k*x can be huge -- ftype is long long here, keep it.
// ============================================================================
typedef long long ftype;
typedef complex<ftype> point;
#define x real
#define y imag

ftype dot(point a, point b) {
    return (conj(a) * b).x();
}

ftype f(point a,  ftype x) {
    return dot(a, {x, 1});
}

const int maxn = 2e5;

point line[4 * maxn];

void add_line(point nw, int v = 1, int l = 0, int r = maxn) {
    int m = (l + r) / 2;
    bool lef = f(nw, l) < f(line[v], l);
    bool mid = f(nw, m) < f(line[v], m);
    if(mid) {
        swap(line[v], nw);
    }
    if(r - l == 1) {
        return;
    } else if(lef != mid) {
        add_line(nw, 2 * v, l, m);
    } else {
        add_line(nw, 2 * v + 1, m, r);
    }
}

ftype get(int x, int v = 1, int l = 0, int r = maxn) {
    int m = (l + r) / 2;
    if(r - l == 1) {
        return f(line[v], x);
    } else if(x < m) {
        return min(f(line[v], x), get(x, 2 * v, l, m));
    } else {
        return min(f(line[v], x), get(x, 2 * v + 1, m, r));
    }
}
// ============================================================================
//  MATRIX (mod p)  -- multiply, power, determinant, rank, inverse, solve
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Matrix exponentiation is the standard trick for "compute the n-th term of
//    a linear recurrence with n up to 1e18" and "count walks of length k in a
//    graph". Everything else here (det/rank/inverse/solve) is Gaussian
//    elimination over Z_p -- see gauss.cpp for the floating-point version.
//
//  API
//    Matrix A(n, m);  A[i][j] = ...;       // 0-based, filled with 0
//    Matrix::id(n)                          // identity
//    A * B    A + B    A.pow(e)    A.transpose()
//    A.det()          determinant mod ANY modulus (Euclid-based, no inverse)
//    A.rank()         rank mod a PRIME
//    A.inverse()      inverse mod a PRIME; empty (0x0) matrix if singular
//    A.solve(b)       one solution of A x = b mod a PRIME; {} if inconsistent
//    A.pw(e)          alias for pow, if `pow` clashes with std::pow
//
//  COMPLEXITY  multiply O(n^3), power O(n^3 log e), det/rank/inverse O(n^3)
//              (det is O(n^3 log MD) for a composite modulus).
//              n = 100 and e = 1e18 is about 6e8 ops -- borderline, n <= 60 is
//              comfortable. The `if (a[i][k])` skip makes sparse matrices fly.
//
//  PITFALLS
//    * rank/inverse/solve need MD PRIME (they divide). det does NOT -- it uses
//      a Euclidean row-reduction that works for any modulus.
//    * pow() requires a SQUARE matrix; nothing asserts it.
//    * Exponent e must be >= 0. pow(0) is the identity, as it should be.
//    * `A * B` needs A.m == B.n; also unchecked.
//    * If you want mint instead: swap `int` for `mint` and delete every `% MD`.
//
//  RECIPES
//    Linear recurrence f(n) = c1 f(n-1) + ... + ck f(n-k):
//      T = [ c1 c2 ... ck ]      f(n) = (T^(n-k+1) * [f(k-1) ... f(0)]^T)[0]
//          [  1  0 ...  0 ]      i.e. put the coefficients in the top row and
//          [  0  1 ...  0 ]      a shifted identity underneath.
//          [  ...          ]
//    Fibonacci: T = [[1,1],[1,0]], T^n = [[F(n+1), F(n)], [F(n), F(n-1)]].
//    Walks: if A is the adjacency matrix, (A^k)[u][v] = # of walks u->v of
//    length exactly k. Swap (+, *) for (min, +) to get shortest walks with
//    exactly k edges (see floydwarshall.cpp for the all-pairs version).
//    Add a constant term to a recurrence by appending a row/column of 1s.
//    MATRIX-TREE (Kirchhoff): # of spanning trees of a graph = any cofactor of
//    (degree matrix - adjacency matrix), i.e. delete one row and column and
//    take det(). For a directed graph's arborescences rooted at r, delete row
//    and column r from (in-degree matrix - adjacency).
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define all(x) x.begin(), x.end()
typedef vector<int> vi;
typedef vector<vi> vvi;

const int MD = 1000000007;      // PRIME for rank/inverse/solve; any value for det

int powmodM(int a, int e, int m) {
  int r = 1 % m; a %= m; if (a < 0) a += m;
  while (e > 0) { if (e & 1) r = r * a % m; a = a * a % m; e >>= 1; }
  return r;
}

struct Matrix {
  int n, m;
  vvi a;

  Matrix(int n = 0, int m = 0) : n(n), m(m), a(n, vi(m, 0)) {}
  Matrix(const vvi &v) : n(v.size()), m(v.empty() ? 0 : v[0].size()), a(v) {}

  static Matrix id(int n) {
    Matrix r(n, n);
    rep(i, 0, n) r.a[i][i] = 1 % MD;
    return r;
  }

  vi &operator[](int i) { return a[i]; }
  const vi &operator[](int i) const { return a[i]; }

  Matrix operator*(const Matrix &o) const {
    Matrix r(n, o.m);
    rep(i, 0, n) rep(k, 0, m) if (a[i][k])            // skip zeros: big win on
      rep(j, 0, o.m)                                  // sparse transition tables
        r.a[i][j] = (r.a[i][j] + a[i][k] * o.a[k][j]) % MD;
    return r;
  }
  Matrix operator+(const Matrix &o) const {
    Matrix r(n, m);
    rep(i, 0, n) rep(j, 0, m) r.a[i][j] = (a[i][j] + o.a[i][j]) % MD;
    return r;
  }
  Matrix operator*(int c) const {
    Matrix r(n, m);
    rep(i, 0, n) rep(j, 0, m) r.a[i][j] = a[i][j] % MD * (c % MD) % MD;
    return r;
  }
  // matrix times column vector
  vi operator*(const vi &v) const {
    vi r(n, 0);
    rep(i, 0, n) rep(j, 0, m) r[i] = (r[i] + a[i][j] * v[j]) % MD;
    return r;
  }
  Matrix transpose() const {
    Matrix r(m, n);
    rep(i, 0, n) rep(j, 0, m) r.a[j][i] = a[i][j];
    return r;
  }

  // A^e, e >= 0. Square matrices only.
  Matrix pow(int e) const {
    Matrix r = id(n), b = *this;
    while (e > 0) { if (e & 1) r = r * b; b = b * b; e >>= 1; }
    return r;
  }
  Matrix pw(int e) const { return pow(e); }

  // Determinant mod ANY modulus. Instead of dividing by the pivot it runs the
  // Euclidean algorithm between rows (repeated subtract-and-swap), so no
  // modular inverse is ever needed. O(n^3 log MD).
  int det() const {
    Matrix b = *this;
    int res = 1;
    rep(i, 0, n) {
      rep(j, i + 1, n) {
        while (b.a[j][i]) {
          int t = b.a[i][i] / b.a[j][i];
          rep(k, i, n) b.a[i][k] = ((b.a[i][k] - t * b.a[j][k]) % MD + MD) % MD;
          swap(b.a[i], b.a[j]);
          res = (MD - res) % MD;                 // a row swap flips the sign
        }
      }
      res = res * b.a[i][i] % MD;
      if (!res) return 0;
    }
    return res;
  }

  // Row-reduce (in place) to reduced row echelon form mod a PRIME, treating
  // the last `extra` columns as an augmented block. Returns the rank and the
  // index of the pivot row for each column (-1 if that column is free).
  int rrefImpl(vi &where, int extra = 0) {
    where.assign(m - extra, -1);
    int r = 0;
    rep(c, 0, m - extra) {
      if (r == n) break;
      int piv = -1;
      rep(i, r, n) if (a[i][c]) { piv = i; break; }
      if (piv < 0) continue;
      swap(a[r], a[piv]);
      int inv = powmodM(a[r][c], MD - 2, MD);
      rep(j, c, m) a[r][j] = a[r][j] * inv % MD;
      rep(i, 0, n) if (i != r && a[i][c]) {
        int f = a[i][c];
        rep(j, c, m) a[i][j] = ((a[i][j] - f * a[r][j]) % MD + MD) % MD;
      }
      where[c] = r;
      r++;
    }
    return r;
  }

  int rank() const {                             // MD must be PRIME
    Matrix b = *this;
    vi where;
    return b.rrefImpl(where);
  }

  // Inverse mod a PRIME. Returns a 0x0 matrix when the matrix is singular.
  Matrix inverse() const {
    Matrix aug(n, 2 * n);
    rep(i, 0, n) {
      rep(j, 0, n) aug.a[i][j] = a[i][j];
      aug.a[i][n + i] = 1;
    }
    vi where;
    if (aug.rrefImpl(where, n) < n) return Matrix(0, 0);
    Matrix r(n, n);
    rep(i, 0, n) rep(j, 0, n) r.a[i][j] = aug.a[i][n + j];
    return r;
  }

  // One solution of A x = b mod a PRIME, or {} if the system is inconsistent.
  // Free variables are set to 0; the solution is unique iff rank == m.
  vi solve(const vi &b) const {
    Matrix aug(n, m + 1);
    rep(i, 0, n) {
      rep(j, 0, m) aug.a[i][j] = a[i][j];
      aug.a[i][m] = ((b[i] % MD) + MD) % MD;
    }
    vi where;
    int r = aug.rrefImpl(where, 1);
    rep(i, r, n) if (aug.a[i][m]) return {};     // 0 = nonzero  ->  no solution
    vi x(m, 0);
    rep(c, 0, m) if (where[c] != -1) x[c] = aug.a[where[c]][m];
    return x;
  }

  void print() const {
    rep(i, 0, n) { rep(j, 0, m) cout << a[i][j] << " \n"[j == m - 1]; }
  }
};

// ---- ready-made linear recurrence ------------------------------------------
// f(i) = sum_j c[j] * f(i-1-j)   for i >= k, with the first k values in `init`
// (init[0] = f(0), ..., init[k-1] = f(k-1)). Returns f(n) mod MD in O(k^3 log n).
// For large k use linrec.cpp (Kitamasa, O(k^2 log n)).
int linearRecMatrix(const vi &c, const vi &init, int n) {
  int k = c.size();
  if (n < k) return ((init[n] % MD) + MD) % MD;
  Matrix T(k, k);
  rep(j, 0, k) T.a[0][j] = ((c[j] % MD) + MD) % MD;
  rep(i, 1, k) T.a[i][i - 1] = 1;
  Matrix P = T.pow(n - k + 1);
  int res = 0;
  rep(j, 0, k) res = (res + P.a[0][j] * (((init[k - 1 - j] % MD) + MD) % MD)) % MD;
  return res;
}

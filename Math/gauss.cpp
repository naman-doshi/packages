// ============================================================================
//  GAUSSIAN ELIMINATION  -- real-valued systems, determinants, and GF(2)
// ----------------------------------------------------------------------------
//  WHAT / WHEN
//    Solving A x = b when the entries are REAL (physics-flavoured problems,
//    probability / expected-value systems with cycles, least squares) or over
//    GF(2) (systems of XOR equations, light-switch puzzles, Gaussian
//    elimination on a bitset). For arithmetic mod a prime use matrix.cpp.
//
//    Expected value with cycles is the classic use: write E[v] = 1 + sum
//    p(v,u) E[u] for every state, move the unknowns to the left, and solve.
//
//  API
//    solveLinear(A, b, x)     -> rank, or -1 if inconsistent. x holds one
//                                solution; rank < m means free variables exist.
//    detReal(A)               -> determinant (destroys A)
//    solveLinearBinary(A, b, x, m)   same contract over GF(2), A is bitsets
//
//  COMPLEXITY  O(n^2 m) real, O(n^2 m / 64) for the bitset version.
//
//  PITFALLS
//    * FLOATING POINT: EPS = 1e-9 is a guess. Scale it to your input, or the
//      "pivot is zero" test will misfire. Full pivoting (used here) is much
//      more stable than partial pivoting -- keep it.
//    * When rank < m the returned x sets free variables to 0. If you need the
//      whole solution space, read off the null space from the reduced matrix.
//    * solveLinear MODIFIES A and b (they're taken by reference).
//    * MAXC in the binary version is a compile-time bitset width -- set it to
//      the largest column count you'll see, and no larger (memory).
//    * "unique solution" means rank == number of variables m, NOT rank == n.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define all(x) x.begin(), x.end()
typedef vector<double> vd;
typedef vector<int> vi;

const double EPS = 1e-9;

// ---- real systems ----------------------------------------------------------
// Solves A x = b for x. n equations, m unknowns (x must already be sized m).
// Returns the rank, or -1 when the system has NO solution.
//   rank == m  -> unique solution
//   rank <  m  -> infinitely many; x is one of them, free variables set to 0
// Uses full pivoting (searches the whole remaining submatrix for the largest
// entry) and remembers the column permutation so x comes back in the right
// order. A and b are destroyed.
int solveLinear(vector<vd> &A, vd &b, vd &x) {
  int n = A.size(), m = x.size(), rank = 0, br = 0, bc = 0;
  vi col(m);
  iota(all(col), 0);
  rep(i, 0, n) {
    double v, bv = 0;
    rep(r, i, n) rep(c, i, m)
      if ((v = fabs(A[r][c])) > bv) br = r, bc = c, bv = v;
    if (bv <= EPS) {                       // nothing left to pivot on
      rep(j, i, n) if (fabs(b[j]) > EPS) return -1;   // 0 = nonzero
      break;
    }
    swap(A[i], A[br]);
    swap(b[i], b[br]);
    swap(col[i], col[bc]);
    rep(j, 0, n) swap(A[j][i], A[j][bc]);
    bv = 1 / A[i][i];
    rep(j, i + 1, n) {
      double fac = A[j][i] * bv;
      b[j] -= fac * b[i];
      rep(k, i + 1, m) A[j][k] -= fac * A[i][k];
    }
    rank++;
  }
  x.assign(m, 0);
  for (int i = rank; i--;) {               // back substitution
    b[i] /= A[i][i];
    x[col[i]] = b[i];
    rep(j, 0, i) b[j] -= A[j][i] * b[i];
  }
  return rank;
}

// Determinant of a square real matrix. Destroys A. O(n^3).
double detReal(vector<vd> &A) {
  int n = A.size();
  double res = 1;
  rep(i, 0, n) {
    int b = i;
    rep(j, i + 1, n) if (fabs(A[j][i]) > fabs(A[b][i])) b = j;
    if (i != b) swap(A[i], A[b]), res *= -1;
    res *= A[i][i];
    if (fabs(res) < EPS) return 0;
    rep(j, i + 1, n) {
      double v = A[j][i] / A[i][i];
      if (v != 0) rep(k, i + 1, n) A[j][k] -= v * A[i][k];
    }
  }
  return res;
}

// ---- GF(2) -----------------------------------------------------------------
// Same contract as solveLinear but every coefficient is a bit: each row of A is
// a bitset of m coefficients, b[i] is 0 or 1, and all arithmetic is XOR.
// Use it for "flip switch i toggles lamps S_i, reach the all-on state".
// Set MAXC to your maximum number of variables.
const int MAXC = 1000;
int solveLinearBinary(vector<bitset<MAXC>> &A, vi &b, bitset<MAXC> &x, int m) {
  int n = A.size(), rank = 0, br;
  vi col(m);
  iota(all(col), 0);
  rep(i, 0, n) {
    for (br = i; br < n; br++) if (A[br].any()) break;
    if (br == n) {
      rep(j, i, n) if (b[j]) return -1;
      break;
    }
    int bc = i;                            // first set column at or after i
    while (bc < m && !A[br][bc]) bc++;     // (GCC's A[br]._Find_next(i-1), but
                                           //  that is libstdc++-only)
    swap(A[i], A[br]);
    swap(b[i], b[br]);
    swap(col[i], col[bc]);
    rep(j, 0, n) if (A[j][i] != A[j][bc]) {
      A[j].flip(i); A[j].flip(bc);
    }
    rep(j, i + 1, n) if (A[j][i]) {
      b[j] ^= b[i];
      A[j] ^= A[i];
    }
    rank++;
  }
  x = bitset<MAXC>();
  for (int i = rank; i--;) {
    if (!b[i]) continue;
    x[col[i]] = 1;
    rep(j, 0, i) b[j] ^= A[j][i];
  }
  return rank;                             // (multiple solutions if rank < m)
}

// ============================================================================
//  HUNGARIAN ALGORITHM  --  ASSIGNMENT PROBLEM  in O(n^2 * m)
// ----------------------------------------------------------------------------
//  Given an n x m cost matrix (n <= m), assign each of the n rows to a DISTINCT
//  column so the total cost is MINIMISED (a min-cost perfect matching of the
//  rows). Much faster than MCMF on dense matrices; the standard choice when
//  n, m are up to a few hundred / low thousand.
//
//  API  (this implementation is 1-indexed internally; pass a 0-indexed matrix)
//    Hungarian hu;
//    auto [cost, asg] = hu.solve(a);   // a is n x m (vector<vector<ll>>), n<=m
//    // cost      = minimum total assignment cost
//    // asg[i]    = column (0-indexed) assigned to row i
//
//  TRICKS
//    * MAXIMISE instead: negate all costs, or subtract each from a big const.
//    * need a square perfect matching: pad the smaller side with zero-cost
//      (or large-cost, if forbidding) dummy rows/columns.
//    * forbidden pairs: set that cell to a huge value (e.g. 1e15), never INF,
//      so potentials don't overflow.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define ll long long

const ll BIG = 1e18;

struct Hungarian {
	// a: n x m cost matrix, n <= m. Returns {min cost, assignment (row->col)}.
	pair<ll, vector<int>> solve(const vector<vector<ll>> &a) {
		int n = a.size(), m = a[0].size();
		// u,v = dual potentials; p[j] = row matched to column j (0 = none)
		vector<ll> u(n + 1, 0), v(m + 1, 0);
		vector<int> p(m + 1, 0), way(m + 1, 0);
		for (int i = 1; i <= n; i++) {
			p[0] = i;
			int j0 = 0;
			vector<ll> minv(m + 1, BIG);
			vector<char> used(m + 1, 0);
			do {
				used[j0] = 1;
				int i0 = p[j0], j1 = -1;
				ll delta = BIG;
				for (int j = 1; j <= m; j++) if (!used[j]) {
					ll cur = a[i0 - 1][j - 1] - u[i0] - v[j];
					if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
					if (minv[j] < delta) { delta = minv[j]; j1 = j; }
				}
				for (int j = 0; j <= m; j++)
					if (used[j]) { u[p[j]] += delta; v[j] -= delta; }
					else minv[j] -= delta;
				j0 = j1;
			} while (p[j0] != 0);
			do { int j1 = way[j0]; p[j0] = p[j1]; j0 = j1; } while (j0);
		}
		vector<int> asg(n);
		for (int j = 1; j <= m; j++) if (p[j] > 0) asg[p[j] - 1] = j - 1;
		ll cost = 0;
		for (int i = 0; i < n; i++) cost += a[i][asg[i]];
		return {cost, asg};
	}
};

int main() {
	// 3 workers x 3 jobs; minimise total cost.
	vector<vector<ll>> a = {
		{4, 1, 3},
		{2, 0, 5},
		{3, 2, 2},
	};
	Hungarian hu;
	auto [cost, asg] = hu.solve(a);
	cout << "min total cost = " << cost << "\n";
	for (int i = 0; i < (int)asg.size(); i++)
		cout << "  worker " << i << " -> job " << asg[i] << "\n";
	return 0;
}

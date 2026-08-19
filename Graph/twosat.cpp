// ============================================================================
//  2-SAT  (n boolean variables, clauses of size <= 2, O(V + E) via Tarjan SCC)
// ----------------------------------------------------------------------------
//  LITERALS are signed ints:  i  == "x_i is true",  ~i  == "x_i is false".
//  (~i is -i-1, so ~0 = -1 is a valid literal for variable 0. Never write -i.)
//
//  Internally literal l is node  2*l   if l >= 0
//                                2*(~l)+1 if l < 0.
//  A clause (a OR b) becomes the two implications  !a -> b  and  !b -> a.
//  Satisfiable iff no variable has x_i and !x_i in the same SCC; the assignment
//  is read off by comparing topological positions of the two components.
//
//  API  (variables 0-based)
//    TwoSAT ts; ts.init(N);       // N variables
//    ts.either(a, b);             // (a OR b)
//    ts.implies(a, b);            // a -> b
//    ts.setValue(a);              // force literal a true
//    ts.notBoth(a, b);            // at most one of a, b
//    ts.atLeastTwo(a, b, c);      // >= 2 of the 3 literals true
//    ts.equal(a, b);              // a <-> b
//    int k = ts.addVar();         // extra variable (allowed before solve())
//    if (ts.solve()) ... ts.val[i] ...   // val[i] in {0,1}
//
//  NOTE: SCC::dfs is recursive -- depth up to 2N. Fine to ~2e5 variables on
//  Codeforces (256MB stack), otherwise raise the stack or use iterative SCC.
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
#define int long long

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define all(x) x.begin(), x.end()
typedef vector<int> vi;
typedef pair<int, int> pii;

struct SCC {
	int N, ti = 0;
	vector<vector<int>> adj;
	vector<int> disc, low, comp, st, comps;
	void init(int _N) {
		N = _N; ti = 0; st.clear(); comps.clear();
		adj.assign(N, {}), disc.assign(N, 0), low.assign(N, 0), comp.assign(N, -1);
	}
	void ae(int x, int y) { adj[x].push_back(y); }
	void dfs(int x) {
		disc[x] = low[x] = ++ti;
		st.push_back(x);  // disc[y] != 0 -> in stack
		for (int y : adj[x])
			if (comp[y] == -1) low[x] = min(low[x], disc[y] ?: (dfs(y), low[y]));
		if (low[x] == disc[x]) {  // make new SCC, pop off stack until you find x
			comps.push_back(x);
			for (int y = -1; y != x;) comp[y = st.back()] = x, st.pop_back();
		}
	}
	void gen() {
		for (int i = 0; i < N; i++)
			if (!disc[i]) dfs(i);
		reverse(begin(comps), end(comps));  // now topological order, sources first
	}
};

struct TwoSAT {
	int n = 0;
	vector<pii> edges;   // implication edges, built lazily so addVar() stays legal
	vi val;              // val[i] = final assignment of x_i (valid after solve())

	void init(int _n) { n = _n; edges.clear(); val.clear(); }
	int addVar() { return n++; }

	int node(int l) { return l >= 0 ? 2 * l : 2 * (~l) + 1; }

	void either(int a, int b) {                    // (a OR b)
		edges.push_back({node(~a), node(b)});
		edges.push_back({node(~b), node(a)});
	}
	void implies(int a, int b)  { either(~a, b); }   // a -> b
	void setValue(int a)        { either(a, a); }    // a must be true
	void notBoth(int a, int b)  { either(~a, ~b); }  // at most one of a, b
	void equal(int a, int b)    { implies(a, b), implies(b, a); }
	void atLeastTwo(int a, int b, int c) {           // >= 2 of 3  <=>  no 2 false
		either(a, b), either(a, c), either(b, c);
	}

	bool solve() {
		SCC scc; scc.init(2 * n);
		for (auto &e : edges) scc.ae(e.first, e.second);
		scc.gen();
		vi pos(2 * n);                               // topological rank of each comp
		rep(i, 0, (int)scc.comps.size()) pos[scc.comps[i]] = i;
		val.assign(n, 0);
		rep(i, 0, n) {
			if (scc.comp[2 * i] == scc.comp[2 * i + 1]) return false;
			val[i] = pos[scc.comp[2 * i]] > pos[scc.comp[2 * i + 1]];
		}
		return true;
	}

	// at most one of the given literals is true; O(k) clauses + O(k) aux variables
	// (the naive "every pair" version is O(k^2) -- use this when k is large)
	void atMostOne(const vi &lits) {
		if ((int)lits.size() <= 1) return;
		int cur = ~lits[0];
		rep(i, 2, (int)lits.size()) {
			int nxt = addVar();
			either(cur, ~lits[i]);
			either(cur, nxt);
			either(~lits[i], nxt);
			cur = ~nxt;
		}
		either(cur, ~lits[1]);
	}
};

signed main() {
	// demo: (x0 OR x1) AND (!x0 OR x2) AND (!x1 OR !x2) AND x0
	TwoSAT ts; ts.init(3);
	ts.either(0, 1);
	ts.either(~0, 2);
	ts.either(~1, ~2);
	ts.setValue(0);

	if (!ts.solve()) { cout << "UNSAT\n"; return 0; }
	cout << "SAT:";
	rep(i, 0, ts.n) cout << " x" << i << "=" << ts.val[i];
	cout << "\n";
	return 0;
}

// ============================================================================
//  GRIDS AND IMPLICIT GRAPHS  --  the patterns that actually appear
// ----------------------------------------------------------------------------
//  A grid is just a graph you never build: the neighbours of (i,j) are the four
//  (or eight) cells around it. Do not construct an adjacency list for a
//  1000x1000 grid -- generate neighbours on the fly.
//
//  CONTENTS
//    1.  direction arrays and the bounds check
//    2.  gridBFS         single source shortest path on a grid
//    3.  multiSourceBFS  distance to the NEAREST of many sources (one BFS!)
//    4.  gridZeroOne     0-1 BFS on a grid (some moves free)
//    5.  floodFill       component labelling / counting
//    6.  escapeEnemy     binary search + BFS (the lecture's Escape From Enemy
//                        Territory) -- the single most reusable pattern here
//    7.  implicitBFS     BFS over states that are not cells at all
//    8.  NOTES
// ============================================================================
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

const int INF = 1e9;

// ============================================================================
//  1. DIRECTIONS
// ----------------------------------------------------------------------------
//  Write these down once and never index them by hand again.
// ============================================================================
const int di4[4] = {-1, 1, 0, 0};
const int dj4[4] = {0, 0, -1, 1};

const int di8[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
const int dj8[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

// knight moves, for the "minimum knight moves" family
const int diK[8] = {-2, -2, -1, -1, 1, 1, 2, 2};
const int djK[8] = {-1, 1, -2, 2, -2, 2, -1, 1};

// ============================================================================
//  2. GRID BFS -- shortest path when every move costs 1
// ----------------------------------------------------------------------------
//  grid[i][j] == '#' is a wall. Returns the distance matrix (INF = unreachable)
//  and, if you pass one, a predecessor matrix for reconstructing the path.
// ============================================================================
vector<vector<int>> gridBFS(const vector<string>& g, int si, int sj,
                            vector<vector<int>>* par = nullptr) {
	int n = (int)g.size(), m = (int)g[0].size();
	vector<vector<int>> dist(n, vector<int>(m, INF));
	if (par) par->assign(n, vector<int>(m, -1));
	queue<pair<int, int>> q;
	dist[si][sj] = 0;
	q.push({si, sj});
	while (!q.empty()) {
		auto [i, j] = q.front();
		q.pop();
		for (int d = 0; d < 4; d++) {
			int ni = i + di4[d], nj = j + dj4[d];
			if (ni < 0 || ni >= n || nj < 0 || nj >= m) continue;  // off board
			if (g[ni][nj] == '#') continue;                        // wall
			if (dist[ni][nj] != INF) continue;                     // seen
			dist[ni][nj] = dist[i][j] + 1;
			if (par) (*par)[ni][nj] = i * m + j;
			q.push({ni, nj});
		}
	}
	return dist;
}

// ============================================================================
//  3. MULTI-SOURCE BFS -- distance to the CLOSEST source
// ----------------------------------------------------------------------------
//  Seed the queue with EVERY source at distance 0 and run one BFS. This is
//  equivalent to adding a virtual node joined to all sources, and it is O(nm)
//  regardless of how many sources there are. Running one BFS per source is the
//  classic TLE.
//
//  Uses: "distance from each cell to the nearest fire / shop / enemy", and as
//  the preprocessing step for the binary-search pattern in (6).
// ============================================================================
vector<vector<int>> multiSourceBFS(const vector<string>& g,
                                   const vector<pair<int, int>>& sources) {
	int n = (int)g.size(), m = (int)g[0].size();
	vector<vector<int>> dist(n, vector<int>(m, INF));
	queue<pair<int, int>> q;
	for (auto [i, j] : sources) {
		if (dist[i][j] != 0) { dist[i][j] = 0; q.push({i, j}); }
	}
	while (!q.empty()) {
		auto [i, j] = q.front();
		q.pop();
		for (int d = 0; d < 4; d++) {
			int ni = i + di4[d], nj = j + dj4[d];
			if (ni < 0 || ni >= n || nj < 0 || nj >= m) continue;
			if (g[ni][nj] == '#') continue;
			if (dist[ni][nj] != INF) continue;
			dist[ni][nj] = dist[i][j] + 1;
			q.push({ni, nj});
		}
	}
	return dist;
}

// ============================================================================
//  4. 0-1 BFS ON A GRID -- some moves are free
// ----------------------------------------------------------------------------
//  cost[ni][nj] is 0 or 1. Deque: free moves to the front, paid moves to the
//  back. "You may break up to k walls" is usually this (cost 1 to enter a wall)
//  followed by reading off dist[target] <= k, or a (cell, walls broken) state.
// ============================================================================
vector<vector<int>> gridZeroOne(const vector<string>& g, int si, int sj) {
	int n = (int)g.size(), m = (int)g[0].size();
	vector<vector<int>> dist(n, vector<int>(m, INF));
	deque<pair<int, int>> dq;
	dist[si][sj] = 0;
	dq.push_back({si, sj});
	while (!dq.empty()) {
		auto [i, j] = dq.front();
		dq.pop_front();
		for (int d = 0; d < 4; d++) {
			int ni = i + di4[d], nj = j + dj4[d];
			if (ni < 0 || ni >= n || nj < 0 || nj >= m) continue;
			int w = (g[ni][nj] == '#') ? 1 : 0;      // entering a wall costs 1
			if (dist[i][j] + w < dist[ni][nj]) {
				dist[ni][nj] = dist[i][j] + w;
				if (w == 0) dq.push_front({ni, nj});
				else        dq.push_back({ni, nj});
			}
		}
	}
	return dist;
}

// ============================================================================
//  5. FLOOD FILL -- label connected components
// ----------------------------------------------------------------------------
//  Iterative (a recursive flood fill on a 1000x1000 grid WILL stack overflow).
//  Returns comp[i][j] = component id, or -1 for walls, and the component count.
// ============================================================================
int floodFill(const vector<string>& g, vector<vector<int>>& comp) {
	int n = (int)g.size(), m = (int)g[0].size();
	comp.assign(n, vector<int>(m, -1));
	int id = 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++) {
			if (g[i][j] == '#' || comp[i][j] != -1) continue;
			vector<pair<int, int>> st = {{i, j}};
			comp[i][j] = id;
			while (!st.empty()) {
				auto [a, b] = st.back();
				st.pop_back();
				for (int d = 0; d < 4; d++) {
					int na = a + di4[d], nb = b + dj4[d];
					if (na < 0 || na >= n || nb < 0 || nb >= m) continue;
					if (g[na][nb] == '#' || comp[na][nb] != -1) continue;
					comp[na][nb] = id;
					st.push_back({na, nb});
				}
			}
			id++;
		}
	return id;
}

// ============================================================================
//  6. BINARY SEARCH + BFS   -- "maximise the minimum clearance"
// ----------------------------------------------------------------------------
//  THE LECTURE'S "Escape From Enemy Territory". Learn this shape; it recurs
//  constantly and it is nearly always worth marks.
//
//  Optimisation problem: get from A to B while staying as FAR as possible from
//  every enemy. Turn it into a DECISION problem: "can I get from A to B never
//  coming within X of an enemy?" That is a plain BFS on the grid restricted to
//  cells with distToEnemy >= X. Feasibility is MONOTONE in X (if X works, so
//  does X-1), so binary search X.
//
//  Step 1 -- ONE multi-source BFS from all enemies gives distToEnemy[][].
//  Step 2 -- binary search X; each check is an O(nm) BFS.
//  Total O(nm log(n+m)). Returns {best X, shortest path length at that X}.
//
//  The general recipe: whenever an optimisation is hard but the corresponding
//  "is value X achievable" is easy AND monotone, binary search the answer.
// ============================================================================
pair<int, int> escapeEnemy(const vector<string>& g,
                           const vector<pair<int, int>>& enemies,
                           int i1, int j1, int i2, int j2) {
	int n = (int)g.size(), m = (int)g[0].size();
	vector<vector<int>> de = multiSourceBFS(g, enemies);

	// can we get from (i1,j1) to (i2,j2) staying >= X away? returns path len or -1
	auto check = [&](int X) -> int {
		if (de[i1][j1] < X || de[i2][j2] < X) return -1;
		vector<vector<int>> dist(n, vector<int>(m, INF));
		queue<pair<int, int>> q;
		dist[i1][j1] = 0;
		q.push({i1, j1});
		while (!q.empty()) {
			auto [i, j] = q.front();
			q.pop();
			for (int d = 0; d < 4; d++) {
				int ni = i + di4[d], nj = j + dj4[d];
				if (ni < 0 || ni >= n || nj < 0 || nj >= m) continue;
				if (g[ni][nj] == '#') continue;
				if (de[ni][nj] < X) continue;            // too close to an enemy
				if (dist[ni][nj] != INF) continue;
				dist[ni][nj] = dist[i][j] + 1;
				q.push({ni, nj});
			}
		}
		return dist[i2][j2] == INF ? -1 : dist[i2][j2];
	};

	int lo = 0, hi = min(de[i1][j1], de[i2][j2]), bestX = -1, bestLen = -1;
	while (lo <= hi) {                       // find the largest feasible X
		int mid = lo + (hi - lo) / 2;
		int len = check(mid);
		if (len != -1) { bestX = mid; bestLen = len; lo = mid + 1; }
		else hi = mid - 1;
	}
	return {bestX, bestLen};
}

// ============================================================================
//  7. IMPLICIT BFS OVER NON-GRID STATES
// ----------------------------------------------------------------------------
//  The lecture's Two Buttons: from n, reach m using x -> 2x and x -> x-1.
//  There is no adjacency list; neighbours are computed when you pop.
//
//  THE IMPORTANT PART IS BOUNDING THE STATE SPACE. Here: after the first
//  doubling you never want two decrements in a row (n -> 2n -> 2n-1 -> 2n-2
//  costs 3 moves, while n -> n-1 -> 2n-2 costs 2), so you never need a value
//  above m+1. Prove a bound like that, then use a flat array over it. A hash
//  map keyed on the state is usually too slow at 1e7 states.
//
//  Same pattern: sliding puzzles, "reach target by these operations", word
//  ladders, "cheapest sequence of button presses". Encode the state as an
//  integer, keep dist[] as an array, generate transitions inline.
// ============================================================================
int twoButtons(int n, int m) {
	if (n >= m) return n - m;                 // doubling can never help
	int LIM = 2 * m + 2;
	vector<int> dist(LIM, INF);
	queue<int> q;
	dist[n] = 0;
	q.push(n);
	while (!q.empty()) {
		int x = q.front();
		q.pop();
		if (x == m) return dist[x];
		if (x - 1 > 0 && dist[x] + 1 < dist[x - 1]) {
			dist[x - 1] = dist[x] + 1;
			q.push(x - 1);
		}
		if (2 * x < LIM && dist[x] + 1 < dist[2 * x]) {
			dist[2 * x] = dist[x] + 1;
			q.push(2 * x);
		}
	}
	return -1;
}

// ============================================================================
//  8. NOTES
// ----------------------------------------------------------------------------
//  INDEXING. Flatten (i,j) to i*m + j when you need it as an int (DSU, a
//  priority queue, a visited bitset). Padding the grid with a border of walls
//  removes every bounds check, at the cost of some confusion -- your call.
//
//  BFS ONLY GIVES SHORTEST PATHS WHEN ALL EDGES COST THE SAME. The moment
//  moves have different costs it is Dijkstra (or 0-1 BFS if the costs are 0/1).
//  A very common wrong answer is BFS on a weighted grid.
//
//  MARK AS SEEN WHEN YOU PUSH, NOT WHEN YOU POP. Marking on pop lets the same
//  cell enter the queue many times; on a big grid that is a TLE or an MLE.
//
//  BIPARTITENESS. A grid with 4-adjacency is bipartite (colour by (i+j)%2 --
//  every move flips the parity, so there are no odd cycles). Knight moves too.
//  This is why "cover the grid with dominoes", "maximum independent set of
//  cells", and "minimum cells to hit all adjacent pairs" are matching/min-cut
//  problems on a bipartite graph -- see flow/modeling.cpp and the rows-and-
//  columns trick (the lecture's Irrigation problem: rows and columns are the
//  two sides, a flower is an edge, the answer is a minimum vertex cover).
//
//  STATE EXPLOSION. (cell, keys collected) with 6 keys is 64 * nm states --
//  fine. With 20 keys it is not. Check the state count before you write it.
//
//  0-1 BFS VS DIJKSTRA: if weights are only 0 and 1, the deque version is
//  O(V+E) with a tiny constant. Weights in {0..k}: Dial's algorithm with k+1
//  buckets, or just use Dijkstra and stop worrying.
// ============================================================================

int main() {
	vector<string> g = {
		"....#.....",
		".##.#.###.",
		".#..#...#.",
		".#.###.##.",
		".........."
	};

	vector<vector<int>> dist = gridBFS(g, 0, 0);
	cout << "BFS dist to (4,9) = " << dist[4][9] << "\n";

	vector<vector<int>> comp;
	cout << "open components = " << floodFill(g, comp) << "\n";

	vector<pair<int, int>> enemies = {{0, 9}, {4, 0}};
	auto [X, len] = escapeEnemy(g, enemies, 0, 0, 4, 9);
	cout << "escape: max clearance = " << X << ", path length = " << len << "\n";

	cout << "twoButtons(4, 6) = " << twoButtons(4, 6) << "\n";   // 4->3->6 = 2
	cout << "twoButtons(15, 16) = " << twoButtons(15, 16) << "\n";
	return 0;
}

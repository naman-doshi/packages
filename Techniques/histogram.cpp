// ============================================================================
//  LARGEST RECTANGLE IN A HISTOGRAM  --  O(n) monotonic stack
// ----------------------------------------------------------------------------
//  PROBLEM
//    Given bar heights h[0..n-1] (each bar width 1), find the area of the
//    largest axis-aligned rectangle that fits under the skyline.
//
//  THE O(n) IDEA
//    Every maximal rectangle is limited by some bar that is its SHORTEST bar.
//    So for each bar i, the biggest rectangle of height h[i] stretches left and
//    right until it hits a STRICTLY SMALLER bar. If L = index of nearest smaller
//    bar to the left and R = nearest smaller to the right, the rectangle spans
//    (R - L - 1) bars, giving area h[i] * (R - L - 1).
//    A monotonic stack of INCREASING heights finds those boundaries for all bars
//    in one pass: when a bar is popped, the current index is its right boundary
//    and the new stack top is its left boundary. Each index is pushed/popped
//    once -> O(n).
//
//  API
//    long long largestRectangle(const vector<long long>& h);   // max area, 0 if empty
//
//  COMPLEXITY  O(n) time, O(n) stack.
//
//  PITFALLS
//    * OVERFLOW: area = height * width can exceed 32 bits -> use long long.
//    * The trailing sentinel (cur = -1 at i == n) is what FLUSHES the stack so
//      bars still on it at the end are counted. Don't drop it.
//    * `>=` when popping treats equal-height bars as poppable; that's fine (a
//      later equal bar accounts for the full width). `>` also works.
//
//  RELATED (same nearest-smaller machinery)
//    * Maximal rectangle of 1s in a 0/1 matrix -> build a histogram per row and
//      call this on each row: O(rows * cols). (maximalRectangle below.)
//    * Sum of minimums over all subarrays; stock-span; "trapping rain water".
// ============================================================================
#include <bits/stdc++.h>
using namespace std;

long long largestRectangle(const vector<long long>& h) {
    int n = (int)h.size();
    vector<int> st;                       // indices; heights strictly increasing
    long long best = 0;
    for (int i = 0; i <= n; i++) {
        long long cur = (i == n) ? -1 : h[i];         // -1 sentinel flushes at end
        while (!st.empty() && h[st.back()] >= cur) {
            long long height = h[st.back()];
            st.pop_back();
            int left = st.empty() ? -1 : st.back();    // nearest smaller to the left
            long long width = i - left - 1;            // right boundary is i
            best = max(best, height * width);
        }
        st.push_back(i);
    }
    return best;
}

// ---- Application: largest all-1s rectangle in a 0/1 grid, O(rows*cols) -------
long long maximalRectangle(const vector<vector<int>>& g) {
    if (g.empty()) return 0;
    int C = (int)g[0].size();
    vector<long long> h(C, 0);
    long long best = 0;
    for (const auto& row : g) {
        for (int c = 0; c < C; c++)
            h[c] = row[c] ? h[c] + 1 : 0;              // height of 1s ending at this row
        best = max(best, largestRectangle(h));
    }
    return best;
}

/*  USAGE
      vector<long long> h = {2, 1, 5, 6, 2, 3};
      cout << largestRectangle(h) << "\n";   // 10  (bars 5,6 -> 5 * 2)

      vector<vector<int>> g = {{1,1,1,0},{1,1,1,0},{0,0,1,1}};
      cout << maximalRectangle(g) << "\n";   // 6  (the 2x3 block of 1s top-left)
*/

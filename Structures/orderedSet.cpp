// ============================================================================
//  GNU PBDS — Order Statistics Tree  (a.k.a. "ordered_set")
// ----------------------------------------------------------------------------
//  A balanced binary search tree (red-black) from libstdc++'s policy-based
//  data structures, augmented so it supports two extra O(log n) queries that
//  std::set does NOT:
//
//      find_by_order(k)  ->  iterator to the k-th smallest element (0-indexed)
//      order_of_key(x)   ->  number of elements STRICTLY LESS THAN x
//                            (== the rank/index where x would sit)
//
//  Everything else behaves like std::set: insert/erase/find/lower_bound, all
//  O(log n), iteration in sorted order.
//
//  COMPILER: GNU g++ only (libstdc++). Does NOT work with MSVC or clang's
//  libc++ (so macOS's default `g++`, which is clang, will fail — use a real
//  g++, or just rely on the online judge / CSE machines which use GNU g++).
// ============================================================================

#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace std;
using namespace __gnu_pbds;

// ---- The workhorse: a set of unique keys with order statistics. ------------
template <class T>
using ordered_set = tree<
    T,                              // key type
    null_type,                      // null_type => it's a SET (use a value type here for a MAP)
    less<T>,                        // comparator (strict weak ordering)
    rb_tree_tag,                    // red-black tree
    tree_order_statistics_node_update>;

// ---- Want duplicates (a multiset)? Use less_equal<T>. ----------------------
//  CAVEAT: with less_equal, `find` is broken and `erase(value)` won't work —
//  erase by iterator instead:  it = ms.upper_bound(x); if (it!=ms.begin()) ms.erase(--it);
//  order_of_key / find_by_order still work and are the usual reason to do this.
template <class T>
using ordered_multiset = tree<
    T, null_type, less_equal<T>, rb_tree_tag, tree_order_statistics_node_update>;

// ---- As a MAP with order statistics: put the mapped type instead of null_type
//      tree<int, string, less<int>, rb_tree_tag, tree_order_statistics_node_update>

int main() {
    ordered_set<int> s;
    for (int x : {10, 20, 30, 40, 50}) s.insert(x);   // {10,20,30,40,50}

    // --- find_by_order(k): the k-th smallest, 0-indexed ---------------------
    cout << *s.find_by_order(0) << "\n";   // 10  (smallest)
    cout << *s.find_by_order(2) << "\n";   // 30  (3rd smallest)
    cout << *s.find_by_order(4) << "\n";   // 50  (largest)
    // find_by_order(k) == end() when k >= size(), so guard it:
    if (s.find_by_order(9) == s.end()) cout << "(no 10th element)\n";

    // --- order_of_key(x): how many elements are STRICTLY LESS than x --------
    cout << s.order_of_key(30)  << "\n";   // 2   -> {10,20}  (strict: 30 itself excluded)
    cout << s.order_of_key(35)  << "\n";   // 3   -> {10,20,30}
    cout << s.order_of_key(10)  << "\n";   // 0
    cout << s.order_of_key(999) << "\n";   // 5   -> all of them

    // --- rank of an element, and the median --------------------------------
    int v = 40;
    cout << "rank of 40 = " << s.order_of_key(v) << "\n";          // 3 (0-indexed)
    cout << "median = " << *s.find_by_order(s.size() / 2) << "\n"; // 30

    // --- erase works just like std::set for the unique-key version ---------
    s.erase(20);                            // {10,30,40,50}
    cout << "size = " << s.size() << "\n";  // 4
    cout << *s.find_by_order(1) << "\n";    // 30

    return 0;
}

// ============================================================================
//  COMMON USES
//    * k-th smallest / k-th largest in a dynamic set                (order stat)
//    * rank of a value / "how many are less than x"                 (order_of_key)
//    * running median                                               (find_by_order(size/2))
//    * counting inversions while inserting left-to-right:
//        for i: inversions += i - st.order_of_key(a[i]); st.insert(a[i]);
//    * "count of smaller numbers after self" (insert from the right)
//
//  HANDLING DUPLICATES WITHOUT less_equal (safer):
//    store pairs {value, uniqueTimestamp} in an ordered_set<pair<int,int>>.
//    Then order_of_key({value, 0}) = count strictly less than `value`.
//
//  GOTCHA: a global `#define int long long` changes the tree's key type too;
//    it still works, but keep that in mind if you mix it with this snippet.
// ============================================================================

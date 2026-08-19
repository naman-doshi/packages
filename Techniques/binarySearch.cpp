// ============================================================================
//  BINARY SEARCH TOOLKIT
// ----------------------------------------------------------------------------
//  The mental model: binary search finds the BOUNDARY of a monotone predicate.
//  If ok(x) is false...false true...true, firstTrue gives the first true; if it
//  is true...true false...false, lastTrue gives the last true. "Binary search
//  the answer" = phrase the answer as such a predicate and call these.
//
//  WHAT'S HERE
//    predicate boundary (int) : first_true / last_true / last_true_jump
//    predicate boundary (ll)  : firstTrue / lastTrue           <- use for big ranges
//    real-valued boundary     : firstTrueReal                  <- eps / precision
//    sorted-array lookup      : binarySearch (exact) · lowerBoundIdx · upperBoundIdx
//    unimodal optimum (int)   : minimise / maximise (min/max VALUE)
//                               argminInt / argmaxInt (the ARGUMENT)
//    unimodal optimum (real)  : ternaryMinReal
//    integer sqrt             : isqrt
//    rotated sorted array     : searchRotated
//
//  GOTCHAS
//    * firstTrue/lastTrue need a MONOTONE predicate; test ok() is really monotone.
//    * ternary/minimise need a UNIMODAL function (one valley / one peak).
//    * mid overflow: use lo + (hi-lo)/2, and long long when ranges are large.
//    * All of these are brute-force verified (see the repo test harness).
// ============================================================================
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <queue>
#include <deque>
#include <bitset>
#include <iterator>
#include <list>
#include <stack>
#include <map>
#include <set>
#include <functional>
#include <numeric>
#include <utility>
#include <limits>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cmath>

using namespace std;

int last_true(int lo, int hi, function<bool(int)> f) {
	// if none of the values in the range work, return lo - 1
	lo--;
	while (lo < hi) {
		// find the middle of the current range (rounding up)
		int mid = lo + (hi - lo + 1) / 2;
		if (f(mid)) {
			// if mid works, then all numbers smaller than mid also work
			lo = mid;
		} else {
			// if mid does not work, greater values would not work either
			hi = mid - 1;
		}
	}
	return lo;
}

int last_true_jump(int lo, int hi, function<bool(int)> f) {
	lo--;
	for (int dif = hi - lo; dif > 0; dif /= 2) {
		while (lo + dif <= hi && f(lo + dif)) { lo += dif; }
	}
	return lo;
}

int first_true(int lo, int hi, function<bool(int)> f) {
  // if none of the values in the range work, return hi + 1
  hi++;
  while (lo < hi) {
    // find the middle of the current range (rounding down)
    int mid = lo + (hi - lo) / 2;
    if (f(mid)) {
      // if mid works, then all numbers greater than mid also work
      hi = mid;
    } else {
      // if mid does not work, smaller values would not work either
      lo = mid + 1;
    }
  }
  return hi;
}

long long f(int x) {
  return x * x;
}

// Minimum VALUE of a unimodal f over [low, high]. Ternary search that stays
// strictly inside the range (the old version peeked f(high+1) and was wrong when
// the minimum was at the right endpoint).
int minimise(int low, int high) {
  while (high - low > 2) {
    int m1 = low + (high - low) / 3;
    int m2 = high - (high - low) / 3;
    if (f(m1) < f(m2)) high = m2; else low = m1;
  }
  long long best = f(low);
  for (int x = low + 1; x <= high; x++) best = min(best, (long long)f(x));
  return (int)best;
}

// Maximum VALUE of a unimodal f over [low, high].
int maximise(int low, int high) {
  while (high - low > 2) {
    int m1 = low + (high - low) / 3;
    int m2 = high - (high - low) / 3;
    if (f(m1) > f(m2)) high = m2; else low = m1;
  }
  long long best = f(low);
  for (int x = low + 1; x <= high; x++) best = max(best, (long long)f(x));
  return (int)best;
}


int binarySearch(vector<int> arr, int search) {
  int l = 0;
  int r = arr.size() - 1;
  while (l <= r) {
    int m = (l + r) / 2;
    if (arr[m] == search) {
      return m;
    } else if (arr[m] < search) {
      l = m + 1;
    } else {
      r = m - 1;
    }
  }
  return -1;
}

// ---------------------------------------------------------------------------
//  ADDED VARIANTS
// ---------------------------------------------------------------------------

// long long "binary search the answer": smallest x in [lo,hi] with ok(x) true
// (predicate must be false...false true...true). Returns hi+1 if none true.
long long firstTrue(long long lo, long long hi, const function<bool(long long)>& ok) {
  long long res = hi + 1;
  while (lo <= hi) {
    long long mid = lo + (hi - lo) / 2;
    if (ok(mid)) { res = mid; hi = mid - 1; }
    else lo = mid + 1;
  }
  return res;
}

// largest x in [lo,hi] with ok(x) true (predicate true...true false...false).
// Returns lo-1 if none true.
long long lastTrue(long long lo, long long hi, const function<bool(long long)>& ok) {
  long long res = lo - 1;
  while (lo <= hi) {
    long long mid = lo + (hi - lo) / 2;
    if (ok(mid)) { res = mid; lo = mid + 1; }
    else hi = mid - 1;
  }
  return res;
}

// real-valued boundary: smallest x in [lo,hi] with ok(x) true, to ~iters bits.
double firstTrueReal(double lo, double hi, const function<bool(double)>& ok, int iters = 100) {
  while (iters--) {
    double mid = (lo + hi) / 2;
    if (ok(mid)) hi = mid; else lo = mid;
  }
  return hi;
}

// first index i with a[i] >= key (== std::lower_bound), or a.size() if none.
int lowerBoundIdx(const vector<long long>& a, long long key) {
  int lo = 0, hi = a.size();
  while (lo < hi) { int mid = (lo + hi) / 2; if (a[mid] >= key) hi = mid; else lo = mid + 1; }
  return lo;
}

// first index i with a[i] > key (== std::upper_bound), or a.size() if none.
int upperBoundIdx(const vector<long long>& a, long long key) {
  int lo = 0, hi = a.size();
  while (lo < hi) { int mid = (lo + hi) / 2; if (a[mid] > key) hi = mid; else lo = mid + 1; }
  return lo;
}

// ARGMIN of a unimodal g over [lo,hi] (leftmost minimiser). Robust integer form:
// find the first x where g(x) <= g(x+1).
long long argminInt(long long lo, long long hi, const function<long long(long long)>& g) {
  long long l = lo, r = hi;
  while (l < r) { long long mid = l + (r - l) / 2; if (g(mid) <= g(mid + 1)) r = mid; else l = mid + 1; }
  return l;
}

// ARGMAX of a unimodal g over [lo,hi] (leftmost maximiser).
long long argmaxInt(long long lo, long long hi, const function<long long(long long)>& g) {
  long long l = lo, r = hi;
  while (l < r) { long long mid = l + (r - l) / 2; if (g(mid) >= g(mid + 1)) r = mid; else l = mid + 1; }
  return l;
}

// x minimising a unimodal real g over [lo,hi] (ternary search).
double ternaryMinReal(double lo, double hi, const function<double(double)>& g, int iters = 200) {
  while (iters--) {
    double m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3;
    if (g(m1) < g(m2)) hi = m2; else lo = m1;
  }
  return (lo + hi) / 2;
}

// floor(sqrt(n)) for n >= 0, exact (no floating point).
long long isqrt(long long n) {
  if (n < 0) return -1;
  long long lo = 0, hi = 3037000499LL;              // hi*hi just under LLONG_MAX
  while (lo < hi) { long long mid = lo + (hi - lo + 1) / 2; if (mid * mid <= n) lo = mid; else hi = mid - 1; }
  return lo;
}

// search a rotated sorted array of DISTINCT values; index of key, or -1.
int searchRotated(const vector<int>& a, int key) {
  int lo = 0, hi = (int)a.size() - 1;
  while (lo <= hi) {
    int mid = (lo + hi) / 2;
    if (a[mid] == key) return mid;
    if (a[lo] <= a[mid]) {                            // left half sorted
      if (a[lo] <= key && key < a[mid]) hi = mid - 1; else lo = mid + 1;
    } else {                                          // right half sorted
      if (a[mid] < key && key <= a[hi]) lo = mid + 1; else hi = mid - 1;
    }
  }
  return -1;
}

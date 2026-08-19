#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
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
#include <iomanip>

using namespace std;

const int INF = 2000000000;

// minimum number of coins to get a given sum
int minCoins(vector<int> coins, int sum) {
  vector<int> dp(sum + 1, INF);
  dp[0] = 0;
  for (int i = 1; i <= sum; i++) {
    for (int j = 0; j < coins.size(); j++) {
      if (coins[j] <= i) {
        dp[i] = min(dp[i], dp[i - coins[j]] + 1);
      }
    }
  }
  return dp[sum];
}

// also printing best solution
int minCoinWithPrint(vector<int> coins, int sum) {
  vector<int> dp(sum + 1, INF);
  vector<int> last(sum + 1, -1);
  dp[0] = 0;
  for (int i = 1; i <= sum; i++) {
    for (int j = 0; j < coins.size(); j++) {
      if (coins[j] <= i) {
        if (dp[i] > dp[i - coins[j]] + 1) {
          dp[i] = dp[i - coins[j]] + 1;
          last[i] = j;
        }
      }
    }
  }
  int cur = sum;
  while (cur > 0) {
    cout << coins[last[cur]] << " ";
    cur -= coins[last[cur]];
  }
  cout << endl;
  return dp[sum];
}

// all sums possible
vector<int> allSumsPossible(vector<int> coins) {
  vector<int> dp(100000, INF);
  dp[0] = 0;
  for (int i = 1; i < dp.size(); i++) {
    for (int j = 0; j < coins.size(); j++) {
      if (coins[j] <= i) {
        dp[i] = min(dp[i], dp[i - coins[j]] + 1);
      }
    }
  }
  vector<int> res;
  for (int i = 0; i < dp.size(); i++) {
    if (dp[i] != INF) {
      res.push_back(i);
    }
  }
  return res;
}

// number of ways to get sum
int numWays(vector<int> coins, int sum) {
  vector<int> dp(sum + 1, 0);
  dp[0] = 1;
  for (int x = 1; x <= sum; x++) {
    for (auto c : coins) {
      if (x - c >= 0) {
        dp[x] += dp[x - c];
      }
    }
  }
  return dp[sum];
}

// LIS
int LIS(vector<int> arr) {
  // note: for nondecreasing, switch to upper bound
  int k = 0;
  vector<int> L(arr.size(), 0);
  for (int i = 0; i < arr.size(); ++i) {
    int pos = lower_bound(L.begin(), L.begin()+k, arr[i]) - L.begin();
    L[pos] = arr[i];
    if (pos == k) {
      k++;
    } 
  }
  return k;
}

vector<int> getLIS(vector<int> arr) {
  vector<int> dp(arr.size(), 1);
  vector<int> pipeline(arr.size(), -1);
  for (int k = 0; k < arr.size(); k++) {
    dp[k] = 1;
    for (int i = 0; i < k; i++) {
      if (arr[i] < arr[k]) {
        if (dp[k] <= dp[i] + 1) {
          dp[k] = dp[i] + 1;
          pipeline[k] = i;
        }
      }
    }
  }
  int maxSpot = max_element(dp.begin(), dp.end()) - dp.begin();
  vector<int> res;
  while (maxSpot != -1) {
    res.push_back(arr[maxSpot]);
    maxSpot = pipeline[maxSpot];
  }
  reverse(res.begin(), res.end());
  return res;
}

// Maximum Subarray Sum
int maxSubarraySum(vector<int> arr) {
  vector<int> dp(arr.size(), 0);
  dp[0] = arr[0];
  for (int i = 1; i < arr.size(); i++) {
    dp[i] = max(dp[i - 1] + arr[i], arr[i]);
  }
  return *max_element(dp.begin(), dp.end());
}

void knapsack(int m, vector<pair<int, int>> items) {
  vector<int> dp(m+1, 0);
  vector<vector<int>> included;
  included.push_back({});
  for (int i = 1; i <= m+1; i++) {
    dp[i] = dp[i-1];
    included.push_back(included[i-1]);
    for (int j = 0; j < items.size(); j++) {
      if (i - items[j].first >= 0) {
        int i1 = dp[i];
        int i2 = dp[i-items[j].first] + items[j].second;
        if (i1 < i2) {
          dp[i] = i2;
          included[i] = included[i-items[j].first];
          included[i].push_back(j+1);
        }
      }
    }
  }
  cout << dp[m] << endl;
  vector<int> ans = included[m];
}

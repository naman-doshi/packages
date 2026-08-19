#include <bits/stdc++.h>
using namespace std;
#define int long long
#define ll long long
 
 
const int INF = 9000372036854775800;
 
#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define trav(a, x) for(auto& a : x)
#define all(x) x.begin(), x.end()
#define print(arr) for (auto i : arr) cout << i << " "; cout << endl;
typedef vector<vector<int>> vvi;
typedef vector<pair<int, int>> vpi;
typedef vector<int> vi;
typedef pair<int, int> pii;

// inversion counting

int mergeAndCount(vector<int>& arr, int l, int m, int r) {
    int count = 0;
    int i = l, j = m, k = 0;
    vector<int> temp(r-l+1);
    while (i < m && j <= r) {
        if (arr[i] <= arr[j]) {
            temp[k++] = arr[i++];
        } else {
            temp[k++] = arr[j++];
            count += m - i; 
        }
    }
    while (i < m) temp[k++] = arr[i++];
    while (j <= r) temp[k++] = arr[j++];
    for (i = l, k = 0; i <= r;) arr[i++] = temp[k++];
    return count;
}

int countInversions(vector<int>& arr, int l, int r) {
    if (l >= r) return 0;
    int m = l + (r - l) / 2;
    int count = 0;
    count += countInversions(arr, l, m);
    count += countInversions(arr, m + 1, r);
    count += mergeAndCount(arr, l, m + 1, r);
    return count;
}

int countInversions(vector<int>& arr) {
    return countInversions(arr, 0, arr.size() - 1);
}
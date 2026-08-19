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

vector<int> succ;
vvi precomp;

// U is max depth
void precompute(int U) {
  int logu = log2(U) + 1;
  precomp.resize(succ.size(), vector<int>(logu, -1));
  rep(i, 0, succ.size()) {
    precomp[i][0] = succ[i];
  }
  rep(j, 1, logu) {
    rep(i, 0, succ.size()) {
      if (precomp[i][j-1] != -1) {
        precomp[i][j] = precomp[precomp[i][j-1]][j-1];
      }
    }
  }
}

int successorK(int cur, int k) {
  int logk = log2(k) + 1;
  for (int j = logk; j >= 0; j--) {
    if (k & (1 << j)) {
      cur = precomp[cur][j];
    }
  }
  return cur;
}

int cycleLength(int x) {
  int a = succ[x];
  int b = succ[succ[x]];
  while (a != b) {
    a = succ[a];
    b = succ[succ[b]];
  }
  a = x;
  while (a != b) {
    a = succ[a];
    b = succ[b];
  }
  int first = a;
  b = succ[a];
  int length = 1;
  while (a != b) {
    b = succ[b];
    length++;
  }
  return length;
}



int main() {
  succ[1] = 3;
  succ[2] = 5;
  succ[3] = 7;
  succ[4] = 6;
  succ[5] = 2;
  succ[6] = 2;
  succ[7] = 1;
  succ[8] = 6;
  succ[9] = 3;
  precompute(15);
  cout << successorK(4, 11) << endl;
}
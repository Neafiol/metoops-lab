#include <bits/stdc++.h>
using namespace std;

vector <double> gauss(vector <double> const& b, vector <vector <double> > const& c) {
    int n = b.size();
    vector <vector <double> > arr = c;
    for (int i = 0; i < n; i++) {
        arr[i].push_back(b[i]);
    }

    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            if (arr[j][i] != 0) {
                swap(arr[i], arr[j]);
            }
        }
        if (arr[i][i] == 0) continue;
        double tmp = 1 / arr[i][i];
        for (int j = n; j >= i; j--)
            arr[i][j] *= tmp;
        for (int j = i + 1; j < n; j++) {
            tmp = arr[j][i];
            for (int k = n; k >= i; k--)
                arr[j][k] -= tmp * arr[i][k];
        }
    }
    vector <double> ans(n);
    ans[n - 1] = arr[n - 1][n];
    for (int i = n - 2; i >= 0; i--) {
        double sum = arr[i][n];
        for (int j = i + 1; j < n; j++)
            sum -= arr[i][j] * ans[j];
        ans[i] = sum;
    }
    return ans;
}

int main() {
    vector <double> ans = gauss({11, 1101}, {{1, 10}, {100, 1001}});
    for (double x : ans)
        cout << x << " ";
    cout << endl;
}

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <math.h>

vector <double> gaus(vector <int> b, vector <vector <int> > c) {
    int n = b.size();
    vector <vector <double> > arr(n, vector <double> (n + 1));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            if (j + 1 <= n) arr[i][j] = c[i][j];
            else arr[i][j] = b[i];
        }
    }
 
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            if (arr[j][i] != 0) {
                swap(arr[i], arr[j]);
            }
        }
        double tmp = arr[i][i];
        if (!tmp) continue;
        for (int j = n; j >= i; j--) arr[i][j] /= tmp;
        for (int j = i + 1; j < n; j++) {
            tmp = arr[j][i];
            for (int k = n; k >= i; k--) arr[j][k] -= tmp * arr[i][k];
        }
    }
    vector <double> ans(n);
    ans[n - 1] = arr[n - 1][n];
    for (int i = n - 2; i >= 0; i--) {
        double sum = arr[i][n];
        for (int j = i + 1; j < n; j++) sum -= arr[i][j] * ans[j];
        ans[i] = sum;
    }
    return ans;
}

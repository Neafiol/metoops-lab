#pragma once

#include <iostream>
#include <cmath>

using namespace std;

void print(vector<vector<double>> const& matrix) {
    cout.precision(15);
    for (auto& v : matrix) {
        for (auto &x : v) {
            cout << x << ' ';
        }
        cout << '\n';
    }
}

vector<double> multiply(vector<vector<double>> const& m, vector<double> const& v) {
    auto size = v.size();
    vector<double> res(size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            res[i] += m[i][j] * v[j];
        }
    }
    return res;
}

double absoluteError(vector<double>& x, vector<double>& y) {
    double norm = 0;
    for (int i = 0; i < x.size(); ++i) {
        double diff = x[i] - y[i];
        norm += diff * diff;
    }
    return sqrt(norm);
}

double relativeError(vector<double>& x, vector<double>& y) {
    double norm1 = absoluteError(x, y);
    double norm2 = 0;
    for (double t : x) {
        norm2 += t * t;
    }
    return norm1 / sqrt(norm2);
}
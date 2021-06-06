#include "LUdecomposition.h"
#include "util.h"
#include <random>

using namespace std;

vector<vector<double>> getHilbertMatrix(int size) {
    vector<vector<double>> matrix(size, vector<double>(size));
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            matrix[i][j] = 1.0 / (i + j + 1);
        }
    }
    return matrix;
}

vector<vector<double>> getRandomMatrix(int size, int k) {
    vector<vector<double>> matrix(size, vector<double>(size));

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i != j) {
                matrix[i][j] = -rand() % 5;
                matrix[i][i] -= matrix[i][j];
            }
        }
    }
    matrix[0][0] += pow(0.1, k);
    return matrix;
}

vector<double> getRange(int from, int count) {
    vector<double> range(count);
    for (int i = 0; i < count; ++i) {
        range[i] = from + i;
    }
    return range;
}

string getName(int size) {
    return "hilbert" + to_string(size) + "_";
}

string getName2(int size, int k) {
    return "random" + to_string(size) + "_" + to_string(k) + "_";
}

void generate(vector<int>& sizes) {
    for (int size : sizes) {
        auto A = getHilbertMatrix(size);
//        print(A);
        auto x = getRange(1, size);
        auto f = multiply(A, x);
        createTest(A, f, getName(size));
    }
}

void generate2(vector<int>& sizes, vector<int>& ks) {
    for (int size : sizes) {
        for (int k : ks) {
            auto A = getRandomMatrix(size, k);
    //        print(A);
            auto x = getRange(1, size);
            auto f = multiply(A, x);
            createTest(A, f, getName2(size, k));
        }
    }
}

void complete(vector<int>& sizes) {
    for (int size : sizes) {
        completeTest(1, getName(size));
    }
}

void complete2(vector<int>& sizes, vector<int>& ks) {
    for (int size : sizes) {
        for (int k : ks) {
            cout << size << " " << k << endl;
            completeTest(1, getName2(size, k));
        }
    }
}

void printResultTable(vector<int>& sizes) {
    for (int size : sizes) {
        auto res = readResult(1, size, getName(size));
        auto x = getRange(1, size);
        printf("%4d %15f %f\n", size, absoluteError(x, res), relativeError(x, res));
    }
}

void printResultTable2(vector<int>& sizes, vector<int>& ks) {
    for (int size : sizes) {
        for (int k : ks) {
            auto res = readResult(1, size, getName2(size, k));
            auto x = getRange(1, size);
            printf("%4d %2d %.15f %.15f\n", size, k, absoluteError(x, res), relativeError(x, res));
        }
    }
}

void test1() {
    vector<int> sizes = {10, 20, 35, 50, 75, 100, 150, 200, 350, 500, 750, 1000};
    generate(sizes);
    complete(sizes);
    printResultTable(sizes);
}

void test2() {
    vector<int> sizes = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
    vector<int> ks = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    generate2(sizes, ks);
    complete2(sizes, ks);
    printResultTable2(sizes, ks);
}

int main() {
    test2();
}

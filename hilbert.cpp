#include "LUdecomposition.h"
#include "util.h"

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

vector<double> getRange(int from, int count) {
    vector<double> range(count);
    for (int i = 0; i < count; ++i) {
        range[i] = from + i;
    }
    return range;
}

void test2() {
    vector<vector<double>> A = {{1, 1}, {2, 2}};
    vector<double> b = {1, 2};
//    createTest(A, b);
    completeTest(2);
}

string getName(int size) {
    return "hilbert" + to_string(size) + "_";
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

void complete(vector<int>& sizes) {
    for (int size : sizes) {
        completeTest(1, getName(size));
    }
}

void printResultTable(vector<int>& sizes) {
    for (int size : sizes) {
        auto res = readResult(1, size, getName(size));
        auto x = getRange(1, size);
        printf("%4d %15f %f\n", size, absoluteError(x, res), relativeError(x, res));
    }
}

int main() {
    vector<int> sizes = {10, 20, 35, 50, 75, 100, 150, 200, 350, 500, 750, 1000};
    generate(sizes);
    complete(sizes);
    printResultTable(sizes);
}

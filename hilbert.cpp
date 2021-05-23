#include "LUdecomposition.h"

using namespace std;

int main() {
    vector<vector<double>> A = {{3, 0,  13, 0, 0},
                                {0, 11, 1,  0, 0},
                                {2, 0,  25, 4, 0},
                                {0, 9,  5,  7, 0},
                                {0, 0,  0,  0, 19}};
    vector<double> b = {1, 2, 3, 4, 5};
//    createTest(A, b, "hilbert");
    completeTest(1, "hilbert");
}

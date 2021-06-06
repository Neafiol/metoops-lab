#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <iomanip>

const long double EPS = 1e-8;

using Vector = std::vector <long double>;
using Matrix = std::vector <Vector>;

// Output overload
std::ostream & operator<<(std::ostream & out, const Vector & v) {
    out << "{";
    for (int i = 0; i < v.size() - 1; i++) {
        out << v[i] << ", ";
    }
    out << v.back() << "}";
    return out;
}

std::ostream & operator<<(std::ostream & out, const Matrix & matrix) {
    for (const auto & rows : matrix) {
        out << rows << "\n";
    }
    return out;
}

Vector operator + (const Vector & a, const Vector & b) {
    Vector c = a;
    for (int i = 0; i < c.size(); i++) {
        c[i] += b[i];
    }
    return c;
}


// Unary minus for vector
Vector operator-(const std::vector<long double> & v) {
    auto y = v;
    for (int i = 0; i < y.size(); i++) {
        y[i] = - y[i];
    }
    return y;
}

// Scalar product
long double operator * (const Vector & a, const Vector & b) {
    long double result = 0;
    for (int i = 0; i < a.size(); i++) {
        result += a[i] * b[i];
    }
    return result;
}

// Gaus xD
Vector gauss(Vector const& b, std::vector <Vector > const& c) {
    int n = b.size();
    std::vector <Vector > arr = c;
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
        long double tmp = 1 / arr[i][i];
        for (int j = n; j >= i; j--)
            arr[i][j] *= tmp;
        for (int j = i + 1; j < n; j++) {
            tmp = arr[j][i];
            for (int k = n; k >= i; k--)
                arr[j][k] -= tmp * arr[i][k];
        }
    }
    Vector ans(n);
    ans[n - 1] = arr[n - 1][n];
    for (int i = n - 2; i >= 0; i--) {
        long double sum = arr[i][n];
        for (int j = i + 1; j < n; j++)
            sum -= arr[i][j] * ans[j];
        ans[i] = sum;
    }
    return ans;
}

// Replace this if you need to. However, make sure your function has non-negative Hessian, ok?
long double f (const Vector & x) {
    return 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1 - x[0]) * (1 - x[0]);
}

long double fxs (const Vector & x, const Vector & s, const long double l) {
    auto x1 = x;
    for (int i = 0; i < s.size(); i++) {
        x1[i] += s[i] * l;
    }
    return f(x1);
}


// Returns both Guessian and Gradient
// Uses + EPS, defined in this file
std::pair<std::vector<long double>, Matrix>
    grad(
        const Vector & x,
        const long double fx) {
    Vector g(x.size(), 0.0);
    Vector F(x.size(), 0.0);
    Matrix H(x.size(), Vector (x.size(), 0.0));
    auto n = x.size();
    auto y = x;
    auto z = x;
    long double a = 1 / (2 * EPS);
    long double b = 1 / (EPS * EPS);
    for (int i = 0; i < n; i++) {
        y[i] = x[i] + EPS;
        auto fy = f(y);
        z[i] = x[i] - EPS;
        auto fz = f(z);
        g[i] = a * (fy - fz);
        H[i][i] = b * (fy - 2 * fx + fz);
        y[i] = x[i];
        z[i] = x[i];
        F[i] = fy;
    }
    for (int i = 0; i < n - 1; i++) {
        y[i] = x[i] + EPS;
        for (int j = i + 1; j < n; j++) {
            y[j] = x[j] + EPS;
            auto fy = f(y);
            H[j][i] = H[i][j] = b * (fy - F[i] - F[j] + fx);
            y[j] = x[j];
        }
        y[i] = x[i];
    }
    // std::cout << "Hessian\n" << H << "\n";
    std::cout << "Gradient\n" << g << "\n";
    return {g, H};
}

// || V ||:
long double norm (const Vector & x) {
    long double accumulate = 0;
    for (const auto e : x) {
        accumulate += e * e;
    }
    return sqrtl(accumulate);
}



void newton_minimization(Vector & x) {
    Vector s;
    do {
        const auto & [g, H] = grad(x, f(x));
        auto negate_g = -g;
        s = gauss(negate_g, H);
        for (int i = 0 ; i < x.size(); i++) {
            x[i] += s[i];
        }
        std::cout << "New X\n" << x << "\n";
    } while (norm(s) > EPS);
}

constexpr long double golden_ratio = (1.0 + sqrt(5.0)) / 2.0;

long double onedim_search(const Vector & x, const Vector & s, long double l = -1000000, long double r = 1000000) {
    auto resphi =  2 - golden_ratio;
    auto x1 = l + resphi * (r - l);
    auto x2 = r - resphi * (r - l);
    auto f1 = fxs(x, s, x1);
    auto f2 = fxs(x, s, x2);
    do {
        if (f1 < f2) {
            r = x2;
            x2 = x1;
            f2 = f1;
            x1 = l + resphi * (r - l);
            f1 = fxs(x, s, x1);
        } else {
            l = x1;
            x1 = x2;
            f1 = f2;
            x2 = r - resphi * (r - l);
            f2 = fxs(x, s, x2);
        }
    } while (fabsl(r - l) > EPS);
    return (r + l) / 2;
}

void newton_minimization_onedim_search(Vector & x) {
    Vector s;
    do {
        const auto & [g, H] = grad(x, f(x));
        auto negate_g = -g;
        s = gauss(negate_g, H);
        auto lmbd = onedim_search(x, s);
        for (int i = 0 ; i < x.size(); i++) {
            x[i] += lmbd * s[i];
        }
        std::cout << "New X\n" << x << "\n";
    } while (norm(s) > EPS);
}

void newton_minimization_direct_downfall(Vector & x) {
    const auto & [ext, _] = grad(x, f(x));
    auto d = -ext;
    auto r = onedim_search(x, d);
    Vector s(x.size());
    Vector cur;
    for (int i = 0; i < d.size(); i++) {
        s[i] = d[i] * r;
    }
    for (int i = 0; i < x.size(); i++) {
        x[i] += s[i];
    }
    std::cout << "New X:\n" << x << "\n";
    do {
        const auto & [g, H] = grad(x, f(x));
        auto negate_g = -g;
        s = gauss(negate_g, H);
        if ((s * g) < 0) {
            d = s;
        } else {
            d = negate_g;
        }
        r = onedim_search(x, d);
        for (int i = 0; i < s.size(); i++) {
            s[i] = r * d[i];
        }
        std::cout << "New X:\n";
        for (int i = 0; i < x.size(); i++) {
            x[i] += s[i];
        }
        std::cout << x << "\n";
    } while(norm(s) > EPS);
}

int main() {
    std::cout << std::fixed << std::setprecision(8);
    Vector x = {-12, 1};
    Vector y = {-12, 1};
    Vector z = {-12, 1};
    newton_minimization(x);
    std::cout << "==============\n";
    newton_minimization_onedim_search(y);
    std::cout << "==============\n";
    newton_minimization_direct_downfall(z);
    std::cout << "==============\n";
    return 0;
}

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <functional>
#include <ctime>

#define MAX_ITER 1000

using namespace std;

using vect_t = vector<double>;
using matr_t = vector<vect_t>;

ostream& operator<<(ostream& out, vect_t const& v) {
    int n = v.size();
    for (int i = 0; i < n; i++) {
        cout << v[i];
        if (i < n - 1) {
            cout << ' ';
        }
    }
    return out;
}

ostream& operator<<(ostream& out, matr_t const& v) {
    int n = v.size();
    for (int i = 0; i < n; i++) {
        cout << v[i] << endl;
    }
    return out;
}

double dotProduct(vect_t const& v, vect_t const& u) {
    auto n = v.size();
    double x = 0;
    for (int i = 0; i < n; i++) {
        x += v[i] * u[i];
    }
    return x;
}

double norm2(vect_t const& v) {
    return dotProduct(v, v);
}

matr_t idenityMatrix(int n) {
    matr_t v(n, vect_t(n));
    for (int i = 0; i < n; i++) {
        v[i][i] = 1;
    }
    return v;
}

vect_t operator*(double x, vect_t const& v) {
    auto n = v.size();
    vect_t u = v;
    for (int i = 0; i < n; i++) {
        u[i] *= x;
    }
    return u;
}

vect_t operator-(vect_t const& v) {
    return -1 * v;
}

vect_t operator/(vect_t const& v, double x) {
    return 1 / x * v;
}

matr_t operator*(vect_t const& v, vect_t const& u) {
    auto n = v.size();
    matr_t g(n, vect_t(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            g[i][j] = v[i] * u[j];
        }
    }
    return g;
}

vect_t operator*(matr_t const& g, vect_t const& u) {
    auto n = g.size();
    vect_t x(u.size());
    for (int i = 0; i < n; i++) {
        x[i] = dotProduct(g[i], u);
    }
    return x;
}

vect_t operator+=(vect_t& v, vect_t const& u) {
    auto n = v.size();
    for (int i = 0; i < n; i++) {
        v[i] += u[i];
    }
    return v;
}

vect_t operator+(vect_t const& v, vect_t const& u) {
    auto w = v;
    return w += u;
}

vect_t operator-(vect_t const& v, vect_t const& u) {
    auto n = v.size();
    vect_t w(n);
    for (int i = 0; i < n; i++) {
        w[i] = v[i] - u[i];
    }
    return w;
}

matr_t operator+(matr_t const& v, matr_t const& u) {
    auto n = v.size();
    matr_t w(n);
    for (int i = 0; i < n; i++) {
        w[i] = v[i] + u[i];
    }
    return w;
}

matr_t operator-(matr_t const& v, matr_t const& u) {
    auto n = v.size();
    matr_t w(n);
    for (int i = 0; i < n; i++) {
        w[i] = v[i] - u[i];
    }
    return w;
}

matr_t transpose(matr_t const& v) {
    auto n = v.size();
    matr_t w(n, vect_t(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            w[i][j] = v[j][i];
        }
    }
    return w;
}

double getParabolaMin(double x1, double y1, double x2, double y2, double x3, double y3) {
    double z1 = x1 * (y3 - y2);
    double z2 = x2 * (y1 - y3);
    double z3 = x3 * (y2 - y1);
    double a = z1 + z2 + z3;
    double b = x1 * z1 + x2 * z2 + x3 * z3;
    return b / (2 * a);
}

template<typename Func>
double brent(double a, double b, double eps, Func&& func) {
    const double K = 0.3819660112501051518;
    double x1, x2, x3;
    double f1, f2, f3;
    x1 = x2 = x3 = (a + b) / 2;
    f1 = f2 = f3 = func(x1);
    double step, prev_step;
    step = prev_step = b - a;

    while (step > eps) {
        double prev2_step = prev_step;
        prev_step = step;
        double u;

        if ((x1 != x2 && x2 != x3 && x3 != x1) &&
            (f1 != f2 && f2 != f3 && f3 != f1) &&
            (u = getParabolaMin(x1, f1, x2, f2, x3, f3),
                (u >= a + eps && u <= b - eps && 2 * abs(u - x1) < prev2_step)))
        {
            step = abs(u - x1);
        } else if (2 * x1 < b - a) {
            step = b - x1;
            u = x1 + K * step;
        } else {
            step = x1 - a;
            u = x1 - K * step;
        }

        double fu = func(u);
        if (fu <= f1) {
            (u >= x1 ? a : b) = x1;
            x3 = x2;
            x2 = x1;
            x1 = u;
            f3 = f2;
            f2 = f1;
            f1 = fu;
        } else {
            (u >= x1 ? b : a) = u;
            if (fu <= f2 || x2 == x1) {
                x3 = x2;
                x2 = u;
                f3 = f2;
                f2 = fu;
            } else if (fu <= f3 || x3 == x1 || x3 == x2) {
                x3 = u;
                f3 = fu;
            }
        }
    }
    return x1;
}

template<typename Vect, typename Matr, typename Func, typename GradFunc, typename NextG>
Vect quasiNewton(Vect&& x, double eps, double sig, Func&& f, GradFunc&& grad, Matr&& G0, NextG&& nextG) {
    const double minAlpha = 0;
    const double maxAlpha = 1e6;

    cout << "init: " << x << endl;
    auto G = G0;
    int n = x.size();
    auto w = -grad(x);
    for (int i = 1; i <= MAX_ITER; i++) {
        auto p = G * w;
        auto alpha = brent(minAlpha, maxAlpha, sig, [&](double alpha){ return f(x + alpha * p); });
        auto dx = alpha * p;
        x += dx;
        // cout << x << endl;
        if (norm2(dx) < eps * eps) {
            cout << "iter: " << i << endl;
            break;
        }
        auto w1 = -grad(x);
        if (i % 50 == 0) {
            G = G0;
        } else {
            G = nextG(G, dx, w1 - w);
        }
        w = move(w1);

    }

    cout << "min: " << x << endl;
    cout << "val: " << f(x) << "\n\n";
    return x;
}

template<typename Vect, typename Matr, typename Func, typename GradFunc>
Vect BFS(Vect x, double eps, double sig, Func&& f, GradFunc&& grad, Matr G) {
    cout << "BFS" << endl;
    return quasiNewton(x, eps, sig, f, grad, G, [] (auto const& G, auto const& dx, auto const& dw) {
        auto z = G * dw;
        auto rho = dotProduct(z, dw);
        auto zp = z / rho;
        auto xxw = dx / dotProduct(dx, dw);
        auto r = zp - xxw;
        return G - xxw * dx - zp * z + rho * r * r;
    });
}

template<typename Vect, typename Matr, typename Func, typename GradFunc>
Vect powell(Vect x, double eps, double sig, Func&& f, GradFunc&& grad, Matr G) {
    cout << "powell" << endl;
    return quasiNewton(x, eps, sig, f, grad, G, [] (auto const& G, auto const& dx, auto const& dw) {
        auto dt = dx + G * dw;
        return G - dt * (dt / dotProduct(dw, dt));
    });
}

// 100*(x2-x1^2)^2+(1-x1)^2, min=(1, 1)
double func1(double x1, double x2) {
    double a1 = x2 - x1 * x1;
    double a2 = 1 - x1;
    return 100 * a1 * a1 + a2 * a2;
}

double func1v(vect_t const& v) {
    return func1(v[0], v[1]);
}

vect_t grad1(double x1, double x2) {
    return {x1*(400*x1*x1-400*x2+2)-2, -200*(x1*x1-x2)};
}

vect_t grad1v(vect_t const& v) {
    return grad1(v[0], v[1]);
}

// (x1^2+x2-11)^2+(x1+x2^2-7)^2, min=
// (3.0000000000000000, 2.0000000000000000)
// (-2.805118086952744, 3.1313125182505729),
// (-3.779310253377746, -3.283185991286169),
// (3.5844283403304917, -1.848126526964403),
double func2(double x1, double x2) {
    double a1 = x1 * x1 + x2 - 11;
    double a2 = x1 + x2 * x2 - 7;
    return a1 * a1 + a2 * a2;
}

double func2v(vect_t const& v) {
    return func2(v[0], v[1]);
}

vect_t grad2(double x1, double x2) {
    return {x1*(4*x1*x1-42)+x2*(4*x1+2*x2)-14,
            x1*(2*x1+4*x2)+x2*(4*x2*x2-26)-22};
}

vect_t grad2v(vect_t const& v) {
    return grad2(v[0], v[1]);
}

// (x1+10x2)^2+5(x3-x4)^2+(x2-2x3)^4+10(x1-x4)^4, min=(0, 0, 0, 0)
double func3(double x1, double x2, double x3, double x4) {
    double a1 = x1 + 10 * x2;
    double a2 = x3 - x4;
    double a3 = x2 - 2 * x3;
    double a4 = x1 - x4;
    a3 *= a3;
    a4 *= a4;
    return a1 * a1 + 5 * a2 * a2 + a3 * a3 + a4 * a4;
}

double func3v(vect_t const& v) {
    return func3(v[0], v[1], v[2], v[3]);
}

// 100-2/(1+((x1-1)/2)^2+((x2-1)/3)^2)-1/(1+((x1-2)/2)^2+((x2-1)/3)^2)
// min=(1.2916430315174929, 1)
double func4(double x1, double x2) {
    double a1 = (x1 - 1) / 2;
    double a2 = (x2 - 1) / 3;
    double a3 = (x1 - 2) / 2;
    a1 *= a1;
    a2 *= a2;
    a3 *= a3;
    return 100 - 2 / (1 + a1 + a2) - 1 / (1 + a3 + a2);
}

double func4v(vect_t const& v) {
    return func4(v[0], v[1]);
}

template<typename Func>
auto grad(Func&& f) {
    return [&] (vect_t const& v) {
        const double eps = 1e-6;
        const double sig = 1e-15;
        int n = v.size();
        auto w = v;
        vect_t u(n);
        for (int i = 0; i < n; i++) {
            auto dv = v[i] * eps + sig;
            w[i] = v[i] + dv;
            auto f1 = f(w);
            w[i] = v[i] - dv;
            auto f2 = f(w);
            w[i] = v[i];
            u[i] = (f1 - f2) / (2 * dv);
        }
        return u;
    };
}

void test0() {
    vect_t x = {14, -15};
    double eps = 1e-12;
    double sig = 1e-10;
    auto f = [] (vect_t const& v) {
        double x = v[0];
        double y = v[1];
        double a = x - 1;
        return a * a * a * a + 100 * y * y;
    };
    auto grad = [] (vect_t const& v) -> vect_t {
        double x = v[0];
        double y = v[1];
        double a = x - 1;
        return {4 * a * a * a, 200 * y};
    };
    matr_t G = idenityMatrix(x.size());
    BFS(x, eps, sig, f, grad, G);
    powell(x, eps, sig, f, grad, G);
}

void test1() {
    vect_t x = {0.9, 1.1};
    double eps = 1e-12;
    double sig = 1e-8;
    matr_t G = idenityMatrix(x.size());
    BFS(x, eps, sig, func1v, grad(func1v), G);
    powell(x, eps, sig, func1v, grad(func1v), G);
}

void test2() {
    vect_t x = {19, 15};
    double eps = 1e-12;
    double sig = 1e-8;
    matr_t G = idenityMatrix(x.size());
    auto run1 = [&] (vect_t const& x) {
        BFS(x, eps, sig, func2v, grad(func2v), G);
    };
    auto run2 = [&] (vect_t const& x) {
        powell(x, eps, sig, func2v, grad(func2v), G);
    };

    run1({-19, 15});
    run1({19, 15});
    run1({19, -15});
    run1({59, -15});

    run2({-19, 15});
    run2({19, 15});
    run2({19, -15});
    run2({59, -15});
}

void test3() {
    vect_t x = {1, -1, 1, -1};
    double eps = 1e-12;
    double sig = 1e-8;
    matr_t G = idenityMatrix(x.size());

    BFS(x, eps, sig, func3v, grad(func3v), G);
    powell(x, eps, sig, func3v, grad(func3v), G);
}

void test4() {
    vect_t x = {1, -1};
    double eps = 1e-12;
    double sig = 1e-8;
    matr_t G = idenityMatrix(x.size());

    BFS(x, eps, sig, func4v, grad(func4v), G);
    powell(x, eps, sig, func4v, grad(func4v), G);
}

void gradTest() {
    cout << grad1(2, 2) << endl;
    cout << grad(func1v)({2, 2}) << endl;
    cout << grad2(-1, 2) << endl;
    cout << grad(func2v)({-1, 2}) << endl;
}

int main() {
    cout.precision(12);
    test4();
    // gradTest();
}

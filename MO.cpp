#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <math.h>
#define PI 3.14159265358979323846
#define PHI 1.61803398874989484820
#define TAU 0.61803398874989484820

double func(double x) {
    return - (3 * x * sin(0.75 * x)) + exp(-2 * x);
}

double fib(int x) {
    return (pow(PHI, x) - pow(1 - PHI, x)) / sqrt(5);
}

double sign(double x) {
    return (x > 0) - (x < 0);
}

struct result {
    double left;
    double right;
    double min;
    double min_value;
};

result getResult(double a, double b, double x, double fx) {
//    std::cout << "Final interval is : [ " << a << " ; " << b << " ]\n";
//    std::cout << "The min value is : " << fx << '\n';
    return {a, b, x, fx};
}

result getResult(double a, double b, double x) {
    return getResult(a, b, x, func(x));
}

result getResult(double a, double b) {
    return getResult(a, b, (a + b) / 2);
}

result bisection(double a, double b, double eps) {
    while (b - a > 2 * eps) {
        double x1 = (a + b - eps) / 2;
        double x2 = (a + b + eps) / 2;
        //std::cout << a << " & " << b << " & " << b - a << " & " << x1 << " & " << x2 << " & " << func(x1) << " & " << func(x2) << "\n";
        if (func(x1) <= func(x2)) {
            b = x2;
        } else {
            a = x1;
        }
    }
    return getResult(a, b);
}

result goldenRatio(double a, double b, double eps) {
    double x1 = a + (1 - TAU) * (b - a);
    double x2 = a + TAU * (b - a);
    while (b - a > 2 * eps) {
        //std::cout << a << " & " << b << " & " << b - a << " & " << x1 << " & " << x2 << " & " << func(x1) << " & " << func(x2) << "\n";
        if (func(x1) > func(x2)) {
            a = x1;
            x1 = x2;
            x2 = a + TAU * (b - a);
        } else {
            b = x2;
            x2 = x1;
            x1 = a + (1 - TAU) * (b - a);
        }
    }
    return getResult(a, b);
}

result fibonacci(double a, double b, double eps, int n) {
    double len = (b - a) * fib(n) / fib(n + 1) + (pow(-1, n) * eps) / fib(n + 1);
    double x2 = a + len;
    for (int k = 1; k < n; k++) {
        double x1 = a + (b - x2);
        //std::cout << a << " & " << b << " & " << b - a << " & " << x1 << " & " << x2 << " & " << func(x1) << " & " << func(x2) << "\n";
        double f1 = func(x1);
        double f2 = func(x2);
        if (x1 > x2 && f1 < f2) {
            a = x2;
            x2 = x1;
        } else if (x1 < x2 && f1 < f2) {
            b = x2;
            x2 = x1;
        } else if (x1 < x2 && f1 >= f2) {
            a = x1;
        } else {
            b = x1;
        }
    }
    return getResult(a, b);
}

result parabola(double x1, double x2, double x3, double eps) {
    double x = x1;
    double prev_x = x;
    double step = x3 - x1;
    while (step > eps) {
        double f1 = func(x1);
        double f2 = func(x2);
        double f3 = func(x3);
//        double a1 = f1;
        double a2 = (f2 - f1) / (x2 - x1);
        double a3 = (((f3 - f1) / (x3 - x1)) - ((f2 - f1) / (x2 - x1))) / (x3 - x2);
        x = 0.5 * (x1 + x2 - (a2 / a3));
        double f = func(x);
        step = abs(x - prev_x);
        //std::cout << x1 << " & " << x2 << " & " << x3 << " & " << step << " & " << a1 << " & " << a2 << " & " << a3 << " & " << x << " & " << f2 << " & " << f << "\n";
        if (x1 < x && x < x2 && f >= f2) {
            x1 = x;
        } else if (x1 < x && x < x2 && f < f2) {
            x3 = x2;
            x2 = x;
        } else if (x2 < x && x < x3 && f2 >= f) {
            x1 = x2;
            x2 = x;
        } else {
            x3 = x;
        }
        prev_x = x;
    }
    return getResult(x, x3, x);
}

double getParabolaMin(double x1, double y1, double x2, double y2, double x3, double y3) {
    double a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2));
    double b = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3));
    return -b / (2 * a);
}

result brent(double a, double b, double eps) {
    double x1, x2, x3;
    double f1, f2, f3;
    x1 = x2 = x3 = (a + b) / 2;
    f1 = f2 = f3 = func(x1);
    double step, prev_step;
    step = prev_step = b - a;

    while (step > eps) {
//        std::cout << a << " & " << b << " & " << x1 << " & " << x2 << " & " << x3 << '\n';
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
            u = x1 + (1 - TAU) * step;
        } else {
            step = x1 - a;
            u = x1 - (1 - TAU) * step;
        }
        if (abs(u - x1) < eps) {
            u = x1 + sign(u - x1) * eps;
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
    return getResult(a, b, x1);
}

int main() {
    std::cout << std::fixed;
    std::cout << std::setprecision(10);

    double eps = 0.00001;
    double a = 0;
    double b = 2 * PI;

    bisection(a, b, eps);
    goldenRatio(a, b, eps);
    fibonacci(a, b, eps, 30);

    double n = (b - a) / 20;
    double i;
    for (i = a + n; i < b - n; i += n) {
        if (func(a) >= func(i) && func(i) <= func(b)) break;
    }

    parabola(a, i, b, eps);
    brent(a, b, eps);
    return 0;
}

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

void clear(double arr[100]) {
    for (int i = 0; i < 100; i++) {
        arr[i] = 0;
    }
}

double * bisection(double eps, int n) {
    static double res[100];
    clear(res);
    int i = 2;
    double a = 0;
    double b = 2 * PI;

    while (b - a > 2 * eps) {
        double x1 = (a + b - eps) / 2;
        double x2 = (a + b + eps) / 2;
        if (func(x1) <= func(x2)) {
            b = x2;
            res[i++] = x2;
        } else {
            a = x1;
            res[i++] = x1;
        }
    }
    res[0] = (b - a) / 2;
    res[1] = i - 2;
    return res;
}

double * goldenRatio(double eps, int n) {
    static double res[100];
    clear(res);
    int i = 2;
    double a = 0;
    double b = 2 * PI;

    double x1 = a + (1 - TAU) * (b - a);
    double x2 = a + TAU * (b - a);

    while (b - a > 2 * eps) {
        if (func(x1) > func(x2)) {
            a = x1;
            x1 = x2;
            x2 = a + TAU * (b - a);
            res[i++] = x2;
        } else {
            b = x2;
            x2 = x1;
            x1 = a + (1 - TAU) * (b - a);
            res[i++] = x1;
        }
    }
    res[0] = (b - a) / 2;
    res[1] = i - 2;
    return res;
}

int fib_index(double x) {
    double fib1 = 1;
    double fib2 = 1;
    int n = 2;

    while (fib2 < x) {
        double fib3 = fib1 + fib2;
        fib1 = fib2;
        fib2 = fib3;
        ++n;
    }
    return n;
}

double * fibonacci(double eps, int n2) {
    static double res[100];
    clear(res);

    int i = 2;
    double a = 0;
    double b = 2 * PI;

    int n = fib_index((b - a) / eps) - 2;
    double x1 = a + (b - a) * fib(n) / fib(n + 2);
    res[i++] = x1;

    for (int k = 1; k < n; k++) {
        double x2 = a + (b - x1);
        res[i++] = x2;
        double f2 = func(x2);
        double f1 = func(x1);
        if (x2 > x1 && f2 < f1) {
            a = x1;
            x1 = x2;
        } else if (x2 < x1 && f2 < f1) {
            b = x1;
            x1 = x2;
        } else if (x2 < x1 && f2 >= f1) {
            a = x2;
        } else {
            b = x2;
        }
    }
    res[0] = (b - a) / 2;
    res[1] = n+1;
    return res;
}

double * parabola(double eps, int n) {
    static double res[100];
    clear(res);
    int i = 2;
    double x1 = 0;
    double x2;
    double x3 = 2 * PI;

    double m = (x3 - x1) / 20;
    for (x2 = x1 + m; x2 < x3 - m; x2 += m) {
        if (func(x1) >= func(x2) && func(x2) <= func(x3)) break;
    }

    double x = x1;
    double prev_x = x;
    double step = x3 - x1;
    while (step > eps) {
        double f1 = func(x1);
        double f2 = func(x2);
        double f3 = func(x3);
        double a2 = (f2 - f1) / (x2 - x1);
        double a3 = (((f3 - f1) / (x3 - x1)) - ((f2 - f1) / (x2 - x1))) / (x3 - x2);
        x = 0.5 * (x1 + x2 - (a2 / a3));
        double f = func(x);
        step = fabs(x - prev_x);
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
        res[i++] = x;
        prev_x = x;
    }
    res[0] = step;
    res[1] = i - 2;
    return res;
}

double getParabolaMin(double x1, double y1, double x2, double y2, double x3, double y3) {
    double a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2));
    double b = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3));
    return -b / (2 * a);
}

double * brent(double eps, int n) {
    static double res[100];
    clear(res);
    int i = 2;
    double a = 0;
    double b = 2 * PI;

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
                (u >= a + eps && u <= b - eps && 2 * fabs(u - x1) < prev2_step)))
        {
            step = fabs(u - x1);
        } else if (2 * x1 < b - a) {
            step = b - x1;
            u = x1 + (1 - TAU) * step;
        } else {
            step = x1 - a;
            u = x1 - (1 - TAU) * step;
        }
        if (fabs(u - x1) < eps) {
            u = x1 + sign(u - x1) * eps;
        }

        double fu = func(u);
        if (fu <= f1) {
            if(u >= x1) {
            a=x1;
            }
            else {
                b=x1;
            }
            x3 = x2;
            x2 = x1;
            x1 = u;
            f3 = f2;
            f2 = f1;
            f1 = fu;
        } else {
            if (u >= x1) {
                b = u;
            }
            else {
                a=u;
            }
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
        res[i++] = u;
    }
    res[0] = step;
    res[1] = i - 2;
    return res;
}

//int main() {
//    std::cout << std::fixed;
//    std::cout << std::setprecision(10);
//    double *arr = fibonacci(0.0001, 30);
//    for (int i = 0; i < 100; i++) {
//        std::cout << arr[i] << '\n';
//    }
//    return 0;
//}
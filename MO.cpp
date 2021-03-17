#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <math.h>
#define N 12
#define PI 3.14159265358979323846
#define PHI 1.61803398874989484820
#define EPSILON 0.1

double func(double x) {
    return - (3 * x * sin(0.75 * x)) + exp(-2 * x);
}

double f(int x) {
    return (pow(PHI, x) - pow(1 - PHI, x)) / sqrt(5);
}

void bisection(double a, double b, double l) {
    double p = a, q = b;
    while ((b - a) >= l) {
        p = (a + b) / 2 - EPSILON;
        q = (a + b) / 2 + EPSILON;
        //std::cout << a << " & " << b << " & " << b - a << " & " << p << " & " << q << " & " << func(p) << " & " << func(q) << "\n";
        if (func(p) < func(q))
            b = q;
        else
            a = p;
    }
    std::cout << "\nFinal interval is : [ " << a << " ; " << b << " ]";
    std::cout << "\nThe min value is : " << func(b);
}

void goldenRatio(double a, double b, double l) {
    double p = a + (1 - 0.618) * (b - a), q = a + 0.618 * (b - a);
    while ((b - a) >= l) {
        //std::cout << a << " & " << b << " & " << b - a << " & " << p << " & " << q << " & " << func(p) << " & " << func(q) << "\n";
        if (func(p) > func(q)) {
            a = p;
            p = q;
            q = a + 0.618 * (b - a);
        }
        else {
            b = q;
            q = p;
            p = a + (1 - 0.618) * (b - a);
        }
    }
    std::cout << "\nFinal interval is : [ " << a << " ; " << b << " ]";
    std::cout << "\nThe min value is : " << func(b);
}

void fibonacci(double a, double b, double l) {
    int k = 1;
    l = (l * f(N)) / f(N + 1) + (pow(-1, N) * EPSILON) / f(N + 1);
    double p, q = a + l;
    while (k < N) {
        p = a + (b - q);
        //std::cout << a << " & " << b << " & " << b - a << " & " << p << " & " << q << " & " << func(p) << " & " << func(q) << "\n";
        if (p > q && func(p) < func(q)) {
            a = q;
            q = p;
        }
        else if (p < q && func(p) < func(q)) {
            b = q;
            q = p;
        }
        else if (p < q && func(p) >= func(q)) {
            a = p;
        }
        else {
            b = p;
        }
        k++;
    }
    std::cout << "\nFinal interval is : [ " << a << " ; " << b << " ]";
    std::cout << "\nThe min value is : " << func(b);
}

void parabola(double x1, double x2, double x3) {
    double a1, a2, a3, prev = 0, x, d = x3 - x1;
    while (d > EPSILON) {
        a1 = func(x1);
        a2 = (func(x2) - func(x1)) / (x2 - x1);
        a3 = (((func(x3) - func(x1)) / (x3 - x1)) - ((func(x2) - func(x1)) / (x2 - x1))) / (x3 - x2);
        x = 0.5 * (x1 + x2 - (a2 / a3));
        d = fabs(x - prev);
        //std::cout << x1 << " & " << x2 << " & " << x3 << " & " << d << " & " << a1 << " & " << a2 << " & " << a3 << " & " << x << " & " << func(x2) << " & " << func(x) << "\n";
        if (x1 < x && x < x2 && func(x) >= func(x2)) {
            x1 = x;
        }
        else if (x1 < x && x < x2 && func(x) < func(x2)) {
            x3 = x2;
            x2 = x;
        }
        else if (x2 < x && x < x3 && func(x2) >= func(x)) {
            x1 = x2;
            x2 = x;
        }
        else {
            x3 = x;
        }
        prev = x;
    }
    std::cout << "\nFinal interval is : [ " << x << " ; " << x3 << " ]";
    std::cout << "\nThe min value is : " << func(x);
}

void brent(double a, double b) {
    /*double k = (3 - sqrt(5) / 2), x = (a + b) / 2, w = x, v = x, d = b - a, prev = d;
    while (d > EPSILON) {
        a1 = func(x1);
        a2 = (func(x2) - func(x1)) / (x2 - x1);
        a3 = (((func(x3) - func(x1)) / (x3 - x1)) - ((func(x2) - func(x1)) / (x2 - x1))) / (x3 - x2);
        x = 0.5 * (x1 + x2 - (a2 / a3));
        d = fabs(x - prev);
        if (x1 < x && x < x2 && func(x) >= func(x2)) {
            x1 = x;
        }
        else if (x1 < x && x < x2 && func(x) < func(x2)) {
            x3 = x2;
            x2 = x;
        }
        else if (x2 < x && x < x3 && func(x2) >= func(x)) {
            x1 = x2;
            x2 = x;
        }
        else {
            x3 = x;
        }
        prev = x;
    }*/
    //std::cout << "\nFinal interval is : [ " << x << " ; " << x3 << " ]";
    //std::cout << "\nThe min value is : " << func(x);
}

int main() {
    std::cout << std::fixed;
    std::cout << std::setprecision(10);
    double a = 0, b = 2 * PI, l = 0.5;
    bisection(a, b, l);
    goldenRatio(a, b, l);
    fibonacci(a, b, b - a);
    double i, n = (b - a) / 20;
    for (i = a + n; i < b - n; i += n) {
        if (func(a) >= func(i) && func(i) <= func(b)) break;
    }
    parabola(a, i, b);
    return 0;
}
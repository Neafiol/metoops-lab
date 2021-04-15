#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846
#define TAU 0.61803398874989484820

typedef struct vector {
    double x, y;
} vector;

typedef double (*function)(double, double);

void printResult(vector* result);
double func1(double x, double y);
double func2(double x, double y);
double func3(double x, double y);
double func4(double x, double y);
double vectorFunc(vector v);
vector findGradient(vector* v);
vector nextValue(vector* v, vector* antiGrad, double lambda);
double lambdaFunc(vector* v, vector* antiGrad, double lambda);
double norm2(vector* v);
double norm(vector* v);
double* gradientDescent(double eps, int mode, int funcIndex);
vector gradientDescent1(vector v, double eps);
vector gradientDescent2(vector current, double eps);
double sign(double x);
double getParabolaMin(double x1, double y1, double x2, double y2, double x3, double y3);
double brent(double a, double b, double eps, vector* v, vector* antiGrad);

function func = func1;
function funcArr[] = {func1, func2, func3, func4};

int main() {
    vector v;
    double eps = 1e-4;
    v.x = 10;
    v.y = -10;
    vector result1 = gradientDescent1(v, eps);
    printResult(&result1);
    printf("\n");
    vector result2 = gradientDescent2(v, eps);
    printResult(&result2);
    return 0;
}

void printResult(vector* result) {
    printf("x = %0.17f, y = %0.17f\n", result->x, result->y);
    fflush(stdout);
}

double func1(double x, double y) {
    return 64 * x * x + 126 * x * y + 64 * y * y - 10 * x + 30 * y + 13;
}

double func2(double x, double y) {
    return (x * x + y * y);
}

double func3(double x, double y) {
    return x * x - y * y;
}

double func4(double x, double y) {
    return 64 * x * x + 64 * y * y - 10 * x + 30 * y + 13;
}

double vectorFunc(vector v) {
    return func(v.x, v.y);
}

vector findGradient(vector* v) {
    vector grad;
    double h = 0.1;
    double dx = (func(v->x + h, v->y) - func(v->x - h, v->y)) / (2 * h);
    double dy = (func(v->x, v->y + h) - func(v->x, v->y - h)) / (2 * h);
    grad.x = dx;
    grad.y = dy;
    return grad;
}

vector nextValue(vector* v, vector* antiGrad, double lambda) {
    vector w;
    w.x = v->x - lambda * antiGrad->x;
    w.y = v->y - lambda * antiGrad->y;
    return w;
}

double lambdaFunc(vector* v, vector* antiGrad, double lambda) {
    return vectorFunc(nextValue(v, antiGrad, lambda));
}

double norm2(vector* v) {
    return v->x * v->x + v->y * v->y;
}

double norm(vector* v) {
    return sqrt(norm2(v));
}

double* gradientDescent(double eps, int mode, int funcIndex) {
    static double data[100];
    return data;
}

vector gradientDescent1(vector current, double eps) {
    vector last;
    double currentValue = vectorFunc(current);
    double lastValue;
    double lambda = 1;
    for (int i = 0; i < 500; i++) {
        vector antiGrad = findGradient(&current);
        if (norm2(&antiGrad) < eps * eps) {
            break;
        }

        last = current;
        lastValue = currentValue;

        for (;;) {
            current = nextValue(&last, &antiGrad, lambda / norm(&antiGrad));
            currentValue = vectorFunc(current);
            if (currentValue <= lastValue) {
                break;
            }
            lambda *= 0.5;
        }
        if (current.x == last.x && current.y == last.y) {
            break;
        }
        // printResult(&current);
    };
    return current;
}

vector gradientDescent2(vector current, double eps) {
    vector last;
    double currentValue = vectorFunc(current);
    double lastValue;
    for (int i = 0; i < 100; i++) {
        vector antiGrad = findGradient(&current);
        if (norm2(&antiGrad) < eps * eps) {
            break;
        }

        last = current;
        lastValue = currentValue;

        double lambda = brent(-100, 100, eps, &current, &antiGrad);
        current = nextValue(&current, &antiGrad, lambda);
        currentValue = vectorFunc(current);
        if (current.x == last.x && current.y == last.y) {
            break;
        }
        // printResult(&current);
    };
    return current;
}

double sign(double x) {
    return (x > 0) - (x < 0);
}

double getParabolaMin(double x1, double y1, double x2, double y2, double x3, double y3) {
    double z1 = x1 * (y3 - y2);
    double z2 = x2 * (y1 - y3);
    double z3 = x3 * (y2 - y1);
    double a = z1 + z2 + z3;
    double b = x1 * z1 + x2 * z2 + x3 * z3;
    return b / (2 * a);
}

double brent(double a, double b, double eps, vector* v, vector* antiGrad) {
    double x1, x2, x3;
    double f1, f2, f3;
    x1 = x2 = x3 = (a + b) / 2;
    f1 = f2 = f3 = lambdaFunc(v, antiGrad, x1);
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
        }
        else if (2 * x1 < b - a) {
            step = b - x1;
            u = x1 + (1 - TAU) * step;
        }
        else {
            step = x1 - a;
            u = x1 - (1 - TAU) * step;
        }

        double fu = lambdaFunc(v, antiGrad, u);
        if (fu <= f1) {
            if (u >= x1) a = x1;
            else b = x1;
            x3 = x2;
            x2 = x1;
            x1 = u;
            f3 = f2;
            f2 = f1;
            f1 = fu;
        }
        else {
            if (u >= x1) b = u;
            else a = u;
            if (fu <= f2 || x2 == x1) {
                x3 = x2;
                x2 = u;
                f3 = f2;
                f2 = fu;
            }
            else if (fu <= f3 || x3 == x1 || x3 == x2) {
                x3 = u;
                f3 = fu;
            }
        }
    }
    return x1;
}
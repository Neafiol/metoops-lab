#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846
#define TAU 0.61803398874989484820

typedef struct vector {
    double x, y;
} vector;

void printResult(vector* result);
double func(double x, double y);
double func0(vector v);
vector findGradient(vector* v);
vector nextValue(vector* v, vector* antiGrad, double lambda);
vector gradientDescent(vector* v, double eps, int mode);
double brent(double a, double b, double eps, vector* v, vector* antiGrad);
double sign(double x);
double lambdaFunc(vector* v, vector* antiGrad, double lambda);
double getParabolaMin(double x1, double y1, double x2, double y2, double x3, double y3);

int main() {
    vector v;
    double eps;
    v.x = 10;
    v.y = -10;
    eps = 1e-4;
    vector result0 = gradientDescent(&v, eps, 0);
    printResult(&result0);
    vector result1 = gradientDescent(&v, eps, 1);
    printResult(&result1);
    return 0;
}

void printResult(vector* result) {
    printf("x = %0.6f, y = %0.6f\n", result->x, result->y);
    fflush(stdout);
}

double func(double x, double y) {
    // return x * x + y * y;
    // return x * x - y * y;
    return 64 * x * x + 126 * x * y + 64 * y * y - 10 * x + 30 * y + 13;
}

double func0(vector v) {
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

vector gradientDescent(vector* v, double eps, int mode) {
    vector current = *v;
    vector last;
    double currentValue = func0(current);
    double lastValue;
    double difference;
    do {
        last = current;
        lastValue = currentValue;
        vector antiGrad = findGradient(&current);
        double lambda = 0.001;
        switch (mode) {
            case 0:
                lambda = 0.001;
                break;
            case 1:
                lambda = brent(0.000001, 10, eps, v, &antiGrad);
                break;
            default:
                break;
        }
        current = nextValue(&current, &antiGrad, lambda);
        currentValue = func0(current);
        difference = fabs(currentValue - lastValue);
    } while (difference > eps);
    return current;
}

double lambdaFunc(vector* v, vector* antiGrad, double lambda) {
    return func0(nextValue(v, antiGrad, lambda));
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
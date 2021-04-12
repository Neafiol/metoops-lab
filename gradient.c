#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846
#define TAU 0.61803398874989484820

typedef struct vector {
    double x, y;
} vector;

double func(struct vector x);
vector findGradient(vector* x);
vector nextValue(vector* x, vector* antiGrad, double* lambda);
vector gradientDescent(vector* x, double* eps, int mode);
double brent(double a, double b, double eps, vector* x, vector* antiGrad);
double sign(double x);
double lambdaFunc(vector* x, vector* antiGrad, double* lambda);
double getParabolaMin(double x1, double y1, double x2, double y2, double x3, double y3);

int main() {
    vector x;
    double eps;
    x.x = 10;
    x.y = -10;
    eps = 0.0001;
    vector result = gradientDescent(&x, &eps, 0);
    printf("x = %0.6f, y = %0.6f\n", result.x, result.y);
    result = gradientDescent(&x, &eps, 1);
    printf("x = %0.6f, y = %0.6f\n", result.x, result.y);
    return 0;
}

double func(struct vector x) {
    //return x.x * x.x + x.y * x.y;
    //return x.x * x.x - x.y * x.y;
    return 64 * x.x * x.x + 126 * x.x * x.y + 64 * x.y * x.y - 10 * x.x + 30 * x.y + 13;
}

vector findGradient(vector* x) {
    vector grad;
    double h = 0.1;
    double dx = (func((vector) { .x = x->x + h, .y = x->y }) - func((vector) { .x = x->x - h, .y = x->y })) / (2 * h);
    double dy = (func((vector) { .x = x->x, .y = x->y + h }) - func((vector) { .x = x->x, .y = x->y - h })) / (2 * h);
    grad.x = dx;
    grad.y = dy;
    return grad;
}

vector nextValue(vector* x, vector* antiGrad, double* lambda) {
    vector v;
    v.x = x->x - *lambda * antiGrad->x;
    v.y = x->y - *lambda * antiGrad->y;
    return v;
}

vector gradientDescent(vector* x, double* eps, int mode) {
    vector current = *x;
    vector last;
    double difference;
    do {
        last = current;
        vector antiGrad = findGradient(&current);
        double lambda = 0.001;
        switch (mode) {
        case 0:
            lambda = 0.001;
            break;
        case 1:
            lambda = brent(0.000001, 10, *eps, &x, &antiGrad);
            break;
        default:
            break;
        }
        current = nextValue(&current, &antiGrad, &lambda);
        difference = fabs(func(current) - func(last));
    } while (difference > *eps);
    return current;
}

double lambdaFunc(vector* x, vector* antiGrad, double* lambda) {
    return func((vector) { .x = x->x - *lambda * antiGrad->x, .y = x->y - *lambda * antiGrad->y });
}

double sign(double x) {
    return (x > 0) - (x < 0);
}

double getParabolaMin(double x1, double y1, double x2, double y2, double x3, double y3) {
    double a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2));
    double b = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3));
    return -b / (2 * a);
}

double brent(double a, double b, double eps, vector* x, vector* antiGrad) {
    double x1, x2, x3;
    double f1, f2, f3;
    x1 = x2 = x3 = (a + b) / 2;
    f1 = f2 = f3 = lambdaFunc(x, antiGrad, &x1);
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
        if (fabs(u - x1) < eps) {
            u = x1 + sign(u - x1) * eps;
        }

        double fu = lambdaFunc(&x, &antiGrad, &u);
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
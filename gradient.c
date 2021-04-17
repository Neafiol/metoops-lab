#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846
#define TAU 0.61803398874989484820
#define MAX_ITER 2000
#define MAX_VECTOR_SIZE 1000

struct vector;
struct function2d;
struct functionNd;
struct baseFunction;

typedef struct vector vector;
typedef struct function2d function2d;
typedef struct functionNd functionNd;
typedef struct baseFunction baseFunction;

struct vector {
    int size;
    double data[MAX_VECTOR_SIZE];
};

typedef enum classType {_2d, Nd} classType;

struct baseFunction {
    classType type;
};

vector getVector2d(double x, double y);
vector getVectorNd(double* data, int size);
function2d getFunction2d(double xx, double xy, double yy, double x, double y, double c);
functionNd getFunctionNd(vector *v);

struct functionNd {
    baseFunction base;
    vector H;
};

double calcNd(functionNd* f, vector* v) {
    double res = 0;
    for (int i = 0; i < f->H.size; i++) {
        res += f->H.data[i] * v->data[i] * v->data[i];
    }
    return res;
}

double diffNd(functionNd* f, int i, double x, double dx) {
    double xh = x + dx;
    return f->H.data[i] * (xh * xh - x * x);
}

vector gradNd(functionNd* f, vector* v, double value0) {
    const double h = 0.01;
    const double invh = 1 / h;
    vector grad = getVector2d(0, 0);
    for (int i = 0; i < f->H.size; i++) {
        grad.data[i] = diffNd(f, i, v->data[i], h) * invh - f->H.data[i] * h;
    }
    return grad;
}

struct function2d {
    baseFunction base;
    const double xx, xy, yy, x, y, c;
};

double calc2d0(function2d* f, double x, double y) {
    return f->xx * x * x + f->xy * x * y + f->yy * y * y + f->x * x + f->y * y + f->c;
}

double calc2d(function2d* f, vector* v) {
    return calc2d0(f, v->data[0], v->data[1]);
}

vector grad2d0(function2d* f, double x0, double y0, double value0) {
    const double h = 0.01;
    const double invh = 1 / h;
    double dx = (calc2d0(f, x0 + h, y0) - value0) * invh - f->xx * h;
    double dy = (calc2d0(f, x0, y0 + h) - value0) * invh - f->yy * h;
    return getVector2d(dx, dy);
}

vector grad2d(function2d* f, vector* v, double value0) {
    return grad2d0(f, v->data[0], v->data[1], value0);
}

double calc(baseFunction* f, vector* v);
vector findGradient(baseFunction* f, vector* v, double value0);
vector nextValue(vector* v, vector* antiGrad, double lambda);
double lambdaFunc(baseFunction* func, vector* v, vector* antiGrad, double lambda);
double norm2(vector* v);
double norm(vector* v);
void normalize(vector* v);
vector gradientDescent1(vector current, double eps, baseFunction* func);
vector gradientDescent2(vector current, double eps, baseFunction* func);
double sign(double x);
double getParabolaMin(double x1, double y1, double x2, double y2, double x3, double y3);
double brent(double a, double b, double eps, vector* v, vector* antiGrad, baseFunction* func);

void printResult(vector* result);
double* gradientDescent(double eps, int mode, int funcIndex);

static function2d func2dArr[2] = {
    {_2d, 64, 126, 64, -10, 30, 13}, // 64x^2 + 126xy + 64y^2 - 10x + 30y + 13
    {_2d, 1, 0, 1, 0, 0, 0}, // x^2 + y^2
};

int main() {
    baseFunction* func2d = (baseFunction*) (&func2dArr[0]);
    vector init2d = getVector2d(0, 0);
    double eps = 1e-8;

    vector result1 = gradientDescent1(init2d, eps, func2d);
    printResult(&result1);
    printf("\n");

    vector result2 = gradientDescent2(init2d, eps, func2d);
    printResult(&result2);
    printf("\n");

    double* res = gradientDescent(eps, 1, 1);
    printf("x0 = %0.17f, x1 = %0.17f\n", res[0], res[1]);
    printf("\n");

    double funcData[] = {1, 2, 5}; // x^2 + 2y^2 + 5z^2
    vector vecFunc = getVectorNd(funcData, 3);
    functionNd funcNd0 = getFunctionNd(&vecFunc);
    baseFunction* funcNd = (baseFunction*) (&funcNd0);
    double initData[] = {0, 0, 0};
    vector initNd = getVectorNd(initData, 3);

    vector result3 = gradientDescent1(initNd, eps, funcNd);
    printResult(&result3);
    printf("\n");

    return 0;
}

void printResult(vector* result) {
    for (int i = 0; i < result->size; i++) {
        printf("x%d = %0.17f", i, result->data[i]);
        if (i < result->size - 1) {
            printf(", ");
        } else {
            printf("\n");
        }
    }
    fflush(stdout);
}

vector nextValue(vector* v, vector* antiGrad, double lambda) {
    vector w = getVector2d(0, 0);
    for (int i = 0; i < v->size; i++) {
        w.data[i] = v->data[i] - lambda * antiGrad->data[i];
    }
    return w;
}

double lambdaFunc(baseFunction* f, vector* v, vector* antiGrad, double lambda) {
    vector val = nextValue(v, antiGrad, lambda);
    return calc(f, &val);
}

double norm2(vector* v) {
    double res = 0;
    for (int i = 0; i < v->size; i++) {
        res += v->data[i] * v->data[i];
    }
    return res;
}

double norm(vector* v) {
    return sqrt(norm2(v));
}

void normalize(vector* v) {
    double invNorm = 1 / norm(v);
    for (int i = 0; i < v->size; i++) {
        v->data[i] *= invNorm;
    }
}

double* gradientDescent(double eps, int mode, int funcIndex) {
    static double data[100];
    baseFunction* func = (baseFunction*)(&func2dArr[funcIndex - 1]);
    vector init = getVector2d(0, 0);
    vector res;
    switch (mode) {
        case 1:
            res = gradientDescent1(init, eps, func);
            break;
        case 2:
            res = gradientDescent2(init, eps, func);
            break;
        default:
            exit(1);
    }
    data[0] = res.data[0];
    data[1] = res.data[1];
    return data;
}

vector gradientDescent1(vector current, double eps, baseFunction* func) {
    vector last = getVector2d(0, 0);
    double currentValue = calc(func, &current);
    double lastValue;
    double lambda = 1;
    for (int i = 0; i < MAX_ITER; i++) {
        vector antiGrad = findGradient(func, &current, currentValue);
        if (norm2(&antiGrad) < eps * eps) {
            break;
        }

        last = current;
        lastValue = currentValue;

        for (;;) {
            current = nextValue(&last, &antiGrad, lambda / norm(&antiGrad));
            currentValue = calc(func, &current);
            if (currentValue <= lastValue) {
                break;
            }
            lambda *= 0.5;
        }
        // printResult(&current);
    };
    return current;
}

vector gradientDescent2(vector current, double eps, baseFunction* func) {
    double lambdaMin = 0;
    double lambdaMax = 1e9;
    vector last = getVector2d(0, 0);
    double currentValue = calc(func, &current);
    double lastValue;
    for (int i = 0; i < MAX_ITER; i++) {
        vector antiGrad = findGradient(func, &current, currentValue);
        if (norm2(&antiGrad) < eps * eps) {
            break;
        }

        last = current;
        lastValue = currentValue;

        normalize(&antiGrad);
        double lambda = brent(lambdaMin, lambdaMax, eps, &current, &antiGrad, func);
        current = nextValue(&current, &antiGrad, lambda);
        currentValue = calc(func, &current);
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

double brent(double a, double b, double eps, vector* v, vector* antiGrad, baseFunction* func) {
    double x1, x2, x3;
    double f1, f2, f3;
    x1 = x2 = x3 = (a + b) / 2;
    f1 = f2 = f3 = lambdaFunc(func, v, antiGrad, x1);
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

        double fu = lambdaFunc(func, v, antiGrad, u);
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

vector findGradient(baseFunction* f, vector* v, double value0) {
    switch (f->type) {
        case _2d: return grad2d((function2d*) f, v, value0);
        case Nd: return gradNd((functionNd*) f, v, value0);
        default: exit(1);
    }
}

double calc(baseFunction* f, vector* v) {
    switch (f->type) {
        case _2d: return calc2d((function2d*) f, v);
        case Nd: return calcNd((functionNd*) f, v);
        default: exit(1);
    }
}

vector getVector2d(double x, double y) {
    vector v;
    v.size = 2;
    v.data[0] = x;
    v.data[1] = y;
    return v;
}

vector getVectorNd(double* data, int size) {
    vector v;
    v.size = size;
    for (int i = 0; i < size; i++) {
        v.data[i] = data[i];
    }
    return v;
}

function2d getFunction2d(double xx, double xy, double yy, double x, double y, double c) {
    function2d f = {_2d, xx, xy, yy, x, y, c};
    return f;
}

functionNd getFunctionNd(vector *v) {
    functionNd f = {Nd, *v};
    return f;
}

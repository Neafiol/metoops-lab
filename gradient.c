#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846
#define TAU 0.61803398874989484820
#define MAX_ITER 2000
#define MAX_VECTOR_SIZE 10000

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

vector getVectorNd(double* data, int size);
vector getZeroVector(int size);

struct functionNd {
    baseFunction base;
    int size;
    double H[MAX_VECTOR_SIZE];
    double b[MAX_VECTOR_SIZE];
    double c;
};

double calcNd(functionNd* f, vector* v) {
    double res = 0;
    for (int i = 0; i < f->size; i++) {
        res += (f->H[i] * v->data[i] + f->b[i]) * v->data[i];
    }
    return res + f->c;
}

vector gradNd(functionNd* f, vector* v) {
    vector grad;
    grad.size = v->size;
    for (int i = 0; i < f->size; i++) {
        grad.data[i] = 2 * f->H[i] * v->data[i] + f->b[i];
    }
    return grad;
}

void multiplyArrayByCoordinates(double* a, double* b, double* c, int size) {
    for (int i = 0; i < size; i++) {
        c[i] = a[i] * b[i];
    }
}

vector multiplyByCoordinates(vector* a, vector* b) {
    vector v;
    v.size = a->size;
    multiplyArrayByCoordinates(a->data, b->data, v.data, v.size);
    return v;
}

double dotProduct(vector* a, vector* b) {
    vector v = multiplyByCoordinates(a, b);
    double res = 0;
    for (int i = 0; i < v.size; i++) {
        res += v.data[i];
    }
    return res;
}

void mulByScalar(vector* v, double k) {
    for (int i = 0; i < v->size; i++) {
        v->data[i] *= k;
    }
}

void addToVector(vector* v, vector* x) {
    for (int i = 0; i < v->size; i++) {
        v->data[i] += x->data[i];
    }
}

void subFromVector(vector* v, vector* x) {
    for (int i = 0; i < v->size; i++) {
        v->data[i] -= x->data[i];
    }
}

vector applyOperatorNd(functionNd* f, vector* p) {
    vector w;
    w.size = p->size;
    multiplyArrayByCoordinates(f->H, p->data, w.data, w.size);
    mulByScalar(&w, 2);
    return w;
}

void printFunctionNd(functionNd* f) {
    printf("%ga^2", f->H[0]);
    for (int i = 1; i < f->size; i++) {
        printf(" + %g%c^2", f->H[i], i + 'a');
    }
    for (int i = 0; i < f->size; i++) {
        printf(" + %g%c", f->b[i], i + 'a');
    }
    printf(" + %g\n", f->c);
}

struct function2d {
    baseFunction base;
    double xx, xy, yy, x, y, c;
};

double calc2d0(function2d* f, double x, double y) {
    return f->xx * x * x + f->xy * x * y + f->yy * y * y + f->x * x + f->y * y + f->c;
}

double calc2d(function2d* f, vector* v) {
    return calc2d0(f, v->data[0], v->data[1]);
}

vector grad2d0(function2d* f, double x, double y) {
    double dx = 2 * f->xx * x + f->xy * y + f->x;
    double dy = 2 * f->yy * y + f->xy * x + f->y;
    return (vector) {2, dx, dy};
}

vector grad2d(function2d* f, vector* v) {
    return grad2d0(f, v->data[0], v->data[1]);
}

vector applyOperator2d(function2d* f, vector* p) {
    double dx = 2 * f->xx * p->data[0] + f->xy * p->data[1];
    double dy = 2 * f->yy * p->data[1] + f->xy * p->data[0];
    return (vector) {2, dx, dy};
}

void printFunction2d(function2d* f) {
    printf("%gx^2 + %gxy + %gy^2 + %gx + %gy + %g\n", f->xx, f->xy, f->yy, f->x, f->y, f->c);
}

double calc(baseFunction* f, vector* v);
vector findGradient(baseFunction* f, vector* v);
vector applyOperator(baseFunction* f, vector* p);
vector nextValue(vector* v, vector* grad, double lambda);
double lambdaFunc(baseFunction* func, vector* v, vector* grad, double lambda);
double norm2(vector* v);
double norm(vector* v);
void normalize(vector* v);
vector gradientDescent1(vector current, double eps, baseFunction* func);
vector gradientDescent2(vector current, double eps, baseFunction* func);
vector gradientDescent3(vector x, double eps, baseFunction* func);
double sign(double x);
double getParabolaMin(double x1, double y1, double x2, double y2, double x3, double y3);
double brent(double a, double b, double eps, vector* v, vector* grad, baseFunction* func);

void printVector(vector* result);
void printTable(vector* n, vector* k);
double* gradientDescent(double eps, int method, int funcIndex); // for ui

double trace[2 * MAX_ITER]; // for ui
int iterCount;

function2d func2dArr[] = {
    {_2d, 64, 126, 64, -10, 30, 13}, // 64x^2 + 126xy + 64y^2 - 10x + 30y + 13
    {_2d, 1, 0, 1, -2, 5, 0}, // x^2 + y^2 - 2x + 5y
    {_2d, 25, -6, 1, -20, 300, 1}, // 25x^2 - 6xy + y^2 - 20x + 300y + 1
};

functionNd funcNdArr[] = {
    {Nd, 4, {1, 2, 5, 2}, {2, 3, -5, 10}, 5}, // x^2 + 2y^2 + 5z^2 + 2t^2 + 2x + 3y - 5z + 10t + 5
    {Nd, 2, {1, 1}, {-2, 5}, 0}, // x^2 + y^2 - 2x + 5y
};

typedef vector (*gradientFunc)(vector, double, baseFunction*);

gradientFunc gradientArr[] = {
    gradientDescent1,
    gradientDescent2,
    gradientDescent3
};

// test for 2d functions
void test1(int funcId) {
    function2d f = func2dArr[funcId];
    baseFunction* func = (baseFunction*) (&f);
    vector init = {2, 0, 0};
    double eps = 1e-16;

    printFunction2d(&f);
    printf("\n");

    for (int i = 0; i < 3; i++) {
        printf("result%d:\n", i);
        vector result = gradientArr[i](init, eps, func);
        printVector(&result);
        printf("\n");
    }
}

// test for n-dimensional functions
void test2(int funcId) {
    functionNd f = funcNdArr[funcId];
    baseFunction* func = (baseFunction*) (&f);
    vector init = getZeroVector(f.size);
    double eps = 1e-16;

    printFunctionNd(&f);
    printf("\n");

    for (int i = 0; i < 3; i++) {
        printf("result%d:\n", i);
        vector result = gradientArr[i](init, eps, func);
        printVector(&result);
        printf("\n");
    }
}

// trace test for 2d functions
void test3(int methodId, int funcId) {
    double* res = gradientDescent(1e-16, methodId + 1, funcId + 1);
    printf("trace:\n");
    for (int i = 0; i < MAX_ITER; ++i) {
        double x = res[2 * i];
        double y = res[2 * i + 1];
        if (x == 0 && y == 0) {
            break;
        }
        printf("x = %0.17f, y = %0.17f\n", x, y);
    }
    printf("\n");
}

int getRand(int from, int to) {
    return rand() % (to - from + 1) + from;
}

//test for random n-dimensional functions
void test4(int n, int k, bool out) {
    double eps = 1e-16;
    int base = 16;

    functionNd f = {Nd, n};
    for (int i = 0; i < n; i++) {
        f.H[i] = getRand(base, base * k);
        f.b[i] = getRand(-base * (n + k), base * (n + k));
    }
    f.H[0] = base;
    f.H[1] = base * k;

    if (n <= 26 && out) {
        printFunctionNd(&f);
        printf("\n");
    }

    baseFunction* func = (baseFunction*) (&f);
    vector init = getZeroVector(n);

    vector result = gradientDescent3(init, eps, func);
    if (out) {
        printf("count of iterations: %d\n", iterCount);
        printf("result:\n");
        printVector(&result);
        printf("\n");
    }
}

void printIterationTable() {
    vector n = {7, 10, 50, 100, 500, 1000, 5000, 10000};
    vector k = {17, 1, 2, 5, 10, 20, 50, 100, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000};
    printTable(&n, &k);
}

int main() {
    srand(3);

    // test1(2);
    // test2(0);
    // test3(2, 2);
    // test4(10, 50, 1);
    printIterationTable();
}

void printVector(vector* result) {
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

vector nextValue(vector* v, vector* grad, double lambda) {
    vector w;
    w.size = v->size;
    for (int i = 0; i < v->size; i++) {
        w.data[i] = v->data[i] - lambda * grad->data[i];
    }
    return w;
}

double lambdaFunc(baseFunction* f, vector* v, vector* grad, double lambda) {
    vector val = nextValue(v, grad, lambda);
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

// for ui
void clearTrace() {
    for (int i = 0; i < 2 * MAX_ITER; i++) {
        trace[i] = 0;
    }
}

// for ui
double* gradientDescent(double eps, int methodId, int funcId) {
    clearTrace();
    baseFunction* func = (baseFunction*) (&func2dArr[funcId - 1]);
    gradientFunc method = gradientArr[methodId - 1];
    vector init = {2, 1, 1};
    method(init, eps, func);
    return trace;
}

vector gradientDescent1(vector current, double eps, baseFunction* func) {
    vector last;
    double currentValue = calc(func, &current);
    double lastValue;
    double lambda = 1;
    for (int i = 0; i < MAX_ITER; i++) {
        trace[2 * i] = current.data[0]; // for ui
        trace[2 * i + 1] = current.data[1]; // for ui

        vector grad = findGradient(func, &current);
        if (norm2(&grad) < eps * eps) {
            break;
        }

        last = current;
        lastValue = currentValue;

        for (int j = 0; j < 32; j++) {
            current = nextValue(&last, &grad, lambda / norm(&grad));
            currentValue = calc(func, &current);
            if (currentValue < lastValue) {
                break;
            }
            lambda *= 0.5;
        }
        // printVector(&current);
    };
    return current;
}

vector gradientDescent2(vector current, double eps, baseFunction* func) {
    double lambdaMin = 0;
    double lambdaMax = 1e9;
    vector last;
    for (int i = 0; i < MAX_ITER; i++) {
        trace[2 * i] = current.data[0]; // for ui
        trace[2 * i + 1] = current.data[1]; // for ui

        vector grad = findGradient(func, &current);
        if (norm2(&grad) < eps * eps) {
            break;
        }

        last = current;

        normalize(&grad);
        double lambda = brent(lambdaMin, lambdaMax, eps, &current, &grad, func);
        current = nextValue(&current, &grad, lambda);
        // printVector(&current);
    };
    return current;
}

vector gradientDescent3(vector x, double eps, baseFunction* func) {
    trace[0] = x.data[0]; // for ui
    trace[1] = x.data[1]; // for ui

    vector grad = findGradient(func, &x); // grad = Hf
    double gradNorm2 = norm2(&grad); // gradNorm2 = ||Hf||^2
    vector p = grad;
    mulByScalar(&p, -1); // p = -Hf

    for (int i = 1; i <= MAX_ITER && gradNorm2 > eps * eps; i++) {
        vector Ap = applyOperator(func, &p); // Ap = A*p
        double Ap_p = dotProduct(&Ap, &p); // Ap_p = (A*p, p)
        double a = gradNorm2 / Ap_p; // a = ||Hf||^2 / (A*p, p)
        vector ap = p;
        mulByScalar(&ap, a); // ap = a*p
        addToVector(&x, &ap); // x = x' + a*p
        // printVector(&x);
        mulByScalar(&Ap, a); // Ap = a*A*p
        addToVector(&grad, &Ap); // grad = Hf1 = Hf + a*A*p
        double grad1Norm2 = norm2(&grad); // grad1Norm2 = ||Hf1||^2
        double b = grad1Norm2 / gradNorm2; // b = ||Hf1||^2 / ||Hf||^2
        mulByScalar(&p, b); // p = b*p'
        subFromVector(&p, &grad); // p = -Hf1 + b*p'
        gradNorm2 = grad1Norm2;

        trace[2 * i] = x.data[0]; // for ui
        trace[2 * i + 1] = x.data[1]; // for ui
        iterCount = i;
    }
    return x;
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

double brent(double a, double b, double eps, vector* v, vector* grad, baseFunction* func) {
    double x1, x2, x3;
    double f1, f2, f3;
    x1 = x2 = x3 = (a + b) / 2;
    f1 = f2 = f3 = lambdaFunc(func, v, grad, x1);
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

        double fu = lambdaFunc(func, v, grad, u);
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

vector findGradient(baseFunction* f, vector* v) {
    switch (f->type) {
        case _2d: return grad2d((function2d*) f, v);
        case Nd: return gradNd((functionNd*) f, v);
        default: exit(1);
    }
}

vector applyOperator(baseFunction* f, vector* p) {
    switch (f->type) {
        case _2d: return applyOperator2d((function2d*) f, p);
        case Nd: return applyOperatorNd((functionNd*) f, p);
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

vector getVectorNd(double* data, int size) {
    vector v;
    v.size = size;
    for (int i = 0; i < size; i++) {
        v.data[i] = data[i];
    }
    return v;
}

vector getZeroVector(int size) {
    vector v;
    v.size = size;
    for (int i = 0; i < v.size; i++) {
        v.data[i] = 0;
    }
    return v;
}

void printTable(vector* n, vector* k) {
    const int iter = 20;

    printf("n, k ");
    for (int j = 0; j < k->size; j++) {
        printf(" & %-4d", (int) k->data[j]);
    }
    printf("\n");

    for (int i = 0; i < n->size; i++) {
        printf("%-5d", (int) n->data[i]);
        for (int j = 0; j < k->size; j++) {
            int iterCnt = 0;
            for (int t = 0; t < iter; t++) {
                test4(n->data[i], k->data[j], 0);
                iterCnt += iterCount;
            }
            iterCnt = (iterCnt + iter / 2) / iter;
            printf(" & %-4d", iterCnt);
        }
        printf("\n");
        fflush(stdout);
    }
}

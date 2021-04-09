#include <stdio.h>
#include <math.h>

typedef struct vector {
    double x, y;
} vector;

double func(struct vector x);
vector findGradient(vector* x);
vector nextValue(vector* x, vector* antiGrad, double* lambda);
vector gradientDescent(vector* x, double* eps);

int main() {
    vector x;
    double eps;
    x.x = -10;
    x.y = -10;
    eps = 0.0001;
    vector result = gradientDescent(&x, &eps);
    printf("x = %0.6f, y = %0.6f\n", result.x, result.y);
    return 0;
}

double func(struct vector x) {
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

vector gradientDescent(vector* x, double* eps) {
    vector current = *x;
    vector last;
    double difference;
    do {
        last = current;
        vector antiGrad = findGradient(&current);
        double lambda = 0.001;
        current = nextValue(&current, &antiGrad, &lambda);
        difference = abs(func(current) - func(last));
    } while (difference > *eps);
    return current;
}
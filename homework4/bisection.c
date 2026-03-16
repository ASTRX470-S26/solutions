#include <stdio.h>
#include <math.h> // for nan, fabs
#include "methods.h"

float bisection(float xmin, float xmax, int nmax, float eps, float (*f)(float)) {
    int i;
    for (i = 0; i < nmax; i++) {
        float x = (xmin + xmax) / 2;
        if (f(xmin) * f(x) > 0) {
            xmin = x;
        } else {
            xmax = x;
        }
        if (fabs(f(x)) < eps) {
            break;
        }
    }
    if (i == nmax) {
        printf("Warning: failed to converge. nmax: %d, eps: %g.\n",nmax,eps);
        return nan("");
    } else {
        return (xmin + xmax) / 2;
    }
}

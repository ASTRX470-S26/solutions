#include <stdio.h>
#include <math.h>
#include "methods.h"

float newton(float x0, float dx, int nmax, float eps, float (*f)(float)) {
  
    float df, dxt;
    int i;
    for (i = 0; i < nmax; ++i) {
        df = (f(x0 + dx/2) - f(x0 - dx/2)) / dx;
        dxt = -f(x0) / df;
        while (fabs(f(x0 + dxt)) > fabs(f(x0))) {
            dxt /= 2;
        }
        x0 = x0 + dxt;
        if (fabs(f(x0)) < eps) {
            break;
        }
    }
    if (i == nmax) {
        printf("Warning: failed to converge. nmax: %d, eps: %g.\n",nmax,eps);
        return nan("");
    } else {
        return x0;
    }

}

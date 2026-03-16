#include <stdio.h>
#include <math.h>
#include "methods.h"

float secant(float x0, float x1, int nmax, float eps, float (*f)(float)) {
    float x;
    int i;
    for (i = 0; i < nmax; ++i) {
        x = x1 - (*f)(x1) * (x1 - x0) / ((*f)(x1) - (*f)(x0));
        
        // Update the guesses for the next iteration
        x0 = x1;
        x1 = x;

	if (fabs((*f)(x1)) < eps) {
            break;
        }
    }
    if (i == nmax) {
        printf("Warning: failed to converge. nmax: %d, eps: %g.\n",nmax,eps);
        return nan("");
    } else {
        return x;
    }
}

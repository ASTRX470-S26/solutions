#include <stdio.h>
#include <stdlib.h>

#define NINT 100

void spline_quadratic(double *x, double *y, double *xint, double *yint, int size_x, int size_xint) {
  
  double z[size_x];
  z[0] = 1;
    
  // compute slopes and zs
  double m[size_x - 1];
  for (int i = 0; i < size_x - 1; i++) {
    m[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]);
  }
  
  for (int i = 1; i < size_x; i++) {
    z[i] = -z[i-1] + 2 * m[i-1];
  }
  
  int i = 0;
  for (int j = 0; j < size_xint; j++) {
    // ensure xint is between x[i] and x[i+1]
    while (xint[j] < x[i] || xint[j] > x[i+1]) {
      i += 1;
    }
    // compute interpolated value
    double delta = xint[j] - x[i];
    yint[j] = y[i] + z[i] * delta + (z[i+1] - z[i]) / (2 * (x[i+1] - x[i])) * delta * delta;
  }
}

int main() {
  
  double x[NINT];
  double y[NINT];
  
  unsigned int seed = 131523;
  srand(seed);   // Initialization
 
  for (int i=0; i<NINT; ++i) {
    x[i] = (double)i*0.1;
    // return random double between 0 and 1
    double r = (double)(rand()/RAND_MAX);
    y[i] = 3*x[i]+5 + 2.*(1.-r);
  }
  
  double xint[] = {0.5, 2.333, 4.3, 8.02};
  double yint[3];
 
  int size_xint = 4;
  int size_x = NINT;
  spline_quadratic(x, y, xint, yint, size_x, size_xint);

  printf("Interpolated values:\n");
  for (int i = 0; i < size_xint; ++i) {
        printf("%.5f ", yint[i]);
  }
  printf("\n");
  
  return 0;
}

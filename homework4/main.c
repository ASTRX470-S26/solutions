#include <stdio.h> // for printf
#include <math.h> // for sin
#include <stdlib.h>  // for atof, atoi
#include <string.h> // for strcomp
#include "methods.h"


// Example function
float example_func(float x) {
  return pow(x,3)-2*x+2;
}

#include <stdlib.h>  // for atof, atoi

int main(int argc, char *argv[]) {
  // Check for minimum number of arguments and print usage message
  if (argc < 2) {
    printf("Usage:\n");
    printf("  %s bisection <xmin> <xmax> [nmax] [eps]\n", argv[0]);
    printf("  %s secant <guess1> <guess2> [nmax] [eps]\n", argv[0]);
    printf("  %s newton <guess> <dx> [nmax] [eps]\n", argv[0]);
    return 1;
  }
  
  char *method = argv[1];
  printf("Using: %s \n", method);
  
  // Default values
  int nmax = 1000;
  float eps = 1.e-7;
    
  float root;
  if (strcmp(method, "bisection") == 0) {
    if (argc < 4) {
      printf("Error: bisection requires <xmin> <xmax>\n");
      return 1;
    }
    float xmin = atof(argv[2]);
    float xmax = atof(argv[3]);
    if (argc > 4) nmax = atoi(argv[4]);
    if (argc > 5) eps = atof(argv[5]);
    
    root = bisection(xmin, xmax, nmax, eps, example_func);
    
  } else if (strcmp(method, "secant") == 0) {
    if (argc < 4) {
      printf("Error: secant requires <guess1> <guess2>\n");
      return 1;
    }
    float guess1 = atof(argv[2]);
    float guess2 = atof(argv[3]);
    if (argc > 4) nmax = atoi(argv[4]);
    if (argc > 5) eps = atof(argv[5]);
    
    root = secant(guess1, guess2, nmax, eps, example_func);
    
  } else if (strcmp(method, "newton") == 0) {
    if (argc < 4) {
      printf("Error: newton requires <guess> <dx>\n");
      return 1;
    }
    float guess = atof(argv[2]);
    float dx = atof(argv[3]);
    if (argc > 4) nmax = atoi(argv[4]);
    if (argc > 5) eps = atof(argv[5]);
    
    root = newton(guess, dx, nmax, eps, example_func);
    
  } else {
    printf("Error: method %s not supported\n", method);
    return 1;
  }
  
  printf("Root: %.6f\n", root);
  return 0;
}

float bisection(float xmin, float xmax, int nmax, float eps, float (*f)(float));
float secant(float x0, float x1, int nmax, float eps, float (*f)(float));
float newton(float x0, float dx, int nmax, float eps, float (*f)(float));

"""Homework 1 solutions"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import leggauss


def read2dict(name):
    """Read a dictionary"""

    arr = np.loadtxt(name)
    data = {}
    data['i'] = arr[:, 0]
    data['x1v'] = arr[:, 1]
    data['rho'] = arr[:, 2]
    data['press'] = arr[:, 3]
    data['vel1'] = arr[:, 4]
    data['vel2'] = arr[:, 5]
    data['vel3'] = arr[:, 6]

    return data


def plot_var(data, var, col='k'):
    """Plot a variable"""

    plt.plot(data['x1v'], data[var], col, lw=2)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('x1v')
    plt.ylabel(var)
    plt.show()


def trapezoid(f, a, b, n):
    """Use trapezoid method to compute integral of functinon f"""
    h = (b-a) / (n-1)
    quad = (f(a)+f(b))/2.*h
    for i in range(1, n-1):
        quad += h * f(h*i+a)
    return quad


def simpson(f, a, b, n):
    """Use Simpson's rule to compute integral of functinon f"""
    if n % 2 == 0:
        print("n needs to be odd for Simpson's rule")
        return None

    h = (b-a)/(n-1)
    quad = (f(a)+f(b))/3.*h+f(b-h)*4./3.*h
    for i in range(1, n-2):
        if i % 2 == 0:
            quad += 2/3*h*f(a+h*i)
        if i % 2 == 1:
            quad += 4/3*h*f(a+h*i)

    return quad


def gaussian(f, a, b, n):
    """Use Gauss-Legendre quad. to compute integral of functinon f"""
    x, w = leggauss(n)

    # transform the interval
    y = (b-a)/2*f((b-a)/2*x+(b+a)/2)

    return np.dot(y, w)

"""Homework 2 solutions"""

# import packages here
from random import randint
from numpy.polynomial.legendre import leggauss
import numpy as np


def joke():
    """ Tell a joke """

    print("There are 10 types of people in the world.  Those that "
          "understand binary and those who don't.")


def simpson_slow(f, a, b, n):
    """Use loop for Simpson's rule to compute integral of function """
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


def simpson_fast(f, a, b, n):
    """ Simpson's rule to compute integral of functinon -- no loop"""

    if n % 2 == 0:
        print("n needs to be odd for Simpson's rule")
        return None

    dx = (b-a) / (n-1)
    x = np.arange(a, b+dx, dx)
    fa = f(x)
    fa[1::2] *= 4
    fa[0::2] *= 2
    fa[0] *= 0.5
    fa[-1] *= 0.5

    return dx/3*fa.sum()


def gaussian(f, a, b, n):
    """ Implement Gauss Legendre quadrature """

    # obtain weights
    x, w = leggauss(n)

    # transform the interval
    y = (b-a)/2*f((b-a)/2*x+(b+a)/2)

    return np.dot(y, w)


def deriv2cd(f, x, dx):
    """ Implement centered difference second derivative """

    return ((f(x+dx)-f(x)) - (f(x)-f(x-dx)))/dx**2


def deriv2cd_5(f, x, dx):
    """ Implement higher order centered difference second derivative """

    return (-f(x+2*dx)-f(x-2*dx)+16*f(x+dx)+16*f(x-dx)-30*f(x))/12/dx**2


class Dice:
    """ Class for simulating dice rolls"""

    __results = []

    def __init__(self, faces=6, number=1):
        self.faces = faces
        self.number = number

    def roll(self):
        """ simulate dice roll """

        self.__results = []
        for _ in range(self.number):
            self.__results.append(randint(1, self.faces))

    def results(self):
        """ return current results """

        return self.__results

    def sort(self):
        """ sort results """

        self.__results.sort(reverse=True)

    def clear(self):
        """ clear roll results """

        self.__results = []

    def change_number(self, number):
        """ Change number of dice """

        self.number = number


def compute_probability(ntrials):
    """ Function for computing probability of survival """

    dustin = Dice(20, 3)

    success = 0
    for _ in range(ntrials):
        dustin.roll()
        dustin.sort()
        if dustin.results()[0] == 20:
            success += 1

    return success/ntrials


class Territory:
    """ Class for simulating risk territory"""

    def __init__(self, armies):
        self.armies = armies

    def add_armies(self, gains):
        """ Add armies to territory """

        self.armies += gains

    def remove_armies(self, losses):
        """ Remove armies from territory """

        self.armies -= losses
        self.armies = max(self.armies, 0)


def risk_battle(attacker, defender):
    """ Function for simulating risk battle """

    if attacker.armies <= 1:
        print("attacker must have at least one army")
        return
    if defender.armies < 1:
        print("defdnder must have at least one army")
        return

    # assume both attacker and defender alway roll
    # the maximum number of dice available
    natt = min(attacker.armies-1, 3)
    ndef = min(defender.armies, 2)

    # create dice for each territory
    attdice = Dice(6, natt)
    defdice = Dice(6, ndef)

    # roll the dice
    attdice.roll()
    defdice.roll()

    # sort the dice
    attdice.sort()
    defdice.sort()

    ntest = min(natt, ndef)

    for i in range(ntest):
        if attdice.results()[i] <= defdice.results()[i]:
            attacker.remove_armies(1)
        else:
            defender.remove_armies(1)

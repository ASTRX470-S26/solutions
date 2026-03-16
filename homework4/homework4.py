import numpy as np
from scipy.optimize import minimize


def least_squares_fit(x, y, model, guess):
    """
    Perform least squares fitting using scipy.optimize.minimize.

    Parameters:
    x (numpy array): Independent variable.
    y (numpy array): Dependent variable.
    model : Model function to fit the data.
    guess (numpy array): Initial guess for the parameters of the model.

    returns the result
    """

    # Define the fit statistic as the function to be minimized
    def statistic(params):
        return np.sum((y - model(*params, x))**2)

    #  Minimize the least squares statistic
    res = minimize(statistic, guess)

    return res.x

def chi_squared_fit(x, y, yerr, model, guess):
    """
    Perform least squares fitting using scipy.optimize.minimize.

    Parameters:
    x (numpy array): Independent variable.
    y (numpy array): Dependent variable.
    model : Model function to fit the data.
    guess (numpy array): Initial guess for the parameters of the model.

    returns the result
    """

    # Define the fit statistic as the function to be minimized
    def statistic(params):
        resid = (y - model(*params, x)) / yerr
        return np.sum(resid**2)

    #  Minimize the least squares statistic
    res = minimize(statistic, guess)

    return res.x, res.fun

import numpy as np
import astropy.constants as cons
from scipy.optimize import curve_fit

def bilinear_interpolation(x, y, z, xint, yint):
    """
    Perform bilinear interpolation.

    Args:
    x (numpy.ndarray): x-coordinates of the known points.
    y (numpy.ndarray): y-coordinates of the known points.
    z (numpy.ndarray): Values at the known points. z[i][j] corresponds to the value at (x[i], y[j]).
    xint (float): x-coordinate of the point to interpolate.
    yint (float): y-coordinate of the point to interpolate.

    Returns:
    float: Interpolated value at the point (xint, yint).
    """

    # Quick way to get x values of neighboring ponts
    xm = np.max(x[x <= xint])
    xp = np.min(x[x >= xint])
    ym = np.max(y[y <= yint])
    yp = np.min(y[y >= yint])

    if xm == xp or ym == yp:
        return None

    im = np.where(x == xm)[0][0]
    ip = np.where(x == xp)[0][0]
    jm = np.where(y == ym)[0][0]
    jp = np.where(y == yp)[0][0]

    # Bilinear interpolation formula
    z00 = z[im, jm]
    z10 = z[ip, jm]
    z01 = z[im, jp]
    z11 = z[ip, jp]
    
    x_delta = xp - xm
    y_delta = yp - ym

    s = (xint - xm) / (xp - xm)
    t = (yint - ym) / (yp - ym)

    interpolated_value = (1 - s) * (1 - t) * z00 + s * (1 - t) * z10 + (1 - s) * t * z01 + s * t * z11

    return interpolated_value

def fit_spectrum(nu, inu, model, guess):
    """Fit spectrum with model """

    params, pcov = curve_fit(model, nu, inu, p0=guess)

    return params

def blackbody(nu, temp, area):
    """ Compute the Planck function """
    
    h = cons.h.cgs.value
    c = cons.c.cgs.value
    kb = cons.k_B.cgs.value
    
    bb = area*2*h/c**2*nu**3/(np.exp(h*nu/kb/temp)-1)

    return bb

def compute_jacobian(t, func, x, h=1e-4, tiny=1e-8):
    J = np.zeros([t.size,x.size])
    for i in range(x.size) :
        delta = h*abs(x[i]) + tiny
        x1, x2 = x.copy(), x.copy()
        x1[i] -= delta
        x2[i] += delta
        J[:,i] =  (func(t,x2) - func(t,x1))/(2*delta)
    return J
    
def nonlinear_least_squares(t, y, func, guess, eps=1e-6, nmax=100, tiny=1e-10, alpha=0.5): 
    """
    t - data points along x-axis
    y - data points in y
    func - function to fit
    guess - parameters to start
    nmax - max number of iterations
    eps - max error
    tiny - minimum number
    alpha - parameter controlling how large of a change to make

    Returns minimum, error, iterations or None, None, None if the interval is violated
    """

    error = None
    x = guess
    for i in range(nmax):
        r = func(t, x) - y 
        J = compute_jacobian(t, func, x)
        JT = np.transpose(J)
        dx = -np.linalg.inv(JT@J)@JT@r
        x += alpha*dx
        error = (np.abs(dx)/(np.abs(x)+tiny)).sum()
        if(error < eps): 
            return x, error, i+1

    print( "nmax exceeded")
    return x, error, nmax
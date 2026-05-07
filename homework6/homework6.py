import numpy as np
from scipy.integrate import solve_ivp

def rk2(derivatives, t, dt, y):
    k1 = derivatives(t, y)*dt
    k2 = derivatives(t+dt/2, y+k1/2)*dt
    return y + k2


# You can generallize the example above for RK2
def rk4(derivatives, t, dt, y, args=()):

    k1 = derivatives(t, y, *args)
    k2 = derivatives(t+0.5*dt, y+k1*0.5*dt, *args)
    k3 = derivatives(t+0.5*dt, y+k2*0.5*dt, *args)
    k4 = derivatives(t, y+k3*dt, *args)

    return y+dt*(k1+2*(k2+k3)+k4)/6

def matrix_multiply(matrix1, matrix2):
    """
    Multiply two matrices together
    """
    l1 = matrix1.shape[1]
    l2 = matrix2.shape[0]

    # Check if matrices can be multiplied
    if l1 != l2:
        print("Matrices cannot be multiplied. Inner dimensions must match.")
        return None

    # Initialize result matrix with zeros
    result = [[0 for _ in range(matrix2.shape[1])] for _ in range(matrix1.shape[0])]

    # Perform matrix multiplication
    for i in range(matrix1.shape[0]):
        for j in range(matrix2.shape[1]):
            for k in range(matrix2.shape[0]):
                result[i][j] += matrix1[i][k] * matrix2[k][j]

    return np.array(result)

def matrix_multiply_new(matrix1, matrix2):
    """
    Multiply two matrices together
    """
    return matrix1 @ matrix2

def invert_matrix(mat):
    """
    Invert a matrix
    """

    return np.linalg.inv(mat)

def double_pendulum(t, y, m, l, g):
    """
    Equations of motion for a double compound pendulum
    """

    theta1 = y[0]
    theta2 = y[1]
    p1 = y[2]
    p2 = y[3]

    # these are the equations of motion one gets from a Lagrangian
    # approach.  Note these assume a compound pendulum (e.g., a rod),
    # and the moment of inertia term is that for a rod.

    # moment of inertia-like term
    I = m*l**2
    
    denom = 16 - 9*np.cos(theta1 - theta2)**2
    t1dot = 6/I * (2*p1 - 3*np.cos(theta1 - theta2)*p2)/denom
    t2dot = 6/I * (8*p2 - 3*np.cos(theta1 - theta2)*p1)/denom
    
    p1dot = -0.5*I*(t1dot*t2dot*np.sin(theta1 - theta2) + 3*g*np.sin(theta1)/l)
    p2dot = -0.5*I*(-t1dot*t2dot*np.sin(theta1 - theta2) + g*np.sin(theta2)/l)

    return np.asarray([t1dot, t2dot, p1dot, p2dot])


def run_integrator(method, th1, th2, max_dt=np.inf, t0=0, t1=100):
    """
    Function for running the integration with scipy solve_ivp
    """
    
    g = 9.81
    l = 4
    m = 1.0

    y = np.empty(4)
    y[0] = np.radians(th1) # theta1 = 80.
    y[1] = np.radians(th2) # theta2 = 10.
    y[2] = 0. # p_theta1 = 0.
    y[3] = 0. # p_theta2 = 0.
    
    solution = solve_ivp(double_pendulum, (t0,t1), y, method, args=(m,l,g),
                         max_step = max_dt)
        
    return solution.t, solution.y

def sound_speed(r, v, bern, gamma):

    return np.sqrt((gamma-1)*(bern-v**2/2+2/r))

def rho(cs, gamma):

    return cs**(2/(gamma-1))

def mdot(r, v, rho):

    return -4.*np.pi*rho*v*r**2

def deriv_fluid(r, v, bern, gamma):
    cs = sound_speed(r, v, bern, gamma)
    return -2*v*cs**2/(cs**2-v**2)/r*(1-1/cs**2/r)

def bc(yl, yr):

    return np.array([yl[0]-2,yr[0]-3])
     
def deriv_poisson(x, y):

    dpdx = y[1]
    dodx = -np.pi**2/4*(y[0]-1)**2

    return np.array([dpdx,dodx])

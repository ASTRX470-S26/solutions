import numpy as np
import random as ran
import math


def initial_conditions(x, L, prof): 
    """
    Initial conditions options: tophat, gaussian, sine
    """
    
    nx = x.size
    f = np.ones(nx)*0.5
  
    # top hat
    if(prof == 'tophat') :
        f[int(0.25*nx):int(0.75*nx)] = 1.5
    elif(prof == 'gaussian') : 
        # gaussian profile
        sigma = 0.125*L
        f += np.exp(-(x-x.mean())**2/sigma**2)
    elif(prof == 'sine'):
        f += 0.25*np.sin(2*np.pi*x)

    return f

def calculate_flux(u, v): 
    return u*v

def periodic_bc(u, ng):
    """ Periodic boundary conditions """

    u[0:ng] = u[-2*ng:-ng]
    u[-ng:] = u[ng:2*ng]

    return u

def upwind(FL, FR, v):
    """ Upwind flux """

    F = FL
    if(v < 0):
        F = FR 
    return F


def minmod(a, b):
    """ Min/mod slope limiter """
 
    slope = np.zeros(a.size)
    slope[a*b <= 0] = 0 

    boolarr = np.logical_and(a*b > 0, np.abs(a) < np.abs(b))
    slope[boolarr] = a[boolarr]

    boolarr = np.logical_and(a*b > 0, np.abs(a) > np.abs(b))
    slope[boolarr] = b[boolarr]

    return slope

def reconstruct(u, first_order, limiter):
    """
    Reconstruct the left and right states of u:
    Two options:
    Use first order or second order integration
    If second_order, use a limiter or not
    """

    if first_order:
        uL = u.copy()[:-1]
        uR = u.copy()[1:]
    else:
        uL = u.copy()[1:-2]
        uR = u.copy()[2:-1]
        if limiter:
            dup = u[2:] - u[1:-1]
            dum = u[1:-1] - u[:-2]
            du = minmod(dup, dum)
        else:
            du = (u[2:] - u[:-2])/2
        uL += 0.5*du[:-1]
        uR -= 0.5*du[1:]
        
    return uL, uR

def divergence(u, v, dx, first_order, limiter):
    """ Evolve on step of the advection algorithm """

    # perform reconstructions
    uL, uR = reconstruct(u, first_order, limiter)

    # calculate left and right state fluxes
    FL = calculate_flux(uL, v)
    FR = calculate_flux(uR, v)
    
    # use upwind to set the flux
    F = upwind(FL, FR, v)
    
    # return divergence
    return (F[1:]-F[:-1])/dx

def run_advection_2ndorder(x, L, cfl=0.4, tend=1, first_order=False, limiter=True,
                  prof='tophat'):
    """ Run the advection problem from start to finish """

    v = 1 
    t = 0
    dx = x[1]-x[0]
    nx = x.size
    if first_order:
        ng = 1
    else:
        ng = 2

    # initialize f with ghost zones on each end
    u = np.zeros(nx+2*ng)
    u[ng:-ng] = initial_conditions(x, L, prof)


    while(t < tend) : 
        dt = min(tend-t, np.abs(cfl*dx/v))

        # apply boundary conditions
        u = periodic_bc(u, ng)

        # compute divergence and update intermediate state
        u1 = u.copy()
        div1 = divergence(u, v, dx, first_order, limiter)
        u1[ng:-ng] -= div1*dt*0.5

        # apply boundary conditions
        u1 = periodic_bc(u1, ng)
        
        div = divergence(u1, v, dx, first_order, limiter)
        u[ng:-ng] -= div*dt

        t += dt
        

    # return u without ghost zones
    return u[ng:-ng]

def iterate_sor(omega, itermax=10000, eps=1.e-3):
    """
    Solve Poisson's equation using SOR
    """

    # initialize grid and rhs on this process
    ny = nx = 101
    x0 = 0
    x1 = 1
    y0 = 0
    y1 = 1
    
    # set up calculation grid
    x = np.linspace(x0, x1, nx)
    y = np.linspace(y0, y1, ny) 
    dx = x[1]-x[0]
    
    phi = np.zeros((nx, ny))
    rho = np.zeros_like(phi)
    rho[50, 50] = 1

    err = 10*eps # ensure while executes at least once
    iterations = 0
    while err > eps:
        iterations += 1
        phinew = np.copy(phi)
        for i in range(1, nx-1):                                                
            for j in range(1, ny-1):

                dphi = 0.25*(phi[i-1, j]+phi[i+1, j]+phi[i, j-1]+phi[i, j+1]-dx*dx*rho[i, j])-phi[i, j]
                phi[i, j] += dphi*omega

        TINY = 1e-10 # avoids potential NaN
        err = np.abs((phi-phinew)[1:-1,1:-1]/(0.5*np.abs(phi+phinew)+TINY)[1:-1,1:-1]).sum()
        
        if iterations == itermax:
            break

    print("Total iteration: {0} {1:.3e}".format(iterations, err))

    return x, y, phi

def analytic_solution(x, y, mu, c):

    X, Y = np.meshgrid(x, y)
    rad = np.sqrt((X-0.5)**2 + (Y-0.5)**2)

    G = 0.25/np.pi
    return np.log(rad)*mu*2*G + c

def run_transfer(nphot, eps, tau0):
    """ Perform Monte Carlo radiation transfer """

    nesc = 0.
    nscat = 0
    for i in range(nphot):

        z = -tau0
        # compute initial direction
        mu = np.sqrt(ran.uniform(0, 1))
        
        evolving = True
        while (evolving):
            
            z += mu*(-np.log(ran.uniform(0, 1)))
            # check if photon goes through boundary
            if (z >= 0):
                evolving = False
                break
            elif (z <= -tau0):
                mu = -mu
                z = -2*tau0-z
            if (ran.uniform(0, 1) < eps):
                # scatter photon isotropically
                mu = ran.uniform(-1, 1)
                nscat += 1
            else:
                # absorb photon
                evolving = False
                break

        if z >= 0:
            nesc += 1

    return nscat, nesc

# import packages here
import numpy as np
import astropy.constants as cons
import scipy.sparse as sparse
from scipy import integrate

def problem1():
    """ implements problem 1"""
    a1 = np.arange(16).reshape((4, 4)).transpose()
    print(a1)
    a2 = np.arange(20).reshape((4, 5))[:,:-1]
    print(a2)
    print(a1 * a2)
    print(a1 @ a2)


def my_gauss_seidel(amat, bvec, xinit=None, eps=1.e-10, itermax=10000):
    """ My implementation of Gauss-Seidel """

    # check if input matrices are properly formatted
    if amat.shape[0] != amat.shape[1]:
        print("Error: amat is not a square matrix")
        return
    if amat.shape[0] != bvec.shape[0]:
        print("Error: amat and bvec do not have same dimensions")
        return
    if (xinit is not None) and (xinit.shape[0] != bvec.shape[0]):
        print("Error: xinit and bvec do not have same dimensions")
        
    n = amat.shape[0]

    if xinit is None:
        x = bvec / amat.diagonal()
    else:
        x = xinit

    delta = np.empty(n)
    for iter in range(itermax):
        for i in range(n):
            delta[i] = (bvec[i] - np.dot(amat[i,:], x)) / amat[i,i]
            x[i] += delta[i]

        # compute relative change in x
        change = abs(delta)
        inds = np.where(abs(x) > 1.e-14)
        change[inds] /= abs(x[inds])

        # break if error below eps
        if np.max(change) < eps:
            break
            
    if (iter == itermax-1):
        print("Warning: maximum iterations {:d} exceeded with ".format(itermax))
        print("maximum change {:e}.".format(np.max(change)))

    return x


def absorption_extinction(temp, rho, nu):
    """ Computes free-free extinction coefficient """

    h = cons.h.cgs.value
    kb = cons.k_B.cgs.value
    mp =  cons.m_p.cgs.value
    me = cons.m_e.cgs.value
    ec = cons.e.esu.value
    c = cons.c.cgs.value
    
    norm = 4*ec**6/(3*mp**2*me*h*c)*(2*np.pi/(3*kb*me))**0.5
    return norm*rho**2/temp**0.5/nu**3*(1-np.exp(-h*nu/(kb*temp)))


def scattering_extinction(rho):
    """ Computes Thomson scattering extinction coefficient """

    sigt = cons.sigma_T.cgs.value
    mp =  cons.m_p.cgs.value
    
    return sigt/mp*rho


def planck(temp, nu):
    """ Compute the Planck function """
    
    h = cons.h.cgs.value
    c = cons.c.cgs.value
    kb = cons.k_B.cgs.value
    
    return 2*h/c**2*nu**3/(np.exp(h*nu/kb/temp)-1)


def transfer_matrices(tau, bnu, enu):
    """ Construct the transfer matrix and vector """

    n = tau.shape[0]
    
    # construct b vector
    bvec = np.copy(bnu)
    bvec[0] = 0.

    # construct a matrix
    amat = np.zeros((n,n))
    
    # top boundary i=0 
    amat[0,0] = (tau[1]-tau[0]+1/3**0.5)
    amat[0,1] = -1/3**0.5
    
    # lower boundary i=n-1
    amat[-1,-1] = 1.
    
    # interior points
    for i in range(1,n-1):
        dtaup = tau[i+1]-tau[i]
        dtaum = tau[i]-tau[i-1]
        dtaut = dtaup+dtaum
        fac = 2/3/enu[i]
        
        amat[i,i-1] = -fac/dtaut/dtaum
        amat[i,i] = fac/dtaup/dtaum + 1
        amat[i,i+1] = -fac/dtaut/dtaup

    return amat, bvec


def formal_solution(mu, tau, snu, ibase):
    """Integrate to perfor the formal solution"""
    
    integ = integrate.simpson(snu*np.exp((tau[0]-tau)/mu)/mu,tau)
 
    return ibase*np.exp((tau[0]-tau[-1])/mu) + integ


def solve_tridiagonal(amat, bvec):
    """Solve matrix with sparse solver"""
    
    n = amat.shape[0]
    data = np.zeros((3,n))
    data[0,:-1] = np.diagonal(amat, offset=-1)
    data[1,:] = np.diagonal(amat, offset=0)
    data[2,1:] = np.diagonal(amat, offset=1)
    
    diags = [-1, 0, 1]
    asp = sparse.spdiags(data,diags,n,n,format='csc')
    return sparse.linalg.spsolve(asp,bvec)


def compute_spectrum(z, rho, temp, nu, mu, taumin = 0.001, method='npsolve'):
    """ Compute the emission at the surface of the atmosphere """

    # initialize some arrays
    inu = np.empty_like(nu)
    dz = np.gradient(z)
    n = z.shape[0]
    tau = np.empty(n)
    
    for j,freq in enumerate(nu):

        # evaluate and store functions of opacities, etc.
        abnu = absorption_extinction(temp, rho, freq)
        scnu = scattering_extinction(rho)
        bnu = planck(temp, freq)
        enu = abnu / (scnu + abnu)

        # compute optical depth
        tau[0] = taumin
        for i in range(1,n):
            tau[i] = tau[i-1] + (abnu[i-1]+scnu[i-1]+abnu[i]+scnu[i])*dz[i]
        
        amat, bvec = transfer_matrices(tau, bnu, enu)
        if method == 'npsolve':
            jnu = np.linalg.solve(amat, bvec)
        elif method == 'tridiag':
            jnu = solve_tridiagonal(amat, bvec)
        elif method == 'lufact':
            p, l, u = lu_factorization(amat)
            jnu = lu_solve(p, l, u, bvec)
        else:
            print("Error: method: "+method+" not supported.")
            return None

        snu = enu*bnu+(1-enu)*jnu
        inu[j] = formal_solution(mu, tau, snu, bnu[-1])

    return inu

def lu_factorization(a):
    """Perform LU factorization with partial pivoting"""
    
    n = a.shape[0]

    # Initialize p, l, and  matrices
    p = np.eye(n) # identity matrix
    l = np.eye(n)
    u = np.copy(a)

    for j in range(n):
        # Update lower triangular matrix
        for i in range(j):
            u[i, j] = u[i, j] - u[i, :i] @ u[:i, j]
            #print("u:", j, i, u[i, j], a[i, j])
        # Update upper triangular matrix
        for i in range(j,n):
            u[i, j] = (u[i, j] - u[i, :j] @ u[:j, j])
            #print("l:",j, i, u[i, j], a[i, j])
        # Find pivot row 
        k = np.argmax(np.abs(u[j:, j])) + j

        # permute rows in u and p matrices
        u[[j, k], :] = u[[k, j], :]      
        p[[j, k], :] = p[[k, j], :]
        
        for i in range(j+1,n):
            u[i, j] /= u[j, j]

    # output lower and upper diagonals in seperate matrices
    for j in range(n):
        for i in range(j+1,n):
            l[i,j] = u[i,j]
            u[i,j] = 0.
            
    return p, l, u

def lu_solve(p, l, u, b):
    """Solve matrix equation used lu decomposition"""

    bp = p @ b
    #l y = b

    n = b.shape[0]
    y = np.empty(n)

    for i in range(n):
        y[i] = bp[i]
        for j in range(i):
            y[i] -= l[i,j]*y[j]
        #y[i] /= l[i,i]

    x = np.empty(n)
    for i in range(n-1, -1, -1):
        x[i] = y[i]
        for j in range(i+1, n):
            x[i] -= u[i,j]*x[j]
        x[i] /= u[i,i]

    return x

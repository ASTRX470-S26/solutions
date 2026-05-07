
# import packages here
import numpy as np
import astropy.constants as cons
import matplotlib.pyplot as plt

def planck(temp, nu):
    """ Compute the Planck functino """
    
    h = cons.h.cgs.value
    c = cons.c.cgs.value
    kb = cons.k_B.cgs.value
    
    return 2*h/c**2*nu**3/(np.exp(h*nu/kb/temp)-1)

def main():

    nu = np.logspace(14., 16.5, 64)

    area = 1.e10
    
    temp = 1.e5
    inu = planck(temp, nu)

    # generate random noise with a normal distribution
    amp = 0.2
    rng = np.random.default_rng()
    delta = amp*rng.normal(size=(nu.shape[0],))
    inu *= (1.+delta)

    plt.plot(nu, inu)

    out = np.empty((2,64))
    out[0,:] = nu
    out[1,:] = inu*area
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig('bb.png')
    np.savetxt("spectrum.txt",out)



if __name__=="__main__": 
    main()

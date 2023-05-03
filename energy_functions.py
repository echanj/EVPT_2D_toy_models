import numpy as np


def twoD_Gaussian(x,y, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g= offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    return g 



def coupled_oscillator_energy_force(x,y,d0=5.0,a=1.0,k=1.0,lamb=2.278):
    #calculate the energy on force on the right hand side of the equal signs
    energy = d0*(x**2.0-a**2.0)**2.0 + 0.5*k*y**2.0 + lamb*x*y 
    force_x = -(2.0*d0*(x**2.0-a**2.0)*2.0*x + lamb*y) 
    force_y = -(k*y + lamb*x)
 
    return energy, force_x, force_y


def doubleWellPot(x,a=0.5,b=0.25):

# The functional for of a symetric double well is -a*X^2 + b*X^4
# u = -0.5*x**2 + 0.25*x**4
 u = -a*x**2 + b*x**4
 return u 

def bin_centers(bin_edges):
    return (bin_edges[1:]+bin_edges[:-1])/2.

import numpy as np
from scipy.optimize import minimize
from scipy.spatial.distance import cdist

def rosen_function(x):
    """
    The energy function that we want to minimize.
    In this example, we'll use the Rosenbrock function.
    """
    return np.sum(100 * (x[1:] - x[:-1] ** 2) ** 2 + (1 - x[:-1]) ** 2,axis=0) # important that we need to specifiy the sum is along axis=0 if we want to broadcast 

def local_minimization(x0,energy_function):
    """
    Perform a local minimization starting from x0.
    In this example, we'll use the L-BFGS-B algorithm.
    """
#    res = minimize(energy_function, x0, method='L-BFGS-B')
    res = minimize(energy_function, x0,  method='CG', options={'maxiter':10, 'gtol':0.1})
    return res.x

def local_minimization_EV(x0,s,kspr,energy_function):
    """
    Perform a local minimization starting from x0.
    In this example, we'll use the L-BFGS-B algorithm.
    """
#    res = minimize(energy_function, x0, method='L-BFGS-B')
    res = minimize(energy_function, x0, args=(s,kspr),  method='CG', options={'maxiter':10, 'gtol':0.1})
    return res.x

def perturb_coordinates(x, stepsize):
    """
    Perturb the coordinates of x by a random displacement.
    The magnitude of the displacement is given by stepsize.
    """
    return x + stepsize * np.random.randn(*x.shape)

def acceptance_probability(delta_energy, temperature):
    """
    Compute the acceptance probability for a move with
    energy difference delta_energy and temperature T.
    """
    if delta_energy < 0:
        return 1.0
    else:
        return np.exp(-delta_energy / temperature)

def basin_hopping(x0, energy_function, niter=100, stepsize=0.5, temperature=1.0, **kwargs):
    """
    Perform a basin hopping optimization starting from x0.
    """
    x = x0.copy()
    minima = []
    for i in range(niter):
        # Perturb the coordinates
        x_new = perturb_coordinates(x, stepsize)

        # Perform a local minimization
        x_new = local_minimization(x_new,energy_function)

        # Compute the energy difference between the new and old configurations
        delta_energy = energy_function(x_new,  **kwargs) - energy_function(x, **kwargs)

        # Accept or reject the move based on the acceptance probability
        if acceptance_probability(delta_energy, temperature) > np.random.rand():
            x = x_new

        # Save the minimum energy configuration
        minima.append(x)

    # Return the minimum energy configuration found
    minima = np.array(minima)
    distances = cdist(minima, minima)
    index = np.argmin(distances.sum(axis=1))
    return minima[index],minima

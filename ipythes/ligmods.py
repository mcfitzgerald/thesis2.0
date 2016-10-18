#Useful functions for ligand binding simulation and fitting

#Load Dependencies
import numpy as np


#Simulating data

#Dilution series

def dilser(low=0.001, limit=100.0, dilfactor = 2.0):
    """Returns a numpy array containing a dilution series that ranges from
    "low" to "limit" by "dilfactor".
    """
    a = [low]
    while a[-1] <= limit:
        a.append(a[len(a)-1]*dilfactor)
    
    return np.asarray(a)

#Noise -- add noise to simulated data

def noiser(data, percent=0.05):
    """Takes a numpy array and randomly moves data point +/- a given "percent"
    of that data point's value. Moves are drawn from a normal distribution with
    mean = 1 and standard deviation of "percent". Returns a numpy array of the
    noisy data.
    """

    return np.random.normal(1,percent,len(data))*data

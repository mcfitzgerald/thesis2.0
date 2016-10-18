import numpy as np
import lmfit
import matplotlib.pyplot as plt

def klotz1(flig,k1):
    return (flig*k1)/(1+(flig*k1))

def objfun(params, x, data, eps=None):
    k1 = params['k1']
    model = klotz1(x,k1)
    
    if eps is None:
        return (data - model)
    else:
        weights = 1/(np.square(eps))
        return (data - model)*weights
        

    

#model functions for ligbind

import numpy as np

def klotz1_obj(params, x, data, eps=None):
    k1 = params['k1']
    model = (k1 * x)/(1 + (k1 * x))
    weights = 1/np.square(eps)
    
    return weights*(model - data)
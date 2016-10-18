def objfun(params,x,data,eps=None):
    nexp, npts = data.shape #number of experiments and number of data points per exp
    resid = np.zeros(data.shape)
    k1 = params['k1']
    model = klotz1(x,k1)
    
    if eps is None:
        for i in range(nexp):
            resid[i, :] = (data[i, :] - model)
            return resid.flatten()
    else:
        for i in range(nexp):
            weights = 1/(np.square(eps[i, :]))
            resid[i, :] = (data[i, :] - model)*weights
            return resid.flatten()
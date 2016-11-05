def dilser(low=0.001, limit=100.0, dilfactor = 2.0):
    """Returns a list containing a dilution series that ranges from
    "low" to "limit" by "dilfactor".
    """
    a = [low]
    while a[-1] <= limit:
        a.append(a[len(a)-1]*dilfactor)
    
    return np.array(a)

def wyman_sim_parmvec(parm,lig,rtot):
    """ Generates bound fraction for Wyman model given parameters and returns as 
    nested array simulated dataset  
    parm is a 1-D numpy array that must describe [k11,k21,k22,l20]
    lig is a 1-D array or 1-D array of 1-D arrays of ligand concentrations
    rtot is a 1-D numpy array of total receptor concentrations for each data set
    size of rtot and lig must match 
    """
    k11 = parm[0]
    k21 = parm[1]
    k22 = parm[2]
    l20 = parm[3]
    
    holder = []
    
    for i in range(lig.size):
        rfree = (((-1 - k11*lig[i]) + \
        (np.sqrt((1 + k11*lig[i])**2 + 8*l20*rtot[i]*(1 + k21*lig[i] + \
        k21*k22*(lig[i]**2)))))/(4*l20*(1 + k21*lig[i] + k21*k22*(lig[i]**2))))
        
        bfrac = (k11*lig[i] + l20*k21*rfree*lig[i] + \
        2*l20*k21*k22*rfree*(lig[i]**2))/(1 + 2*l20*rfree + k11*lig[i] + \
        2*l20*k21*rfree*lig[i] + 2*l20*k21*k22*rfree*(lig[i]**2))
        
        holder.append(bfrac)
        
    return np.array(holder)
    
def wyman_obj_lmfit(parm,lig,data,rtot,eps=None):
    """maximizes vector ops -- needs lmfit parameter class"""
    #define an objective function for use with lmfit (must mesh with lmfit.Parameter class)
    k11 = parm['k11']
    k21 = parm['k21']
    k22 = parm['k22']
    l20 = parm['l20']
    
    ligc = np.concatenate(lig)
    datac = np.concatenate(data)
    
    rfree = (((-1 - k11*ligc)) + \
    (np.sqrt(np.concatenate((np.square((1 + k11*lig)) + \
    8.*l20*rtot*(1 + k21*lig + k21*k22*(np.square(lig)))))))) \
    / (4*l20*(1 + k21*ligc + k21*k22*(np.square(ligc))))    
       
    bfrac = (k11*ligc + l20*k21*rfree*ligc + \
    2*l20*k21*k22*rfree*(np.square(ligc))) \
    / (1 + 2*l20*rfree + k11*ligc + \
    2*l20*k21*rfree*ligc + 2*l20*k21*k22*rfree*(np.square(ligc)))
    
    residual = (bfrac - datac)
    
    if eps is None:
        return residual
    else:
        weights = 1/np.square(np.concatenate(eps))
        return (residual*weights)

def wyman_obj_vecmin(parm,lig,data,rtot,eps=None):        
    #define an objective function for use with scipy optimize (obviates parameter class)
    k11 = parm[0]
    k21 = parm[1]
    k22 = parm[2]
    l20 = parm[3]
    
    ligc = np.concatenate(lig)
    datac = np.concatenate(data)
    
    rfree = (((-1 - k11*ligc)) + \
    (np.sqrt(np.concatenate((np.square((1 + k11*lig)) + \
    8.*l20*rtot*(1 + k21*lig + k21*k22*(np.square(lig)))))))) \
    / (4*l20*(1 + k21*ligc + k21*k22*(np.square(ligc))))    
       
    bfrac = (k11*ligc + l20*k21*rfree*ligc + \
    2*l20*k21*k22*rfree*(np.square(ligc))) \
    / (1 + 2*l20*rfree + k11*ligc + \
    2*l20*k21*rfree*ligc + 2*l20*k21*k22*rfree*(np.square(ligc)))
    
    residual = (bfrac - datac)

    if eps is None:
        return residual
    else:
        weights = 1/np.square(np.concatenate(eps))
        return (residual*weights)
        
def lnprob_lmfit(parm,lig,data,rtot,eps):
    return -0.5 * np.sum((wyman_obj_lmfit(parm,lig,data,rtot,eps))**2 + np.log(np.concatenate(2 * np.pi * eps**2)))
    
def lnprob_vecmin(parm,lig,data,rtot,eps):
    return -0.5 * np.sum((wyman_obj_vecmin(parm,lig,data,rtot,eps))**2 + np.log(np.concatenate(2 * np.pi * eps**2)))
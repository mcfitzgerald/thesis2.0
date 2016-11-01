import numpy as np

def dilser(low=0.001, limit=100.0, dilfactor = 2.0):
    """Returns a list containing a dilution series that ranges from
    "low" to "limit" by "dilfactor".
    """
    a = [low]
    while a[-1] <= limit:
        a.append(a[len(a)-1]*dilfactor)
    
    return np.array(a)

lig = dilser()

a = np.delete(lig,[15,16])

b = np.delete(lig,[14,17])

c = np.delete(lig,[13,15])

ligset = np.array([lig,a,lig,b,c])

noiseset = np.array([np.random.normal(1,0.05,ligset[i].size) for i in range(ligset.shape[0])])

rtots = np.array([0.001,0.005, 0.01, 0.02, 0.05])

sim_parms = dict({'k11':3.7, 'k21':1.8, 'k22':0.12, 'l20':293.0})

@profile
def wyman_sim(lig,rtot,parm):
    k11 = parm['k11']
    k21 = parm['k21']
    k22 = parm['k22']
    l20 = parm['l20']
    
    
    ligf = ligset
    ligc = np.concatenate(ligset)
    rtotf = rtots
    
    rfreetop = (np.concatenate((-1 - k11*ligf)) + \
    (np.sqrt(np.concatenate((np.square((1 + k11*ligf)) + \
    8.*l20*rtotf*(1 + k21*ligf + k21*k22*(np.square(ligf))))))))
    
    rfreebottom = np.concatenate(4*l20*(1 + k21*ligf + k21*k22*(np.square(ligf))))
    
    rfree = rfreetop/rfreebottom
       
    
    bfrac = (k11*ligc + l20*k21*rfree*ligc + \
    2*l20*k21*k22*rfree*(np.square(ligc)))/(1 + 2*l20*rfree + k11*ligc + \
    2*l20*k21*rfree*ligc + 2*l20*k21*k22*rfree*(np.square(ligc)))
    
    return bfrac





wyman_sim(ligset,rtots,sim_parms)

# coding: utf-8

# In[1]:

import numpy as np


# In[2]:

def dilser(low=0.001, limit=100.0, dilfactor = 2.0):
    """Returns a list containing a dilution series that ranges from
    "low" to "limit" by "dilfactor".
    """
    a = [low]
    while a[-1] <= limit:
        a.append(a[len(a)-1]*dilfactor)
    
    return np.array(a)


# In[3]:

lig = dilser()


# In[4]:

a = np.delete(lig,[15,16])


# In[5]:

b = np.delete(lig,[14,17])


# In[6]:

c = np.delete(lig,[13,15])


# In[7]:

ligset = np.array([lig,a,lig,b,c])


# In[8]:

noiseset = np.array([np.random.normal(1,0.05,ligset[i].size) for i in range(ligset.shape[0])])


# In[9]:

rtots = np.array([0.001,0.005, 0.01, 0.02, 0.05])


# In[10]:

sim_parms = dict({'k11':3.7, 'k21':1.8, 'k22':0.12, 'l20':293.0})


# In[11]:
@profile
def wyman_sim(lig,rtot,parm):    
    k11 = parm['k11']
    k21 = parm['k21']
    k22 = parm['k22']
    l20 = parm['l20']
    
    holder = []
    
    for i in range(lig.size):
        
        ligf = ligset
        rtotf = rtots
        
        rfreetop = ((-1 - k11*ligf) + \
        (np.sqrt(np.square((1 + k11*ligf)) + \
        8.*l20*rtotf*(1 + k21*ligf + k21*k22*(np.square(ligf))))))
        
        rfreebottom = (4*l20*(1 + k21*ligf + k21*k22*(np.square(ligf))))
        
        rfree = rfreetop/rfreebottom
        
        bfrac = (k11*ligf + l20*k21*rfree*ligf + \
        2*l20*k21*k22*rfree*(ligf**2))/(1 + 2*l20*rfree + k11*ligf + \
        2*l20*k21*rfree*ligf + 2*l20*k21*k22*rfree*(ligf**2))
        
        holder.append(bfrac)
        
    return np.array(holder)
wyman_sim(ligset,rtots,sim_parms)



import lmfit


# In[23]:

parms = lmfit.Parameters()


# In[24]:

parms.add('k11', value=10., min=0.01)
parms.add('k21', value=10., min=0.01)
parms.add('k22', value=10., min=0.01)
parms.add('l20', value=100., min=0.01)

def wym_obj2(parm,lig,rtot,data,eps=None):
    """
    parms
    --------
    lmfit.Parameters() object containin definitions for k11, k21, k22, and l20
    
    x
    --------
    array of arrays (probably assymetric elements without "formal" dimensions) containing 
    free ligand concentrations, must be of same size (no. of elements) as rtots, and same "shape"
    as data (same number of arrays within the array and the same number of elements in each subarray)
    
    Continue documentation
    """
    
    k11 = parm['k11']
    k21 = parm['k21']
    k22 = parm['k22']
    l20 = parm['l20']
    
    model = []
    
    for i in range(lig.size):
        rfree = (((-1 - k11*lig[i]) +         (np.sqrt((1 + k11*lig[i])**2 + 8*l20*rtot[i]*(1 + k21*lig[i] +         k21*k22*(lig[i]**2)))))/(4*l20*(1 + k21*lig[i] + k21*k22*(lig[i]**2))))
        
        bfrac = (k11*lig[i] + l20*k21*rfree*lig[i] +         2*l20*k21*k22*rfree*(lig[i]**2))/(1 + 2*l20*rfree + k11*lig[i] +         2*l20*k21*rfree*lig[i] + 2*l20*k21*k22*rfree*(lig[i]**2))
        
        model.append(bfrac)
        
    if eps is None:
        return np.concatenate((np.array(model)) - data)
    else:
        weights = 1/(np.square(eps))
        return np.concatenate(((np.array(model)) - data)*weights)



import numpy as np
from scipy import optimize

def dilser(low=0.001, limit=100., dilfactor=2.):
    '''returns a numpy array dilution series from low to limit'''
    #replace this with a generator/iterator someday
    a = [low]
    while a[-1] <= limit:
        a.append(a[len(a)-1]*dilfactor)
    return np.array(a)


def lfunc(model,data,err=None):
    if err is None:
        return 

    
    
def fit(self,err=None,**kwargs):
    '''Wraps scipy optimize routine and uses data/parameters from model object.
    If weights are included, a weighted fit is performed'''
    if err is None:
        return [optimize.least_squares(self.fitfunc,self.guess,bounds=self.bounds, \
                args=(self.model.ligs,self.model.meanset[i],self.model.rtot),**kwargs) \
                for i in range(self.model.meanset.shape[0])]
    else:
        return [optimize.least_squares(self.fitfunc,self.guess,bounds=self.bounds, \
                args=(self.model.ligs,self.model.meanset[i],self.model.rtot,err[i]),**kwargs) \
                for i in range(self.model.meanset.shape[0])]
    
def fitwrap(self,weight,**kwargs):
    if weight == 1:
        return optwrap(self,err=self.model.stdset,**kwargs)
    elif weight == 2: 
        return optwrap(self,err=(self.model.noise*self.model.meanset),**kwargs)
    elif weight == 3:
        return optwrap(self,err=(self.model.meanset),**kwargs)
    else:
        return optwrap(self,**kwargs)



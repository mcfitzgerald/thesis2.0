# -*- coding: utf-8 -*-
from . import helpers
from . import models
import numpy as np

class WymSim:
    
    #default init values
    __dparms = np.array([3.,2.,0.1,100.])
    __drtot = np.array([0.001,0.0025,0.005,0.01,0.025,0.05])
    __dlrange = np.array([0.001,100.,2.])
    __drp = 3
    __dst = 3333
    __dns = 0.05
    __dmodfunc = models.wymfunc
    
    def __init__(self,modfunc=__dmodfunc,parms=__dparms,rtot=__drtot,ligrange=__dlrange,reps=__drp,sets=__dst,noise=__dns):
        '''Returns a Wyman model object'''
        self.__sets = sets
        self.__reps = reps
        self.model = model
        self.parms = parms
        self.rtot = rtot
        self.noise = noise
        self.ligs = np.array([helpers.dilser(ligrange[0],ligrange[1],ligrange[2]) for i in range(len(self.rtot))])
        self.bfrac = np.array([self.modfunc(self.parms,self.ligs[i],self.rtot[i]) for i in range(len(self.rtot))])
        self.__noised = np.array([[np.random.normal(self.bfrac,(self.noise*self.bfrac)) for i in range(self.__reps)] for i in range(self.__sets)])
        self.meanset = np.array([self.__noised[i].mean(axis=0) for i in range(len(self.__noised))])
        self.stdset = np.array([self.__noised[i].std(axis=0) for i in range(len(self.__noised))])
        
        
class WymFit:
    dguess = [10.,10.,10.,100.]
    dbounds = ((0.,0.,0.,0.,),(100.,100.,100.,10000.))
    
    def __init__(self, model, guess=dguess, bounds=dbounds, weight=0, **kwargs):
        '''Takes Wyman Model Object'''
        self.model = model
        self.guess = np.array(guess)
        self.bounds = bounds
        self.fits = WymanSimFit.fitwrap(self,weight,**kwargs)
        self.ests = np.array([self.fits[i].x for i in range(len(self.fits))])
        self.k11 = self.ests[:,0]
        self.k21 = self.ests[:,1]
        self.k22 = self.ests[:,2]
        self.l20 = self.ests[:,3]
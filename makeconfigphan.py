#!/usr/bin/env python
"""
Created on Fri Apr 24 13:13:39 2015

@author: John Swoboda
"""

import numpy as np
from SimISR.makeConfigFiles import makepicklefile

def main():

    beamlist = np.loadtxt('spcorbco.txt')
    beamlist.astype(np.int)
    radarname = 'pfisr'

    Tint=4.0*60.0
    time_lim = 900+Tint
    pulse = np.ones(14)
    rng_lims = [150,500]
    IPP = .0087
    NNs = 28
    NNp = 100
    simparams =   {'IPP':IPP,
                   'TimeLim':time_lim,
                   'RangeLims':rng_lims,
                   'Pulse':pulse,
                   'Pulsetype':'long',
                   'Tint':Tint,
                   'Fitinter':Tint,
                   'NNs': NNs,
                   'NNp':NNp,
                   'dtype':np.complex128,
                   'ambupsamp':30,
                   'species':['O+','NO+','O2+','e-'],
                   'numpoints':128,
                   'SUMRULE': np.array([[-2,-3,-3,-4,-4,-5,-5,-6,-6,-7,-7,-8,-8,-9]
                       ,[1,1,2,2,3,3,4,4,5,5,6,6,7,7]])}

    fname = 'PFISRphantomprocspcor'


    makepicklefile(fname+'.pickle',beamlist,radarname,simparams)
if __name__ == "__main__":
    main()
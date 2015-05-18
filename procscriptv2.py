#!/usr/bin/env python


import os, inspect, glob,pdb
import scipy as sp
import numpy as np
from RadarDataSim.utilFunctions import makepicklefile,GenBarker
import RadarDataSim.runsim as runsim
def makepicklelongpulse(filepath):
    beamlist = np.loadtxt('spcorbco.txt')
    beamlist.astype(np.int)
    radarname = 'pfisr'

    Tint=4.0*60.0
    time_lim = 900+Tint
    pulse = GenBarker(13)
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
                   'startfile':os.path.join(filepath,'avedata.h5'),
                   'SUMRULE': np.array([[-2,-3,-3,-4,-4,-5,-5,-6,-6,-7,-7,-8,-8,-9]
                       ,[1,1,2,2,3,3,4,4,5,5,6,6,7,7]])}

    fname = os.path.join(filepath,'PFISRphantomprocspcor')


    makepicklefile(fname+'.pickle',beamlist,radarname,simparams)
def makepicklebarkercode(filepath):
    beamlist = np.loadtxt('spcorbco.txt')
    beamlist.astype(np.int)
    radarname = 'pfisr'

    Tint=4.0*60.0
    time_lim = 900+Tint
    pulse = np.ones(14)
    rng_lims = [80,250]
    IPP = .0087
    NNs = 28
    NNp = 100
    simparams =   {'IPP':IPP,
                   'TimeLim':time_lim,
                   'RangeLims':rng_lims,
                   'Pulse':pulse,
                   'Pulsetype':'barker',
                   'Tint':Tint,
                   'Fitinter':Tint,
                   'NNs': NNs,
                   'NNp':NNp,
                   'dtype':np.complex128,
                   'ambupsamp':30,
                   'startfile':os.path.join(filepath,'avedata.h5'),
                   'species':['O+','NO+','O2+','e-'],
                   'numpoints':128,
                   'SUMRULE': np.array([[0],[0]])}

    fname = os.path.join(filepath,'PFISRphantombarker')
    makepicklefile(fname+'.pickle',beamlist,radarname,simparams)
    
def makepickleline(filepath):
    beamlist = np.loadtxt('linelist.txt')
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
                   'startfile':os.path.join(filepath,'avedata.h5'),
                   'SUMRULE': np.array([[-2,-3,-3,-4,-4,-5,-5,-6,-6,-7,-7,-8,-8,-9]
                       ,[1,1,2,2,3,3,4,4,5,5,6,6,7,7]])}

    fname = os.path.join(filepath,'PFISRphantomprocline')
    makepicklefile(fname+'.pickle',beamlist,radarname,simparams)
    
def runstuff(datapath,picklefilename):
    funcnamelist=['spectrums','radardata','fitting']
    runsim.main(funcnamelist,datapath,os.path.join(datapath,picklefilename),True)
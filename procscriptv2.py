#!/usr/bin/env python


import os, inspect, glob,pdb
import scipy as sp
import numpy as np
from RadarDataSim.utilFunctions import makeconfigfile
from RadarDataSim.IonoContainer import IonoContainer
import RadarDataSim.runsim as runsim
from turn2geodata import fit2geodata
def makepicklelongpulse(filepath):
    beamlist = np.loadtxt('spcorbco.txt')
    beamlist.astype(np.int)
    radarname = 'pfisr'

    Tint=5.25*60.0# will give around 300 pulses
    fittedtimeint = 60.0 # Each incriment will be

    time_lim = 900+Tint
    rng_lims = [100,600]
    IPP = .0087
    NNs = 28
    NNp = 100
    simparams =   {'IPP':IPP,
                   'TimeLim':time_lim,
                   'RangeLims':rng_lims,
                   'Pulselength':280e-6,
                   'Pulsetype':'long', # type of pulse can be long or barker,
                   't_s': 20e-6,
                   'Tint':Tint,
                   'Fitinter':fittedtimeint,
                   'NNs': NNs,
                   'NNp':NNp,
                   'dtype':np.complex128,
                   'ambupsamp':30,
                   'species':['O+','NO+','O2+','e-'],
                   'numpoints':128,
                   'startfile':os.path.join(filepath,'avedata.h5')}

    fname = os.path.join(filepath,'PFISRphantomprocspcor.ini')


    makeconfigfile(fname+'.pickle',beamlist,radarname,simparams)
def makepicklebarkercode(filepath):
    beamlist = np.loadtxt('spcorbco.txt')
    beamlist.astype(np.int)
    radarname = 'pfisr'

    Tint=120 # 113 pulses
    time_lim = 900+Tint
    pulse = np.ones(14)
    rng_lims = [80,250]
    IPP = .0087
    NNs = 28
    NNp = 100
    simparams =   {'IPP':IPP,
                   'TimeLim':time_lim,
                   'RangeLims':rng_lims,
                   'Pulselength':140e-6,
                   't_s': 10e-6,
                   'Pulsetype':'barker',
                   'Tint':Tint,
                   'Fitinter':Tint,
                   'NNs': NNs,
                   'NNp':NNp,
                   'dtype':np.complex128,
                   'ambupsamp':30,
                   'startfile':os.path.join(filepath,'avedata.h5'),
                   'species':['O+','NO+','O2+','e-'],
                   'numpoints':128}

    fname = os.path.join(filepath,'PFISRphantombarker')
    makeconfigfile(fname+'.pickle',beamlist,radarname,simparams)

def makepickleline(filepath):
    IPP = .0087
    NNs = 28
    NNp = 100

    beamlist = np.loadtxt('linelist.txt')
    beamlist.astype(np.int)
    radarname = 'pfisr'
    nbeams = len(beamlist)

    Tint=30.0 # integrates right around 200 pulses
    fittedtimeint = 30.0 # Each incriment will be
    time_lim = 900+Tint
    rng_lims = [100,600]

    simparams =   {'IPP':IPP,
                   'TimeLim':time_lim,
                   'RangeLims':rng_lims,
                   'Pulselength':280e-6,
                   't_s': 20e-6,
                   'Pulsetype':'long',
                   'Tint':Tint,
                   'Fitinter':fittedtimeint,
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
    makeconfigfile(fname+'.pickle',beamlist,radarname,simparams)

def reducedata(inputdir,newcoords):

    dirlist = glob.glob(os.path.join(inputdir,'*.mat'))
    numlist = [os.path.splitext(os.path.split(x)[-1])[0] for x in dirlist]
    numdict = {numlist[i]:dirlist[i] for i in range(len(dirlist))}
    slist = sorted(numlist,key=ke)

    #coordlims = {'x':[0,300],'y':[0,400],'z':[0,700]}
    for inum in slist:

        curfile = numdict[inum]
        curiono = IonoContainer.readmat(curfile)
        if curiono.Time_Vector[0]==1e-6:
            curiono.Time_Vector[0] = 0.0

        if type(newcoords)==dict:
            curiono.coordreduce(newcoords)
        elif type(newcoords)==np.ndarray:
            curiono.interp(newcoords)
        curiono.saveh5(os.path.join(inputdir,inum+' red.h5'))
#        curiono.makespectruminstanceopen(specfuncs.ISRSspecmake,sensdict,npts).saveh5(outfile)
#%% For stuff
def ke(item):
    if item[0].isdigit():
        return int(item.partition(' ')[0])
    else:
        return float('inf')

def runstuff(datapath,picklefilename):
    funcnamelist=['spectrums','radardata','fitting']
    runsim.main(funcnamelist,datapath,os.path.join(datapath,picklefilename),True)
    fit2geodata(os.path.join(datapath,'Fitted','fitteddata.h5'))
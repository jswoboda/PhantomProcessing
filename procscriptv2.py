#!/usr/bin/env python


import os, inspect, glob,pdb
import scipy as sp
import numpy as np
from RadarDataSim.utilFunctions import makepicklefile,GenBarker
from RadarDataSim.IonoContainer import IonoContainer
import RadarDataSim.runsim as runsim
from turn2geodata import fit2geodata
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
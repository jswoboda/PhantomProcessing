#!/usr/bin/env python
"""
Created on Mon Mar 16 12:14:42 2015

@author: John Swoboda
"""
from __future__ import print_function

#imported basic modules
import os, inspect, time, sys, getopt, glob
from datetime import datetime
import traceback
import pdb
# Imported scipy and matplotlib modules
import scipy as sp
import scipy.io as sio
from matplotlib import rc
import matplotlib.pylab as plt
import tables
# My modules
from RadarDataSim.IonoContainer import IonoContainer
from RadarDataSim.radarData import RadarData, RadarDataFile
import RadarDataSim.specfunctions as specfuncs
import RadarDataSim.const.sensorConstants as sensconst
from RadarDataSim.const.physConstants import v_C_0, v_Boltz
from beamtools.bcotools import getangles
from RadarDataSim.utilFunctions import make_amb
from RadarDataSim.specfunctions import ISRSfitfunction
from RadarDataSim.fitterMethodGen import Fitterionoconainer

def makespectrums(inputdir,outputdir,optinputs):

    dirlist = glob.glob(os.path.join(inputdir,'*.mat'))
    numlist = [os.path.splitext(os.path.split(x)[-1])[0] for x in dirlist]
    numdict = {numlist[i]:dirlist[i] for i in range(len(dirlist))}
    slist = sorted(numlist,key=ke)

    sensdict = optinputs[0]

    if len(optinputs)<2:
        npts = 128
    else:
        npts = int(optinputs[1])
    coordlims = {'x':[-300,300],'y':[-300,300],'z':[0,700]}
    for inum in slist:

        outfile = os.path.join(outputdir,inum+' spectrum.h5')
        curfile = numdict[inum]
        print('Processing file {} starting at {}\n'.format(os.path.split(curfile)[1],datetime.now()))
        curiono = IonoContainer.readmat(curfile)
        curiono.coordreduce(coordlims)
        curiono.makespectruminstanceopen(specfuncs.ISRSspecmake,sensdict,npts).saveh5(outfile)
        print('Finished file {} starting at {}\n'.format(os.path.split(curfile)[1],datetime.now()))
def makeradardata(inputdir,outputdir,optinputs):

    pulse = sp.ones(14)
    rng_lims = [150,500]
    IPP = .0087
    angles = getangles('spcorbco.txt')
    ang_data = sp.array([[iout[0],iout[1]] for iout in angles])
    sensdict = sensconst.getConst('pfisr',ang_data)
    Tint=3.0*60.0
    time_lim = 900.0+Tint

    NNs = 28
    NNp = 100
    rng_gates = sp.arange(rng_lims[0],rng_lims[1],sensdict['t_s']*v_C_0*1e-3)
    simparams =   {'IPP':IPP,'angles':angles,'TimeLim':time_lim,'Pulse':pulse,\
    'Timevec':sp.arange(0,time_lim,Tint),'Tint':Tint,'Rangegates':rng_gates,\
    'Noisesamples': NNs,'Noisepulses':NNp,'dtype':sp.complex128}
    simparams['SUMRULE'] = sp.array([[-2,-3,-3,-4,-4,-5,-5,-6,-6,-7,-7,-8,-8,-9],[1,1,2,2,3,3,4,4,5,5,6,6,7,7]])
    simparams['amb_dict'] = make_amb(sensdict['fs'],30,sensdict['t_s']*len(pulse),len(pulse))

    dirlist = glob.glob(os.path.join(inputdir,'*.h5'))
    filelist = [os.path.split(item)[1] for item in dirlist]
    timelist = [int(item.partition(' ')[0]) for item in filelist]
    Ionodict = {timelist[it]:dirlist[it] for it in range(len(dirlist))}

    dirlist2 = glob.glob(os.path.join(outputdir,'*.h5'))
    if dirlist2:
        numlist2 = [os.path.splitext(os.path.split(x)[-1])[0] for x in dirlist2]
        numdict2 = {numlist2[i]:dirlist2[i] for i in range(len(dirlist2))}
        slist2 = sorted(numlist2,key=ke)
        outlist2 = [numdict2[ikey] for ikey in slist2]
    else:
        outlist2 = None

    rdata = RadarDataFile(Ionodict,sensdict,simparams,outputdir,outfilelist=outlist2)
    timearr = sp.linspace(0.0,time_lim,num=220)
    ionoout = rdata.processdataiono(timearr,Tint)
    ionoout.readh5(os.path.join(outputdir,'00lags.h5'))
    (DataLags,NoiseLags) = rdata.processdata(timearr,Tint)
    sio.savemat(os.path.join(outputdir,'ACFdata.mat'),mdict=DataLags)
    sio.savemat(os.path.join(outputdir,'Noisedata.mat'),mdict=NoiseLags)
    return ()
def fitdata(inputdir,outputdir,optinputs):
    dirlist = glob.glob(os.path.join(inputdir,'*lags.h5'))

    # Species in matlab format 1=O+,2=NO+,3=N2+,4=O2+,5=N+, 6=H+,7=e-
    # kept 1,2,4,6,7
    pulse = sp.ones(14)
    rng_lims = [150,500]
    IPP = .0087
    angles = getangles('spcorbco.txt')
    ang_data = sp.array([[iout[0],iout[1]] for iout in angles])
    sensdict = sensconst.getConst('pfisr',ang_data)
    Tint=3.0*60.0
    time_lim = 900.0+Tint

    NNs = 28
    NNp = 100
    rng_gates = sp.arange(rng_lims[0],rng_lims[1],sensdict['t_s']*v_C_0*1e-3)
    simparams =   {'IPP':IPP,'angles':angles,'TimeLim':time_lim,'Pulse':pulse,\
    'Timevec':sp.arange(0,time_lim,Tint),'Tint':Tint,'Rangegates':rng_gates,\
    'Noisesamples': NNs,'Noisepulses':NNp,'dtype':sp.complex128}
    simparams['SUMRULE'] = sp.array([[-2,-3,-3,-4,-4,-5,-5,-6,-6,-7,-7,-8,-8,-9],[1,1,2,2,3,3,4,4,5,5,6,6,7,7]])
    simparams['amb_dict'] = make_amb(sensdict['fs'],30,sensdict['t_s']*len(pulse),len(pulse))

    species = ['O+','NO+','O2+','e-']
    sensdict['species'] = species
    Ionoin=IonoContainer.readh5(dirlist[0])
    fitterone = Fitterionoconainer(Ionoin,sensdict,simparams)
    (fitteddata,fittederror) = fitterone.fitdata(ISRSfitfunction,startvalfunc)
    (Nloc,Ntimes,nparams)=fitteddata.shape
    fittederronly = fittederror[:,:,range(nparams),range(nparams)]
    paramlist = sp.concatenate((fitteddata,fittederronly),axis=2)
    paramnames = []
    for isp in species[:-1]:
        paramnames.append('Ni_'+isp)
        paramnames.append('Ti_'+isp)
    paramnames = paramnames+['Ne','Te','Vi']
    paramnamese = ['n'+ip for ip in paramnames]
    paranamsf = sp.array(paramnames+paramnamese)


    Ionoout=IonoContainer(Ionoin.Sphere_Coords,paramlist,Ionoin.Time_Vector,coordvecs = Ionoin.Coord_Vecs, paramnames=paranamsf,species=species)

    Ionoout.saveh5(os.path.join(outputdir,'fitteddatasphere.h5'))

    return ()
def startvalfunc(Ne_init, loc,time,exinputs):
    """ """
    h5file = tables.openFile('avedata.h5')
    zdata = h5file.root.zkm.read()
    datast = h5file.root.stdata.read()
    vel = h5file.root.velocity.read()
    numel =sp.prod(datast[-2:]) +1

    xarray = sp.zeros((loc.shape[0],numel))
    for ilocn, iloc in enumerate(loc):
        indx = sp.argmin(sp.absolute(zdata-iloc[2]))
        xarray[ilocn,:-1]=sp.reshape(datast[indx,0],numel-1)
        locmag = sp.sqrt(sp.sum(iloc*iloc))
        xarray[ilocn,-1] = vel[indx,0]*iloc/locmag

    xarray = sp.repeat(xarray[:,sp.newaxis,:],(1,len(time),1))


    return xarray

#%% For stuff
def ke(item):
    if item[0].isdigit():
        return int(item.partition(' ')[0])
    else:
        return float('inf')


if __name__ == "__main__":
    dfilename = 'diary.txt'
    inputsep = '***************************************************************\n'
    argv = sys.argv[1:]
    outstr = 'procscript.py -f <function: spectrums radardata or fitting> -i <inputdir> -o <outputdir> -n <number of samples>'
    funcdict = {'spectrums':makespectrums, 'radardata':makeradardata, 'fitting':fitdata}
    #pdb.set_trace()
    try:
        opts, args = getopt.getopt(argv,"hf:i:o:n:")
    except getopt.GetoptError:
        print(outstr)
        sys.exit(2)


    outdirexist = False
    for opt, arg in opts:
        if opt == '-h':
            print(outstr)
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputdir = arg
        elif opt in ("-o", "--ofile"):
            outdirexist = True
            outdir = arg
        elif opt in ("-f", "--func"):
            curfunc = funcdict[arg]
	elif opt in ('-n', "--num"):
	    npts = int(arg)

    if not outdirexist:
        outdir = inputdir
    sensdict = sensconst.getConst('risr')
    full_path = os.path.realpath(__file__)
    path, file = os.path.split(full_path)

    dfullfilestr = os.path.join(path,dfilename)
    f= open(dfullfilestr,'a')
    f.write(inputsep)
    f.write(curfunc.__name__+'\n')
    f.write(time.asctime()+'\n')

    try:
        stime = datetime.now()
        curfunc(inputdir,outdir,[sensdict,npts])
        ftime = datetime.now()
        ptime = ftime-stime
        f.write('Success!\n')
        f.write('Duration: {}\n'.format(ptime))
        f.write('Input directory: {}\n'.format(inputdir))
        f.write('Output directory: {}\n'.format(outdir))
    except Exception, e:
        f.write('Failed!\n')
        ftime = datetime.now()
        ptime = ftime-stime
        f.write('Duration: {}\n'.format(ptime))
        f.write('Input directory: {}\n'.format(inputdir))
        f.write('Output directory: {}\n'.format(outdir))
        traceback.print_exc(file=sys.stdout)
        traceback.print_exc(file = f)
        #pdb.set_trace()
    f.write(inputsep)
    f.close()



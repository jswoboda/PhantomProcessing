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
# My modules
from RadarDataSim.IonoContainer import IonoContainer
from RadarDataSim.radarData import RadarData, RadarDataFile
import RadarDataSim.specfunctions as specfuncs
import RadarDataSim.const.sensorConstants as sensconst
from RadarDataSim.const.physConstants import v_C_0, v_Boltz
from beamtools.bcotools import getangles
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
    sensdict = optinputs[0]
    pulse = sp.ones(14)
    rng_lims = [150,500]
    IPP = .0087
    angles = getangles('spcorbco.txt')
    ang_data = sp.array([[iout[0],iout[1]] for iout in angles])
    sensdict = sensconst.getConst('risr',ang_data)
    Tint=3.0*60.0
    time_lim = 900.0+Tint

    NNs = 28
    NNp = 100
    rng_gates = sp.arange(rng_lims[0],rng_lims[1],sensdict['t_s']*v_C_0*1e-3)
    simparams =   {'IPP':IPP,'angles':angles,'TimeLim':time_lim,'Pulse':pulse,\
    'Timevec':sp.arange(0,time_lim,Tint),'Tint':Tint,'Rangegates':rng_gates,\
    'Noisesamples': NNs,'Noisepulses':NNp,'dtype':sp.complex128}
    dirlist = glob.glob(os.path.join(inputdir,'*.h5'))
    filelist = [os.path.split(item)[1] for item in dirlist]
    timelist = [int(item.partition(' ')[0]) for item in filelist]
    Ionodict = {timelist[it]:dirlist[it] for it in range(len(dirlist))}

    dirlist2 = glob.glob(os.path.join(outputdir,'*.h5'))
    numlist2 = [os.path.splitext(os.path.split(x)[-1])[0] for x in dirlist2]
    numdict2 = {numlist2[i]:dirlist2[i] for i in range(len(dirlist2))}
    slist2 = sorted(numlist2,key=ke)
    outlist2 = [numdict2[ikey] for ikey in slist2]

    rdata = RadarDataFile(Ionodict,sensdict,simparams,outputdir,outfilelist=None)
    timearr = sp.linspace(0.0,time_lim,num=220)
    (DataLags,NoiseLags) = rdata.processdata(timearr,Tint)
    sio.savemat(os.path.join(outputdir,'ACFdata.mat'),mdict=DataLags)
    sio.savemat(os.path.join(outputdir,'Noisedata.mat'),mdict=NoiseLags)
    return ()
def fitdata(inputdir,outputdir,optinputs):
    return ()
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



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
from matplotlib import rc
import matplotlib.pylab as plt
# My modules
from RadarDataSim.IonoContainer import IonoContainer
from RadarDataSim.radarData import RadarData
import RadarDataSim.specfunctions as specfuncs
import RadarDataSim.const.sensorConstants as sensconst

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

    for inum in slist:
        outfile = os.path.join(outputdir,inum+'  spectrum.h5')
        curfile = numdict[inum]
        curiono = IonoContainer.readmat(curfile)
        curiono.makespectruminstanceopen(specfuncs.ISRSspecmake,sensdict,npts).saveh5(outfile)

def makeradardata(inputdir,outputdir,optinputs):
    return () 
#def fittdata(inputdir,outputdir,optinputs):

def ke(item):
    if item[0].isdigit():
        return int(item.partition(' ')[0])
    else:
        return float('inf')


if __name__ == "__main__":
    dfilename = 'diary.txt'
    inputsep = '***************************************************************\n'
    argv = sys.argv[1:]
    outstr = 'procscript.py -f <function: spectrums radardata or fitting> -i <inputdir> -o <outputdir>'
    funcdict = {'spectrums':makespectrums, 'radardata':makeradardata}
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
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

    if not outdirexist:
        outdir = inputdir
    sensdict = sensconst.getConst('risr')
    full_path = os.path.realpath(__file__)
    path, file = os.path.split(full_path)

    dfullfile = os.path.join(path,dfilename)
    print(inputsep,file=dfullfile)
    print(curfunc.__name__+'\n',file=dfullfile)
    print(time.asctime()+'\n',file=dfullfile)

    try:
        stime = datetime.now()
        curfunc(inputdir,outdir,[sensdict] +args)
        ftime = datetime.now()
        ptime = ftime-stime
        print('Success!\n',file=dfullfile)
        print('Duration: {}'.format(ptime),file=dfullfile)
    except Exception, e:
        print('Failed!\n',file=dfullfile)
        traceback.print_exc(file=sys.stdout)
        traceback.print_exc(file = dfullfile)
    print(inputsep,file=dfullfile)




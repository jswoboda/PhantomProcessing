#!/usr/bin/env python
"""
Created on Mon Mar 16 12:14:42 2015

@author: John Swoboda
"""
from __future__ import print_function

#imported basic modules
import os, time, sys, getopt, glob
from datetime import datetime
import traceback
import pdb
# Imported scipy and matplotlib modules
import scipy as sp
import tables
# My modules
from RadarDataSim.IonoContainer import IonoContainer
from RadarDataSim.radarData import RadarDataFile
import RadarDataSim.specfunctions as specfuncs
from RadarDataSim.specfunctions import ISRSfitfunction
from RadarDataSim.fitterMethodGen import Fitterionoconainer
from RadarDataSim.makeConfigFiles import readconfigfile
from turn2geodata import fit2geodata

#%% Make spectrums
def makespectrums(inputdir,outputdir,configfile,remakealldata):

    dirlist = glob.glob(os.path.join(inputdir,'* red.h5'))
    numlist = [os.path.splitext(os.path.split(x)[-1])[0] for x in dirlist]
    numdict = {numlist[i]:dirlist[i] for i in range(len(dirlist))}
    slist = sorted(numlist,key=ke)

    (sensdict,simparams) = readconfigfile(configfile)

#    coordlims = {'x':[-300,300],'y':[-300,300],'z':[0,700]}
    for inum in slist:

        outfile = os.path.join(outputdir,inum+' spectrum.h5')
        curfile = numdict[inum]
        print('Processing file {} starting at {}\n'.format(os.path.split(curfile)[1]
            ,datetime.now()))
        curiono = IonoContainer.readh5(curfile)
        if curiono.Time_Vector[0]==1e-6:
            curiono.Time_Vector[0] = 0.0
#        curiono.coordreduce(coordlims)
#        curiono.saveh5(os.path.join(inputdir,inum+' red.h5'))
        curiono.makespectruminstanceopen(specfuncs.ISRSspecmake,sensdict,
                                     simparams['numpoints']).saveh5(outfile)
        print('Finished file {} starting at {}\n'.format(os.path.split(curfile)[1],datetime.now()))

#%% Make Radar Data
def makeradardata(inputdir,outputdir,configfile,remakealldata):

    dirlist = glob.glob(os.path.join(inputdir,'*.h5'))
    filelist = [os.path.split(item)[1] for item in dirlist]
    timelist = [int(item.partition(' ')[0]) for item in filelist]
    Ionodict = {timelist[it]:dirlist[it] for it in range(len(dirlist))}

    radardatalist = glob.glob(os.path.join(outputdir,'*RawData.h5'))
    if radardatalist and (not remakealldata):
        numlist2 = [os.path.splitext(os.path.split(x)[-1])[0] for x in radardatalist]
        numdict2 = {numlist2[i]:radardatalist[i] for i in range(len(radardatalist))}
        slist2 = sorted(numlist2,key=ke)
        outlist2 = [numdict2[ikey] for ikey in slist2]
    else:
        outlist2 = None

    rdata = RadarDataFile(Ionodict,configfile,outputdir,outfilelist=outlist2)
    ionoout = rdata.processdataiono()
    ionoout.saveh5(os.path.join(outputdir,'00lags.h5'))

    return ()
#%% Fitt data
def fitdata(inputdir,outputdir,configfile,optintputs):
    dirlist = glob.glob(os.path.join(inputdir,'*lags.h5'))


    Ionoin=IonoContainer.readh5(dirlist[0])
    fitterone = Fitterionoconainer(Ionoin,configfile)
    (fitteddata,fittederror) = fitterone.fitdata(ISRSfitfunction,startvalfunc)
    (Nloc,Ntimes,nparams)=fitteddata.shape
    fittederronly = fittederror[:,:,range(nparams),range(nparams)]
    paramlist = sp.concatenate((fitteddata,fittederronly),axis=2)

    paramnames = []
    species = readconfigfile(configfile)[1]['species']
    for isp in species[:-1]:
        paramnames.append('Ni_'+isp)
        paramnames.append('Ti_'+isp)
    paramnames = paramnames+['Ne','Te','Vi']
    paramnamese = ['n'+ip for ip in paramnames]
    paranamsf = sp.array(paramnames+paramnamese)


    Ionoout=IonoContainer(Ionoin.Sphere_Coords,paramlist,Ionoin.Time_Vector,ver =1,coordvecs = Ionoin.Coord_Vecs, paramnames=paranamsf,species=species)
    outfile = os.path.join(outputdir,'fitteddataspherepenalty.h5')
    Ionoout.saveh5(outfile)
    fit2geodata(outfile)

#%% fit function stuff
def startvalfunc(Ne_init, loc,time,exinputs):
    """ """
    h5file = tables.openFile('avedata.h5')
    zdata = h5file.root.zkm.read()
    datast = h5file.root.stdata.read()
    vel = h5file.root.velocity.read()
    numel =sp.prod(datast.shape[-2:]) +1

    xarray = sp.zeros((loc.shape[0],numel))
    for ilocn, iloc in enumerate(loc):
        indx = sp.argmin(sp.absolute(zdata-iloc[2]))
        xarray[ilocn,:-1]=sp.reshape(datast[indx,0],numel-1)
        locmag = sp.sqrt(sp.sum(iloc*iloc))
        xarray[ilocn,-1] = sp.sum(vel[indx,0]*iloc)/locmag
    xarray = sp.repeat(xarray[:,sp.newaxis,:],len(time),axis=1)
    h5file.close()

    return xarray

#%% For sorting
def ke(item):
    if item[0].isdigit():
        return int(item.partition(' ')[0])
    else:
        return float('inf')
#%% Main function
def main(argv):

    configfile = 'PFISRphantomprocspcor.pickle'

    inputsep = '***************************************************************\n'

    outstr = 'procscript.py -f <function: spectrums radardata or fitting> -i <inputdir> -o <outputdir> -r <type y to remake data>'
    funcdict = {'spectrums':makespectrums, 'radardata':makeradardata, 'fitting':fitdata}
    #pdb.set_trace()
    try:
        opts, args = getopt.getopt(argv,"hf:i:o:r:")
    except getopt.GetoptError:
        print(outstr)
        sys.exit(2)


    outdirexist = False
    remakealldata = False
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
        elif opt in ('-r', "--re"):
            if arg.lower() == 'y':
                remakealldata = True



    if not outdirexist:
        outdir = inputdir

    full_path = os.path.realpath(__file__)
    path, file = os.path.split(full_path)
    dfilename = 'diary'+curfunc.__name__+'.txt'
    dfullfilestr = os.path.join(path,dfilename)
    f= open(dfullfilestr,'a')
    f.write(inputsep)
    f.write(curfunc.__name__+'\n')
    f.write(time.asctime()+'\n')

    try:
        stime = datetime.now()
        curfunc(inputdir,outdir,configfile,remakealldata)
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
def runall(basedir):
    paramdir = os.path.join(basedir,'Origparams')
    specdir = os.path.join(basedir,'Spectrums')
    radardir = os.path.join(basedir,'Radardata')
    fitdir = os.path.join(basedir,'Fitted')
    #create spectrums
    input1 = ['-f','spectrums','-i',paramdir,'-o',specdir,'-r','y']
    main(input1)
    #create radardata
    input2 = ['-f','radardata','-i',specdir,'-o',radardir,'-r','y']
    main(input2)
    #create fitted data
    input3 = ['-f','fitting','-i',radardir,'-o',fitdir,'-r','y']
    main(input3)

    return()
if __name__ == "__main__":
    argv = sys.argv[1:]
    if argv[0] == '-b':
        runall(argv[1])
    else:
        main(argv)

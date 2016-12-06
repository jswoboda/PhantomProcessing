#!/usr/bin/env python
"""
Created on Tue Apr 14 18:06:56 2015

@author: John Swoboda
"""
import os, glob
from SimISR.IonoContainer import IonoContainer
from GeoData.GeoData import GeoData
from GeoData.utilityfuncs import readIono
if __name__ == "__main__":
    inputdir = '/Users/Bodangles/Documents/Python/PhantomProc/Origparams'
    outdir = '/Users/Bodangles/Documents/MATLAB/phantoms/Origparams'
    dirlist = glob.glob(os.path.join(inputdir,'*.mat'))
    keep = ['Ne','Te','Ti']
    for ifile in dirlist:
        Ionoin=IonoContainer.readmat(ifile)
        GD = GeoData(readIono,[Ionoin])
        for ikey in GD.data.keys():
            if ikey not in keep:
                del GD.data[ikey]
        fname = os.path.splitext(os.path.split(ifile)[1])[0]
        GD.write_h5(os.path.join(outdir,fname+'.h5'))
        print('Saved '+os.path.join(outdir,fname+'.h5'))
def fit2geodata(inputfile,outfile=None):
    Ionoin=IonoContainer.readh5(inputfile)
    keep = ['Ne','Te','Ti']
    GD=GeoData(readIono,[Ionoin])
    for ikey in GD.data.keys():
        if ikey not in keep:
            del GD.data[ikey]
    indir,fnameall=os.path.split(inputfile)
    fname = os.path.splitext(fnameall)[0]
    if outfile is None:
        outdir = indir
        outname = fname+'GEOD.h5'
    else:
        outdir, outname = os.path.split(outfile)
        
    GD.write_h5(os.path.join(outdir,outname))
    print('Saved '+os.path.join(outdir,outname))
    
def fit2geodatalist(inputfilelist,outfilelist=None):
    if outfilelist is None:
        outfilelist=[None]*len(inputfilelist)
    for ifiles in zip(inputfilelist,outfilelist):
        fit2geodata(ifiles[0],ifiles[1])
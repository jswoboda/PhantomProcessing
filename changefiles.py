# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 12:13:27 2015

@author: swoboj
"""

import os, glob
import pdb
from RadarDataSim.IonoContainer import IonoContainer
from GeoData.GeoData import GeoData
from GeoData.utilityfuncs import readIono

def makegeodata(filename,infiledir,outfiledir,outputtype='h5'):
    
    filelist = glob.glob(os.path.join(infiledir,filename))
    
    for ifile in filelist:
        (fbase,fext) = os.path.splitext(os.path.split(ifile)[-1])
        print('Reading '+ifile)
        if fext=='.h5':
            Ionoin = IonoContainer.readh5(ifile)
        elif fext=='.mat':
            Ionoin = IonoContainer.readmat(ifile)
        GD = GeoData(readIono,[Ionoin])
        
        datakeys =GD.data.keys()
        
        for ikey in datakeys:
            if ('Ti_' in ikey) or ('Ni_' in ikey):
                del GD.data[ikey]
        
        outfilename = os.path.join(outfiledir,fbase+'GEOD.'+outputtype)
        
        GD.write_h5(outfilename)
        print(outfilename+' saved')
        
        
        
            
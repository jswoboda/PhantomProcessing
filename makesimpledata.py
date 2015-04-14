#!/usr/bin/env python
"""
Created on Sun Apr 12 15:06:00 2015

@author: John Swoboda
"""

#!/usr/bin/env python
"""
Created on Mon Apr  6 15:07:16 2015

@author: John Swoboda
"""
import os, inspect
import pdb
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib import cm
import matplotlib as mpl
from matplotlib import ticker
from GeoData.GeoData import GeoData
from GeoData.utilityfuncs import readIono
from RadarDataSim.IonoContainer import IonoContainer
import tables
import scipy as sp

if __name__ == "__main__":
    curpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

    file_name = os.path.join(curpath,'Origparams','0.mat')

    Iono = IonoContainer.readmat(file_name)
    Param_List = Iono.Param_List
    dataloc = Iono.Cart_Coords
    velocity = Iono.Velocity
    zlist,idx = sp.unique(dataloc[:,2],return_inverse=True)
    siz = list(Param_List.shape[1:])
    vsiz = list(velocity.shape[1:])

    outdata = sp.zeros([len(zlist)]+siz)
    outvel = sp.zeros([len(zlist)]+vsiz)

    for izn,iz in enumerate(zlist):
        arr = sp.argwhere(idx==izn)
        outdata[izn]=sp.mean(Param_List[arr],axis=0)
        outvel[izn]=sp.mean(velocity[arr],axis=0)

    h5file=tables.openFile('avedata.h5',mode = "w", title = "Ave data")

    h5file.createArray('/', 'stdata',outdata,'Static array')
    h5file.createArray('/', 'zkm',zlist,'Static array')
    h5file.createArray('/', 'names',Iono.Param_Names,'Static array')
    h5file.createArray('/','velocity',outvel,'Static array')
    h5file.close()





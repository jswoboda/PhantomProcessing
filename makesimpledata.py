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
from RadarDataSim.IonoContainer import IonoContainer
import scipy as sp

def main(curpath):

    file_name = os.path.join(curpath,'Origparams','0.mat')

    Iono = IonoContainer.readmat(file_name)
    Param_List = Iono.Param_List
    dataloc = Iono.Cart_Coords
    velocity = Iono.Velocity
    zlist,idx = sp.unique(dataloc[:,2],return_inverse=True)
    siz = list(Param_List.shape[1:])
    vsiz = list(velocity.shape[1:])

    datalocsave = sp.column_stack((sp.zeros_like(zlist),sp.zeros_like(zlist),zlist))
    outdata = sp.zeros([len(zlist)]+siz)
    outvel = sp.zeros([len(zlist)]+vsiz)

    for izn,iz in enumerate(zlist):
        arr = sp.argwhere(idx==izn)
        outdata[izn]=sp.mean(Param_List[arr],axis=0)
        outvel[izn]=sp.mean(velocity[arr],axis=0)

    Ionoout = IonoContainer(datalocsave,outdata,Iono.Time_Vector,Iono.Sensor_loc,ver=0,paramnames=Iono.Param_Names,
                            species=Iono.Species,velocity=outvel)
    Ionoout.saveh5('avedata.h5')
if __name__ == "__main__":
    main()





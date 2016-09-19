# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 19:42:13 2016

@author: swoboj
"""
import os,glob
import scipy as sp
import GeoData.GeoData as GeoData
from GeoData.plotting import plot3Dslice
import GeoData.GeoData as GeoData

indir =  '/home/swoboj/DATA/Phantoms/longpulse/Fitted'
outdir = '~/Documents/MATLAB/phantoms/Fittedfigsaug'
curfilename = os.path.join(indir,'fitteddataGEODinterp.h5')
indir2 = '/home/swoboj/DATA/Phantoms/longpulse/OrigparamsGEOD/'
fullfprot = os.path.join(indir2,'*.h5')

fdir = glob.glob(fullfprot)
ford = sp.zeros(len(fdir))
for k,curfname in enumerate(fdir):
    
    
    GDtemp=GeoData.read_h5(curfname)
    if k==0:
        GDin=GDtemp
    else:
        GDin.add_times(GDtemp)
    



GDfit = GeoData.read_h5(curfilename)
xlist = [];
ylist = [];
zlist = [200,325,450]
figpos = [57,400,1250,400]
bounds = [[5e10,5e11],[1000,4000]]
key = ['Ne','Te']

keystr = [[key[0][0],'_',key[0][1]],[key[1][0],'_',key[1][1]]]
units = ['m$^{-3}$','$^{\circ}$K']
mvfilename = 'fitted.avi';
if not os.path.isdir(outdir):
    os.mkdir(outdir);
end

time_vec = GDfit.times

curfname = fdir[0]
fmain=os.path.split(curfname)[-1]
reffile = os.path.join(indir2,fmain)
GD2 = GeoData.read_h5(reffile)

coords = GD2.dataloc
xvec = sp.unique(coords[:,0])
xvec = xvec[(xvec>-150)&(xvec<150)]
yvec = sp.unique(coords[:,1])
yvec = yvec[(yvec>-150)&(yvec<150)]
zvec = sp.unique(coords[:,2]);
zvec = zvec[(zvec>0)&(zvec<500)]
[X,Y,Z] = sp.meshgrid(xvec,yvec,zvec)

new_coords = sp.column_stack([X.flatten(),Y.flatten(),Z.flatten()])
newcoordname= 'Cartisian';

[Xplane,Yplane] = sp.meshgrid(xvec,sp.linspace(-63,63,len(xvec)))
Zplane = -2*Yplane+325

timeregcell=GDin.timeregister(GDfit)
for i_in,curlist in enumerate(timeregcell):
    for i_out in curlist:
        
        coords=GDin.dataloc;
        [C,ia,ib] = intersect(coords,new_coords,'rows');
        GDin.dataloc=C;
        curfields = fieldnames(GDin.data);
        for m in range(len(curfields)):
            GDin.data[curfields[m]] = GDin.data[curfields[m]][ia,:]
        
        
        time=floor(GDin.times(i_in,1));
        time2=floor(GDfit.times(i_out,1));
        
        fmain = [ num2str(time,'%d'),'_',num2str(time2,'%d')];

        hfig = figure('Position',figpos,'Color',[1,1,1]);
        hax1 = subplot(1,2,1);
        
        hslice=sliceGD(GDin,xlist,ylist,zlist,'Fig',hfig,'key',key{1},'axh',hax1,...
            'time',i_in,'bounds',bounds{1});
        set(hslice,'FaceColor','interp','EdgeColor','none');axis tight
        hold all
        hslice=sliceGD(GDin,Xplane,Yplane,Zplane,'Fig',hfig,'key',key{1},'axh',hax1,...
            'time',i_in,'bounds',bounds{1});
        set(hslice,'FaceColor','interp','EdgeColor','none');axis tight
        view(100,20)

        c = colorbar;
        c.Label.String=[keystr{1},' in ' units{1}];
        c.Label.FontSize = 14;
        title(['Input ',keystr{1},' at t=',num2str(time),' s'],'FontSize',16)

        hax2 = subplot(1,2,2);
        hslice=sliceGD(GDfit,xlist,ylist,zlist,'Fig',hfig,'key',key{1},'axh',hax2,...
            'time',i_out,'bounds',bounds{1});
        hold all
        set(hslice,'FaceColor','interp','EdgeColor','none');axis tight
         hslice=sliceGD(GDfit,Xplane,Yplane,Zplane,'Fig',hfig,'key',key{1},'axh',hax2,...
            'time',i_out,'bounds',bounds{1});
        view(100,20)

        c = colorbar;
        c.Label.String=[keystr{1},' in ' units{1}];
        c.Label.FontSize = 14;
        title(['Output ',keystr{1},' at t=',num2str(time2),' s'],'FontSize',16)
        axis([hax1,hax2],[hax1.XLim,hax1.YLim,hax1.ZLim])

        drawnow
        saveas(hfig,fullfile(outdir,fmain),'fig');
        pause(3)
        saveas(hfig,fullfile(outdir,fmain),'png');
        close(hfig)
    end
end

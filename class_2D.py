# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 13:47:04 2023

@author: Alexander.Kurapov
"""

class timeSer_2D:
    
    def clipTime(self,dateSTR,dateEND):
        
        import numpy as np
        
        nt=self.dt64.size
        
        xmask = (self.dt64 >= dateSTR) & (self.dt64 <=dateEND)
        
        for attr in self.__dict__.keys():
            
            a = getattr(self,attr)
            
            if isinstance(a,np.ndarray):
                if a.size == nt:
                   setattr(self,attr,a[xmask])
                   
                elif len(a.shape) == 3:
                    if a.shape[2] == nt:
                        # then a 2D field + time
                        # for when 2D time series fields are already read
                        # and additional clipping is needed
                        setattr(self,attr,a[xmask][:][:])
            
            elif isinstance(a,list):    
                if len(a)==nt:
                   b = [a_k for a_k, xmask_k in zip(a,xmask) if xmask_k]
                   setattr(self,attr,b) 
                   
    def clipXY(self,xlims,ylims):
        import numpy as np
        
        ii = np.argwhere( (self.xGrid >= xlims[0]) & (self.xGrid <= xlims[1]) )
        jj = np.argwhere( (self.yGrid >= ylims[0]) & (self.yGrid <= ylims[1]) )        
        # - note: ii, jj are returned as vector columns, size (nx,1), (ny,1)
        # reshape them to (nx,), (ny,)
        ii=np.reshape(ii,(ii.size,))
        jj=np.reshape(jj,(jj.size,))

        self.xGrid = self.xGrid[ii]
        self.yGrid = self.yGrid[jj]

        self.i0 = ii[0]
        self.j0 = jj[0]
                 
#        if hasattr(self,'values'):
#            self.values = self.values[:][jj][ii]

    def dailyFields(self,varName):
        # read all fields from the file, assuming 1 snapshot per file
        import netCDF4 as n4
        import numpy as np
        
        nt = self.dt64.size
        nx = self.xGrid.size
        ny = self.yGrid.size
        
        ii=np.arange(self.i0,self.i0+nx)
        jj=np.arange(self.j0,self.j0+ny)
        
        values = np.zeros((nt,ny,nx))
        
        for k in range(nt):
            print(self.fileList[k])
            nc = n4.Dataset(self.fileList[k])
            if len(nc.variables[varName].shape) == 2:
                # if no time dimension included, like in VIIRS
                values[k,:,:] = nc.variables[varName][jj,ii]
            elif len(nc.variables[varName].shape) == 3:
                # if time dimension included like in OSTIA or ROMS:
                values[k,:,:] = nc.variables[varName][:,jj,ii]
            nc.close()   

        return values

    def dailyField(self,varName,it):
        # read one field from the file, assuming 1 snapshot per file
        import netCDF4 as n4
        import numpy as np
        
        nx = self.xGrid.size
        ny = self.yGrid.size
        
        ii=np.arange(self.i0,self.i0+nx)
        jj=np.arange(self.j0,self.j0+ny)
        
        print(self.fileList[it])
        
        nc = n4.Dataset(self.fileList[it])
        if len(nc.variables[varName].shape) == 2:
            # if no time dimension included, like in VIIRS
            value = nc.variables[varName][jj,ii]
        elif len(nc.variables[varName].shape) == 3:
            # if time dimension included like in OSTIA or ROMS:
            value = nc.variables[varName][:,jj,ii]
        nc.close()   
        
        return value
        
    def sampleDaily(self,x,y,dt64):
        # Sample a 2D field daily time series
        from scipy.interpolate import RegularGridInterpolator
        import numpy as np
        
        tGrid=self.dt64.astype('float') # gridded field times
        t=dt64.astype('float')          # sample times, corresponding to x, y
        
        if self.dt64.dtype == '<M8[us]':
            tGrid = tGrid*1.e-6
        if dt64.dtype == '<M8[us]':
            t=t*1.e-6
        
        nt=dt64.size
        pts=np.hstack([t.reshape(nt,1),y.reshape(nt,1),x.reshape(nt,1)])
        
        vfun = RegularGridInterpolator((tGrid,self.yGrid,self.xGrid), \
                                       self.values,                   \
                                       method='nearest',bounds_error=False)
        
        v = vfun(pts)    
        
        return v
    
    def years(self):
        self.years = self.dt64.astype('datetime64[Y]').astype(int) + 1970
    
    def months(self):
        self.months = self.dt64.astype('datetime64[M]').astype(int) % 12 + 1
   
    def daysSinceRefDate(self,refDate):

        self.refDate = refDate
        self.tdays=(self.dt64-self.refDate).astype('double')/(24*3600)

    def timeAve(self,varName,dSTR,dEND):
        # compute the time average field varName between dSTR,dEND (both are dt64)

        import numpy as np
        import copy
        
        
        cut = copy.copy(self)
        cut.clipTime(dSTR, dEND)
        values = cut.dailyFields(varName)
        
        ave = np.nanmean(values,0)
        
        return ave

class L3_VIIRS(timeSer_2D):
    
    def __init__(self,satID,waveLength):

        import netCDF4 as n4        

        # notes:
        # satID: string like "SNPP" or "JPSS1"
        # waveLength: integer (e.g. 410)
        # xGridName, yGridName: coordinate variable names xGrid, yGrid
        #                       that can be utilized for 2D interpolation
        
        import glob
        import re
        import numpy as np
        
        ######
        # Extract fileList:
        
        channelNo=str(waveLength)
        
        pdir='../../Dat/OceanColor/' + satID + '_VIIRS/WC/*'
        fname=satID+'_VIIRS.'+'*'+'.L3m.DAY.RRS.Rrs_'+channelNo+'*.nc'
        p= pdir + '/' + fname
            
        self.fileList=glob.glob(p)
        
        ######
        # time as datetime64
        # no time variable in the daily files, extract the date from 
        # the file names
        
        YMD=list()

        for fname in self.fileList:
            ymd=re.findall('\d{8}',fname)[0] # in each file name, find 8 digits
            YMD.append(ymd[0:4] + '-' + ymd[4:6] + '-' + ymd[6:8])
            
        self.dt64 = np.array(YMD,dtype='datetime64') + \
                    + np.timedelta64(12,'h') + \
                    + np.timedelta64(0,'m') + np.timedelta64(0,'s')
        ######
        # xGrid, yGrid (as 1D variables):
        
        xGridName='lon'
        yGridName='lat'
            
        nc = n4.Dataset(self.fileList[0])
        self.xGrid = nc.variables[xGridName][:]
        self.yGrid = nc.variables[yGridName][:]
        nc.close()
        
        self.xGridName=xGridName
        self.yGridName=yGridName
        
        #- indices of the first element of the 2d fields in each direction
        # (the length will always be consistent with size of xGrid, yGrid)
        self.i0 = 0
        self.j0 = 0


class OSTIA(timeSer_2D):
    
    def __init__(self,region):

        import netCDF4 as n4        
        import glob
        import re
        import numpy as np
        
        ######
        # Extract fileList:
        
        pdir='../../Dat/OSTIA/' + region + '/*'
        fname='*.ostia.sst.' + region + '.nc'
        p= pdir + '/' + fname
            
        self.fileList=glob.glob(p)
        
        ######
        # time as datetime64
        # extract the date from the file names
        
        YMD=list()

        for fname in self.fileList:
            ymd=re.findall('\d{8}',fname)[0] # in each file name, find 8 digits
            YMD.append(ymd[0:4] + '-' + ymd[4:6] + '-' + ymd[6:8])
            
        self.dt64 = np.array(YMD,dtype='datetime64') + \
                    + np.timedelta64(12,'h') + \
                    + np.timedelta64(0,'m') + np.timedelta64(0,'s')
         
        ######
        # xGrid, yGrid (as 1D variables):
            
        xGridName='lon'
        yGridName='lat'
            
        nc = n4.Dataset(self.fileList[0])
        self.xGrid = nc.variables[xGridName][:]
        self.yGrid = nc.variables[yGridName][:]
        nc.close()
        
        self.xGridName=xGridName
        self.yGridName=yGridName
        
        #- indices of the first element of the 2d fields in each direction
        # (the length will always be consistent with size of xGrid, yGrid)
        self.i0 = 0
        self.j0 = 0

    # def readValues(self,varName):
    #     import netCDF4 as n4
    #     import numpy as np
        
    #     nt = self.dt64.size
    #     nx = self.xGrid.size
    #     ny = self.yGrid.size
        
    #     ii=np.arange(self.i0,self.i0+nx)
    #     jj=np.arange(self.j0,self.j0+ny)
        
    #     self.values = np.zeros((nt,ny,nx))
        
    #     print(ii)
    #     print(jj)
        
    #     for k in range(nt):
    #         print(self.fileList[k])
    #         nc = n4.Dataset(self.fileList[k])
    #         self.values[k,:,:] = nc.variables[varName][:,jj,ii]
    #         nc.close()   
            
class HF_WC_DAILY(timeSer_2D):
    
    def __init__(self):

        import netCDF4 as n4        
        import glob
        import re
        import numpy as np
        
        ######
        # Extract fileList:
        # Note: two different directories
        # 2008-2011: daily obtained from ASCII files (origin?)
        # 2012-2023: daily computed from hourly files loaded from an 
        # opendap server (see load_hf_hourly_from_opendap.m, 
        # hourly_to_daily)
        
        fname='hf_daily_rtv_*.nc'
        
        # 2008-2011:
        pdir='../../Dat/HF/WC_RTV_DAILY/NCfromASCII_DAILY/*'
        p= pdir + '/' + fname       
        flist1=glob.glob(p)
        
        # 2021-present:
        pdir='../../Dat/HF/WC_RTV_DAILY/DOP05_NT12_ECC01/*'
        p= pdir + '/' + fname       
        flist2=glob.glob(p)
        
        self.fileList=flist1 + flist2
        
        ######
        # time as datetime64
        # extract the date from the file names
        
        YMD=list()

        for fname in self.fileList:
            ymd=re.findall('\d{8}',fname)[0] # in each file name, find 8 digits
            YMD.append(ymd[0:4] + '-' + ymd[4:6] + '-' + ymd[6:8])
            
        self.dt64 = np.array(YMD,dtype='datetime64') + \
                    + np.timedelta64(12,'h') + \
                    + np.timedelta64(0,'m') + np.timedelta64(0,'s')

        ######
        # xGrid, yGrid (as 1D variables):
            
        xGridName='lon'
        yGridName='lat'
            
        nc = n4.Dataset(self.fileList[0])
        self.xGrid = nc.variables[xGridName][:]
        self.yGrid = nc.variables[yGridName][:]
        nc.close()
        
        self.xGridName=xGridName
        self.yGridName=yGridName
        
        #- indices of the first element of the 2d fields in each direction
        # (the length will always be consistent with size of xGrid, yGrid)
        self.i0 = 0
        self.j0 = 0             
        
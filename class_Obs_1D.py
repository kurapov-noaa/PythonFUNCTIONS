# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 14:16:05 2023

@author: Alexander.Kurapov
"""

class Obs_1D:
    
    def clip(self,dateSTR,dateEND):
        import numpy as np
        
        nt=self.time.size
        
        xmask = (self.dt64 >= dateSTR) & (self.dt64 <=dateEND)
        
        for attr in self.__dict__.keys():
            a = getattr(self,attr)
            if isinstance(a,np.ndarray):
               if a.size == nt:
                   setattr(self,attr,a[xmask])
            
class CUGN(Obs_1D):
        
    def __init__(self,lineID,vlev=0):
        # note: lineID is a string variable like '90'
        
        import netCDF4 as n4
        import numpy as np
        import datetime as dt
        import dateutil.parser as dparser
        
        gliderDir='C:/Users/Alexander.Kurapov/Documents/Workspace/Dat/Glider_Rudnick'
        fname=gliderDir + '/' + 'CUGN_line_' + lineID + '.nc'
        
        nc  = n4.Dataset(fname)
        
        self.time = nc.variables['time'][:]
        self.timeUnits = nc.variables['time'].units
        
        self.lon = nc.variables['lon'][:]
        self.lat = nc.variables['lat'][:]
        self.temp = nc.variables['temperature'][vlev]
        self.salt = nc.variables['salinity'][vlev]
        
        nc.close()

        xmask=self.time.mask & self.lon.mask & self.lat.mask & \
              self.temp.mask * self.salt.mask

        self.time = self.time[~xmask]
        self.lon  = self.lon[~xmask]
        self.lat  = self.lat[~xmask]
        self.temp = self.temp[~xmask]
        self.salt = self.salt[~xmask]
        
        self.timeRefDate = dparser.parse(self.timeUnits,fuzzy=True)
        
        # make the time naive list to immediately convert it to a datetime64 
        # array where UTC is assumed
        dtimeList = [self.timeRefDate.replace(tzinfo=None) + \
                      dt.timedelta(seconds=self.time[k]) \
                      for k in range(len(self.time))]
            
        # - convert datetime list to np.datetime64
        self.dt64=np.array(dtimeList,dtype='datetime64')

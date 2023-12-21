# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 13:56:52 2023

@author: Alexander.Kurapov
"""



class timeSer():
    
    # CONTENTS:
    # - uniquely defined for each dataset:
    #    self.fileList: a set of files to read from 
    #    self.dt64: time variable, datetime64
    #    self.
    
    def daysSinceRefDate(self,refDate):
        from timeFunctionsAK import daysSinceRefDate as daysFunction
        
        self.refDate = refDate
        self.tdays=daysFunction(self.dt64,refDate)
    
class ndbc(timeSer):
    
    def __init__(self,buoyID):

        #import netCDF4 as n4        
        import glob
        import pandas as pd
        import numpy as np

        sID=str(buoyID)
        
        
        # collect file names for yearly collections:
        pdir='../../Dat/NDBC_buoys/' + sID
        fname= sID + 'h' + '*.txt'
        p = pdir + '/' + fname         
        self.fileList=glob.glob(p)
        
        # add file names for the last year (currently 2023), month by month:
        yr='2023'
        p='../../Dat/NDBC_buoys/' + sID + '/' + yr + '/' +sID + '_*' + yr + '.txt'
        self.fileList.extend(glob.glob(p))    
        
        self.dt64 = []
        self.var= {}
        self.var['Uwind'] = []
        self.var['Vwind'] = []
        self.var['WSPD'] = []
        self.var['angle'] = []
        
        for fname in self.fileList:
            print(fname)
            a = pd.read_csv(fname,sep='\s+',header=0,skiprows=[1])
            yr=a['#YY']
            mo=a['MM']
            dd=a['DD']
            hh=a['hh']
            mm=a['mm']

            # alph: wind direction wrt true east (vector in the dir of the wind)            
            alph=np.pi*(270-np.array(a['WDIR']))/180
            WSPD=np.array(a['WSPD'])
            WSPD[np.argwhere(WSPD==99.0)]=np.nan
            Uwind=WSPD*np.cos(alph)
            Vwind=WSPD*np.sin(alph)

            dt64 = [np.datetime64(str(yr[k]) + '-' \
                               + str(mo[k]).zfill(2) + '-' \
                               + str(dd[k]).zfill(2) + 'T' \
                               + str(hh[k]).zfill(2) + ':' \
                               + str(mm[k]).zfill(2) + ':00') for k in range(yr.size)]   
              
            self.dt64.extend(dt64) 
            self.var['Uwind'].extend(Uwind)
            self.var['Vwind'].extend(Vwind)
            self.var['WSPD'].extend(WSPD)
            self.var['angle'].extend(alph)
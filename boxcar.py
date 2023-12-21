# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 11:14:14 2023

@author: Alexander.Kurapov
"""

def boxcar(t,f,tout,L):
    # Given a time series f(t) (possibly with gaps and of variable time resolution)
    # apply the box car filter averaging f values in the neighborhood (-L/2, L/2) 
    # of each point t[k]
 
    import numpy as np

    halfL=0.5*L
    ntout = tout.size

    fout=np.nan*np.empty_like(tout)
    for it in range(ntout):
    
        ii = np.argwhere( (t>=tout[it]-halfL) & (t<=tout[it]+halfL) )
        if len(ii) > 0:
            if np.any(~np.isnan(f[ii])):
                fout[it] = np.nanmean(f[ii])
    
    return fout


                
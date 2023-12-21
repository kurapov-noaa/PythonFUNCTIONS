# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 09:41:21 2023

@author: Alexander.Kurapov
"""

def interp2(Y,X,F,y,x):
    from scipy.interpolate import RegularGridInterpolator
    
    # X, Y: meshgrid like, the same shape as F
    
    hfun = RegularGridInterpolator( (Y[:,0],X[0,:]), F, method='linear', bounds_error=False)
    f = hfun((y,x))
    
    return f
        
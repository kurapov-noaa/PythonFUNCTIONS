# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 13:41:03 2023

@author: Alexander.Kurapov
"""

def pseudoHeat(h0,mask,Kh,nsteps):
    
    import copy
    import numpy as np
    
    # Smooth a 2D field iterating the solver for the pseudoheat eqn.
    # h0: field to smooth
    # mask: mask (land)
    # Kh: diffusion coeff
    # nsteps: number of steps
    #
    # Boundary conditions: clamping along the open (not masked) edges of
    # the rectangular domain, no gradient along the mask
    
    ny = h0.shape[0]
    nx = h0.shape[1]
    
    # all points:
    ii=np.arange(nx).reshape((1,nx)) 
    jj=np.arange(ny).reshape((ny,1))
    
    # internal points 
    ii_int = ii[:,1:nx-1]
    jj_int = jj[1:ny-1,:]
    
    Mx = np.minimum(mask[1:ny-1,1:nx],mask[1:ny-1,0:nx-1])   
    My = np.minimum(mask[1:ny,1:nx-1],mask[0:ny-1,1:nx-1])
    
    h = copy.copy(h0)
    
    for it in range(nsteps):
        
        # compute gradients:
        Dx = np.diff(h[1:ny-1,:], axis = 1) 
        Dy = np.diff(h[:,1:nx-1], axis = 0)
        
        # apply no-gradient conditions along the masked points
        # "as natural conditions":
        Dx[ (Mx==0) ] = 0    
        Dy[ (My==0) ] = 0    
        
        h[jj_int,ii_int] += Kh * ( np.diff(Dx,axis=1) + np.diff(Dy,axis=0) )
        
    return h


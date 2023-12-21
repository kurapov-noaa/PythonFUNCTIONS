# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 15:46:20 2023

@author: Alexander.Kurapov
"""

def harm_const_multiFile(timeSer,varName,Tp):
    
    import numpy as np
    
    # timeSer: class_2D object for now
    # fields will be read on a fly, one by one, assuming 1 field per file
    
    nx = timeSer.xGrid.size
    ny = timeSer.yGrid.size
    N = nx*ny                # 1-d representation of the spatial arrays
    nt = timeSer.tdays.size
    nper = Tp.size
    
    amp=0 #np.zeros((nper,ny,nx))
    phase=0 #np.zeros((nper,ny,nx))
    
    omega=2*np.pi/Tp
    omega_x_t = np.matmul(timeSer.tdays.reshape(nt,1),omega.reshape(1,nper))
    
    co = np.cos(omega_x_t)
    si = np.sin(omega_x_t)
    
    A = np.concatenate((co,-si),axis=1)
    ApA=np.matmul(A.transpose(),A)
   
    B=np.linalg.solve(ApA,A.transpose())

    c = np.zeros((2*nper,N))
    a0 = np.zeros((1,N))
    ont = 1./nt
    for it in range(nt):
        print(timeSer.fileList[it])
        a = timeSer.dailyField(varName,it).reshape(1,N)
        a0 += ont*a
        c += np.matmul(B[:,it].reshape(2*nper,1),a)
     
    # Correct c to represent the demeaned time series:
    b=np.sum(B,axis=1).reshape(2*nper,1) 
    c -= np.matmul(b,a0)

    zR = c[range(nper),:]
    zI = c[range(nper,2*nper),:]

    amp = np.sqrt( zR*zR + zI * zI ).reshape(nper,ny,nx)
    phase = np.arctan2(zI,zR).reshape(nper,ny,nx) * 180./np.pi

    a0=a0.reshape(ny,nx)
    
    return a0,amp,phase

def harm_const_array(t,f,Tp):
    
    import numpy as np
  
    nt = f.shape[0]
    ny = f.shape[1]
    nx = f.shape[2]
    N = nx * ny               # 1-d representation of the spatial arrays
    nper = Tp.size
    
    omega=2*np.pi/Tp
    omega_x_t = np.matmul(t.reshape(nt,1),omega.reshape(1,nper))
    
    co = np.cos(omega_x_t)
    si = np.sin(omega_x_t)
    
    A = np.concatenate((co,-si),axis=1)
    ApA=np.matmul(A.transpose(),A)
   
    B=np.linalg.solve(ApA,A.transpose())

    c = np.zeros((2*nper,N))
    a0 = np.zeros((1,N))
    ont = 1./nt
    for it in range(nt):
        a = f[it,:,:].reshape(1,N)
        a0 += ont*a
        c += np.matmul(B[:,it].reshape(2*nper,1),a)
     
    # Correct c to represent the demeaned time series:
    b=np.sum(B,axis=1).reshape(2*nper,1) 
    c -= np.matmul(b,a0)

    zR = c[range(nper),:]
    zI = c[range(nper,2*nper),:]

    amp = np.sqrt( zR*zR + zI * zI ).reshape(nper,ny,nx)
    phase = np.arctan2(zI,zR).reshape(nper,ny,nx) * 180./np.pi

    a0=a0.reshape(ny,nx)

    return a0,amp,phase

def write_harm_const(a0,amp,phase,varUnits,Tp,timeUnits,refDate,dateSTR,dateEND,outfile):
    import netCDF4 as n4
    import numpy as np
 
    nx = a0.shape[1]
    ny = a0.shape[0]
    ncons = Tp.size
    
    ymdSTR=dateSTR.item().strftime('%Y%m%d') # item() is used to go from dt64 to dt
    ymdEND=dateEND.item().strftime('%Y%m%d')
    
    ymdREF=refDate.item().strftime('%Y%m%d')
    
    nc = n4.Dataset(outfile,mode='w')

    nc.createDimension('nx', nx)
    nc.createDimension('ny', ny) 
    nc.createDimension('ncons', ncons)
    
    v=nc.createVariable('dateSTR',np.int32)
    v.long_name = 'the harmonic analysis interval start'
    v[:] = int(ymdSTR)

    v=nc.createVariable('dateEND',np.int32)
    v.long_name = 'the harmonic analysis interval ending'
    v[:] = int(ymdEND)
     
    v=nc.createVariable('refDate',np.int32)
    v.long_name = 'reference date'
    v[:] = int(ymdREF)
      
    v=nc.createVariable('Tp',np.float64,('ncons'))
    v.long_name = 'periods'
    v.units = timeUnits
    v[:] = Tp
    
    v=nc.createVariable('time_mean',np.float64,('ny','nx'))
    v.long_name = 'time mean'
    v.units = varUnits
    v[:] = a0
    
    v=nc.createVariable('amp',np.float64,('ncons','ny','nx'))
    v.long_name = 'amplitude'
    v.units = varUnits
    v[:] = amp
    
    v=nc.createVariable('phase',np.float64,('ncons','ny','nx'))
    v.long_name = 'phase, degrees'
    v.units = 'degrees'
    v[:] = phase
    
    nc.close()    
    
def read_harm_const(annualCycleFile):
    import netCDF4 as n4
    import numpy as np
    import datetime as dt
    
    nc = n4.Dataset(annualCycleFile,mode='r')
    
    a0 = nc.variables['time_mean'][:]
    amp = nc.variables['amp'][:]
    phase = nc.variables['phase'][:]
    Tp = nc.variables['Tp'][:]
    
    ymd = nc.variables['refDate'][:]
    ymd = str(ymd)
    ymd = dt.datetime.strptime(ymd,'%Y%m%d')
    refDate=np.datetime64(ymd,'s')
    
    return a0,amp,phase,Tp,refDate
        
    nc.close()      
    
def anom(timeSer,varName,annualCycleFile):
    
    a0,amp,phase,Tp,refDate = read_harm_const(annualCycleFile) 
    timeSer.daysSinceRefDate(refDate) # => add refDate, tdays to timeSer
    F = timeSer.dailyFields(varName)

    t = timeSer.tdays
    
    F_ann = ann(a0,amp,phase,Tp,t)
    F_anom = F - F_ann 
    
    return F_anom

def ann_forLoop(a0,amp,phase,Tp,t):
    import numpy as np
    
    omega = 2*np.pi / Tp
    phase = phase*np.pi/180    

    nx = a0.shape[1]
    ny = a0.shape[0]
    nt = len(t)
    nper = len(Tp)

    a = np.zeros((nt,ny,nx))
    for it in range(nt):
    
        # annual:
        a[it,:,:]=a0
        for k in range(nper):
            a[it,:,:] += amp[k,:,:]*np.cos(omega[k]*t[it] + phase[k,:,:])
    
    return a
    
def ann(a0,amp,phase,Tp,t):
    import numpy as np
    
    omega = 2*np.pi / Tp
    
    phase = phase*np.pi/180    

    nx = a0.shape[1]
    ny = a0.shape[0]
    nt = len(t)
    nper = len(Tp)

    N = nx * ny
    a0=a0.reshape(1,N)
    amp=amp.reshape(nper,N)
    phase=phase.reshape(nper,N)
    t=t.reshape(nt,1)
    omega = omega.reshape(1,3)
    
    
    a = np.zeros((nt,N))
    a+= a0 # adds a0 as a row vector to each row of a matrix
    
    aR = amp * np.cos(phase)
    aI = amp * np.sin(phase)
    
    tom = np.matmul(t,omega)
    cos_tom = np.cos(tom)
    sin_tom = np.sin(tom)
    
    a += np.dot(cos_tom,aR) - np.dot(sin_tom,aI)

    a = a.reshape(nt,ny,nx)

    return a
    
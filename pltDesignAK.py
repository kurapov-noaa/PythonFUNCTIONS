# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 12:52:41 2021

@author: Alexander.Kurapov
"""

def add_year_lines(ax,refDate,clr,lnStyle):
    
    import datetime as dt
    import matplotlib.pyplot as plt
    
    xlims=ax.get_xlim()
    ylims=ax.get_ylim()
    dSTR=refDate+dt.timedelta(xlims[0])
    dEND=refDate+dt.timedelta(xlims[1])
    for yr in range(dSTR.year+1,dEND.year+1):
        print(yr)
        t0=(dt.datetime(yr,1,1,0,0,0)-refDate).total_seconds()/(24*3600)
        plt.plot((t0,t0),ylims,color=clr,linestyle=lnStyle)
        plt.text(t0,ylims[0]+0.9*(ylims[1]-ylims[0]),str(yr))
                
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    
def add_year_lines_64(ax,clr,lnStyle,labelPosition=0.9):
    
    # use when the x axis is already datetime64
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    xlims=ax.get_xlim()
    ylims=ax.get_ylim()
    
    dSTR=np.datetime64(int(xlims[0]),'D')
    dEND=np.datetime64(int(xlims[1]),'D')
    
    yearSTR=dSTR.astype(object).year
    yearEND=dEND.astype(object).year
    yrs=np.arange(yearSTR,yearEND+1)
    syrs=[ str(yr) for yr in yrs]
    jan1=[ np.datetime64(syr,'D') for syr in syrs]
    
    t2d,y2d = np.meshgrid(jan1,ylims)
    
    ax.plot(t2d,y2d,color=clr,linestyle=lnStyle)
    
def summerMonths(ax,color=(255/255,243/255,109/255)):
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    xlims=ax.get_xlim()
    ylims=ax.get_ylim()
    
    dSTR=np.datetime64(int(xlims[0]),'D')
    dEND=np.datetime64(int(xlims[1]),'D')
    
    yearSTR=dSTR.astype(object).year
    yearEND=dEND.astype(object).year
    yrs=np.arange(yearSTR,yearEND+1)
    syrs=[ str(yr) for yr in yrs]
    
    for syr in syrs:
        may1=np.datetime64(syr + '-05-01')
        sep1=np.datetime64(syr + '-09-01')
        
        xy=(may1,ylims[0])
        pp=plt.Rectangle(xy,sep1-may1,np.diff(ylims),color=color)
        ax.add_patch(pp)
        

    
    
    
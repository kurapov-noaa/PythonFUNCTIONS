# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 10:29:57 2023

@author: Alexander.Kurapov
"""

def daysSinceRefDate(dt64,refDate):
    import numpy as np
    
    tdays=(dt64 - refDate)/np.timedelta64(1,'s')/24/3600
    
    return tdays

def days2dt64(tdays,refDate):
    # given an array of fractional days from the refDate, compute
    # a datetime64 array

    import numpy as np
           
    seconds_array = (tdays * 86400).astype(int)
    dt64 = np.array([ (refDate + np.timedelta64(s,'s')) for s in seconds_array]) 
    
    return dt64

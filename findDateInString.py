# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 12:33:53 2020

@author: Alexander.Kurapov
"""

def findDateInString(a):
    # Finds a date in string "a" (such as units attribute)
    # returns datetime object
    import re
    import datetime as dt
    
    match = re.search('\d{4}-\d{2}-\d{2}', a)
    date1 = dt.datetime.strptime(match.group(), '%Y-%m-%d')
    
    return date1
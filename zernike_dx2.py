#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 12:50:22 2020

@author: esthersoriahernandez
"""

def zernike_dx(n,m,x,y):
    import math 
    import numpy as np 
    import div
    from scipy.special import comb

    if (m < 0):
        m_abs = -m
    else:
        m_abs = m
    
    
    if (m > 0):
        dzx = 0
        for s in range(int(((n-m_abs)/2)+1)):
            for j in range(int(((n-m_abs)/2-s)+1)):
                for k in range(int((m_abs)/2+1)):
                    dzx =  dzx + ((pow(-1,(s+k))*math.factorial(n-s))/ \
                                  (math.factorial(s)*math.factorial((n+m_abs)/2-s)*\
                                   math.factorial((n-m_abs)/2-s)))*comb((n-m_abs)/2-s,j)*\
                                  comb(m_abs,2*k)*(n-2*(s+j+k))\
                                  *x**(n-2*(s+j+k)-1)*y**(2*(j+k))
                    
        dzx = math.sqrt(2*(n+1))*dzx
        
    elif (m<0):
        dzx = 0
        for s in range(int(((n-m_abs)/2)+1)):
                for j in range(int(((n-m_abs)/2-s)+1)):
                    for k in range(int(((m_abs-1)/2)+1)):
                        dzx =  dzx + ((pow(-1,(s+k))*math.factorial(n-s))/ \
                              (math.factorial(s)*math.factorial((n+m_abs)/2-s)*math.factorial((n-m_abs)/2-s)))*\
                              comb((n-m_abs)/2-s,j)*comb(m_abs,2*k+1)*(n-2*(s+j+k)-1)*\
                              x**(n-2*(s+j+k)-2)*y**(2*(j+k)+1)
                        
        dzx = math.sqrt(2*(n+1))*dzx
        
    elif(m == 0):        
        dzx = 0  
        for s in range(int(n/2+1)):
            for j in range(int(n/2-s+1)) :
                dzx = dzx +\
                    (((pow(-1,s))*math.factorial(n-s))\
                    / (math.factorial(s)*math.factorial(n/2-s)*math.factorial(n/2-s)))\
                    * comb(n/2-s,j)*(n-2*(s+j))\
                    * x**(n-2*(s+j)-1)*y**(2*j)
 
        dzx = math.sqrt(n+1)*dzx
        
        
    return dzx

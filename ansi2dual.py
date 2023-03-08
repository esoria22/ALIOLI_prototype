#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 11:32:39 2020

@author: esthersoriahernandez

This function changes the notation from  ANSI to dual.

Inputs:
	*j:ANSI index  

Outputs:
	*n: radial order
	*m: angular frenquency
"""

def ansi2dual(j):
    import math
    n=math.ceil((-3+math.sqrt(9+8*j))/2)
    m=2*j-n*(n+2)
    return n,m
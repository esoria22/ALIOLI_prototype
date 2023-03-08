#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 12:35:21 2020

@author: esthersoriahernandez

This function calculate the center of gravity of the input image in different regions
Inputs:
	*Image of interest
	*Thershold factor(between 0 and 1)
	*regx: vector with the starting points x coordinate
	*regy: vector with the starting points y coordinate
	*regh: vector with lenght of each region x coordinate
	*regw: vector with lenght of each region y coordinate
Outputs:
	*vector with COG coordinate x in each region
	*vector with COG coordinate y in each region

"""

def cog_multiple(input, threshold_factor, regx, regy, regh, regw):
    import numpy as np
    import cog
    num_elements = np.size(regx)
    
    cx = np.zeros(num_elements)
    cy = np.zeros(num_elements)
    
    for i in range (num_elements):
        
        curx = int(regx[i])
        cury = int(regy[i])
        curh = int(regh[i])
        curw = int(regw[i])
        cur_reg = input[cury:cury+curh-1, curx:curx+curw-1]
        
        threshold = np.max(cur_reg) * threshold_factor
        
        cx_temp,cy_temp = cog.cog(cur_reg, threshold)
    
        
        cx[i] = cx_temp
        cy[i] = cy_temp
    return cx,cy
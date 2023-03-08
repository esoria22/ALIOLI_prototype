# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 09:19:55 2020

@author: esoria
"""

def g_wfs_precalc(modes, angles, width, height):
    import numpy as np
    import sympy as sp
    import math
    import zernike_dx as zdx
    import zernike_dy as zdy
    import ansi2dual as a2d
    import radon as rd
    import pandas as pd
#    width = np.snmodes = np.size(modes)
    
    # Get number of angles from the input vector of angles.
    nangles = len(angles);
    
    # Prepare the coordinate system
    x = np.linspace(-1, 1, height);
    y = np.linspace(-1, 1, width);
    X, Y = np.meshgrid(x, y);
    r=(X**2+Y**2)**(0.5)
    theta=np.arctan2(Y,X)
    
    # Define the region where the Zernike derivates will be calculated
    idx = r<=1;
    
    # Obtain the vector of valid Radon coordinates
    i1_dummy = np.zeros((width, height))
    rtx,coords = rd.radon(i1_dummy, angles, circle= True)
    
    ncoords = len(coords)
    nmodes = len(modes)
    # The functions "zernike_dx" and "zernike_dy" will be used to calculate
    # the Zernike derivatives all over the input grid. The value "R" sets the
    # radius in pixels where the calculated derivates are valid.
    R = np.floor(min(width, height)/2)
    
    # This avoids a negative sqrt later on. The values outside the radius have
    # to be ignored anyway.
    coords[abs(coords) > R] = R
    
    # Allocate memory for the Zernike derivatives
    z_dx = np.zeros((width, height))
    z_dy = np.zeros((width, height))
    
    # Allocate memory for the output matrix
    H = np.zeros((ncoords, nangles, nmodes))
#    height = np.size(i1,0)
    # This function uses the "zernike" library for calculating the
    # Zernike derivatives.
   
    
    # Get number of modes from the input vector of modes.
    
    
#   Loop through all the specified modes
    for mode_idx in np.arange(0,nmodes):
        
        # Select current mode
        mode = modes[mode_idx]
              
        n, m = a2d.ansi2dual(mode+2)
        z_dx[idx] = zdx.zernike_dx(n, m, X[idx], Y[idx])
        z_dy[idx] = zdy.zernike_dy(n, m, X[idx], Y[idx])
        where_are_NaNs = np.isnan(z_dy)
        z_dy[where_are_NaNs] = 0
        where_are_NaNs = np.isnan(z_dx)
        z_dx[where_are_NaNs] = 0
        
        # Loop through all the specified angles
        for angle_idx in range(nangles):
        
            # Select current angle
            angle = [angles[angle_idx]]
            angled=  np.deg2rad(angle) 
            """
            % Calculate the derivate of the wavefront in the direction
            % perpendicular to the current angle.
            % In this case, the angle has been negated so that it is coherent
            % with the way the "radon" function understand angles. It is easy
            % to check whether the result is correct by inspecting the H matrix
            % elements that correspond to the "defocus" mode. Being the shape
            % of this mode circular, the resulting mean slopes are necessarely
            % equal for all angles.
            """
            z_d = z_dx * np.cos(-angled) + z_dy * np.sin(-angled)
            """
            % This constant enables obtaining the mean slope from the integration
            % that the Radon transform performs."""
            L = 2*np.sqrt(R**2 - coords**2)
           
    	
            # Divisions by 0 need to be avoided. Values outside the radius "R"
            # need to be ignored anyway.
            L[L == 0] = float('inf')
            #hasta aqui bien
        
    	
            # Calculate the mean wavefront slope along perpendicular cuts for
            # the current angle and Zernike mode.
            r = rd.radon(z_d, angle,True)[0]
            H[:, angle_idx, mode_idx] = (r[:,0])/ L
    return H
       

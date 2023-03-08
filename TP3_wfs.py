# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 14:33:38 2021

@author: esoria
"""

def g_wfs(i1, i2, nangles, distance,H_pinv):
    import numpy as np
    import math 
    from scipy.interpolate import InterpolatedUnivariateSpline
    import norm_accum as na
    import findcrossing as fc
    from skimage.transform import radon, rescale
    import find
    import radon as rd
    
    
    """
    % This function executes the wavefront recontruction algorithm presented in:
    % M. van Dam and R. Lane, "Wave-front sensing from defocused images by use of
    % wave-front slopes," Appl. Opt.  41, 5497-5502 (2002).
    %
    %  i1:		Image of the first defocused pupil.
    %  i2:		Image of the second defocused pupil.
    %  angles:	Vector of the angles where the Radon transform will be applied.
    %  z:		Distance between the pupil plane and the two defocused planes.
    %  H_pinv:	Inverse matrix that will be used for the final LMS fit.
    
    % Calculate Radon transform of each input image for all the specified angles.
    % MATLAB supports performing the Radon transform with "gpuArray" data types.
    % It may be positive to discard the pixels outside the area of the pupil, just
    % to discard the noise.
    """
    angles_full = np.linspace(0, 180, nangles + 1)
    angles = angles_full[0:- 1]# 0 degrees is equal to 180 degrees
    '''
    j = int(np.sqrt(np.size(i1)))
    i1=np.reshape(i1,(j,j),order='C')
    i2=np.reshape(i2,(j,j),order='C')
    '''
   
    P1,xp = rd.radon(i1, angles,circle= True)
    P2,xp = rd.radon(i2, angles, circle = True)
   
    
    '''b = np.ceil (math.sqrt (sum (np.array([np.size(i1,0),np.size(i1,1)])**2))/2 + 1)
    xp = np.arange(-b,b+1)
    xp = np.transpose(xp)'''
    ncoords = len(xp)
    """
    % Get a vector of equidistant values in the range (0, 1). These are the values
    % of the CDFs (Cummulative Distribution Functions) whose coordinates will be
    % obtained. The CDF values 0 and 1 always give a zero slope, so they can be
    % safely ignored.
    % Van Dam's paper does not specify the adequate number of values for an accurate
    % fit. This needs further research.
    """
    ncrossings = np.size(i1,0)
    s_vector_full = np.linspace(0.0, 1.0, ncrossings + 2)
    s_vector = s_vector_full[1: ]
    
    #Get number of specified angles
    nangles = np.size(angles)
    z = distance
    #Allocate memory for the estimated mean slopes at perpendicular directions of
    # each radon angle.
    slopes = np.zeros((ncoords, nangles))
    for i in range(nangles):
        	#Get the Radon transform at the current angle
            P1_i = P1[:, i]
            P2_i = P2[:, i]
            
    	
    	#Accumulate and normalize the Radon vectors to get each CDF.Ec24
            c1 = na.norm_accum(P1_i)
            c2 = na.norm_accum(P2_i)
            c1=c1[1:]
            c2=c2[1:]
    	#Find the coordinates where each CDF has the values defined in "s_vector".
    	#Given that the CDFs are discrete, an interpolation needs to be performed.
            u1 = fc.matching(c1, xp, s_vector)
        #u1=u1[1:72]
            u2 = fc.matching(c2, xp, s_vector)
            #u2=u2[1:72]
        	# Estimate the slopes of the wavefront at the pupil plane (z = 0).
            slopes_nointerp= (u2 - u1) / (2*z)
    #        
        	# Get the coordinates where the slope estimates are valid.
            slopes_pos = (u1 + u2) / 2
            s= InterpolatedUnivariateSpline(slopes_pos,slopes_nointerp,k=3)
            slopes[:,i]= s(xp)
            #slopes[:, i] = interp(xp,slopes_pos, slopes_nointerp,kind = 'linear' )
        	# Interpolate the slopes so that they fall in the same discrete coordinates olate.interp1d
        	# of the original Radon vector. This is necessary because the mean slopes of
        	# each Zernike mode have been calculated along those positions.
            # kind = 'linear'  fill_value="extrapolate"
    if H_pinv is None:
        slope = slopes.ravel(order='F')
        return slopes
    else:
        slope = slopes.ravel(order='F')
        W = H_pinv.dot(slope)
        return slopes, W

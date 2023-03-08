 # -*- coding: utf-8 -*-
"""
Created on Thu May  6 12:33:03 2021

@author: esoria


This function simulates the behavior of a SHwfs

Inputs:
	*pupil diameter (in meters)
    *wavelenght_wfs
The parameter lenslet diameter, array_size, and focal lenght are fixed
Outputs:
	*shwfs: my wfs as a surface
	*Tsh_aperture: the aperture which I am working 
"""

from hcipy import *
def SH2(pupil_diameter,wavelength_wfs):
    '''
    import numpy as np
    #fixed parameters (Thorlabs MLA150C)
    
    f_number = 48
    num_lenslets = 12# 40 lenslets along one diameter
  

    magnification = pupil_diameter / 1.52
    magnifier = Magnifier(magnification)
  

    num_lenslets = 49
    focal_length = 0.0146
    lenslet_diameter = 300*10**-6
    tam_array=3.8*10**-3
    
    #select the area on studing, defining the proyection of my pupil
    pup_array_diameter=pupil_diameter
    #generate the pupil array
    pup_array_grid=make_pupil_grid(256,pup_array_diameter)
    swfs = SquareShackHartmannWavefrontSensorOptics(pup_array_grid, f_number, \
                                                 num_lenslets, pupil_diameter)
    #create the aperture
    aperture_sh=circular_aperture(pup_array_diameter)
    #Generate the field
    Tsh_aperture = evaluate_supersampled(aperture_sh, pup_array_grid, 1)

    #Generate microlenses grid
    #first I create a vector with the position of each microlens
    x = np.arange(-pup_array_diameter, pup_array_diameter, lenslet_diameter)
    mla_grid = CartesianGrid(SeparatedCoords((x, x)))
    
    #Definition of the microlenses Array
    micro_lens_array = MicroLensArray(pup_array_grid,mla_grid, focal_length)
    #Defining the wfs, in function of my pupil poryecter over the microleses array
    swfs=ShackHartmannWavefrontSensorOptics(pup_array_grid,micro_lens_array)
    return swfs, Tsh_aperture
    '''
    import numpy as np
    #fixed parameters (Thorlabs MLA150C)
    num_lenslets = 49
    focal_length = 0.0146
    lenslet_diameter = 300*10**-6
    tam_array=9*10**-3
    #select the area on studing, defining the proyection of my pupil
    pup_array_diameter=pupil_diameter
    #generate the pupil array
    pup_array_grid=make_pupil_grid(256,pup_array_diameter)

    #create the aperture
    aperture_sh=circular_aperture(pup_array_diameter)
    #Generate the field
    Tsh_aperture = evaluate_supersampled(aperture_sh, pup_array_grid, 1)

    #Generate microlenses grid
    #first I create a vector with the position of each microlens
    x = np.arange(-pup_array_diameter, pup_array_diameter, lenslet_diameter)
    mla_grid = CartesianGrid(SeparatedCoords((x, x)))
    
    #Definition of the microlenses Array
    micro_lens_array = MicroLensArray(pup_array_grid,mla_grid, focal_length)
    #Defining the wfs, in function of my pupil poryecter over the microleses array
    swfs=ShackHartmannWavefrontSensorOptics(pup_array_grid,micro_lens_array)
    return swfs, Tsh_aperture
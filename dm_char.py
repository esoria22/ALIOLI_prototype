# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 11:07:10 2021

@author: esoria
"""
from hcipy import *
import numpy as np
def DM(N,diameter,pupil_grid,char):
    import numpy as np
    if  char == 1:
        num_modes=N
        dm_modes = make_zernike_basis(num_modes, diameter,pupil_grid, starting_mode=0, ansi=True)
        dm_modes = ModeBasis([mode / np.ptp(mode) for mode in dm_modes], pupil_grid)
#
    #defino el DM
        deformable_mirror = DeformableMirror(dm_modes)
        num_act=0
    else:
        num_actuators_across_pupil = int(np.sqrt(N))
        actuator_spacing = diameter / num_actuators_across_pupil
        influence_functions = make_gaussian_influence_functions(pupil_grid, num_actuators_across_pupil, actuator_spacing)
        deformable_mirror = DeformableMirror(influence_functions)
        num_act = deformable_mirror.num_actuators
        
    return deformable_mirror, num_modes, num_act
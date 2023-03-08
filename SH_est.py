# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 10:24:37 2021

@author: esoria
"""
from hcipy import *
def SH_estimator(detector,star_regx, star_regy, star_regh, star_regw, star_refx,star_refy):
    import numpy as  np
    import cog_multiple
    #normalizo la señal leida en la cámara
    image=detector/ detector.sum()
    j = int(np.sqrt(image.shape))
    imagen=np.reshape(image,(j,j),order='C')
    
    star_threshold_factor = 0.2
    cx, cy = cog_multiple.cog_multiple(imagen, star_threshold_factor, star_regx, star_regy, star_regh, star_regw)
    star_dx = star_regx + cx - star_refx
    star_dy = star_regy + cy - star_refy
    slopes = np.append(star_dx,star_dy)
    
    return slopes
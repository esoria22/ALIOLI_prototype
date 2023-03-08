#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 11:23:42 2020

@author: esthersoriahernandez
"""
#calculo de centroides imagen de referencia
def setref(ubicacion,tf,bias):
    import pyfits as ft
    import imagesc as imagesc
    import numpy as np
    import cog
    import math 
    import draw_regs
    import cog_multiple
    import matplotlib
    import hcipy as hp
    import matplotlib.pyplot as plt
    
    if ubicacion=='*.fits':
        ref=hp.read_fits(ubicacion)
    else:
        referencia=ubicacion
    #'/Users/esthersoriahernandez/Desktop/Funciones Python/ref_aoli.fits'
    # ref=ref[0].data
    # ref=z=np.reshape(ref,(300,300))
    pupil = referencia
    
    star_threshold_factor = tf
    
    naccum = 1
    bias_level = bias * naccum
    
    # #Subtract bias
    ref = referencia-bias_level
    pupil = pupil - bias_level
    ref[ref<0] = 0
    pupil[pupil<0]= 0
    
    #Get sensor dimension
    sensor_height = np.size(ref,0)
    sensor_width = np.size(ref,1)
    plt.figure()
    plt.imshow(ref,cmap='gray')
    print('Locate 3 spots positioned as in an "L"')
    p1,p2,p3= plt.ginput(3)
    plt.show()
    x_ref=np.array([p1[0],p2[0],p3[0]])
    y_ref=np.array([p1[1],p2[1],p3[1]])
    n_ref_v = float(input('Number of spots in the first line: '))
    n_ref_h = float(input('Number of spots in the second line: '))
    
    #Calculate step between each spot
    step_xh = (x_ref[2] - x_ref[1]) / (n_ref_h - 1)
    step_yh = (y_ref[2] - y_ref[1]) / (n_ref_h - 1)
    step_xv = (x_ref[1] - x_ref[0]) / (n_ref_v - 1)
    step_yv = (y_ref[1] - y_ref[0]) / (n_ref_v - 1)
    
    #Calculate distances between spots
    dist_h = np.sqrt(step_xh**2 + step_yh**2)
    dist_v = np.sqrt(step_xv**2 + step_yv**2)
    
    #Get the minimum distance between spots (for defining search regions)
    min_dist = min(dist_h, dist_v)
    
    # choice = questdlg('Would you like to auto-adjust the 3 marked spots?', ...
    # 	'Auto-adjust', ...
    # 	'Yes','No','Yes');
    
    # if (strcmp(choice, 'Yes') == true)
    #Perform auto-adjust of the three spots based on CoG
    for i in range(np.size(x_ref)):
        search_reg_x1 = np.floor(x_ref[i] - min_dist/2)
        search_reg_x2 = np.floor(x_ref[i] + min_dist/2)
        search_reg_y1 = np.floor(y_ref[i] - min_dist/2)
        search_reg_y2 = np.floor(y_ref[i] + min_dist/2)
        search_reg = ref[int(search_reg_y1-1):int(search_reg_y2-1), int(search_reg_x1-1):int(search_reg_x2-1)]
        cx,cy = cog.cog(search_reg, np.amax(search_reg) * star_threshold_factor)
        x_ref[i] = search_reg_x1 + cx
        y_ref[i] = search_reg_y1 + cy
    
    
    #  Re-calculate steps, needed after auto-adjustment
    step_xh = (x_ref[2] - x_ref[1]) / (n_ref_h - 1)
    step_yh = (y_ref[2] - y_ref[1]) / (n_ref_h - 1)
    step_xv = (x_ref[1] - x_ref[0]) / (n_ref_v - 1)
    step_yv = (y_ref[1] - y_ref[0]) / (n_ref_v - 1)
    
    #      Re-calculate distances between spots
    dist_h = np.sqrt(step_xh**2 + step_yh**2)
    dist_v = np.sqrt(step_xv**2 + step_yv**2)
    
    #  Get the minimum distance between spots (for defining search regions)
    min_dist = min(dist_h, dist_v)
    
    
    # Number of times that the area defined by the "L" region will be expaneded
    # on each side.
    n_reps_h = math.ceil(min(sensor_width/abs(step_xh), sensor_height/abs(step_yh)))
    n_reps_v = math.ceil(min(sensor_width/abs(step_xv), sensor_height/abs(step_yv)))
    
    # % Let the user set the geometric location of the spots
    cx, cy = cog.cog(pupil, np.amax(pupil)/10)
    plt.close()
    plt.figure()
    plt.imshow(pupil)
    
    plt.plot(cx+1, cy+1, 'rx')
    
    #
    plt.colormap='gray'
    
    plt.title='Set the center of the spots, the inner radius and the outer radius'
    p1,p2,p3 = plt.ginput(3)
    x_cir=np.array([p1[0],p2[0],p3[0]])
    y_cir=np.array([p1[1],p2[1],p3[1]])
     # Calculate radius based on the three points that have just been specified
    r_int = np.sqrt((x_cir[1]-x_cir[0])**2 + (y_cir[1]-y_cir[0])**2)
    r_ext = np.sqrt((x_cir[2]-x_cir[0])**2 + (y_cir[2]-y_cir[0])**2)
    r_int_n = r_int / r_ext
    
    #  We start with no spots
    star_refx = np.array([])
    star_refy = np.array([])
    num_stars = 0
    
    #  Sweep over the entire sensor, and determine if we fall within the defined
    # region. If we do, add a spot.
    i_end = n_reps_h*(n_ref_h-1);
    i_start = - i_end;
    j_end = n_reps_v*(n_ref_v-1);
    j_start = - j_end;
    for i in range(int(i_start),int(i_end+1)):
        
        x_init = x_ref[0] + step_xh * i 
        y_init = y_ref[0] + step_yh * i    
        
        for j in range(int(j_start),int(j_end+1)):
            
            star_refx_tmp = x_init + step_xv * j
            star_refy_tmp = y_init + step_yv * j
            
            if star_refx_tmp >= 1 and star_refx_tmp <= sensor_width and star_refy_tmp >= 1 and star_refy_tmp <= sensor_height:                
                center_dist = np.sqrt((star_refx_tmp-x_cir[0])**2 + (star_refy_tmp-y_cir[0])**2)
                if center_dist >= r_int and center_dist <= r_ext:
                    star_refx =np.append(star_refx, star_refx_tmp)
                    star_refy = np.append(star_refy, star_refy_tmp) 
                    num_stars = num_stars + 1
 
    
    star_regx = np.ceil(star_refx - min_dist/2)
    star_regy = np.ceil(star_refy - min_dist/2)
    star_regw = np.tile(int(np.floor(min_dist)), (num_stars,1))
    star_regw=star_regw.astype(int)
    star_regh = star_regw
    
    cx, cy = cog_multiple.cog_multiple(ref, star_threshold_factor, star_regx, star_regy, star_regh, star_regw)
    star_refx = star_regx + cx
    star_refy = star_regy + cy
    
    #Re-calculate regions
    star_regx = np.ceil(star_refx - min_dist/2) # should be integer
    star_regy = np.ceil(star_refy - min_dist/2) #should be integer
    star_regx=star_regx.astype(int)
    star_regy=star_regy.astype(int)
    # #Calculate the normalized coordinates
    star_refx_n = (star_refx - x_cir[0]) / r_ext
    star_refy_n = (star_refy - y_cir[0]) / r_ext
    plt.close()
    return star_refx, star_refy,star_regx,star_regy,star_regw, star_regh,star_refx_n,star_refy_n, x_cir,y_cir, r_ext

    # #  Display the resul
    # #  Display the resul
    # # %figure;
    plt.imshow(ref)
    plt.colormap='gray'
    plt.title=('Centroid positions') 
    plt.plot(star_refx, star_refy, 'g+');
    draw_regs.draw_regs(star_regx, star_regy, star_regw, star_regh);
    if (r_int > 0):
        matplotlib.patches.Rectangle((x_cir[0],y_cir[0]), -r_int, 2*r_int*[1,1],EdgeColor='red')
    
   
    plt.show()
    
    # %save('mat/star_ref_aoli.mat', 'num_stars', 'star_refx', 'star_refy', 'star_regx', 'star_regy', 'star_regw', 'star_regh', 'star_refx_n', 'star_refy_n', 'r_int_n', 'x_cir', 'y_cir');
    
    

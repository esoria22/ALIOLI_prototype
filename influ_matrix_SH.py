# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 10:05:56 2022

@author: esoria
"""
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 16:16:39 2021

@author: esoria
"""

from hcipy import *
def influ_M_SH(star_refx,star_refy,star_regx,star_regy,star_regw, star_regh,num_modes,num_act,influ,wf,deformable_mirror,swfs,binn,dt,probe_amp,TCS_aperture,magnifier,star_zd_pinv):
    import numpy as np
    import proc_pupil
    import matplotlib.pyplot as plt
    import leer_det
    import cog_multiple
    import SH_est
    #lectura= np.zeros((2*np.size(star_regx), 2*num_act))
    i=0
  
    if influ==1:
        num_stars = np.size(star_regx)
        influ_M = np.zeros((num_modes,num_modes))
        star_threshold_factor = 0.1
        #wf.total_power = 1
        deformable_mirror.flatten()
        plt.figure(figsize=(20, 10))
    #genero la matriz de influencia Modos sobre el DM a respuesta en la pirámide
        for ind in range(num_modes):
            slopes = 0
            amps=[-probe_amp,probe_amp]
            for amp in amps:
     
                deformable_mirror.flatten()

                deformable_mirror.actuators[ind]=amp
                dm_wf = deformable_mirror.forward(wf)
                wf_wfs= swfs.forward(magnifier(dm_wf))
                detector = subsample_field(wf_wfs.power, subsampling=2, statistic='sum') * dt 
                detector=detector/detector.sum()
                slopes1 = SH_est.SH_estimator(detector,star_regx, star_regy, star_regh, star_regw,star_refx,star_refy)
                slopes1=star_zd_pinv.dot(slopes1)
                #normalizo la señal leida en la cámara
                slopes += amp * slopes1 /(2*np.var(amps))
              
            # Plot mode response
                plt.clf()
                #plt.suptitle('Mode %d / %d: DM shape' % (h, num_modes))
            
                plt.subplot(1,2,1)
                #plt.title('DM surface') 
                im1 = imshow_field(detector)
            
                plt.subplot(1,2,2)
               
                #plt.title('Slopes map')
                x=np.arange(np.size(slopes1))
                plt.bar(x,slopes1)
               # plt.quiver(star_refx,star_refy,slopes1[:num_stars],slopes1[num_stars:],scale=10,width=0.004,headwidth=3)
                plt.pause(0.1)
            influ_M[:,ind]=slopes
        plt.close()
   
    else:
        num_stars = np.size(star_regx)
        influ_M = np.zeros((num_stars*2,num_act))
        plt.figure(figsize=(20, 10))       
        plt.figure(figsize=(20, 10))
    
        #probe_amp=0.1
        for ind in np.arange(0,num_act):
            slopes = 0
            
        # Probe the phase response
            for s in [-1,1]:
                deformable_mirror.flatten()
                ampli = np.zeros((num_act,))
                ampli[ind] = s * probe_amp
                deformable_mirror.actuators = ampli
                dm_wf = deformable_mirror.forward(wf)
                wf_wfs= swfs.forward(magnifier(dm_wf))
                detector = subsample_field(wf_wfs.power, subsampling=1, statistic='sum') * dt 
                slopes1 = SH_est.SH_estimator(detector,star_regx, star_regy, star_regh, star_regw,star_refx,star_refy)
                #lectura[:,i]=slopes1
                slopes += s * slopes1 /(2*probe_amp)
              
               # i= i+1
             # Plot mode response
                plt.clf()
                #plt.suptitle('Actuator %d / %d: ALPAO 81' % (ind, num_modes))
            
                plt.subplot(1,2,1)
                #plt.title('DM surface')
                im1 = imshow_field(detector)
                #, cmap='RdBu', mask=TCS_aperture)
            
                plt.subplot(1,2,2)
                #plt.title('Slopes map')
                plt.quiver(star_refx,star_refy,slopes1[:num_stars],slopes1[num_stars:],scale=10,width=0.004,headwidth=3)
                plt.pause(0.1)
            influ_M[:,ind]=slopes
        plt.close()
        
    return influ_M

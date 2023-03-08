# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 16:16:39 2021

@author: esoria
"""

from hcipy import *
def influ_M_SH(star_refx,star_refy,star_regx,star_regy,star_regw, star_regh,num_modes,num_act,num_modes_rec,influ,wf,deformable_mirror,swfs,binn,dt,probe_amp,TCS_aperture,magnifier,star_zd_pinv,x_cir,y_cir, r_ext):
    import numpy as np
    import proc_pupil
    import matplotlib.pyplot as plt
    import leer_det
    import cog_multiple
    import SH_est
    #lectura= np.zeros((2*np.size(star_regx), 2*num_act))
    i=0
    t_factor=0
    if influ==1:
        num_stars = np.size(star_regx)
        #influ_M = np.zeros((num_modes_rec,num_modes))
        influ_M = np.zeros((2*num_stars,num_modes))
        star_threshold_factor = 0.1
        wf.total_power = 1
        deformable_mirror.flatten()
        plt.figure(figsize=(20, 10))
    #genero la matriz de influencia Modos sobre el DM a respuesta en la pirámide
        for ind in range(num_modes):
            modes = 0
            
            for s in [-1,1]:
                deformable_mirror.flatten()
                ampli = np.zeros((num_modes,))
                ampli[ind] = s *probe_amp
                deformable_mirror.actuators= ampli
                dm_wf = deformable_mirror.forward(wf)
                wf_wfs= swfs.forward(magnifier(dm_wf))
                detector = subsample_field(wf_wfs.power, subsampling=1, statistic='sum') * dt 
                slopes1 = SH_est.SH_estimator(detector,star_regx, star_regy, star_regh, star_regw,star_refx,star_refy,x_cir,y_cir, r_ext,t_factor)
                mode=star_zd_pinv.dot(slopes1)
                #mode=slopes1
                #normalizo la señal leida en la cámara
                modes +=s* mode /(probe_amp)
              
            # Plot mode response
                plt.clf()
                #plt.suptitle('Mode %d / %d: DM shape' % (h, num_modes))
            
                plt.subplot(1,3,1)
               # plt.title('DM surface') 
                im1 = imshow_field(deformable_mirror.surface, cmap='RdBu', mask=TCS_aperture)
            
                plt.subplot(1,3,2)
                x=np.arange(1,num_modes_rec+1)
                #plt.bar(x,mode)
                imshow_field(detector)
                plt.xlabel('Num mode ANSI')
                plt.ylabel('Amplitude')
                plt.show()
                #plt.title('Slopes map')
               # plt.quiver(star_refx,star_refy,slopes1[:num_stars],slopes1[num_stars:],scale=10,width=0.004,headwidth=3)
                plt.subplot(1,3,3)
                dx=slopes1[: int(np.size(slopes1)/2)]
                dy=slopes1[int(np.size(slopes1)/2):]
                plt.plot(dx,dy,'.')
                plt.xlim([-0.1, 0.1])
                plt.ylim([-0.1, 0.1])
                plt.show()
                plt.pause(0.1)
            influ_M[:,ind]=modes
        plt.close()
   
    else:
        num_stars = np.size(star_regx)
        influ_M = np.zeros((num_act,num_modes_rec))
        plt.figure(figsize=(20, 10))       
        plt.figure(figsize=(20, 10))
        amps = [-probe_amp, probe_amp]
        #probe_amp=0.1
        for ind in np.arange(0,num_act):
            slopes = 0
            
        # Probe the phase response
            for s in [-1,1]:
                deformable_mirror.flatten()
                ampli = np.zeros((num_act,))
                ampli[ind] = s *probe_amp
                deformable_mirror.actuators = ampli
                dm_wf = deformable_mirror.forward(wf)
                wf_wfs= swfs.forward(magnifier(dm_wf))
                detector = subsample_field(wf_wfs.power, subsampling=4, statistic='sum') * dt 
                slopes1 = SH_est.SH_estimator(detector,star_regx, star_regy, star_regh, star_regw,star_refx,star_refy,t_factor)
                #lectura[:,i]=slopes1
                slopes += s * slopes1 / probe_amp
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
            influ_M[:,ind]=star_zd_pinv.dot(slopes)
        plt.close()
        
    return influ_M
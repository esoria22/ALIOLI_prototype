# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 12:43:34 2021

@author: esoria
"""
from hcipy import *
def influ_M(x,imagen_ref,num_modes,influ,wf,deformable_mirror,pwfs,binn,dt,probe_amp,pupilas,mask,TCS_aperture,magnifier_1,th):
    import numpy as np
    import proc_pupil
    import matplotlib.pyplot as plt
    import leer_det
    if influ==1:
        influ_M= np.zeros((2*np.size(x),num_modes))
        plt.figure(figsize=(20, 10))
    #genero la matriz de influencia Modos sobre el DM a respuesta en la pirámide
        for h in np.arange(0,num_modes):
            amps = [-probe_amp, probe_amp]
            influ_M_temp=np.zeros(2*np.size(x))
            slopes=0
            for amp in amps:
                deformable_mirror.flatten()
                deformable_mirror.actuators[h] = amp
                wf_dm_a = deformable_mirror.forward(wf)
                wf_p = pwfs.forward(magnifier_1(wf_dm_a))
                imagen=subsample_field(wf_p.power, subsampling=1, statistic='sum') 
                imagen /= np.sum(imagen)
                influ_M_temp=imagen-imagen-ref
                slopes += amp * influ_M_temp /(np.var(amps))
            # Plot mode response
                plt.clf()
                plt.suptitle('Mode %d / %d: DM shape' % (h, num_modes))
            
                plt.subplot(1,2,1)
                #plt.title('DM surface')
                im1 = imshow_field(deformable_mirror.surface, cmap='RdBu', mask=TCS_aperture)
            
                plt.subplot(1,2,2)
                #plt.title('Slopes map')
                j=int(np.size(influ_M_temp)/2)
                dx=influ_M_temp[:j]
                dy=influ_M_temp[j:]
                x = np.linspace(1,-1,np.size(mask,0))
                y = np.linspace(1,-1,np.size(mask,1))
                x,y = np.meshgrid(x, y)
                x = x[mask>0].ravel()
                y = y[mask>0].ravel()
                plt.quiver(x,y,dx,dy,scale=10,width=0.004,headwidth=2)
                plt.pause(0.1)
            influ_M[:,h]=slopes
        plt.close()
    else:
        influ_M= np.zeros((2*np.size(x),num_modes))
        plt.figure(figsize=(20, 10))
        num_act=num_modes
        amps = [-probe_amp, probe_amp]
        #probe_amp=0.1
        for ind in np.arange(0,num_act):
            slope = 0  
        # Probe the phase response
            for s in [1, -1]:
                amp = np.zeros((num_modes,))
                amp[ind] = s * probe_amp
                deformable_mirror.actuators = amp
        
                wf_dm = deformable_mirror.forward(wf)
                wf_p = pwfs.forward(magnifier_1(wf_dm))     
                imagen = leer_det.leer_det(wf_p,binn,1e-3)
                slope1=proc_pupil.proc_pupil(imagen,pupilas,th)
                slope += s* slope1/(2*np.var(amps))
                
             # Plot mode response
                plt.clf()
                plt.suptitle('Actuator %d / %d: ALPAO 81' % (ind, num_modes))
            
                plt.subplot(1,2,1)
                #plt.title('DM surface')
                im1 = imshow_field(deformable_mirror.surface, cmap='RdBu', mask=TCS_aperture)
            
                plt.subplot(1,2,2)
               # plt.title('Slopes map')
                j=int(np.size(slope1)/2)
                dx=slope1[:j]
                dy=slope1[j:]
                x = np.linspace(1,-1,np.size(mask,0))
                y = np.linspace(1,-1,np.size(mask,1))
                x,y = np.meshgrid(x, y)
                x = x[mask>0].ravel()
                y = y[mask>0].ravel()
                plt.quiver(x,y,dx,dy,scale=10,width=0.002,headwidth=3)
                plt.pause(0.1)
            influ_M[:,ind]=slope
    plt.close()
        
    return influ_M
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 14:06:56 2021

@author: esoria
"""
from proc_pupilTP3 import *
from hcipy import *
def influ_M(angles,coords,num_modes,modes,influ,wf,deformable_mirror,F1,F2,probe_amp,TCS_aperture,binn,dt,distance,H_pinv,magnifier,centro):
    import numpy as np
    import proc_pupil
    import matplotlib.pyplot as plt
    import leer_det
    import TP3_wfs
   
    if influ==1:
        ncoords = np.size(coords)
        nangles=np.size(angles)
        nmodes=np.size(modes)
        influ_M= np.zeros((nmodes, num_modes))
        plt.figure(figsize=(20, 10))
    #genero la matriz de influencia Modos sobre el DM a respuesta en la pirámide
        for h in range(num_modes):
            amps = [-probe_amp, probe_amp]
            influ_M_temp=np.zeros((nmodes))
            slopes=np.zeros((nmodes))
            for amp in amps:
                deformable_mirror.flatten()
                if h==1 or h==0:
                    deformable_mirror.actuators[h] = 4 *amp
                    wf_dm_a = deformable_mirror.forward(wf)
                    wf_wfs= F1.forward(magnifier(wf_dm_a))
                    wf_wfs2= F2.forward(magnifier(wf_dm_a))
                    i1 = subsample_field(wf_wfs.power, subsampling=4, statistic='sum') * dt   
                    i1_=i1/i1.sum()
                    i2 = subsample_field(wf_wfs2.power, subsampling=4, statistic='sum') * dt  
                    i2_=i2/i2.sum()
                    j = int(np.sqrt(np.size(i1_)))
                    i1=np.reshape(i1_,(j,j),order='F')
                    i2=np.reshape(i2_,(j,j),order='F')
                    i1= i1[1:,1:]
                    i2= i2[1:,1:]
                    #i1=proc_pupilTP3(i1,centro)
                    #i2=proc_pupilTP3(i2,centro)
                    influ_M_temp=TP3_wfs.g_wfs(i1, i2, nangles, distance,H_pinv)[1]
                    slopes+= amp*influ_M_temp/(2*4*np.var(amps))
                else:
                    deformable_mirror.actuators[h] = amp
                    wf_dm_a = deformable_mirror.forward(wf)
                    wf_wfs= F1.forward(magnifier(wf_dm_a))
                    wf_wfs2= F2.forward(magnifier(wf_dm_a))
                    i1 = subsample_field(wf_wfs.power, subsampling=4, statistic='sum') * dt   
                    i1_=i1/i1.sum()
                    i2 = subsample_field(wf_wfs2.power, subsampling=4, statistic='sum') * dt  
                    i2_=i2/i2.sum()
                    j = int(np.sqrt(np.size(i1_)))
                    i1=np.reshape(i1_,(j,j),order='F')
                    i2=np.reshape(i2_,(j,j),order='F')
                    i1= i1[1:,1:]
                    i2= i2[1:,1:]
                    #i1=proc_pupilTP3(i1,centro)
                    #i2=proc_pupilTP3(i2,centro)
                    influ_M_temp=TP3_wfs.g_wfs(i1, i2, nangles, distance,H_pinv)[1]
                    slopes+= amp*influ_M_temp/(2*np.var(amps))
               
                
                #Plot mode response
                plt.clf()
                #plt.suptitle('Mode %d / %d: DM shape' % (h, nmodes))
            
                plt.subplot(1,3,1)
                #plt.title('DM surface')
                #x=np.arange(np.size(influ_M_temp))
                #plt.bar(x,influ_M_temp)
                im1 = imshow_field(deformable_mirror.surface, cmap='RdBu', mask=TCS_aperture)       
                plt.subplot(1,3,2)
                imshow_field(i1_,cmap='gray')
                plt.subplot(1,3,3)
                imshow_field(i2_,cmap='gray')
                plt.pause(0.1)
            
            influ_M[:,h]=slopes
        plt.close()
    else:
        ncoords = np.size(coords)
        nangles=np.size(angles)
        nmodes=np.size(modes)
        num_act=num_modes
        influ_M= np.zeros((nmodes,num_act))
        plt.figure(figsize=(20, 10))
        amps = [-probe_amp, probe_amp]
        #probe_amp=0.1
        for ind in np.arange(0,num_act):
            influ_M_temp=np.zeros((nmodes,))
            slopes=np.zeros((nmodes,)) 
        # Probe the phase response
            for s in [1, -1]:
                deformable_mirror.flatten()
                amp = np.zeros((num_act,))
                amp[ind] = s * probe_amp            
                deformable_mirror.actuators = amp
                wf_dm_a = deformable_mirror.forward(wf)
                wf_wfs= F1.forward(magnifier(wf_dm_a))
                wf_wfs2= F2.forward(magnifier(wf_dm_a))
                i1 = subsample_field(wf_wfs.power, subsampling=4, statistic='sum') * dt 
                i1_ = i1/i1.sum()
                i2 = subsample_field(wf_wfs2.power, subsampling=4, statistic='sum') * dt 
                i2_ = i2/i2.sum() 
                j = int(np.sqrt(np.size(i1_)))
                i1=np.reshape(i1_,(j,j),order='F')
                i1= i1[1:,1:]
                i2=np.reshape(i2_,(j,j),order='F')
                i2= i2[1:,1:]
                influ_M_temp=TP3_wfs.g_wfs(i1, i2, nangles, distance,H_pinv)[1]
                slopes+= (s*probe_amp*influ_M_temp)/(2*np.var(amps))
                
                #Plot mode response
                plt.clf()
                #plt.suptitle('Act %d / %d: DM shape' % (ind+1, num_act))
            
                plt.subplot(1,3,1)
                #plt.title('DM surface')
                im1 = imshow_field(deformable_mirror.surface, cmap='RdBu', mask=TCS_aperture)       
                plt.subplot(1,3,2)
                imshow_field(i1_,cmap='gray')
                plt.subplot(1,3,3)
                imshow_field(i2_,cmap='gray')
                plt.pause(0.1)
            influ_M[:,ind]=slopes
        plt.close()
        
    return influ_M
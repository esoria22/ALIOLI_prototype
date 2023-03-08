#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 12:19:30 2020

@author: esthersoriahernandez
 return swfs, Tsh_aperture,pup_array_grid
    shwfs,Tsh_aperture,pupil_array_grid=SH_wfs.SH(telescope_diameter, wavelength_wfs)
"""

from hcipy import *
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
import proc_pupil
import detect_circle
import influ_matrix_Pyr
import influ_matrix_SH
import influ_matrix_TP3
import SH_wfs2
import leer_det
import dm_char
import TP3_wfs
import SH_wfs
import setref2
import cog_multiple
import SH_est
import radon as rd
import g_wfs_precal as gprec
import ansi2dual
import Zernike
import TP3_wfs
import math 
import test_influ
import zernike_dx2
import zernike_dy2
from pup_det_TP3 import *
from proc_pupilTP3 import *
import time
#Defino las caracteristicas de mi sistema, wvl, apertura..
for m in np.arange(2,3):
    plt.close('all')
    type_wfs=m
    if m==1:
        name="Pyramid WFS"
    if m==2:
        name="Shack-Hartmann WFS"
    if m==3:
        name="TP3 WFS"
    influ=1 # 1 si es modos, 2 si quiero mover actuador a actuador
    char=1#1 si es modos, 2 si es actuadores
     #1 wfs pyramide, 2 SH, 3 TP3
    N =66 #num modes or actuators
    pup_TP3=2*10**-3
    #caracteristicas de el detector de WFS
    dt= 1 #tiempo de exposición
    binn =1 # binning
    delta_t=1e-3
    wavelength_wfs = 700.0E-9
    #telescope_diameter = 1.5*(0.25/(1.5*13.8))
    telescope_diameter = 1.52
    pup_resz= 3.6*10**-3
    zero_magnitude_flux = 3.9E10
    if type_wfs == 1:
        pupil_grid = make_pupil_grid(256, telescope_diameter)
    if type_wfs ==2: 
        pupil_grid=make_pupil_grid(256, telescope_diameter)
    if type_wfs ==3: 
        pupil_grid=make_pupil_grid(256, telescope_diameter)
    p1=[-1,-1.2]
    p2=[-2,0]
    spy=0
    if spy==1:       
        spider_width=0.01
        spider_offset = [0, 0]  # meter
        spider_offset = np.array(spider_offset)
        mirror_edge1 = (telescope_diameter / (2 * np.sqrt(2)), telescope_diameter / (2 * np.sqrt(2)))
        mirror_edge2 = (-telescope_diameter / (2 * np.sqrt(2)), telescope_diameter / (2 * np.sqrt(2)))
        mirror_edge3 = (telescope_diameter / (2 * np.sqrt(2)), -telescope_diameter / (2 * np.sqrt(2)))
        mirror_edge4 = (-telescope_diameter / (2 * np.sqrt(2)), -telescope_diameter / (2 * np.sqrt(2)))
        spyder1=make_spider(spider_offset, mirror_edge1, spider_width)
        spyder2=make_spider(spider_offset, mirror_edge2, spider_width)
        spyder3=make_spider(spider_offset, mirror_edge3, spider_width)
        spyder4=make_spider(spider_offset, mirror_edge4, spider_width)
        aperture=make_obstructed_circular_aperture(telescope_diameter, 0.35)
        TCS_aperture = Field(aperture(pupil_grid)*spyder1(pupil_grid)*spyder2(pupil_grid)*spyder3(pupil_grid)*spyder4(pupil_grid),pupil_grid)
        imshow_field(TCS_aperture)
    else:
         aperture=circular_aperture(telescope_diameter)
         TCS_aperture = Field(aperture(pupil_grid),pupil_grid)
    #defino el numero de modos con el que quiero trabajar
    wf= Wavefront(TCS_aperture,wavelength_wfs)
    wf.total_power = zero_magnitude_flux
    
    if type_wfs==1:
        pup_grid=make_pupil_grid(20,2.5*10**-3)
        pwfs=PyramidWavefrontSensorOptics(pup_grid, separation =0.0045, wavelength_0=wavelength_wfs)
        magnification=2.5*10**-3 /telescope_diameter
        magnifier_1 = Magnifier(magnification)
        camera = NoiselessDetector(pwfs.output_grid)
    if type_wfs ==2:   
        swfs, Tsh_aperture = SH_wfs2.SH2(pup_resz,wavelength_wfs) 
        magnification= pup_resz/telescope_diameter
        magnifier_2 = Magnifier(magnification)
    if type_wfs ==3:
        magnification= pup_TP3/telescope_diameter
        magnifier_3= Magnifier(magnification)
        distance =3/1000
        pup_grid=make_pupil_grid(60,pup_TP3)
        F1 = FresnelPropagator(pup_grid, distance, num_oversampling=4, refractive_index=1)
        F2 = FresnelPropagator(pup_grid, -distance, num_oversampling=4, refractive_index=1)
    deformable_mirror,num_modes,num_act = dm_char.DM(N,telescope_diameter,pupil_grid,char) 
   # aper,segments=make_hexagonal_segmented_aperture(1,1, 0.1,starting_ring=0,return_segments=True)
    deformable_mirror2, num_m,num_a=dm_char.DM(N,telescope_diameter,pupil_grid,1)  
    #imshow_field(deformable_mirror.surface, mask=TCS_aperture,cmap='seismic')
    
    #TOMO IMAGEN DE REFERENCIA
    #aplano el DM y tomo como imagen de referencia la lectura en la cámara
    deformable_mirror.flatten()
    deformable_mirror2.flatten()
    wf_dm2 = deformable_mirror2.forward(wf)
    wf_dm_a = deformable_mirror.forward(wf_dm2)    
  


    if type_wfs==1:
        wf_wfs = pwfs.forward(magnifier_1(wf_dm_a))
    if type_wfs ==2:
        wf_wfs= swfs.forward(magnifier_2(wf_dm_a))
    if type_wfs ==3:   
        wf_wfs= F1.forward(magnifier_3(wf_dm_a))
        wf_wfs2= F2.forward(magnifier_3(wf_dm_a))
    #leo el valor de la intensidad del frente de ondas que incide sobre el detector
    if type_wfs == 3:
        i1 = subsample_field(wf_wfs.power, subsampling=4, statistic='sum') * dt 
        i1 = i1/i1.sum()
        i2 = subsample_field(wf_wfs2.power, subsampling=4, statistic='sum') * dt 
        i2 = i2/i2.sum()
        #centro=pup_det(i1,0.1)
   
    if type_wfs==2:
        detector = subsample_field(wf_wfs.power, subsampling=2, statistic='sum') * dt
        imshow_field(detector)
        #normalizo la señal leida en la cámara
        j = int(np.sqrt(detector.shape))
        imagen_ref=np.reshape(detector,(j,j),order='C')
        imagen_ref /= np.sum(imagen_ref)
        #plt.imshow(imagen_ref)
        #me guardo en un archivo fits la imagen de referencia de mi pirámide
        write_fits(detector,'im_ref1.fits')
        th=imagen_ref[int(j/2),int(j/2)]
        #print(th)
    if type_wfs==1:
         imagen_ref=subsample_field(wf_wfs.power, subsampling=1, statistic='sum') * dt          
         imagen_ref /= np.sum(imagen_ref)
         write_fits(imagen_ref,'im_ref1.fits')
         th=0.1
         #image_ref /= np.sum(image_ref)
    #Utilizo la imagen de referencia para detectar las posiciones de las pupilas sobre el detector
    if type_wfs==1:
        '''
        if binn == 1:
            pupilas = detect_circle.detect_circle('im_ref1.fits',10,14)
        if binn == 2:
            pupilas = detect_circle.detect_circle('im_ref1.fits',9,15)
        radios =pupilas[:,:,2]
        rad=min(min(radios))
        mask = np.zeros((2*rad,2*rad))
        for i in np.arange(0,2*rad,1):
            for j in np.arange(0,2*rad,1):
               if np.sqrt((i-rad)**2+(j-rad)**2)<rad:
                    mask[i,j]=1
        #utilizando la máscara anterior me quedó unicamente con los píxeles que caen dentro de la pupila
        x = np.linspace(1,-1,np.size(mask,0))
        y = np.linspace(1,-1,np.size(mask,1))
        x,y = np.meshgrid(x, y)
        x = x[mask>0].ravel()
        y= y[mask>0].ravel()
        probe_amp=0.045*wavelength_wfs
        #HAGO LA CARACTERIZACION ESTÁTICA
        #plt.title('Actuators DM')
        #ref=proc_pupil.proc_pupil(imagen_ref,pupilas,th)
        influ_M=influ_matrix_Pyr.influ_M(x,ref,num_modes,influ,wf,deformable_mirror,pwfs,binn,dt,probe_amp,pupilas,mask,TCS_aperture,magnifier_1,th)
        '''
        probe_amp=0.01*wavelength_wfs
        influ_M=test_influ.influ(num_modes, probe_amp, deformable_mirror, camera, pwfs, imagen_ref,magnifier_1,wf)
        
    if type_wfs==2:
        star_refx, star_refy,star_regx,star_regy,star_regw, star_regh,star_refx_n,star_refy_n, x_cir,y_cir, r_ext=setref2.setref(imagen_ref,0.05,0)
        probe_amp=0.1*wavelength_wfs
        num_modes_rec=66
        num_stars=np.size(star_refx)
        star_zdx=np.zeros((num_stars,num_modes_rec))
        star_zdy=np.zeros((num_stars,num_modes_rec))
        #calcular coeficiente
        for j in range(1,num_modes_rec+1):
            n,m = ansi2dual.ansi2dual(j) 
            star_zdx[:,j-1]=zernike_dx2.zernike_dx(n,m,star_refx_n,star_refy_n)
            star_zdy[:,j-1]=zernike_dy2.zernike_dy(n,m,star_refx_n,star_refy_n)
        star_zd=np.zeros((2*num_stars,num_modes_rec))
        star_zd[0:num_stars,:]=star_zdx
        star_zd[num_stars:,:]=star_zdy 
        star_zd_pinv = inverse_tikhonov(star_zd, rcond=1e-2, svd=None)
        influ_M = influ_matrix_SH.influ_M_SH(star_refx,star_refy,star_regx,star_regy,star_regw, star_regh,num_modes,num_act,influ,wf,deformable_mirror,swfs,binn,dt,probe_amp,TCS_aperture,magnifier_2,star_zd_pinv)
        m=2
    if type_wfs==3:
        centro=0
        j = int(np.sqrt(np.size(i1)))
        i1=np.reshape(i1,(j,j),order='F')
        i2=np.reshape(i2,(j,j),order='F')
        i1= i1[1:,1:]
        i2= i2[1:,1:]
        i1_=i1/i1.sum()
        i2_=i2/i2.sum()
        #i1=proc_pupilTP3(i1,centro)
       # i2=proc_pupilTP3(i2,centro) 
        nangles = 20
        angles_full = np.linspace(0, 180, nangles + 1)
        angles = angles_full[0:- 1]# 0 degrees is equal to 180 degrees
        i1_dummy = np.zeros((j, j))
        i1_dummy= i1_dummy[1:,1:]
        nmodes = 66
        modes = np.arange(1,nmodes+1)
        probe_amp=0.1*wavelength_wfs
        H = gprec.g_wfs_precalc(modes, angles, j-1, j-1)
        H_reshape = np.reshape(H, (np.size(H,0)*np.size(H,1), np.size(H,2)),order="F")
        H_pinv = inverse_tikhonov(H_reshape, rcond=1e-3, svd=None)
        rtx,coords = rd.radon(i1_dummy, angles,circle= True)   
        influ_M= influ_matrix_TP3.influ_M(angles,coords,num_modes,modes,influ,wf,deformable_mirror,F1,F2,probe_amp,TCS_aperture, binn, dt, distance,H_pinv,magnifier_3,centro)
  
    rcond = 0.05
    reconstruction_matrix= inverse_tikhonov(influ_M, rcond=rcond, svd=None)
   
              
   
    rec="rec"+str(m)
    globals()[rec]=reconstruction_matrix
'''
    modes=np.arange(0,66,1)
    mean = 0.0   # some constant  # some constant (standard deviation)
    amp=np.arange(0.05,1.5,0.05)
    for mode in modes:
        rms="RMS_"+str(m)+"_"+str(mode)+"g"
        rms_rec="RMS_"+str(m)+"_"+str(mode)+"r"
        globals()[rms]=np.zeros((np.size(amp),))
        globals()[rms_rec]=np.zeros((np.size(amp),))
        for i in range(np.size(amp)):
           
            #zero_magnitude_flux*10**(-i/2.5)
            #plt.figure(2)
            deformable_mirror.flatten()
            deformable_mirror2.flatten()
            
            deformable_mirror2.actuators[mode] =amp[i]*wavelength_wfs
            
           
            
            #deformable_mirror2.random(0.07*wavelength_wfs)
            wf_dm2 = deformable_mirror2.forward(wf)
            wf_dm_a = deformable_mirror.forward(wf_dm2)
            s=deformable_mirror2.surface
            h=deformable_mirror2.actuators
           # plt.subplot(1,2,1)
            #imshow_field(deformable_mirror2.surface, mask=TCS_aperture,cmap='seismic')
            #plt.colorbar()
           
            y_predicted =np.zeros([np.size(s),])
            MSE = mean_squared_error(s, y_predicted)
            RMSE = math.sqrt(MSE)
            globals()[rms][i]=RMSE
            print('RMS:')
            print(RMSE)
            if type_wfs==1:           
                wf_wfs = pwfs.forward(magnifier_1(wf_dm_a))
                imagen=subsample_field(wf_wfs.power, subsampling=1, statistic='sum')*dt
                imagen /= np.sum(imagen)
                analy=imagen
                #ruido=0*np.random.normal(mean,  10**-6*i, image.shape)
                #ruido2=np.sum(abs(ruido))
                #j = int(np.sqrt(np.size(image)))
                #detector=np.reshape(image,(j,j),order='F')
                #noise=np.reshape(ruido,(j,j),order='F')
                #camara=detector+ noise
                #analy=image+ruido
                #señal=np.sum(abs(camara))
                #print(señal/ruido2)
                #plt.figure(2)
                #imshow_field(analy,cmap='gray')
                slopes = analy-imagen_ref
               
                mod="mod"+str(m)
                globals()[mod]=globals()[rec].dot(slopes)
             
           
              
            if type_wfs ==2:
                
                wf_wfs= swfs.forward(magnifier_2(wf_dm_a))
                
                detector = subsample_field(wf_wfs.power, subsampling=2, statistic='sum') * dt
                detector=detector/detector.sum()
        
                
                slopes = SH_est.SH_estimator(detector,star_regx, star_regy, star_regh, star_regw,star_refx,star_refy)
                slopes=star_zd_pinv.dot(slopes)
                mod="mod"+str(m)
                globals()[mod]=globals()[rec].dot(slopes)
                
            if type_wfs ==3: 
               
                wf_wfs= F1.forward(magnifier_3(wf_dm_a))
                wf_wfs2= F2.forward(magnifier_3(wf_dm_a))
                i1 = subsample_field(wf_wfs.power, subsampling=4, statistic='sum') * dt 
                i1=i1/i1.sum()
                i2 = subsample_field(wf_wfs2.power, subsampling=4, statistic='sum') * dt 
                i2=i2/i2.sum()
                j = int(np.sqrt(np.size(i1)))       
                i1=np.reshape(i1,(j,j),order='F')
                i2=np.reshape(i2,(j,j),order='F')
                i1= i1[1:,1:]
                i2= i2[1:,1:]
                
                #ruido=0*np.random.normal(mean,  10**-6*i, i1.shape)
                #ruido2=np.sum(abs(ruido))
                #rcamara1= i1+ ruido
                #rseñal=np.sum(abs(camara1))
                #rplt.imshow(camara1)
                #rprint(señal/ruido2)
                #rcamara2= i2+ ruido
                
                #i1=proc_pupilTP3(i1,centro)
                #i2=proc_pupilTP3(i2,centro)
                #plt.figure()
                #plt.imshow(camara1,cmap='gray')
                #plt.figure()
                #plt.imshow(camara2,cmap='gray')
                #mod1=rec1.dot(slopes)
                rec="rec"+str(m)
                mod="mod"+str(m)
                slopes = TP3_wfs.g_wfs(i1, i2, nangles, distance,H_pinv)[1]
                
                globals()[mod]=globals()[rec].dot(slopes)
           
            mod="mod"+str(m)
            deformable_mirror.actuators=globals()[mod]
            j2=deformable_mirror.actuators
            sr=deformable_mirror.surface
            #plt.figure()
            #plt.subplot(1,2,2)
            #imshow_field(deformable_mirror.surface, mask=TCS_aperture,cmap='seismic')
            #plt.colorbar()
            MSE_r = mean_squared_error(sr, y_predicted)
            RMSEr = math.sqrt(MSE_r)
            globals()[rms_rec][i]=RMSEr
            print(RMSEr)
            
modes=np.arange(0,66)   
values= np.zeros((66,))     
for mode in modes:
    rms="RMS_"+str(m)+"_"+str(mode)+"g"
    rms_rec="RMS_"+str(m)+"_"+str(mode)+"r"
    fact="fact"+str(mode)
    globals()[fact]=np.mean(globals()[rms])/np.mean(globals()[rms_rec])
    values[mode]=globals()[fact]
plt.figure()
fig, ax = plt.subplots(figsize=(12, 12))
ax.set_xlabel('# Mode ANSI',fontsize=25)
ax.set_ylabel('Gain Factor',fontsize=25)
ax.grid()
plt.plot(modes[1:],values[1:],linewidth=3)
'''   
#####sensibilidaddddd
magnitudes=[0,1,2,4,5,6,7,8,9,10]
amp=[0.001,0.01, 0.1, 0.5, 0.8, 1.2,5]
magnitudes=[0,1,2,4,5,6,7]
mean=0
mode=3
i=0
for magn in magnitudes:
    wf.total_power=zero_magnitude_flux*10**(-magn/2.5)
    rms="RMSs_"+str(m)+"_"+str(mode)+"g"
    rms_rec="RMSs_"+str(m)+"_"+str(mode)+"r"
    globals()[rms]=np.zeros((np.size(magnitudes),))
    globals()[rms_rec]=np.zeros((np.size(magnitudes),))
    
    deformable_mirror.flatten()
    deformable_mirror2.flatten()
    
    
    deformable_mirror2.random(0.1* wavelength_wfs)
    wf_dm2 = deformable_mirror2.forward(wf)
    wf_dm_a = deformable_mirror.forward(wf_dm2)
    s=deformable_mirror2.surface
    h=deformable_mirror2.actuators
    plt.figure(i)
    plt.subplot(1,3,1)
    imshow_field(deformable_mirror2.surface, mask=TCS_aperture,cmap='seismic')
    plt.colorbar()
    
    y_predicted =np.zeros([np.size(s),]) 
    MSE = mean_squared_error(s, y_predicted)  
    RMSE = math.sqrt(MSE)
    globals()[rms][i]=RMSE
    print('RMS:')
    print(RMSE)
    if type_wfs==1:           
        wf_wfs = pwfs.forward(magnifier_1(wf_dm_a))
        imagen=subsample_field(wf_wfs.power, subsampling=1, statistic='sum') * dt             
        ruido=np.random.normal(mean, 10**3, imagen.shape)
        ruido2=np.sum(abs(ruido))
        j = int(np.sqrt(np.size(imagen)))
        analy=imagen+ruido
        analy /= np.sum(analy)
        
       #señal=np.sum(abs(camara))
        #print(señal/ruido2)
        plt.figure(i)
        plt.subplot(1,3,2)
        imshow_field(analy,cmap='gray')
        slope =analy-imagen_ref
        #slopes = proc_pupil.proc_pupil((camara),pupilas,th)
        rec="rec"+str(m)
        mod="mod"+str(m)
        globals()[mod]=globals()[rec].dot(slope)
     
    
       
    if type_wfs ==2:
        wf_wfs= swfs.forward(magnifier_2(wf_dm_a))
        
        detector = subsample_field(wf_wfs.power, subsampling=2, statistic='sum') * dt
        detector=abs(detector)
        j = int(np.sqrt(np.size(detector)))
        detector1=np.reshape(detector,(j,j),order='F')
        noise=np.random.normal(mean, 10**3, detector.shape)
        ruido=np.reshape(noise,(j,j),order='F')
        ruido2=np.sum(abs(ruido))
        camara=detector1+ruido
        señal=np.sum(abs(camara))
        print(señal/ruido2)
        plt.figure(i)
        plt.subplot(1,3,2)
        plt.imshow(camara,cmap='gray')
        t_factor=0
        img=detector+noise
        img=img/img.sum()
        ini=time.time()
        slopes = SH_est.SH_estimator(img,star_regx, star_regy, star_regh, star_regw,star_refx,star_refy)
        mode=star_zd_pinv.dot(slopes)
        mod="mod"+str(m)
        globals()[mod]=globals()[rec].dot(mode)
        fin=time.time()
        print(fin-ini)
         
    if type_wfs ==3: 
       
        wf_wfs= F1.forward(magnifier_3(wf_dm_a))
        wf_wfs2= F2.forward(magnifier_3(wf_dm_a))
        i1 = subsample_field(wf_wfs.power, subsampling=4, statistic='sum') * dt /2
       
        i2 = subsample_field(wf_wfs2.power, subsampling=4, statistic='sum') * dt /2
       
        j = int(np.sqrt(np.size(i1)))       
        i1=np.reshape(i1,(j,j),order='F')
        i2=np.reshape(i2,(j,j),order='F')
        i1= i1[1:,1:]
        i2= i2[1:,1:]
        noise=np.random.normal(mean, 10**3, i1.shape)
        
        ruido2=np.sum(abs(noise))
        i1= i1+ noise
        i1=i1/i1.sum()
        señal1=np.sum(abs(i1))
        señal2=np.sum(abs(i2))
        #plt.imshow(i1)
        print(señal2/ruido2)
        i2= i2+ noise
        i2=i2/i2.sum()
        ratio=señal2/señal1
        #i1=proc_pupilTP3(i1,centro)
        #i2=proc_pupilTP3(i2,centro)
        plt.figure()
        plt.imshow(i1,cmap='gray')
        plt.colorbar()
        
        plt.figure()
        plt.imshow(i2,cmap='gray')
        plt.colorbar()
        #mod1=rec1.dot(slopes)
        rec="rec"+str(m)
        mod="mod"+str(m)
        
        slopes = TP3_wfs.g_wfs(i1, i2, nangles, distance,H_pinv)[1]
       
        globals()[mod]=globals()[rec].dot(slopes)
        
    mod="mod"+str(m)
    deformable_mirror.actuators=globals()[mod]
    plt.figure(i)
    plt.subplot(1,3,2)
    imshow_field(deformable_mirror.surface, mask=TCS_aperture,cmap='seismic')
    plt.colorbar()
    MSE_r = mean_squared_error(sr, y_predicted)
    RMSEr = math.sqrt(MSE_r)
    globals()[rms_rec][i]=RMSEr
    print(RMSEr)
    j2=globals()[mod]
    i=i+1


                 
            
'''        
    plt.close('all')
    plt.figure()
    deformable_mirror.actuators=-j2
    #deformable_mirror2.actuators=h
    
    spatial_resolution = wavelength_wfs / telescope_diameter
    focal_grid = make_focal_grid(q=8, num_airy=20, spatial_resolution=spatial_resolution)
    
    propagator = FraunhoferPropagator(pupil_grid, focal_grid)   
    PSF_in = propagator.forward(deformable_mirror.forward(deformable_mirror2.forward(wf))).power 
    print(np.max(PSF_in))
    
    imshow_field(np.log10(PSF_in)+10)
    plt.colorbar()
    plt.show() 
    plt.figure()
    imshow_field((PSF_in))
    plt.colorbar()
    plt.show() 
   
    
    plt.figure()
    deformable_mirror.flatten()
    #deformable_mirror.actuators=-j2
    #deformable_mirror2.actuators=h
    
    spatial_resolution = wavelength_wfs / telescope_diameter
    focal_grid = make_focal_grid(q=8, num_airy=20, spatial_resolution=spatial_resolution)
    
    propagator = FraunhoferPropagator(pupil_grid, focal_grid)   
    PSF_in = propagator.forward(deformable_mirror.forward(deformable_mirror2.forward(wf))).power 
    print(np.max(PSF_in))
    
    imshow_field(np.log10(PSF_in)+10)
    plt.colorbar()
    plt.show() 
    plt.figure()
    imshow_field((PSF_in))
    plt.colorbar()
    plt.show() 
  

def rebin(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)


coc=RMS_3_21g/RMS_3_21r
###PLOTING###
plt.figure()
fig, ax = plt.subplots(figsize=(12, 12))
#ax.plot(RMS_0_3g*1E9,RMS_0_3r*1E9,'-b',linewidth=3,label="Shack Hartmann")
ax.plot(RMS_3_5g*1E9,RMS_3_5g*1E9,'--k',label="Ideal response")
#ax.plot(RMS_1_3g*1E9,RMS_1_3r*1E9,'-r',linewidth=3,label="Pyramid")
ax.plot(RMS_3_5g*1E9,values[5]*RMS_3_5r*1E9,'-g',linewidth=3,label="TP3")

ax.set_xlabel('RMS of applied aberration (nm)',fontsize=25)
ax.set_ylabel('RMS of Rec aberration (nm)',fontsize=25)
ax.grid()
ax.legend(fontsize=24)


plt.figure()
fig, ax = plt.subplots(figsize=(12, 12))
ax.plot(RMS_0_14g*1E9,RMS_0_14r*1E9,'-b',linewidth=3,label="Shack Hartmann")
ax.plot(RMS_1_14g*1E9,RMS_1_14g*1E9,'--k',label="Ideal response")
ax.plot(RMS_1_14g*1E9,RMS_1_14r*1E9,'-r',linewidth=3,label="Pyramid")
ax.plot(RMS_3_14g*1E9,RMS_3_14r*1E9,'-g',linewidth=3,label="TP3")

ax.set_xlabel('RMS of applied aberration (nm)',fontsize=25)
ax.set_ylabel('RMS of Rec aberration (nm)',fontsize=25)
ax.grid()
ax.legend(fontsize=24)


plt.figure()
fig, ax = plt.subplots(figsize=(12, 12))
ax.plot(RMS_2_43g*1E9,RMS_2_43r*1E9,'-b',linewidth=3,label="Shack Hartmann")
ax.plot(RMS_1_43g*1E9,RMS_1_43g*1E9,'--k',label="Ideal response")
ax.plot(RMS_1_43g*1E9,RMS_1_43r*1E9,'-r',linewidth=3,label="Pyramid")
ax.plot(RMS_3_43g*1E9,RMS_3_43r*1E9,'-g',linewidth=3,label="TP3")

ax.set_xlabel('RMS of applied aberration (nm)',fontsize=25)
ax.set_ylabel('RMS of Rec aberration (nm)',fontsize=25)
ax.grid()
ax.legend(fontsize=24)




fig, ax = plt.subplots(figsize=(12, 12))

ax.plot(RMS_2_4g*1E9,RMS_2_4g*1E9,'--k',label="Ideal response")

ax.plot(RMS_2_4g*1E9,RMS_2_4r*1E9,'-g',linewidth=3,label="SH")

ax.set_xlabel('RMS of applied aberration (nm)',fontsize=25)
ax.set_ylabel('RMS of Rec aberration (nm)',fontsize=25)
ax.grid()
ax.legend(fontsize=24)


plt.figure()
fig, ax = plt.subplots(figsize=(12, 12))

ax.plot(RMS_3_7g*1E9,RMS_3_7g*1E9,'--k',label="Ideal response")

ax.plot(RMS_3_7g*1E9,RMS_3_7r*1E9*values[7],'-g',linewidth=3,label="TP3")

ax.set_xlabel('RMS of applied aberration (nm)',fontsize=25)
ax.set_ylabel('RMS of Rec aberration (nm)',fontsize=25)
ax.grid()
ax.legend(fontsize=24)
'''
 ####ANALIZAR LA MATRIZ DE INFLUENCIA ####
    
'''
     for i in range(np.size(influ_M,0)):
        for j in range(np.size(influ_M,1)):
            if influ_M[i,j]<0:
                influ_M[i,j]=0
    plt.figure()
    influ_Mp=influ_M[5:,:]
    plt.imshow(influ_M,origin='lower')
    plt.colorbar()
    plt.xlabel('Zernike mode applied (ANSI index)')
    plt.ylabel('Zernike mode measured (ANSI index)')
    rcond = 0.05
    reconstruction_matrix= inverse_tikhonov(influ_M, rcond=rcond, svd=None)
   
    plt.imshow(reconstruction_matrix,origin='lower')
    u, s, vh = np.linalg.svd(influ_M, full_matrices=True)
    x=np.arange(1,np.size(influ_M,0)+1)
    plt.figure()
    plt.plot(x,np.log(s))
    plt.xlim
    plt.grid()
    plt.ylabel('Eigenvalues')
    plt.xlabel('Number of eigenmode')
    plt.show()
'''
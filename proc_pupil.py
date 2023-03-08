# -*- coding: utf-8 -*-t
"""
Created on Thu Feb  4 13:51:53 2021

@author: esoria
"""

def proc_pupil(aberr,pupilas,th):
    import numpy as  np
    centers_y=pupilas[:,:,1]
    centers_x=pupilas[:,:,0]
    radios =pupilas[:,:,2]
#me quedo con el valor mínimo del vector radios para que las cuatro matrices correspondientes a cada pupila sean iguales
    rad=min(min(radios))
    #me guardo una máscara con las posiciones de la pupila en cada matriz
    mask = np.zeros((2*rad,2*rad))
    aberr=aberr
    for i in np.arange(0,2*rad,1):
        for j in np.arange(0,2*rad,1):
           if np.sqrt((i-rad)**2+(j-rad)**2)<rad:
                mask[i,j]=1

    #genero las 4 matrices cuadradas sobre las que cae cada pupila, el procedimiento es el siguiente:
        #- genero la matriz cuadrada de 0 sobre la que mi pupila está circunscrita.
        #- recorro esta matriz  y si ek valor de la coodenada r(en polares) cae dentro de mi pupila
        #le asigno a ese punto el valor de la intensidad del punto cooresponiente en la imagen grande, la del detector.
    

    I_d = np.zeros((2*rad,2*rad))
    for i in np.arange(0,2*rad,1):
        for j in np.arange(0,2*rad,1):
           if np.sqrt(abs(i-rad)**2+(j-rad)**2)<rad:
                I_d[i,j]=aberr[i-rad+centers_x[0,1], j - rad+ centers_y[0,1]]
  
    I_b = np.zeros((2*rad,2*rad))
    for i in np.arange(0,2*rad,1):
        for j in np.arange(1,2*rad,1):
           if np.sqrt(abs(i-rad)**2+(j-rad)**2)<rad:
                I_b[i,j]=aberr[i-rad+centers_x[0,2], j - rad+ centers_y[0,2]]
  
    I_a = np.zeros((2*rad,2*rad))
    for i in np.arange(0,2*rad,1):
        for j in np.arange(1,2*rad,1):
           if np.sqrt(abs(i-rad)**2+(j-rad)**2)<rad:
                I_a[i,j]=aberr[i-rad+centers_x[0,0], j - rad+ centers_y[0,0]]
 
    I_c = np.zeros((2*rad,2*rad))
    for i in np.arange(0,2*rad,1):
        for j in np.arange(1,2*rad,1):
           if np.sqrt(abs(i-rad)**2+(j-rad)**2)<rad:
                I_c[i,j]=aberr[i-rad+centers_x[0,3], j - rad+ centers_y[0,3]]

       #aplico el algoritmo 
    
    norm = I_a + I_b + I_c + I_d
    #rm = np.zeros_like(norm)
    #rm[norm != 0] = 1 / norm[norm != 0]

    I_x = (I_a + I_b - I_c - I_d) / norm
    I_y = (I_a - I_b - I_c + I_d) / norm
    #me quedo con los valores de las pendientes de los puntos que caen dentro de la pupila
    
    I_x = I_x[mask>0].ravel()
    
    I_y = I_y[mask>0].ravel()
    
    #doy como salida los valores de las pendientes medidas en un vector
    dx=np.append(I_x,I_y)
    return dx
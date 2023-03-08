import cv2
import numpy as np
from astropy.io import fits
import os
import pyfits as ft
'''

#This function reads the reference image and gives back the position and ratio of the pupils founded.

Inputs:
	*direction_image: the path of the reference image taken (fits file)
	*minR: minimum ratio expected (in pixels)
	*maxR: maximum ratio expected (in pixels)  

Outputs:
	*circles: matrix with the circles data
		** first column: centers x
		** second column: centers y
		** third column: ratios
'''


def detect_circle(direction_image, minR,maxR):
    ff=ft.open(direction_image,uint=True)     #reading the image
    I=ff[0].data
    h=np.shape(I)
    if np.shape(h)==3:
        I=I[0,:,:]
    #I = fits.getdata(os.path.join(inputdir,d[0]))
    #u
    #change the data type to uint8
    img=I.astype(np.float32)  #I.astype(np.uint8)
    #normalizing
    cv2.normalize(img,img,0,255, cv2.NORM_MINMAX, cv2.CV_32F)
    img= img.astype(np.uint8)
    
    #cimg = cv2.cvtColor(img,cv2.COLOR_GRAY2BGR)
    #aplying median filter
    cimg = cv2.medianBlur(img,5)
    #aplying the Hought filter
    circles = cv2.HoughCircles(img,cv2.HOUGH_GRADIENT,1,20, param1=10,param2=10, minRadius=minR,maxRadius=maxR)
    circles = np.uint16(np.around(circles))
    #plotting the detectes circles over the image 
    for i in circles[0,:]:
        # draw the outer circle
        cv2.circle(cimg,(i[0],i[1]),i[2],(0,255,0),1)
        # draw the center of the circle
        cv2.circle(cimg,(i[0],i[1]),2,(0,0,255),1)
    
    cv2.imshow('detected circles',cimg)
    cv2.waitKey(0)
    cv2.destroyAllWindows()
    #saving the image in png format with the detected circles
    outputdir="C:/Users/esoria/Documents/Sin_pupilas_TP3/result_img.png"
    
    cv2.imwrite(outputdir,cimg)
    return circles

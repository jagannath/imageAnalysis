#! /home/jaggu/anaconda/bin/python2.7

from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
import os
import sys
import cv2
import collections


class BeadExpt(object):
    def __init__(self,subDir,pathDir):
        self.subDir = subDir
        self.pathDir = pathDir
        self.pathSubDir = os.path.join(pathDir,self.subDir)
        self.imageDir = os.path.join(self.pathSubDir,"Graphs")
        if not os.path.exists(self.imageDir): os.makedirs(self.imageDir)
        self.channel = {1:'DIC',2:'DAPI',3:'FITC',4:'TRITC',5:'CY5'}

    def _getImageName(self,c):                           
        name = self.subDir + 'c' + str(c) + '.tif'
        return os.path.join(self.pathDir,self.subDir,name)

    def _openImages(self,fname,flag=0):
        orig = cv2.imread(fname,flag)
        corig = cv2.imread(fname,0)
        cimg = cv2.cvtColor(corig,cv2.COLOR_GRAY2BGR)
        return (orig, cimg)

    def _show(self,img):
        plt.imshow(img,cmap='gray')
        plt.show()
        return True

    def _saveImage(self,fname,img):
        plt.imshow(img,cmap='gray')
        plt.savefig(fname,bbox_inches='tight')
        return True

    def findHoughCircles(self,orig):
        def _blur(img):
            blur = cv2.medianBlur(orig,5)
            return blur
        
        blur = _blur(orig)
        circles = (
        cv2.HoughCircles(blur,cv2.cv.CV_HOUGH_GRADIENT,dp=1,minDist=30,param1=90,param2=60,minRadius=25,maxRadius=75)
        )
        if circles is None: circles = []
        return circles

    def _saveHoughCircles(self,circles,cimg,idx=1,putNbr=False):
        if circles is not None:
            print "Number of Circle found : %d ..."%len(circles[0])                                                                        
            nbr = 0
            #if len(circles[0])>20: return circMask, ringMask #Something is really wrong if there are more than 20 circles detected
            for c in circles[0][0:20]:
                nbr+=1
            # c[0],c[1] = (x,y) center; c[1] = r (radius)
                rad = int(c[2])
                cv2.circle(cimg,(c[0],c[1]),rad,( 71,144,48),thickness=1)               
                if putNbr: cv2.putText(cimg,str(nbr),(c[0],c[1]),cv2.FONT_HERSHEY_SIMPLEX,1,255)
        else: print "No Circles Found"
        name = os.path.join(self.imageDir,self.subDir+'_'+str(idx)+'_HoughCircles.png')
        self._saveImage(name,cimg)
        return True
    
    def drawRingMask(self,orig,x,y,rad,thickness=1):
        ringMask = np.zeros(orig.shape,dtype=np.uint16)
        cv2.circle(ringMask,(x,y),rad,1, thickness)
        return ringMask

    def applyMask(self,mask,img):
        intensity = np.sum(np.multiply(img,mask))
        area = np.count_nonzero(mask)
        return intensity,area,intensity/area
    

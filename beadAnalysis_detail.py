#! /home/jaggu/anaconda/bin/python2.7

"""
This programme does an extensive image processing of bead data to try
understanding a different way of extracting intensity information under the
fluorescent beads. Hough circles are used to identify the beads (approx) and
then different calculation/graphing methods used.
"""

from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import scipy
from PIL import Image
import os
import sys
from subprocess import call
import cv2
from commonFunctions import locate
import re
import cPickle as pickle
import collections
import shelve
import time

class BeadExpt(object):
    def __init__(self,subDir):
        self.subDir = subDir
        self.pathSubDir = os.path.join(pathDir,self.subDir)
        self.imageDir = os.path.join(self.pathSubDir,"Graphs")
        if not os.path.exists(self.imageDir): os.makedirs(self.imageDir)
        self.channel = {1:'DIC',2:'DAPI',3:'FITC',4:'TRITC',5:'CY5'}

    def _getImageName(self,c):                           
        name = self.subDir + 'c' + str(c) + '.tif'
        return os.path.join(pathDir,self.subDir,name)

    def _openImages(self,fname,flag=0):
        orig = cv2.imread(fname,flag)
        cimg = cv2.cvtColor(orig,cv2.COLOR_GRAY2BGR)
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
    
    def drawRingMask(self,orig,x,y,rad,thickness):
        ringMask = np.zeros(orig.shape,dtype=np.uint16)
        cv2.circle(ringMask,(x,y),rad,1, thickness)
        return ringMask

    def applyMask(self,mask,img):
        intensity = np.sum(np.multiply(img,mask))
        area = np.count_nonzero(mask)
        return intensity,area,intensity/area


class PlotRadial(object):
    """
    This gathers the information about mean intensity,area under different ring
    masks for each circle as a dictionary. The idea is to be able to plot it to
    see how to make sense of the different graphs.
    """
    def __init__(self,circInfo,imageDir,subDir):
        self.circInfo = circInfo
        self.imageDir = imageDir
        self.subDir = subDir
        #self.fig = plt.Figure(figsize=( 30,30),dpi=100,tight_layout=True)
        self.fig, self.axTupTup = plt.subplots(nrows=5,ncols=4)
        self.axTup = [ax for tup in self.axTupTup for ax in tup]
        font = {'size':6}
        matplotlib.rc('font',**font)
    
    def _plotRadial(self,ax,x,y):
        ax.plot(x,y)
        return ax
    
    def _plotLabels(self,ax,title,xlabel,ylabel):
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.text( 0.5,0.7,title,
                horizontalalignment='center',
                transform = ax.transAxes)
        return ax
    
    def _setTickIntervals(self,ax,x,y,xstep,ystep):
        ax.xaxis.set_ticks(np.arange(0,np.max(x),xstep))
        ax.yaxis.set_ticks(np.arange(0,np.amax(y),ystep))
        return ax

    def _setLimits(self,ax,xmax,ymax):
        ax.set_xlim([0,xmax])
        ax.set_ylim([0,ymax])
        return ax

    def _showFig(self):
        plt.tight_layout(pad=0.8,w_pad=0.8,h_pad=1.2)
        plt.show()

    def _saveFig(self,fname):
        plt.tight_layout(pad=0.8,w_pad=0.8,h_pad=2.0)
        plt.savefig(fname)
        plt.close(self.fig)


def getTextParams(typeIntensity):
    labelDict = {'meanIntensity':( 0.5,10000,
                                  'Normalized Radial Distance',
                                  'Mean Intensity'),
                 'minMaxNormIntensity':( 0.5,0.5,
                                        'Normalized Radial Distance',
                                        'Normalized Intensity ')
                }
                 
    return labelDict[typeIntensity]



def plotIntensityList(circles,ax,rrange,intList,nbr,typeIntensity):
    (xstep,ystep,xlabel,ylabel) = getTextParams(typeIntensity)
    xmax, ymax = 1.2, 1.1 * np.amax(intList)
    title = 'Circle # '+str(nbr+1)

    circles._plotRadial(ax,rrange,intList)
    circles._setTickIntervals(ax,rrange,intList,xstep,ystep)
    circles._setLimits(ax,xmax,ymax)
    circles._plotLabels(ax,title,xlabel,ylabel)

def plot_RadialIntensity(circInfo,imageDir,subDir):
    for typeIntensity in ['meanIntensity','minMaxNormIntensity']:
        print "Type of Intensity: %s"%(typeIntensity)
        circles = PlotRadial(circInfo,imageDir,subDir)
        for nbr in range(len(circInfo)):
            ax = circles.axTup[nbr]
            cInfo = circles.circInfo[nbr]
            rrange = [c[1] for c in cInfo]
            intList = [c[4] for c in cInfo]
            if typeIntensity == 'minMaxNormIntensity':
                minInt,maxInt = np.amin(intList),np.amax(intList)
                minMaxNormInt = [((val-minInt)/(maxInt-minInt)) for val in intList]
                plotIntensityList(circles,ax,rrange,minMaxNormInt,nbr, typeIntensity)           
            else:
                plotIntensityList(circles,ax,rrange,intList,nbr,typeIntensity)
        name = os.path.join(imageDir,subDir+'_'+typeIntensity+'.png')
        circles._saveFig(name)


def get_intensity(subDir):
    expt = BeadExpt(subDir)
    imageDir = expt.imageDir
    dicImage = expt._getImageName(1)
    orig,cimg = expt._openImages(dicImage)
    circles = expt.findHoughCircles(orig)
    circInfo = collections.defaultdict(list)
    if len(circles)>1: 
        expt._saveHoughCircles(circles,cimg,putNbr=True)
        dyeImage = expt._getImageName(4)
        dyeImg,cDyeImg = expt._openImages(dyeImage,flag=-1)#Dont use cDyeImg now
        expt._saveHoughCircles(circles,dyeImg,idx=4,putNbr=True)
        for i,c in enumerate(circles[0]):
            x,y,rad = map(int,c)
            thickness = int( 0.1*rad)
            #print "Circle # %d : %d,%d,%d"%(i,x,y,rad)
            for r in np.arange( 0.1,1.2,0.1):
                newrad = int(r*rad)
                ringMask = expt.drawRingMask(orig,x,y,newrad,thickness)
                intensity,area,meanIntensity = expt.applyMask(ringMask,dyeImg)
                circInfo[i].append([i, r, intensity,area,meanIntensity])
    return (circInfo,imageDir)





#sourceDir = "/project/marcotte/jagannath/projectfiles/EpiMicroscopy/rawFiles/2014-Aug"
#dateStamp = "2014-08-19"
#pathDir = os.path.join(sourceDir,dateStamp,"rawImages")
#subDir = "TentagelNH2_EDC_JSPR001_Before_1fps_flds014"

def makeRadialPlot(pathDir):
    print "Generating Masks for folders in %s ..."%(pathDir)
    for subDir in next(os.walk(pathDir))[1]:
        print "Processing %s ..."%(subDir)
        circInfo,imageDir = get_intensity(subDir)
        plot_RadialIntensity(circInfo,imageDir,subDir)

def test_case(subDir):
    print "Testing %s"%(os.path.join(pathDir,subDir))
    (circInfo,imageDir) = get_intensity(subDir)
    if circInfo>0: plot_RadialIntensity(circInfo,imageDir,subDir)

if __name__ == '__main__':
    monthIs = {'05':'May','06':'June','07':'July','08':'Aug','09':'Sept'}

    [ARG,dateStamp] = sys.argv[1:]
    epiDir = "/project/marcotte/jagannath/projectfiles/EpiMicroscopy"
    rawDataDir = os.path.join(epiDir,"rawFiles")
    sourceDir ="/project/marcotte/jagannath/projectfiles/EpiMicroscopy/rawFiles"
    month = monthIs[dateStamp.split('-')[1]]
    exptDir = os.path.join(sourceDir,"2014-"+month,dateStamp)
    pathDir = os.path.join(exptDir,"rawImages")
    t0 = time.clock()
    if ARG == 'RADIALPLOT':
        makeRadialPlot(pathDir)
    elif ARG == 'TEST':
        print "Testing"
        month = monthIs[dateStamp.split('-')[1]] 
        pathDir = os.path.join(sourceDir,"2014-"+month,dateStamp,"rawImages")
        subDir = "TentagelNH2_CR1NHS_200uM_Before_1s_flds005"
        test_case(subDir)
    else:
        raise SystemExit("Incorrect argument")
    t1 = time.clock()
    print ("Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
          )

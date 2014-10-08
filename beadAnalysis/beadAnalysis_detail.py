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
sys.path.append('/project/marcotte/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis')
from subprocess import call
import cv2
from commonFunctions import locate
import re
import cPickle as pickle
import collections
import shelve
import time
from beadExperiment import BeadExpt

class PlotRadial(object):
    """
    This gathers the information about mean intensity,area under different ring
    masks for each circle as a dictionary. The idea is to be able to plot it to
    see how to make sense of the different graphs.
    """
    def __init__(self,circInfo,imageDir,subDir):
        #plt.close() #Ensuring it is clean
        self.circInfo = circInfo
        self.imageDir = imageDir
        self.subDir = subDir
        #self.fig = plt.Figure(figsize=( 30,30),dpi=100,tight_layout=True)
        self.fig, self.axTupTup = plt.subplots(nrows=2,ncols=2)
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
    
    def _settickintervals(self,ax,x,y,xstep,ystep):
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


class PlotLineProfile(object):
    """
    This class plots a simple line profile, given the coordinates of the line.
    It will also overlay this profile on the actual image and above the
    coordinates of the line.
    """

    def __init__(self,subDir,imageDir,circles):
        #plt.close() # Closes older plots
        self.imageDir = imageDir
        self.subDir = subDir
        self.nbrCircles = len(circles)
        self.lineImageDir = os.path.join(imageDir,"LineProfiles")
        if not os.path.exists(self.lineImageDir): os.makedirs(self.lineImageDir)

    def defineFigure(self,nrows=13,ncols=1,figsize=( 10,10),fontsize=6,dpi=72):
        self.fig, self.axTupTup = plt.subplots(nrows=nrows,ncols=ncols,squeeze=False,figsize=figsize)
        self.axTup = [ax for tup in self.axTupTup for ax in tup]
        font = {'size':fontsize}
        matplotlib.rc('font',**font)
        return self.fig,self.axTup
    
    def _getImageName(self,c):                           
        name = self.subDir + 'c' + str(c) + '.tif'
        return os.path.join(pathDir,self.subDir,name)

    def _openImages(self,fname,flag=0):
        orig = cv2.imread(fname,flag)
        corig = cv2.imread(fname,0)#Always to 8bit color
        cimg = cv2.cvtColor(corig,cv2.COLOR_GRAY2BGR)
        return (orig, cimg)

    def _overLayLine(self,cimg,start,end,col=( 0,0,255)):
        cv2.line(cimg,start,end,col,thickness=2)
        return cimg

    def _plotLine(self,ax,x,y,col='#88CCEE'):
        ax.plot(x,y,linewidth=2.0,color=col)
        return ax

    def _plotLabels(self,ax,title,xlabel,ylabel,size='small'):
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.text( 0.5,1.2,title,
                horizontalalignment='center', transform = ax.transAxes,
                size='small',weight='bold')
        return ax

    def _showFig(self):
        plt.tight_layout(pad=0.8,w_pad=0.8,h_pad=1.2)
        plt.show()

    def _saveFig(self,fname):
        plt.tight_layout(pad=0.8,w_pad=0.8,h_pad=2.0)
        plt.savefig(fname)
        plt.close(self.fig)

    def getIntensityList(self,img,orientation,(x0,y0,r)):
        extfactor = 1.2
        extRad = r*extfactor
        def __checkEdges(x):
            if x > 512: x = 512
            elif x < 0: x = 0
            else: x
            return x

        if orientation is 'Xaxis':
            x1 = x0 - extRad 
            x2 = x0 + extRad
            x1,x2 = map(__checkEdges,[x1,x2])
            yx = [(int(y0),int(i)) for i in np.arange(x1,x2)]
            xaxis = np.arange(x1,x2)
            start = (int(x1),int(y0))
            end = (int(x2),int(y0))
        elif orientation is 'Yaxis':
            y1 = y0 - extRad
            y2 = y0 + extRad
            y1,y2 = map(__checkEdges,[y1,y2])
            yx = [(j,x0) for j in np.arange(y1,y2)]
            xaxis = np.arange(y1,y2)
            start = (int(x0),int(y1))
            end = (int(x0),int(y2))
        elif orientation is 'Slant':
            x1 = x0 - extRad*np.sin(np.pi/4)
            y1 = y0 - extRad*np.sin(np.pi/4)
            x2 = x0 + extRad*np.sin(np.pi/4)
            y2 = y0 + extRad*np.sin(np.pi/4)
            xN = [__checkEdges(x) for x in np.arange(x1,x2)]
            yN = [__checkEdges(y) for y in np.arange(y1,y2)]
            yx = zip(yN,xN)
            start = (int(x1),int(y1))
            end = (int(x2),int(y2))
            xaxis = [item[1] for item in yx]
        intList = [img[i] for i in yx]
        
        return (xaxis,intList,yx,start,end)

    def normalizeXY(self,xVal,yVal,binSize=2,method='mean'):
        xDict = collections.defaultdict(list)
        med = np.median(yVal)
        ynorm = [((i-med)/med)*100 for i in yVal]
        xnorm = [((i-np.amin(xVal))/(np.ptp(xVal))) for i in xVal]
        zipNorm = zip(xnorm,ynorm)
        for (x,y) in zipNorm:
            key = np.around(x,decimals=binSize)
            xDict[key].append(y)
        xyList = list()
        if method is 'mean':
            for xF,yList in xDict.items():
                yF = np.mean(yList)
                xyList.append((xF,yF))
        if method is 'max':
            for xF,yList in xDict.items():
                yF = np.amax(yList)
                xyList.append((xF,yF))

        xyList = sorted(xyList,key=lambda xy:xy[0])
        return zip(*xyList)



########### COMMON FUNCTIONS #################
def getCircles(subDir):
    expt = BeadExpt(subDir,pathDir)
    imageDir = expt.imageDir
    dyeImage = expt._getImageName(4)
    orig,cimg = expt._openImages(dyeImage)
    circles = expt.findHoughCircles(orig)
    orig = cv2.imread(dyeImage,-1)
    return circles,orig,cimg

######### LINE PLOTTING FUNCTION ######################

def sumProfiles(allXYDict,method='mean'):
    allXY = list()
    for xVal,yVals in sorted(allXYDict.items()):
        if method is 'mean': y = np.mean(yVals)
        if method is 'max': y = np.amax(yVals)
        allXY.append((xVal,y))
    x,y = zip(*allXY)
    return (x,y)


def plot_LineIntensity(circles,orig,cimg,imageDir,subDir,orientation='Xaxis'):
    # Note = For accesing the images, the coordinates are read as (y,x) [rows
    # and cols]

    colBlue = [( 51,34,136),'#332288']
    colGreen = [( 153,153,51),'#999933']
    colRed = [( 136,34,85),'#882255']
    colors = [colBlue,colGreen,colRed]
#    orientationList = ['Slant','Xaxis','Yaxis']
    allXYDict = collections.defaultdict(list)
    method = 'mean'
    
    nbrCircles = len(circles[0])
    line = PlotLineProfile(imageDir,subDir,circles[0])
    fig, axTup = line.defineFigure(nrows=nbrCircles,ncols=1,figsize=( 2,10),dpi=300)
     
    for i,cCoords in enumerate(circles[0]):
        x0,y0,r =  cCoords
        print "Circle #%d at (%f,%f); radius = %f"%(i,x0,y0,r)
        xVal,yVal,xyCoords,start,end =  line.getIntensityList(orig,orientation,(x0,y0,r))
        cimg = line._overLayLine(cimg,start,end)
        cv2.putText(cimg,str(i),(cCoords[0],int(cCoords[1]-10)),cv2.FONT_HERSHEY_SIMPLEX,1,color=colRed[0][0])
        xVal,yVal = line.normalizeXY(xVal,yVal)
        for j,key in enumerate(xVal):allXYDict[key].append(yVal[j])
        
        ax= axTup[i]
        ax = line._plotLine(ax,xVal,yVal,colBlue[1])
        ax = line._plotLabels(ax,'Circle #'+str(i),'Range (pixels)','')

        xstep,ystep = 0.25,50
        #ax.xaxis.set_ticks(np.arange(int(min(xVal)),int(max(xVal)),xstep))
        ax.xaxis.set_ticks(np.arange(np.amin(xVal),np.max(xVal),xstep))
        ax.yaxis.set_ticks(np.arange(0,np.amax(yVal),ystep))    
        
    fname = os.path.join(imageDir,subDir+'_'+orientation+'_LINE_NORMALIZE.png')
    print fname
    line._saveFig(fname)
    plt.imshow(cimg,cmap='gray')
    name = fname[:-4]+'_'+method+'_'+orientation+'.CIRCLES.png'
    plt.savefig(name)
    plt.close()
    
    x,y = sumProfiles(allXYDict,method='mean')

    fig,axTup = line.defineFigure(nrows=1,ncols=1,figsize=(5,5),dpi=300,fontsize=12)
    ax = axTup[0]
    ax = line._plotLine(ax,x,y,colRed[1])
    ax = line._plotLabels(ax,'Mean of Profiles','Normalized Distance','Normalized Median Intensity')
    ax = ax.text( 0.8,0.9,'n='+str(nbrCircles),fontsize=12,transform=ax.transAxes)   
    plt.savefig(fname[:-4]+'_Sum'+method+'_'+orientation+'.MeanOfProfiles.png')
    plt.close()
     
    return allXYDict



"""
    def __plotCircle(cChoice,cCoords):
        cCoords = circles[0][cChoice]
        line = PlotLineProfile(imageDir,subDir,nbrCircles)
        #axLeft,axRight = profile.axTup[cChoice],profile.axTup[cChoice+1]
        x0,y0,r =  cCoords
        print "Circle at (%f,%f); radius = %f"%(x0,y0,r)
        
        xVal,yVal,xyCoords,start,end =  line.getIntensityList(cimg,'Xaxis',(x0,y0,r))
#        profile._overLayLine(cdyeImg,start,end,col=colors[0][0])        
#        axLeft = axLeft.imshow(cdyeImg,cmap='gray')
#        profile._plotLine(axRight,xVal,yVal,col=colors[0][1])

        #for i,ax in enumerate(profile.axTup[:-1]):
         #   xVal,yVal,xyCoords,start,end = profile.getIntensityList(dyeImg,orientation[i],(x0,y0,r))
          #  profile._overLayLine(cdyeImg,start,end,col=colors[i][0])
           # profile._plotLine(ax,xVal,yVal,col=colors[i][1])
           # profile._plotLabels(ax,orientation[i] + ' Profile','Range (Pixels)', 'Raw Intensity')
        return (start,end,(xVal,yVal))
    # End of function
"""

######## RADIAL PLOTTING FUNCTION ###################

def plot_RadialIntensity(allCircles,circInfo,imageDir,subDir):
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
#        circles._setLimits(ax,xmax,ymax)
        circles._plotLabels(ax,title,xlabel,ylabel)


    for typeIntensity in ['meanIntensity','minMaxNormIntensity']:
        print "Type of Intensity: %s"%(typeIntensity)
        circles = PlotRadial(circInfo,imageDir,subDir)
#       for nbr in range(len(circInfo)):
        for nbr in range(4):
            nbr = 3
            print allCircles[0][nbr]
            ax = circles.axTup[nbr]
            cInfo = circles.circInfo[nbr]
            rrange = [c[1] for c in cInfo]
            intList = [c[4] for c in cInfo]
            sys.exit(1)
            if typeIntensity == 'minMaxNormIntensity':
                minInt,maxInt = np.amin(intList),np.amax(intList)
                minMaxNormInt = [((val-minInt)/(maxInt-minInt)) for val in intList]
                plotIntensityList(circles,ax,rrange,minMaxNormInt,nbr, typeIntensity)
            else:
                print rrange
                print intList
                plotIntensityList(circles,ax,rrange,intList,nbr,typeIntensity)
        #circles._showFig()                
        name = os.path.join(imageDir,subDir+'_bin_0.001_'+typeIntensity+'.png')
        circles._saveFig(name)
        #circles._showFig()
        sys.exit(1)


def get_intensity(subDir):
    circles,orig,cimg = getCircles(subDir)

    circInfo = collections.defaultdict(list)
    if len(circles)>0: 
        expt._saveHoughCircles(circles,cimg,putNbr=True)
        dyeImage = expt._getImageName(4)
        dyeImg,cDyeImg = expt._openImages(dyeImage,flag=-1)#Dont use cDyeImg now
        expt._saveHoughCircles(circles,dyeImg,idx=4,putNbr=True)
        
        for i,c in enumerate(circles[0]):
            extFactor = 1.2
            i = 3 
            c = circles[0][i]
            x0,y0,rad = map(int,c)
            extrad = int(rad*extFactor)
            # Estimating the baseline for the circle
            ringMask = expt.drawRingMask(orig,x0,y0,int(rad/2))
            intCIRC,areaCIRC,meanINTCIRC = expt.applyMask(ringMask,dyeImg)
            print "Mean Baseline Intensity : %f"%(meanINTCIRC)
            for r in range(1,extrad):
                ringMask = expt.drawRingMask(orig,x0,y0,r)
                cv2.circle(cDyeImg,(x0,y0),r,(255,0,0),1)
                aboveBGList = list()
                aboveRatio = list()
                for (i,j),val in np.ndenumerate(ringMask):
                    INT = (val * dyeImg.item((i,j)))
                    if INT>meanINTCIRC: 
                        aboveBGList.append(INT)
                        aboveRatio.append(INT/meanINTCIRC)
                if not aboveBGList: 
                    aboveBGList.append(meanINTCIRC)
                    aboveRatio.append(1)
                intensity,area,meanIntensity = expt.applyMask(ringMask,dyeImg)
                print r, meanIntensity, area, np.mean(aboveBGList),len(aboveBGList),np.mean(aboveRatio)
                circInfo[i].append([i, r, intensity,area,meanIntensity])
        return (circles,circInfo,imageDir)
    else:
        return [] , 0, imageDir



#sourceDir = "/project/marcotte/jagannath/projectfiles/EpiMicroscopy/rawFiles/2014-Aug"
#dateStamp = "2014-08-19"
#pathDir = os.path.join(sourceDir,dateStamp,"rawImages")
#subDir = "TentagelNH2_EDC_JSPR001_Before_1fps_flds014"

def makeRadialPlot(pathDir):
    print "Generating Masks for folders in %s ..."%(pathDir)
    for subDir in next(os.walk(pathDir))[1]:
        print "Processing %s ..."%(subDir)
        circles,circInfo,imageDir = get_intensity(subDir)
        plot_RadialIntensity(circInfo,imageDir,subDir)


def combineLinePlots(pathDir):
    print "Making Line Plots for folders in %s ..."%(pathDir)
    #fnameXYDict = collections.defaultdict(dict)
    globalXYDict = collections.defaultdict(list)
    orientation = 'Xaxis'
    method = 'mean'
    names = ["TentagelNH2_JSPR003D_"+i for i in ["Before","Mock","Edman","4hTFA"]]
    circleCount = 0
    for subDir in next(os.walk(pathDir))[1]:
        name = names[0] 
        if subDir.startswith(name):
            print "Processing %s ..."%(subDir)
            imageDir = os.path.join(pathDir,subDir,"LineProfiles")
            if not os.path.exists(imageDir): os.makedirs(imageDir)
            print "Image Dir : %s ...."%(imageDir)
            circles,orig,cimg = getCircles(subDir)

            if len(circles)>0:
                circleCount+=len(circles[0])                           
                allXYDict = plot_LineIntensity(circles,orig,cimg,imageDir,subDir,orientation=orientation)
                for k in allXYDict.keys():
                    globalXYDict[k].extend(allXYDict[k])
    
    print "Number of Circles Counted : %d"%(circleCount)
    plt.close()
    fig,ax = plt.subplots(figsize=( 5,5),dpi=300)
    x,y = sumProfiles(globalXYDict,method='mean')
    ax.plot(x,y,linewidth=3.0,color='#582A72')
    ax.set_xlabel('Normalized distance')
    ax.set_ylabel('Normalized Median Intensity')
    ax = ax.text( 0.8,0.9,'n='+str(circleCount),fontsize=12,transform=ax.transAxes)   
    fname = os.path.join(sourceImageDir,name+'_Sum_'+method+'_'+orientation+'.combinedNormalized.LineProfile.png')
    plt.savefig(fname)




def test_case(subDir):
    print "Testing %s"%(os.path.join(pathDir,subDir))
    allXYDict = collections.defaultdict(list)
    #(circles,circInfo,imageDir) = get_intensity(subDir)
    #if circInfo>0: plot_RadialIntensity(circles,circInfo,imageDir,subDir)
    circles,orig,cimg = getCircles(subDir)
    if len(circles)>0: plot_LineIntensity(circles,orig,cimg,imageDir,subDir)

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
    elif ARG == 'LINEPLOTS':
        sourceImageDir = os.path.join(exptDir,"LinePlotsAll")
        if not os.path.exists(sourceImageDir): os.makedirs(sourceImageDir)
        combineLinePlots(pathDir)
    elif ARG == 'TEST':
        print "Testing"
        month = monthIs[dateStamp.split('-')[1]] 
        pathDir = os.path.join(sourceDir,"2014-"+month,dateStamp,"rawImages")
        subDir = "TentagelNH2_JSPR003C_Before_5s_flds008"
        imageDir = os.path.join(pathDir,subDir,"Graphs")
        test_case(subDir)
    else:
        raise SystemExit("Incorrect argument")
    t1 = time.clock()
    print ("Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
          )

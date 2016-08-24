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
        #if not os.path.exists(self.lineImageDir): os.makedirs(self.lineImageDir)

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

    def getIntensityList(self,img,orientation,(x0,y0,r),theta=(np.pi/4)):
        extfactor = 1.2
        extRad = r*extfactor
        def __checkEdges(x):
            if x > 511: x = 511
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
            x1 = int(x0 - extRad*np.sin(theta))
            y1 = int(y0 - extRad*np.cos(theta))
            x2 = int(x0 + extRad*np.sin(theta))
            y2 = int(y0 + extRad*np.cos(theta))
            xList = np.arange(x1,x2,1)
            yList = np.arange(y1,y2,1)
            xN = [__checkEdges(x) for x in xList]
            yN = [__checkEdges(y) for y in yList]
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
def getCircles(subDir,maskSource='DIC'):
    expt = BeadExpt(subDir,pathDir)
    imageDir = expt.imageDir
    if maskSource is 'DIC':
        dyeImage = expt._getImageName(1)
    else:
        dyeImage = expt._getImageName(4)
    orig,cimg = expt._openImages(dyeImage)
    circles = expt.findHoughCircles(orig)
    orig,cimg = expt._openImages(expt._getImageName(4),-1)
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


def plot_LineIntensity(circles,orig,cimg,imageDir,subDir,method='mean',orientation='Xaxis',maskSource='DIC'):
    # Note = For accesing the images, the coordinates are read as (y,x) [rows
    # and cols]

    colBlue = [( 51,34,136),'#332288']
    colGreen = [( 153,153,51),'#999933']
    colRed = [( 136,34,85),'#882255']
    colors = [colBlue,colGreen,colRed]
    orientationList = ['Slant','Xaxis','Yaxis']
    allXYDict = collections.defaultdict(list)
    
    nbrCircles = len(circles[0])
    line = PlotLineProfile(imageDir,subDir,circles[0])
    fig, axTup = line.defineFigure(nrows=nbrCircles,ncols=1,figsize=( 2,10),dpi=300)
     
    for i,cCoords in enumerate(circles[0]):
        x0,y0,r =  cCoords
        print "Circle #%d at (%f,%f); radius = %f"%(i,x0,y0,r)
        xValList = list()
        yValList = list()
        for orientation in orientationList: 
            xVal,yVal,xyCoords,start,end =  line.getIntensityList(orig,orientation,(x0,y0,r))
            xVal,yVal = line.normalizeXY(xVal,yVal)
            for j,key in enumerate(xVal):allXYDict[key].append(yVal[j])       
            cimg = line._overLayLine(cimg,start,end)

        cv2.putText(cimg,str(i),(cCoords[0],int(cCoords[1]-10)),cv2.FONT_HERSHEY_SIMPLEX,1,color=colRed[0][0])
        
        ax= axTup[i]
        ax = line._plotLine(ax,xVal,yVal,colBlue[1])
        ax = line._plotLabels(ax,'Circle #'+str(i),'Range (pixels)','')

        xstep,ystep = 0.25,50
        ax.xaxis.set_ticks(np.arange(np.amin(xVal),np.amax(xVal),xstep))
        ax.yaxis.set_ticks(np.arange(0,np.amax(yVal),ystep))    
    
    fname = os.path.join(imageDir,subDir+'_'+maskSource+'_'+orientation+'_ALL_LINE_NORM.png')
    line._saveFig(fname)
    plt.imshow(cimg,cmap='gray')
    name = fname[:-4]+'_'+maskSource+'_'+method+'_'+orientation+'.CIRCLES.png'
    plt.savefig(name)
    plt.close()
    
    x,y = sumProfiles(allXYDict,method='mean')

    fig,axTup = line.defineFigure(nrows=1,ncols=1,figsize=(5,5),dpi=300,fontsize=12)
    ax = axTup[0]
    ax = line._plotLine(ax,x,y,colRed[1])
    ax = line._plotLabels(ax,'Mean of Profiles','Normalized Distance','Normalized Median Intensity')
    ax = ax.text( 0.8,0.9,'n='+str(nbrCircles),fontsize=12,transform=ax.transAxes)   
    plt.savefig(fname[:-4]+'_'+maskSource+'_'+method+'_'+orientation+'.All.MeanOfProfiles.png')
    
    plt.close()
     
    return allXYDict

def skip(name):
    flag = False
    skipList = ['zStep','Zstep']
    for i in skipList:
        if i in name:flag=True
    return flag

def makeSymmetry(x,y):
    # Half value is at 0.5
    ylist = list()
    xlen = int(len(x)/2)
    for i in xrange(xlen):
        yF = y[xlen-i]+y[xlen+i]
        ylist.append(yF)
    xlist = np.arange(0,1.,0.02)
    return (xlist,ylist)

def calcPosArea(x,y):
    yPos = [i if i>0 else 0 for i in y ]
    area = np.trapz(yPos,dx=0.02)
    return area 


def storeNegRef(db,fname):
    ofile = open(fname,'w')
    pickle.dump(db,ofile)
    ofile.close()

def makeLineProfilePlot(db,name,circleCount='NA'):
    plt.close()
    print "Number of Circles Counted : %d"%(circleCount)
    fig,ax = plt.subplots(figsize=( 5,5),dpi=300)
    x,y = sumProfiles(db,method='mean')
    ax.plot(x,y,linewidth=3.0,color='#582A72')
    ax.set_xlabel('Normalized Radial distance (with 10% outside bead region)')
    ax.set_ylabel('Normalized Intensity deviation from median')
    ax.set_ylim([-50,100])
    plt.axhline(ls='-.',linewidth=1.0,color="#113B52")
    ax = ax.text( 0.6,0.9,'n = '+str(circleCount)+' beads',fontsize=12,transform=ax.transAxes)   
    fname = os.path.join(sourceImageDir,name+'_all.LineProfile')
    fname1 = fname+'.svg'
    fname2 = fname+'.png'
    print "File Summed : %s "%(fname1)
    plt.plot(x,y,linewidth=3.0,color="#47036D")
    plt.savefig(fname1)
    plt.savefig(fname2)
    print "Done"
    #sys.exit(1) #Quiting before the negative subtraction


def combineLinePlots(pathDir):
    print "Making Line Plots for folders in %s ..."%(pathDir)
    #fnameXYDict = collections.defaultdict(dict)
    globalXYDict = collections.defaultdict(list)
    orientation = 'Xaxis'
    method = 'mean'
    maskSource = 'DIC'
    names = ["TentagelNH2_JSPR011_2_EDC_"+i for i in
             ["Before_MeOH","Before_TFA1h","deFmoc_MeOH","deFmoc_TFA1h","Mock1_MeOH","Mock1_TFA1h","Edman1_MeOH","Edman1_TFA1h","Edman2_MeOH","Edman2_TFA1h"]]
    circleCount = 0
    for subDir in next(os.walk(pathDir))[1]:
        name = names[0] 
        if subDir.startswith(name):
            if not skip(subDir):
                print "Processing %s ..."%(subDir)
                imageDir = os.path.join(pathDir,subDir,"LineProfiles")
                if not os.path.exists(imageDir): os.makedirs(imageDir)
                print "Image Dir : %s ...."%(imageDir)
                circles,orig,cimg = getCircles(subDir,maskSource)

                if len(circles)>0:
                    circleCount+=len(circles[0])                           
                    allXYDict = plot_LineIntensity(circles,orig,cimg,imageDir,subDir,orientation=orientation,maskSource=maskSource)
                    for k in allXYDict.keys():
                        globalXYDict[k].extend(allXYDict[k])

    fnamePkl = os.path.join(sourceImageDir,name+'_all.LineProfile.pkl')   
    #storeNegRef(globalXYDict,fnamePkl)
    
    makeLineProfilePlot(globalXYDict,name,circleCount)
    #sys.exit(1)
    plt.close()
    print "Number of Circles Counted : %d"%(circleCount)
    fig,ax = plt.subplots(figsize=( 5,5),dpi=300)
    x,y = sumProfiles(globalXYDict,method='mean')
    ax.plot(x,y,linewidth=3.0,color='#582A72')
    ax.set_xlabel('Normalized Radial distance (with 10% outside bead region)')
    ax.set_ylabel('Normalized Intensity deviation from median')
    ax.set_ylim([-50,100])
    plt.axhline(ls='-.',linewidth=1.0,color="#113B52")
    ax = ax.text( 0.6,0.9,'n = '+str(circleCount)+' beads',fontsize=12,transform=ax.transAxes)   
    fname = os.path.join(sourceImageDir,name+'_all.LineProfile')
    fname1 = fname+'.svg'
    fname2 = fname+'.png'
    print "File Summed : %s "%(fname1)
    plt.plot(x,y,linewidth=3.0,color="#47036D")
    plt.savefig(fname1)
    plt.savefig(fname2)
    print "Done"
    sys.exit(1) #Quiting before the negative subtraction


    
    negRefPlot = "/project/marcotte/jagannath/projectfiles/EpiMicroscopy/rawFiles/2014-Oct/2014-10-17/LinePlotsAll/TentagelNH2_RhodBCOOH_Zero_all.LineProfile.xyList.NEGREF.pkl" 
    negRefxyList = pickle.load(open(negRefPlot))
    xref,yref = zip(*negRefxyList)

    yDiff = [y[j]-yref[j] for j in xrange(len(y))]
    xF,yF = makeSymmetry(x,yDiff)
    plt.plot(xF,yF,linewidth=3.0,color="#690000")
    
    print fname
    plt.savefig(fname)

    area = calcPosArea(xF,yF)
    print "Positive Area under curve : %f"%(area)
    
    return True

def test_case(subDir):
    print "Testing %s"%(os.path.join(pathDir,subDir))
    allXYDict = collections.defaultdict(list)
    #(circles,circInfo,imageDir) = get_intensity(subDir)
    #if circInfo>0: plot_RadialIntensity(circles,circInfo,imageDir,subDir)
    circles,orig,cimg = getCircles(subDir)
    if len(circles)>0: plot_LineIntensity(circles,orig,cimg,imageDir,subDir)

if __name__ == '__main__':
    monthIs = {'05':'May','06':'June','07':'July','08':'Aug','09':'Sept','10':'Oct','12':'Dec'}

    [ARG,dateStamp] = sys.argv[1:]
    epiDir = "/project/marcotte/jagannath/projectfiles/EpiMicroscopy"
    rawDataDir = os.path.join(epiDir,"rawFiles")
    sourceDir ="/project/marcotte/jagannath/projectfiles/EpiMicroscopy/rawFiles"
    month = monthIs[dateStamp.split('-')[1]]
    year = dateStamp.split('-')[0]
    exptDir = os.path.join(sourceDir,"2015-"+month,dateStamp)
    pathDir = os.path.join(exptDir,"rawImages")

    t0 = time.clock()
    if ARG == 'LINEPLOTS':
        sourceImageDir = os.path.join(exptDir,"LinePlotsAll")
        if not os.path.exists(sourceImageDir): os.makedirs(sourceImageDir)
        combineLinePlots(pathDir)
    elif ARG == 'TEST':
        print "Testing"
        month = monthIs[dateStamp.split('-')[1]] 
        pathDir = os.path.join(sourceDir,"2014-"+month,dateStamp,"rawImages")
        subDir = "TentagelNH2_JSPR003C_Before_5s_flds008"
        #imageDir = os.path.join(pathDir,subDir,"Graphs")
        test_case(subDir)
    else:
        raise SystemExit("Incorrect argument")
    t1 = time.clock()
    print ("Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
          )


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



def makeRadialPlot(pathDir):
    print "Generating Masks for folders in %s ..."%(pathDir)
    for subDir in next(os.walk(pathDir))[1]:
        print "Processing %s ..."%(subDir)
        circles,circInfo,imageDir = get_intensity(subDir)
        plot_RadialIntensity(circInfo,imageDir,subDir)






























"""

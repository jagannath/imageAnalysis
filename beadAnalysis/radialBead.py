#! /home/jaggu/anaconda/bin/python2.7 

from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
import os
import sys
sys.path.append('/project/marcotte/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis')
from commonFunctions import locate,simpleShow
import cPickle as pickle
import collections
import time
from beadExperiment import BeadExpt
from radialProfileFigure import RadialProfile


class Bead:
    """
    This class performs operations on a single bead. It creates a lineplot for
    one bead across a few lines (xaxis,yaxis,slant45,slant135). It then
    normalizes the distance and creates a summed line plot. 
    """
    def __init__(self,circle,orig,cimg):
        self.circle = circle #Tuple of [x0,y0,r]
        self.x0, self.y0,self.r = map(int,circle)
        self.orig = orig
        self.cimg = cimg #Needed for drawing lines
        self.approxRad = 50 
        # This determines the granularity of the binning; Essential when dealing
        # with comparing beads of differing radius.
        self.lineExt = self.r * 1.2
        self.xyNorm_dict = collections.defaultdict(list)
        
    def bounds(self,x):
        if x>=511: x = 511
        if x<=0: x = 0
        return int(x)
    
    def makeLine(self,theta=0):
        angle = np.radians(theta)
        startx = self.x0 
        stopx = self.x0 + self.lineExt*np.cos(angle)
        starty = self.y0 
        stopy = self.y0 + self.lineExt*np.sin(angle)
        (x0,y0) = map(self.bounds,(startx,starty))
        (x1,y1) = map(self.bounds,(stopx,stopy))
        return [(x0,y0), (x1,y1)]
    
    def calcMeanVal_dict(self,vList_d):
        tmpDict = dict()
        for k, v in vList_d.items():
            tmpDict[k] = np.mean(v)
        return tmpDict

    def makeLineProfile(self):
        def _calcLineEqn((x0,y0),(x1,y1)):
            m = (y0-y1)/(x0-x1)
            c = (x1*y0-x0*y1)/(x1-x0)
            return m,c

        def _lineFn(m,c,xList):
            yList = m*xList + c
            return yList
        
        def _getPxlList(theta):
            [(x0,y0),(x1,y1)] = self.makeLine(theta)
            if x0 == x1:
                xList = np.linspace(x0,x1,self.lineExt)
                yList = np.linspace(y0,y1,self.lineExt)
            else:
                m,c = _calcLineEqn((x0,y0),(x1,y1))
                xList = np.linspace(x0,x1,self.lineExt)
                yList = _lineFn(m,c,xList)
        
            hw_list = zip(yList.astype(int),xList.astype(int))
        
            allPxlVals = list()
            for h,w in hw_list:
                pxlVal = self.orig.item((h,w))
                allPxlVals.append(pxlVal)
            
            return (allPxlVals)
        
        def _normRadList(pxL):
            binDict = collections.defaultdict(list)
            binsize = int(self.r)
            allX = list()
            for i,val in enumerate(pxL):
                normR = int(i*self.approxRad/binsize)
                allX.append(normR)
                binDict[normR].append(val)
            
            for k, v in binDict.items():
                self.xyNorm_dict[k].append(np.mean(v))
            return (self.xyNorm_dict)

        allLineProfiles = list()
        allNormProfiles = list()
        thetaList = np.arange( 0,361,15)
        for theta in thetaList:
            allPxlVals =  _getPxlList(theta)
            _normRadList(allPxlVals)
                    
        return self.xyNorm_dict

class BeadsImage:
    def __init__(self,subDir,pathDir,isNeg=False):
        self.subDir = subDir
        self.pathDir = pathDir
        self.isNeg = isNeg
        if not self.isNeg: 
            negF = '/project/marcotte/jagannath/projectfiles/EpiMicroscopy/rawFiles/2014-Oct/2014-10-17/rawImages/negShape.pkl'
            self.negShapeProfile = pickle.load(open(negF))
        self.beadsImg = BeadExpt(self.subDir,self.pathDir)
        self.dicFname = self.beadsImg._getImageName(c=1)
        self.flFname = self.beadsImg._getImageName(c=4) #Fluorescent channel

    def subDict(self,adict,negDict):
        tmpDict = dict()
        for k,v in adict.items():
            tmpDict[k] = adict[k] - negDict[k]
        return tmpDict

    def makeHoughCircles(self):
        self.dicFname = self.beadsImg._getImageName(c=1)
        (dic_orig, dic_cimg) = self.beadsImg._openImages(self.dicFname)
        circles = self.beadsImg.findHoughCircles(dic_orig)
        self.beadsImg.overLayHoughCircles_SVG(circles,dic_cimg)#SVG FILE MADE
        overLayFname, overCimg = self.beadsImg.overLayHoughCircles(circles,dic_cimg,idx=1,putNbr= True)
        self.beadsImg._saveImage(overLayFname+'.png', overCimg)
        return circles 
    
    def makeAllLinePlots(self,allCircles,minRad=50):
        (fl_orig,fl_cimg) = self.beadsImg._openImages(self.flFname,flag=-1)
        allNormDeviation_dict = collections.defaultdict(list) 
        for i,circle in enumerate(allCircles[0]):
            if circle[2]>minRad:
                oneBead = Bead(circle,fl_orig,fl_cimg)
                allNormLineProfiles_dict = oneBead.makeLineProfile()
                meanProfile = combineVal_dict(allNormLineProfiles_dict)
                (med,factor), devMedian, normDevMedian = calcDevMed_dict(meanProfile)
                shapeProfile = normDevMedian
                if not self.isNeg:
                    diffShape = self.subDict(shapeProfile,self.negShapeProfile)
                    reNorm = dict([(k,(v*factor)+med) for k, v in diffShape.items()])
                    allNormDeviation_dict = appendDict(allNormDeviation_dict,reNorm)
                    #allNormDeviation_dict = appendDict(allNormDeviation_dict,devMedian)
                else:
                    allNormDeviation_dict = appendDict(allNormDeviation_dict,shapeProfile)

        imageLineProfile = combineVal_dict(allNormDeviation_dict)
        return imageLineProfile

## END OF CLASS ##

def combineVal_dict(adict,type='mean'):
    """
    A dictionary is passed. The values are typically a list of list; Need is to
    combine them; default is mean; Keys have to be maintained
    """
    tmpDict = dict()
    for k,v in adict.items():
        tmpDict[k]=np.mean(v)
    return tmpDict

def calcDevMed_dict(adict):
    """
    Calculates the median of the list of the dictionary. Then returns a dict
    with the same keys, but val be val - med(vals) 
    and normalized (val - med(vals))/max(vals)-min(vals); Some kind of intensity
    normalization. 
    """
    tmpDict = dict()
    normDict = dict()
    med = np.median(adict.values())
    factor = np.amax(adict.values())
    for k,v in adict.items():
        tmpDict[k] = v - med
        normDict[k] = (v - med)/factor

    return (med,factor), tmpDict, normDict

def appendDict(adict,tmpDict):
    """
    Appends the values in the tmpDict to the actual Dict; If there is no
    corresponding key in adict, it adds it
    """
    for k,v in tmpDict.items():
        adict[k].append(v)
    return adict

def test_case(subDir):
    imgStack = BeadsImage(subDir,pathDir)
    allCircles = imgStack.makeHoughCircles()
    allNormDev = imgStack.makeAllLinePlots(allCircles,minRad=50)
    profile = RadialProfile(pathDir)
    fname = os.path.join(pathDir,subDir,'test')
    profile.drawProfile(allNormDev,fname,nbrCircles = len(allCircles[0]))
    #plt.plot(allNormDev.values())
    #plt.show()


def negControl(pathDir,subDir):
    typeDye = "RhodBCOOH_Before"
    allNegTypeProfile_dict = collections.defaultdict(list)
    for subDir in next(os.walk(pathDir))[1]:
        if typeDye in subDir and not 'zS' in subDir:
            negCntrl = BeadsImage(subDir,pathDir,isNeg=True)
            allCircles = negCntrl.makeHoughCircles()
            negShapeProfile = negCntrl.makeAllLinePlots(allCircles,minRad=50)
            allNegTypeProfile_dict = appendDict(allNegTypeProfile_dict,negShapeProfile)
    meanNegProfile = combineVal_dict(allNegTypeProfile_dict)
    fname = os.path.join(pathDir,'negShape.pkl')
    print "Neg Shape Profile .pkl file saved : %s"%(fname)
    pickle.dump(meanNegProfile,open(fname,'w'))
    plt.plot(meanNegProfile.keys(),meanNegProfile.values())
    plt.show()

def drawFigures(pathDir):
    pklDir = os.path.join(os.path.split(pathDir)[0],'pklFiles')
    radPlotDir = os.path.join(os.path.split(pathDir)[0],'allRadialPlots')
    if not os.path.exists(radPlotDir): os.makedirs(radPlotDir)
    pklFiles = [os.path.join(pklDir,f) for f in os.listdir(pklDir)]
    combineExptPlots(pathDir,pklDir,radPlotDir)
    for pklF in pklFiles:
        exptType = os.path.split(pklF)[1].split('.')[0]
        pf = RadialProfile(pathDir)
        fname = os.path.join(radPlotDir,exptType+'.radProfile')
        [nbrCircles,lineprofile] = pickle.load(open(pklF,'r'))
        pf.drawProfile(lineprofile,fname,nbrCircles=nbrCircles)
    

def combineExptPlots(pathDir,pklDir,radPlotDir):
    pep = 'TentagelNH2_JSPR011_2_EDC_'
    
    packedList = (
    [(os.path.join(pklDir,pep+'Before_MeOH_50ms.meanLineProfile.dict.pkl'),'-','blue'),
     (os.path.join(pklDir,pep+'deFmoc_MeOH_50ms.meanLineProfile.dict.pkl'),'--','red'),
     (os.path.join(pklDir,pep+'Mock1_MeOH_50ms.meanLineProfile.dict.pkl'),'-.','orange'),
     (os.path.join(pklDir,pep+'Edman1_MeOH_50ms.meanLineProfile.dict.pkl'),':','green'),
     (os.path.join(pklDir,pep+'Edman2_MeOH_50ms.meanLineProfile.dict.pkl'),'-','purple')
    ])
    fname = os.path.join(radPlotDir,'All_'+pep+'radialProfile')
    pf = RadialProfile(pathDir)
    pf.combinePlots(packedList,fname)
    


def combineProfiles(pathDir,exptType):
    allExptType_dict = collections.defaultdict(list)
    totalBeads = int()
    for subDir in next(os.walk(pathDir))[1]:
        if exptType in subDir and not 'zS' in subDir:
            expt = BeadsImage(subDir,pathDir)
            allCircles = expt.makeHoughCircles()
            totalBeads+=len(allCircles[0])
            print "Number of Circles : %d"%len(allCircles[0])
            allNormDev = expt.makeAllLinePlots(allCircles,minRad=10)
            allExptType_dict = appendDict(allExptType_dict,allNormDev)
    meanExptProfile = combineVal_dict(allExptType_dict)
    return meanExptProfile,totalBeads


def main(pathDir):
    # Makes PKL file : exptType.meanLineProfile.dict.pkl: [nbrBeads,dict]
    pklDir = os.path.join(os.path.split(pathDir)[0],'pklFiles')
    if not os.path.exists(pklDir): os.makedirs(pklDir)
    allFnames = [f for f in os.listdir(pathDir) if
                 os.path.isfile(os.path.join(pathDir,f))]
    allExptTypes = set(['_'.join(f.split('_')[:-1])  for f in allFnames])
    for exptType in allExptTypes:
        print "Expt type: %s"%(exptType)
        meanExptProfile,totalBeads = combineProfiles(pathDir,exptType)
        exptDetails = [totalBeads,meanExptProfile]
        fpkl = os.path.join(pklDir,exptType+'.meanLineProfile.dict.pkl')
        pickle.dump(exptDetails,open(fpkl,'w'))
    # Plotting the radial plots
    drawFigures(pathDir)
    


if __name__ == '__main__':
    monthIs = {'05':'May','06':'June','07':'July','08':'Aug','09':'Sept','10':'Oct','12':'Dec'}

    [ARG,dateStamp] = sys.argv[1:]
    epiDir = "/project/marcotte/jagannath/projectfiles/EpiMicroscopy"
    rawDataDir = os.path.join(epiDir,"rawFiles")
    sourceDir ="/project/marcotte/jagannath/projectfiles/EpiMicroscopy/rawFiles"
    month = monthIs[dateStamp.split('-')[1]]
    exptDir = os.path.join(sourceDir,"2014-"+month,dateStamp)
    pathDir = os.path.join(exptDir,"rawImages")

    t0 = time.clock()
    if ARG == 'RADIALPLOT':
        sourceImageDir = os.path.join(exptDir,"LinePlotsAll")
        if not os.path.exists(sourceImageDir): os.makedirs(sourceImageDir)
        main(pathDir)
    elif ARG == 'PLOT':
        drawFigures(pathDir)
    elif ARG == 'NEGCONTROL':
        print "Processing a Negative Control ..."
        pathDir = os.path.join(sourceDir,"2014-Oct","2014-10-17","rawImages")
        #pathDir = "project/marcotte/jagannath/projectfiles/EpiMicroscopy/rawFiles/2014-Oct/2014-10-17/rawImages"
        subDir = "TentagelNH2_RhodBCOOH_Before_5s_flds"
        negControl(pathDir,subDir)
    elif ARG == 'TEST':
        print "Testing"
        month = monthIs[dateStamp.split('-')[1]] 
        pathDir = os.path.join(sourceDir,"2014-"+month,dateStamp,"rawImages")
        subDir = "TentagelNH2_JSPR010_2_EDC_Before_MeOH_50ms_flds"
        #subDir = "TentagelNH2_JSPR010_2_EDC_Edman1_TFA1h_50ms_flds022"
        #subDir = "TentagelNH2_JSPR010_2_EDC_Edman2_TFA1h_50ms_flds036"
        #imageDir = os.path.join(pathDir,subDir,"Graphs")
        test_case(subDir)
    else:
        raise SystemExit("Incorrect argument")
    t1 = time.clock()
    print ("Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
          )

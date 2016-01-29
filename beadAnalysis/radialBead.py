#! /home/jaggu/anaconda/bin/python2.7 

from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
import os
import sys
sys.path.append('/project/current/project/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis')
from commonFunctions import locate,simpleShow
import cPickle as pickle
import collections
import time
from subprocess import call
import fnmatch
from beadExperiment import BeadExpt
from radialProfileFigure import RadialProfile

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and
    below supplied root directory.'''
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            allFiles.append(os.path.join(path,filename))
    return allFiles

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
        self.approxRad = 30  
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
    def __init__(self,subDir,pathDir,flChannel=4,isNeg=False):
        self.subDir = subDir
        self.pathDir = pathDir
        self.isNeg = isNeg
        if not self.isNeg: 
            negF = '/project/current/project/jagannath/projectfiles/EpiMicroscopy/rawFiles/2014-Oct/2014-10-17/rawImages/negShape.pkl'
            self.negShapeProfile = pickle.load(open(negF))
        self.beadsImg = BeadExpt(self.subDir,self.pathDir)
        self.dicFname = self.beadsImg._getImageName(c=1)
        self.flFname = self.beadsImg._getImageName(c=flChannel) #Fluorescent channel

    def subDict(self,adict,negDict):
        tmpDict = dict()
        for k,v in adict.items():
            tmpDict[k] = adict[k] - negDict[k]
        return tmpDict

    def makeHoughCircles(self):
        self.dicFname = self.beadsImg._getImageName(c=1)
        (dic_orig, dic_cimg) = self.beadsImg._openImages(self.dicFname)
        circles = self.beadsImg.findHoughCircles(dic_orig)
        if not circles is None:
            self.beadsImg.overLayHoughCircles_SVG(circles,dic_cimg)#SVG FILE MADE
            overLayFname, overCimg = self.beadsImg.overLayHoughCircles(circles,dic_cimg,idx=1,putNbr= True)
            self.beadsImg._saveImage(overLayFname+'.png', overCimg)
        return circles 
    
    def makeAllLinePlots(self,allCircles,minRad=30):
        (fl_orig,fl_cimg) = self.beadsImg._openImages(self.flFname,flag=-1)
        allNormDeviation_dict = collections.defaultdict(list) 
        for i,circle in enumerate(allCircles[0]):
            if circle[2]>minRad:
                oneBead = Bead(circle,fl_orig,fl_cimg)
                allNormLineProfiles_dict = oneBead.makeLineProfile()
                meanProfile,stdevProfile = combineVal_dict(allNormLineProfiles_dict)
                (med,factor), devMedian, normDevMedian = calcDevMed_dict(meanProfile)
                shapeProfile = normDevMedian
                if not self.isNeg:
                    diffShape = self.subDict(shapeProfile,self.negShapeProfile)
                    reNorm = dict([(k,(v*factor)+med) for k, v in diffShape.items()])
                    allNormDeviation_dict = appendDict(allNormDeviation_dict,reNorm)
                else:
                    allNormDeviation_dict = appendDict(allNormDeviation_dict,shapeProfile)

        imageLineProfile,stdevLineProfile = combineVal_dict(allNormDeviation_dict)
        return imageLineProfile,stdevLineProfile

class ExperimentType:
    """
    This class performs a number of functions useful for plotting and
    calculating efficiency of Edman. It also makes pklfiles for the general
    rawProfile and also processed profiles
    """
    def __init__(self,pathDir,exptType,flChannel=4):
        self.pathDir = pathDir
        self.exptType = exptType
        self.flChannel = flChannel
        self.pklDir = os.path.join(os.path.split(pathDir)[0],'pklFiles')
        if not os.path.exists(self.pklDir): os.makedirs(self.pklDir)

    def pickleIt(self,data,tag):
        fname = self.exptType + tag
        pklF = os.path.join(self.pklDir,fname)
        pickle.dump(data,open(pklF,'w'))
        return pklF 

    def combineProfiles(self):
        allExptType_dict = collections.defaultdict(list)
        totalBeads = int()
        for subDir in next(os.walk(self.pathDir))[1]:
            if self.exptType in subDir and not 'zS' in subDir:
                expt = BeadsImage(subDir,self.pathDir,flChannel=self.flChannel)
                allCircles = expt.makeHoughCircles()
                if not allCircles is None:
                    totalBeads+=len(allCircles[0])
                    print "Number of Circles : %d"%len(allCircles[0])
                    allNormDev,stdevLineProfile = expt.makeAllLinePlots(allCircles,minRad=30)

                    allExptType_dict = appendDict(allExptType_dict,allNormDev)
        meanExptProfile,stdevExptProfile = combineVal_dict(allExptType_dict)
        return meanExptProfile,stdevExptProfile,totalBeads

    def reduceLineProfile(self,lineProfile,stdevProfile):
        reduceDict = dict()
        reduceStdevDict = dict()
        min = np.amin(lineProfile.values())
        for k,v in lineProfile.items():
            _v = (v - min)
            reduceDict[k] = _v
            _vStdev = stdevProfile[k]#STDEV DOESNT CHANGE 
            reduceStdevDict[k] = _vStdev
        return reduceDict, reduceStdevDict


    def rescaleLineProfile(self,exptProfile,stdevProfile):
        rescaleDict = dict()
        rescaleStdevDict = dict()
        min, max = np.amin(exptProfile.values()),np.amax(exptProfile.values())
        newMax = max
        range = max - min
        for k,v in exptProfile.items():
            _v = (v - min)*(newMax/range)
            rescaleDict[k] = _v
            vStdev = stdevProfile[k]
            _vStdev = vStdev/range
            rescaleStdevDict[k] = _vStdev
        return rescaleDict,rescaleStdevDict

    def calcAUC(self,xyDict,stdevDict):
        x,y = np.array(xyDict.keys()),np.array(xyDict.values())
        area = np.trapz(y,x=x)

        yStdev = np.array(stdevDict.values())
        yPlus,yMinus = y+yStdev,y-yStdev
        areaPlus = np.trapz(yPlus,x=x)
        areaMinus = np.trapz(yMinus,x=x)
        areaStdev = np.std([areaPlus,area,areaMinus])

        return area,areaStdev



class RadialPlotFigure:
    """
    This class handles the plotting of the different function. It mainly relies
    on the pkl Files generated.
    """
    def __init__(self,pathDir,exptType):
        self.pklDir = os.path.join(os.path.split(pathDir)[0],'pklFiles')
        self.pathDir = pathDir
        self.exptType = exptType
        self.xBin = 50
        self.radPlotDir = os.path.join(os.path.split(pathDir)[0],'allRadialPlots')
        if not os.path.exists(self.radPlotDir): os.makedirs(self.radPlotDir)
        print "In RadialPlotFigure"

    def getYStdev(self,lineProfile,stdevProfile):
        yVal = np.array(lineProfile.values())
        yStd = np.array(stdevProfile.values())
        yPlus,yMinus = yVal+yStd,yVal-yStd
        return yPlus,yMinus

    def drawRadialProfile(self,pklTag,fnameTag):
        pklF = os.path.join(self.pklDir,self.exptType+pklTag)
        fname = os.path.join(self.radPlotDir,self.exptType+'.'+fnameTag)
        nbrBeads,lineProfile,stdevProfile = pickle.load(open(pklF))
        yPlus,yMinus = self.getYStdev(lineProfile,stdevProfile)
        rawXList = lineProfile.keys()
        xList = [x/self.xBin for x in rawXList]
        yList = lineProfile.values()
        pf = RadialProfile(self.pathDir)
        ymaxLimit = 1.2*np.amax(yPlus)
        ymaxredLimit = 30000
        if 'reduced' in pklTag: ylim = ( 0,ymaxredLimit)
        else: ylim = (0, ymaxLimit)
        pf.drawProfileWstdev(xList,yList,yPlus,yMinus,fname,nbrCircles=nbrBeads,ylim=ylim)
        return True


## END OF CLASS ##

def combineVal_dict(adict,type='mean'):
    """
    A dictionary is passed. The values are typically a list of list; Need is to
    combine them; default is mean; Keys have to be maintained
    """
    tmpDict = dict()
    tmpDictStdev = dict()
    for k,v in adict.items():
        tmpDict[k]=np.mean(v)
        tmpDictStdev[k] = np.std(v)
    return tmpDict,tmpDictStdev

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
    allNormDev = imgStack.makeAllLinePlots(allCircles,minRad=30)
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
            negShapeProfile = negCntrl.makeAllLinePlots(allCircles,minRad=30)
            allNegTypeProfile_dict = appendDict(allNegTypeProfile_dict,negShapeProfile)
    meanNegProfile,stdevNegProfile = combineVal_dict(allNegTypeProfile_dict)
    fname = os.path.join(pathDir,'negShape.pkl')
    print "Neg Shape Profile .pkl file saved : %s"%(fname)
    pickle.dump(meanNegProfile,open(fname,'w'))
    plt.plot(meanNegProfile.keys(),meanNegProfile.values())
    plt.show()

## MAIN FUNCTIONS

def drawFigures(pathDir):
    pklDir = os.path.join(os.path.split(pathDir)[0],'pklFiles')
    pklFiles = locate('*.pkl',pklDir)
    #combineExptPlots(pathDir,pklDir,radPlotDir)
    for pklF in pklFiles:
        print pklF
        exptType = os.path.split(pklF)[1].split('.')[0]
        
        plotType1 = RadialPlotFigure(pathDir,exptType)
        fnameTag1 = 'radProfile.stdev'
        plotType1.drawRadialProfile('.meanLineProfile.dict.pkl',fnameTag1)

        plotType2 = RadialPlotFigure(pathDir,exptType)
        fnameTag2 = 'radProfile.reduced.stdev'
        plotType2.drawRadialProfile('.meanLineProfile.reduced.stdev.dict.pkl',fnameTag2)
        
        
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
    
def main(pathDir,dateStamp,exptDIr):
    # Makes PKL file : exptType.meanLineProfile.dict.pkl: [nbrBeads,dict]
    channel_dict = {'405': 2, '488': 3,'555': 4,'647': 5}
    ofname = os.path.join(exptDir,dateStamp+'_RADIALSUMMARY.csv')
    ofile = open(ofname,'w')
    header = '\t'.join(['pickleFname','Experimental Step','MEAN AUC','STDEV AUC','\n'])
    ofile.write(header)

    pklDir = os.path.join(exptDir,'pklFiles')
    if not os.path.exists(pklDir): os.makedirs(pklDir)
    allFnames = [f for f in os.listdir(pathDir) if
                 os.path.isfile(os.path.join(pathDir,f))]

    allExptTypes = set(['_'.join(f.split('_')[:-1])  for f in allFnames])
    for exptType in allExptTypes:
            
        #channel_int = channel_dict[exptType.split('_')[-2]]
        channel_int = 2 
        print "Expt type: %s"%(exptType)
        exptStep = ExperimentType(pathDir,exptType,flChannel=channel_int)
        meanExptProfile,stdevExptProfile,totalBeads = exptStep.combineProfiles()
        
        exptDetails = [totalBeads,meanExptProfile,stdevExptProfile]
        exptStep.pickleIt(exptDetails,'.meanLineProfile.dict.pkl')
        
        reduceLine,reduceStdev = exptStep.reduceLineProfile(meanExptProfile,stdevExptProfile)
        pklPath = exptStep.pickleIt([totalBeads,reduceLine,reduceStdev],'.meanLineProfile.reduced.stdev.dict.pkl')
        area,areaStdev = exptStep.calcAUC(reduceLine,reduceStdev)
        
        info = [pklPath,exptType,str(area),str(areaStdev),'\n']
        ofile.write('\t'.join(info))
    return True        

def convertTifs(pathDir):
    allTifs = locate("*c*.tif",pathDir)
    for f1 in allTifs:
        f2 = f1[:-4]+'.jpg'
        cmd = "convert " + f1 + " " + f2
        call(cmd.split(),shell=False)

    return True





if __name__ == '__main__':
    monthIs = {'01':'Jan','02':'Feb','03':'Mar','04':'Apr','05':'May','06':'June','07':'July','08':'Aug','09':'Sept','10':'Oct','11':'Nov','12':'Dec'}

    [ARG,dateStamp] = sys.argv[1:]
    epiDir = "/project/current/project/jagannath/projectfiles/EpiMicroscopy"
    rawDataDir = os.path.join(epiDir,"rawFiles")
    sourceDir ="/project/current/project/jagannath/projectfiles/EpiMicroscopy/rawFiles"
    #sourceDir = "/project/current/project/jagannath/projectfiles/EpiMicroscopy/amber"
    month = monthIs[dateStamp.split('-')[1]]
    year = dateStamp.split('-')[0]
    exptDir = os.path.join(sourceDir,year+'-'+month,dateStamp)
    pathDir = os.path.join(exptDir,"rawImages")

    t0 = time.clock()
    if ARG == 'RADIALPLOT':
        #sourceImageDir = os.path.join(exptDir,"LinePlotsAll")
        #if not os.path.exists(sourceImageDir): os.makedirs(sourceImageDir)
        main(pathDir,dateStamp,exptDir)
    elif ARG == 'PLOT':
        drawFigures(pathDir)
    elif ARG == 'CONVERT':
        convertTifs(pathDir)
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
        #imageDir = os.path.join(pathDir,subDir,"Graphs")
        test_case(subDir)
    else:
        raise SystemExit("Incorrect argument")
    t1 = time.clock()
    print ("Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
          )

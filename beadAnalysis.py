#! /home/jaggu/anaconda/bin/python2.7

"""
This programme does the image analysis of beads (Tentagel beads) with dyes
immobilized on them. The image is a tiff stack containing the image of the beads
in different channels. I had earlier done the mask generation etc using FIJI.
Now I am trying to implement the hough circle algorithm in opencv2 to be able to
better the image analysis and processing. 
"""

from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
import tiffcapture as tc #This was a package I had to install 
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


class BeadStacks(object):
    """
    This class performs the operation of running through the Epimicroscopy bead
    experiment and generate a mask on the '*c1.tif' image (DIC channel). The
    method used is the 'Hough circle detection'. You can see the documentation
    of cv2.HoughCircles(...) on the website -
    http://docs.opencv.org/modules/imgproc/doc/feature_detection.html
    The method then runs through the rest of the c* images with the following
    order - {1:DIC,2:DAPI,3:FITC,4:TRITC,5:CY5}. The mask generated on the DIC
    channel is translated across the other fields. It then makes a little csv
    file in the same directory with the information of the subDirectory,
    Mean Intensity density and STDEV for each of the channels
    """
    def __init__(self,subDir):
        self.subDir = subDir
        self.pathSubDir = os.path.join(pathDir,self.subDir)
        self.channel = {1:'DIC',2:'DAPI',3:'FITC',4:'TRITC',5:'CY5'}

    def _getImageName(self,c):
        name = self.subDir + 'c' + str(c) + '.tif'
        return os.path.join(pathDir,self.subDir,name)
    
    def _saveImage(self,fname,img,show=False):
        plt.imshow(img,cmap='gray')
        if show: plt.show()
        else: plt.savefig(fname,bbox_inches='tight')
        return None

    def _getIntensityUnderMask(self,f,maskArray):
        img = cv2.imread(f,cv2.IMREAD_UNCHANGED)#CRUCIAL
        intensity = np.sum(np.multiply(img,maskArray))
        return intensity

    def computeBeadIntensity(self):
        fnameDIC = self._getImageName(c=1)
        self.circularMask,self.ringMask = self.applyHoughCircles(fnameDIC)
        circularAreaMask,ringAreaMask = map(np.count_nonzero,(self.circularMask,self.ringMask))
        circData = [self.applyMask(c,self.circularMask,circularAreaMask) for c in xrange( 1,6)]
        ringData = [self.applyMask(c,self.ringMask,ringAreaMask) for c in xrange( 1,6)]
        self.writeData(circData,ringData)

    def applyHoughCircles(self,fname):
        def _saveMaskImages():
            circMaskImgName = os.path.join(self.pathSubDir,self.subDir+'_CircularMask.png')
            ringMaskImgName = os.path.join(self.pathSubDir,self.subDir+'_RingMask.png')
            self._saveImage(circMaskImgName,circMask)
            self._saveImage(ringMaskImgName,cimg)
            return None
        
        print fname
        orig = cv2.imread(fname,0)
        img = cv2.medianBlur(orig,5)
        circMask = np.zeros(orig.shape,dtype=np.uint16)
        ringMask = np.zeros(orig.shape,dtype=np.uint16)
        cimg = cv2.cvtColor(orig,cv2.COLOR_GRAY2BGR)
        # These are the parameters that appear to work. So edit with caution.
        #circles = cv2.HoughCircles(img,cv2.cv.CV_HOUGH_GRADIENT,dp=1,
        #                       minDist=30,param1=90,param2=60,minRadius=30,maxRadius=150)
        

        circles = (
        cv2.HoughCircles(img,cv2.cv.CV_HOUGH_GRADIENT,dp=1,minDist=30,param1=90,param2=60,minRadius=20,maxRadius=75)
        )

        if circles is not None:
            print "Number of Circle found : %d ..."%len(circles[0])
            #if len(circles[0])>20: return circMask, ringMask #Something is really wrong if there are more than 20 circles detected
            for c in circles[0][0:20]:
            # c[0],c[1] = (x,y) center; c[1] = r (radius)
                cv2.circle(circMask,(c[0],c[1]),c[2],1,-1)
                cv2.circle(ringMask,(c[0],c[1]),c[2],1, 3)
                cv2.circle(cimg,(c[0],c[1]),c[2],( 71,144,48),3)
            _saveMaskImages()
        return circMask, ringMask
            
    def applyMask(self,c,maskArray,area):
        f = self._getImageName(c)
        intMaskArea = self._getIntensityUnderMask(f,maskArray)
        if area == 0: intensityDensity = 0
        else: intensityDensity = intMaskArea/area
        return (c,intensityDensity)

    def writeData(self,circData,ringData):
        oname = os.path.join(self.pathSubDir,self.subDir+'_EXPT_RESULTS.csv')
        ofile = open(oname,'w')
        ofile.write("CIRCULAR MASK MEAN INTENSITY \n")
        ofile.write("\t".join(["PATH","FNAME","DIC","DAPI","FITC","TRITC","CY5","\n"]))
        ofile.write("\t".join([self.pathSubDir,self.subDir,str(circData[0][1]),str(circData[1][1]),str(circData[2][1]),str(circData[3][1]),str(circData[4][1]),"\n"]))
        ofile.write("RING MASK MEAN INTENSITY \n")
        ofile.write("\t".join([self.pathSubDir,self.subDir,str(ringData[0][1]),str(ringData[1][1]),str(ringData[2][1]),str(ringData[3][1]),str(ringData[4][1]),"\n"]))
        ofile.close()


class BeadExpt(object):
    """
    This essentially looks for the csv file in the subdirectory on this date and
    makes two summary files. In the first summary file it is just the raw data
    compiled together after doing something like a sorting of the file. In a
    different summary, it combines the sorted file (after trimming the 'flds'
    word in the file name. I would then calculate the mean and stdev of the
    intensity for all the files with this name. 
    """

    def __init__(self,exptDir):
        self.exptDir = exptDir
        self.dateStamp = os.path.split(self.exptDir)[1]
        self.allcsv = locate("*_EXPT_RESULTS.csv",self.exptDir)
        self.allcsv.sort()
        self.destDir = exptDir
        self.outFullFile, self.outShortFile = self.makeOutFiles()
        self.nameCircInt_dict = collections.defaultdict(list)
        self.nameRingInt_dict = collections.defaultdict(list)

    def makeOutFiles(self):
        outFull = os.path.join(self.destDir,self.dateStamp+"_FULLSUMMARY.csv")
        outShort = os.path.join(self.destDir,self.dateStamp+"_SHORTSUMMARY.csv")
        return open(outFull,'w'),open(outShort,'w')

    def getData(self,fname):
        with open(fname) as f:
            lines = f.readlines()
        circData = lines[2].split('\t')[:-1]
        ringData = lines[4].split('\t')[:-1]
        return circData,ringData

    def writeFullSummary(self):
        
        headerLabel = "\t".join(["PATH","FNAME","DIC","DAPI","FITC","TRITC","CY5","\n"])    
        self.outFullFile.write(headerLabel)

        def _writeList(lines):
            for l in lines: 
                self.outFullFile.write("\t".join(l)+"\n")
            self.outFullFile.write("\n")
            
        #ofile.write("\t".join(["PATH","FNAME","DIC","DAPI","FITC","TRITC","CY5","\n"]))
        self.outFullFile.write("CIRCULAR MASK MEAN INTENSITY \n")       
        for k, v in sorted(self.nameCircInt_dict.items()):  _writeList(v)
        self.outFullFile.write("\n RING MASK MEAN INTENSITY \n")
        for k, v in sorted(self.nameRingInt_dict.items()): _writeList(v)

    def fullSummary(self):
        for csv in self.allcsv:
            fpath = csv
            fullname = "_".join(fpath.split('_')[:-3])
            name = os.path.split(fullname)[1]
            circData,ringData = self.getData(fpath)
            self.nameCircInt_dict[name].append(circData)
            self.nameRingInt_dict[name].append(ringData)
        self.writeFullSummary()


def test_summary(exptDir):
    test = BeadExpt(exptDir)
    print test.destDir
    f = test.allcsv[0]
    test.fullSummary()
    test.outFullFile.close()
    test.outShortFile.close()

def test_case(subDir):
    test = BeadStacks(subDir)
    print "Processing %s ..."%(subDir)
    print test.pathSubDir
    test.computeBeadIntensity()


def generateMask(pathDir):
    print "Generating Masks for folders in %s ..."%(pathDir)
    for subDir in next(os.walk(pathDir))[1]:
        image = BeadStacks(subDir)
        print "Processing %s ..."%(subDir)
        image.computeBeadIntensity()


if __name__ == '__main__':
    monthIs = {'05':'May','06':'Jun','07':'July','08':'Aug','09':'Sept'}

    #sourceDir ="/project/marcotte/jagannath/projectfiles/EpiMicroscopy/rawFiles/2014-July"
    #dateStamp = "2014-07-31b"
    #pathDir = os.path.join(sourceDir,dateStamp,"rawImages")
    #subDir = "TentagelCOOH_Pept_NHS_Washed_100ms_flds011"
    #test_case()

    [ARG,dateStamp] = sys.argv[1:]
    epiDir = "/project/marcotte/jagannath/projectfiles/EpiMicroscopy"
    rawDataDir = os.path.join(epiDir,"rawFiles")
    sourceDir ="/project/marcotte/jagannath/projectfiles/EpiMicroscopy/rawFiles"
    month = monthIs[dateStamp.split('-')[1]]
    exptDir = os.path.join(sourceDir,"2014-"+month,dateStamp)
    pathDir = os.path.join(exptDir,"rawImages")
    t0 = time.clock()
    if ARG == 'MASK':
        generateMask(pathDir)
    elif ARG == 'ANALYSE':
        print "Analysing"
        test_summary(exptDir)
    elif ARG == 'TEST':
        print "Testing"
        month = monthIs[dateStamp.split('-')[1]] 
        pathDir = os.path.join(sourceDir,"2014-"+month,dateStamp,"rawImages")
        subDir = "TentagelNH2_BODIPYFLMal_10ms_Before_flds002"
        test_case(subDir)
    else:
        raise SystemExit("Incorrect argument")
    t1 = time.clock()
    print ("Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
          )

















#test_mask(fname)


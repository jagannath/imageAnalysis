#! /home/jaggu/anaconda/bin/python2.7

# AUTHOR : Jagannath S; 
# Date Created : 2014-08-24
# Date Last Modified : 2014-08-24

"""
This script tries to complete and summarize the analysis of photolithography
experiments. There are two important parts - (a) Generating random images,
scaled etc and (b) Summarizing the results of the image processing on the
patterned slide.
NOTE: SOURCE and DESTINATION DIR is in project2 and ensure it can be opened
"""
import os.path
import os
import fnmatch
import csv
import collections
import time
import re
import shutil
import sys
from random import sample, choice
from subprocess import call
import numpy as np

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below supplied root directory.'''
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            allFiles.append(os.path.join(path, filename))
    return allFiles

class ConvertImages(object):
    """
    This class, handles all the images in the directory and converts and scales
    them appropriately. At the moment it handles the 'flds' and 'stitch' images
    """
    def __init__(self,imgDir):
        self.imgDir = imgDir
        self.destImgDir = os.path.join(destDir,'images')
        if not os.path.exists(self.destImgDir): os.makedirs(self.destImgDir)

    def removePngs(self,pattern,sourceDir):
        # To be used sparingly. If other pngs were generated, they disapper 
        pattern = '*.png'
        allPngs = locate(pattern, sourceDir)
        for pngFile in allPngs:
            try:
                os.remove(pngFile)
            except OSError:
                pass
        print "All Png files removed"

    def convertBatchImages(self,sourceDir):
        pattern = '*.tif'
        # convert -contrast-stretch 0.15x0.05% [input] [output]
        allTiffs = locate(pattern, sourceDir)
        for inputTif in allTiffs:
            outputTif1 = inputTif + '%d_scaled.png'
            outputTif2 = inputTif + '%d_noScale.png'
            cmd1 = "convert -contrast-stretch 0.015x0.05% " + inputTif + " " + outputTif1
            cmd2 = "convert " + inputTif + " " + outputTif2
            call(cmd1.split(),shell=False)
            call(cmd2.split(),shell=False)

    def convertOneRandomImage(self):
        tifFiles = locate("*.tif",self.imgDir)
        randomFile = choice(tifFiles).split('/')[-1]
        fpath = os.path.join(self.imgDir,randomFile)
        print "Flds: %s" %(fpath)
        outputPNG1 = os.path.join(self.destImgDir,randomFile + '_scaled.png')
        outputPNG2 = os.path.join(self.destImgDir,randomFile + '_noScale.png')
        # This is unnecessary sometimes if the file comes down as an RGB with
        # Red color instead of gray scale :(
        temp = os.path.join(self.imgDir,randomFile+'_temp.png')
        cmd0 = "convert " + fpath + " -colorspace gray " + temp
        cmd1 = "convert -contrast-stretch 0.015x0.05% " + temp + " " + outputPNG1
        cmd2 = "convert " + temp + " " + outputPNG2
        call(cmd0.split(),shell=False)
        call(cmd1.split(),shell=False)
        call(cmd2.split(),shell=False)

    def rescaleStitchImage(self):
        print os.listdir(self.imgDir)
        for tmp in os.listdir(self.imgDir):
            if tmp.endswith('.tif'):f = tmp

        fpath = os.path.join(self.imgDir,f)
        print "Stitch : %s"%(fpath)
        scale = str(5)
        outputSTITCH1 = os.path.join(self.destImgDir,f + '_noScale.png')
        outputSTITCH2 = os.path.join(self.destImgDir,f + '_scaled.png')
        # Same annoying format change to grayscale
        temp = os.path.join(self.imgDir,f+'_temp.png')
        cmd0 = "convert " + fpath + " -colorspace gray " + temp
        cmd1 = "convert " + temp + " -resize " + scale + "% " + outputSTITCH1
        cmd2 = "convert -contrast-stretch 0.015x0.05% " + outputSTITCH1 + " " + outputSTITCH2

        call(cmd0.split(),shell=False)
        call(cmd1.split(),shell=False)
        call(cmd2.split(),shell=False)


class ProcessedImages(object):
    """
    This class performs calculations on the images processed by FIJI. They work
    on the csv file which contains the information about the patterns and the
    intensity. The simplest calculation would be intensity density
    (Counts/pixel) under the patterns. 
    The output is a csv file in the destination directory. It will have the mean
    and stdev calculations for each of the csv file in the directory.
    NOTE: The calculation is done only for the flds images
    """

    def __init__(self,imgDir):
        self.imgDir = imgDir
    
    def calcIntensity(self,csv):
        ifile = open(csv)
        lines = ifile.readlines()[1:]
        areaList, intList = list(), list()
        for line in lines:
            area, intensity = float(line.split(',')[1]),float(line.split(',')[2])
            areaList.append(area)
            intList.append(intensity)
        return (np.mean(intList), np.std(intList))
    
    def calcBGIntensity(self,csv):
        """
        ["IMAGE FILE","AREA OF BACKGROUND", "MEAN BACKGROUND INTENSITY DENSITY"
        """
        fnameBG = csv[:-12]+'_BG_RESULTS.csv'
        try:
            ifileBG = open(fnameBG)
        except IOError:
            return ( 0, 0)
        line = ifileBG.readlines()[1]
        areaBG, bgIntensity = line.split(',')[1],line.split(',')[2]
        return (bgIntensity,areaBG)


    def analyseCSV(self):
        csvFiles = locate("*tif_RESULTS.csv",self.imgDir)
        for csv in csvFiles:
            (mean,stdev) = self.calcIntensity(csv)
            (bgMean,areaBG) = self.calcBGIntensity(csv)
            csvName = csv.split('/')[-1][:-12] 
            ofile.write("\t".join([csv, csvName, str(mean),str(stdev),str(bgMean), str(areaBG),"\n"]))
            print mean, stdev, csvName
            

dateStamp = "2014-09-03"
#sourceDir ="/project2/marcotte/boulgakov/microscope2/jagannath/rawFiles/2014-Sept/" +dateStamp 
sourceDir ="/project2/marcotte/boulgakov/microscope2/jagannath/rawFiles/2014-Sept/" +dateStamp + "/export" 

destDir = "/project2/marcotte/jaggu/dataAnalysis/microscope2/2014-Sept/" + dateStamp
if not os.path.exists(destDir): os.makedirs(destDir)
resultCSV = os.path.join(destDir,dateStamp + "_SUMMARY_RESULT.csv")
ofile = open(resultCSV,'w')
ofile.write("\t".join(["CSV FILE", "NAME", "MEAN OF MEAN INTENSITY/PIXEL","STDEV OF MEAN INTENSITY/PIXEL","INTENSITY DENSITY OF BACKGROUND","AREA OF BACKGROUND","\n"]))


def test():
    DIRNAMES = 1
    for dirname in os.walk(sourceDir).next()[DIRNAMES]:
        imgDir = os.path.join(sourceDir,dirname)
        imgs = ConvertImages(imgDir)
        csvFiles = ProcessedImages(imgDir)
        if 'flds' in imgDir:
            print "Processing Flds "
            imgs.convertOneRandomImage()
            ofile.write("\n")
            csvFiles.analyseCSV()
        if 'stitch' in imgDir:
            print "Processing Stitch"
            imgs.rescaleStitchImage()

test()
#sourceDir = "/project/marcotte/jagannath/projectfiles/EpiMicroscopy/2014-June/"+ dateStamp + "/rawImages"
#print sourceDir
#convertBatchImages(sourceDir)

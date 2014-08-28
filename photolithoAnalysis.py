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

    def convertOneRandomImage(self,destDir):
        tifFiles = locate("*.tif",self.imgDir)
        randomFile = choice(tifFiles).split('/')[-1]
        fpath = os.path.join(self.imgDir,randomFile)
        print "Flds: %s" %(fpath)
        outputPNG1 = destDir + '/images/' + randomFile + '_scaled.png'
        outputPNG2 = destDir + '/images/' + randomFile + '_noScale.png'
        cmd1 = "convert -contrast-stretch 0.015x0.05% " + fpath + " " + outputPNG1
        cmd2 = "convert " + fpath + " " + outputPNG2
        call(cmd1.split(),shell=False)
        call(cmd2.split(),shell=False)

    def rescaleStitchImage(self,destDir):
        [f] = os.listdir(self.imgDir)
        fpath = os.path.join(self.imgDir,f)
        print "Stitch : %s"%(fpath)
        scale = str(5)
        outputSTITCH = destDir + '/images/' + f + '.png'
        cmd = "convert " + fpath + " -resize " + scale + "% " + outputSTITCH
        call(cmd.split(),shell=False)


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
    
    def analyseCSV(self):
        csvFiles = locate("*tif_RESULTS.csv",self.imgDir)
        for csv in csvFiles:
            (mean,stdev) = self.calcIntensity(csv)
            csvName = csv.split('/')[-1][:-12] 
            ofile.write("\t".join([csv, csvName, str(mean), str(stdev),"\n"]))
            print mean, stdev, csvName
            


dateStamp = "2014-08-26"
sourceDir = "/project2/marcotte/boulgakov/microscope2/jagannath/rawFiles/2014-Aug/" + dateStamp 
destDir = "/project2/marcotte/jaggu/dataAnalysis/microscope2/2014-Aug/" + dateStamp
resultCSV = os.path.join(destDir,dateStamp + "_SUMMARY_RESULT.csv")
ofile = open(resultCSV,'w')
ofile.write("\t".join(["CSV FILE", "NAME", "MEAN OF MEAN INTENSITY/PIXEL","STDEV OF MEAN INTENSITY/PIXEL","\n"]))



def test():
    DIRNAMES = 1
    for dirname in os.walk(sourceDir).next()[DIRNAMES]:
        imgDir = os.path.join(sourceDir,dirname)
        imgs = ConvertImages(imgDir)
        csvFiles = ProcessedImages(imgDir)
        if 'flds' in imgDir:
            print "Processing Flds "
            imgs.convertOneRandomImage(destDir)
            ofile.write("\n")
            csvFiles.analyseCSV()
        if 'stitch' in imgDir:
            print "Processing Stitch"
            imgs.rescaleStitchImage(destDir)

test()
#sourceDir = "/project/marcotte/jagannath/projectfiles/EpiMicroscopy/2014-June/"+ dateStamp + "/rawImages"
#print sourceDir
#convertBatchImages(sourceDir)

#! /home/jaggu/anaconda/bin/python2.7

# AUTHOR : Jagannath S; 
# Date Created : 2014-05-14
# Date Last Modified : 2014-05-14

"""
The simple programme aims to convert files in a directory to another format
using appropriate convert options
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

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below supplied root directory.'''
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            allFiles.append(os.path.join(path, filename))
    return allFiles

def removePngs(dateStamp):
    # To be used sparingly. If other pngs were generated, they disapper 
    pattern = '*.png'
    allPngs = locate(pattern, sourceDir)
    for pngFile in allPngs:
        try:
            os.remove(pngFile)
        except OSError:
            pass
    print "All Png files removed"

def convertBatchImages(dateStamp):
    pattern = '*.tif'
    # convert -contrast-stretch 0.15x0.05% [input] [output]
    allTiffs = locate(pattern, sourceDir)
    for inputTif in allTiffs:
        outputTif1 = inputTif + '%d_scaled.png'
        outputTif2 = inputTif + '%d_noScale.png'
        cmd1 = "convert -contrast-stretch 0.15x0.05% " + inputTif + " " + outputTif1
        cmd2 = "convert " + inputTif + " " + outputTif2
        call(cmd1.split(),shell=False)
        call(cmd2.split(),shell=False)

def convertOneRandomImage(imgDir,destDir):
    randomFile = choice(os.listdir(imgDir))
    if randomFile.endswith('pkl'): randomFile = choice(os.listdir(imgDir))
    fpath = os.path.join(imgDir,randomFile)
    outputPNG1 = os.path.join(destDir,randomFile + '_scaled.png')
    outputPNG2 = os.path.join(destDir,randomFile + '_noScale.png')
    cmd1 = "convert -contrast-stretch 0.015x0.05% " + fpath + " " + outputPNG1
    cmd2 = "convert " + fpath + " " + outputPNG2
    call(cmd1.split(),shell=False)
    call(cmd2.split(),shell=False)

def rescaleStitchImage(imgDir,destDir):
    [f] = os.listdir(imgDir)
    fpath = os.path.join(imgDir,f)
    scale = str(5)
    outputSTITCH = destDir + '/images/' + f + '.png'
    cmd = "convert " + fpath + " -resize " + scale + "% " + outputSTITCH
    call(cmd.split(),shell=False)
    print outputSTITCH 


dateStamp = "2014-11-17"
sourceDir = os.path.join("/project2/marcotte/boulgakov/microscope/2014-Nov",dateStamp)
destDir = os.path.join("/project2/marcotte/jaggu/dataAnalysis/microscope1/2014-Nov",dateStamp,"images")
if not os.path.exists(destDir): os.makedirs(destDir)
DIRNAMES = 1
for dirname in os.walk(sourceDir).next()[DIRNAMES]:
    print sourceDir
    imgDir = os.path.join(sourceDir,dirname)
    if 'flds' in imgDir:
        convertOneRandomImage(imgDir,destDir)
        print "flds"
    if 'stitch' in imgDir:
        rescaleStitchImage(imgDir,destDir)


#sourceDir = "/project/marcotte/jagannath/projectfiles/EpiMicroscopy/2014-June/"+ dateStamp + "/rawImages"
#print sourceDir
#convertBatchImages(sourceDir)

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
from commonFunctions import makeDir

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
    allPngs = locate(pattern, pathDir)
    for pngFile in allPngs:
        try:
            os.remove(pngFile)
        except OSError:
            pass
    print "All Png files removed"

def convertBatchImages(dateStamp):
    pattern = '*.tif'
    # convert -contrast-stretch 0.15x0.05% [input] [output]
    allTiffs = locate(pattern, pathDir)
    for inputTif in allTiffs:
        outputTif1 = inputTif + '%d_scaled.png'
        outputTif2 = inputTif + '%d_noScale.png'
        cmd1 = "convert -contrast-stretch 0.15x0.05% " + inputTif + " " + outputTif1
        cmd2 = "convert " + inputTif + " " + outputTif2
        call(cmd1.split(),shell=False)
        call(cmd2.split(),shell=False)

def convertSameFldImages(fldLetter):
    pattern = '*'+fldLetter+'_*_fld*.tif'
    allSameFldTiffs = locate(pattern, pathDir)
    destSameDir = os.path.join(destDir,'sameFlds',sameLetter)
    if not os.path.exists(destSameDir): os.makedirs(destSameDir)
    for inputTif in allSameFldTiffs:
        dirname = os.path.split(os.path.split(inputTif)[0])[1]
        destfname = os.path.split(inputTif)[1] + '_scaled.png'
        newDestDir = os.path.join(destSameDir,dirname)
        if not os.path.exists(newDestDir):os.makedirs(newDestDir)
        outputPNG = os.path.join(newDestDir,destfname)
        cmd = "convert -contrast-stretch 0.015x0.05% " + inputTif + " " + outputPNG
        call(cmd.split(),shell=False)
        print outputPNG
        
def convertDir(newSourceDir): 
    dname = os.path.split(newSourceDir)[1]
    allInputTiffs = locate('*.tif',newSourceDir)
    newDestDir = os.path.join(destDir,dname)
    if not os.path.exists(newDestDir): os.makedirs(newDestDir)
    for inputTif in allInputTiffs:
        destfname = os.path.split(inputTif)[1] + '_scaled.png'
        outputPNG = os.path.join(newDestDir,destfname)
        cmd = "convert -contrast-stretch 0.015x0.05% " + inputTif + " " + outputPNG
        call(cmd.split(),shell=False)
        print outputPNG


def convertOneRandomImage(imgDir,destDir):
    def __makePNG(fname):
        fpath = os.path.join(imgDir,fname)
        outputPNG1 = os.path.join(destDir,fname + '_scaled.png')
        outputPNG2 = os.path.join(destDir,fname + '_noScale.png')
        cmd1 = "convert -contrast-stretch 0.015x0.05% " + fpath + " " + outputPNG1
        cmd2 = "convert " + fpath + " " + outputPNG2
        call(cmd1.split(),shell=False)
        call(cmd2.split(),shell=False)

    randomFile = choice([f for f in os.listdir(imgDir) if f.endswith('.tif')])
    if randomFile.endswith('c1.tif'): #Temp hack 
        f1 = randomFile
        f2 = randomFile[:-5]+'2.tif'
        __makePNG(f1)
        __makePNG(f2)
    elif randomFile.endswith('c2.tif'):
        f1 = randomFile[:-5]+'1.tif'
        f2 = randomFile
        __makePNG(f1)
        __makePNG(f2)
    else:
        fpath = os.path.join(imgDir,randomFile)
        __makePNG(randomFile)
    return True

def rescaleStitchImage(imgDir,destDir):
    [f] = os.listdir(imgDir)
    fpath = os.path.join(imgDir,f)
    scale = str(5)
    outputSTITCH = destDir + '/images/' + f + '.png'
    cmd = "convert " + fpath + " -resize " + scale + "% " + outputSTITCH
    call(cmd.split(),shell=False)
    print outputSTITCH 
    return True


if __name__ == '__main__':
    month = {'01':'Jan','02':'Feb','10':'Oct','11':'Nov','12':'Dec'}
    [ARG, dateStamp] = sys.argv[1:]
    sourceDir = "/project2/marcotte/boulgakov/microscope"

    sameLetter = None
    dirConvert = []
    sameLetter = 'A'
    #dirConvert = ["/project2/marcotte/boulgakov/microscope/2014-Dec/2014-12-11/141211_DiAS1_Blank_OvMeOH_Both_RF_flds005"]

    yearStamp = dateStamp.split('-')[0]
    monthStamp = yearStamp+"-"+month[dateStamp.split('-')[1]]
    pathDir = os.path.join(sourceDir,monthStamp,dateStamp)
    
    destDir = os.path.join("/project2/marcotte/jaggu/dataAnalysis/microscope1",monthStamp,dateStamp,"images","randomImages")
    makeDir(destDir)
    t0 = time.clock()   

    DIRNAMES = 1
    
    if sameLetter: convertSameFldImages(sameLetter)

    if dirConvert: 
        #for dirname in dirConvert: convertDir(dirname)
        dirname = dirConvert[0]
        convertOneRandomImage(dirname,destDir)

    for dirname in os.walk(pathDir).next()[DIRNAMES]:
        print pathDir
        imgDir = os.path.join(pathDir,dirname)
        if 'flds' in imgDir and not 'tflds' in imgDir:
            convertOneRandomImage(imgDir,destDir)
            print "flds"
        if 'stitch' in imgDir:
            rescaleStitchImage(imgDir,destDir)

    
    t1 = time.clock()
    print "Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0,
                                                       time.strftime("%d %b %Y  %H:%M:%S",time.localtime()))

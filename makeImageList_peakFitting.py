#! /home/jaggu/anaconda/bin/python2.7
from __future__ import division
import sys
import os
import time
import cv2
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
from scipy import ndimage,signal,misc,optimize
import pickle
sys.path.append('/project/marcotte/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis')
from commonFunctions import locate,simpleShow,colorImage
from ImageStack import ImageStack
#from rigidRegistration_xy import RigidRegistrationXY_guess


def makeImageList_flds(dateStamp,pathDir):
    ofname = os.path.join(pathDir,'targetImages_list.flds.txt')
    ofile = open(ofname,'w')
    img_pathList = list()

    for subDir in next(os.walk(pathDir))[1]:
        if not 'tflds' in subDir:
            subDirPath = os.path.join(pathDir,subDir)
            fList = [os.path.join(subDirPath,f) for f in os.listdir(subDirPath)
                    if f.endswith('.tif')]
            img_pathList.extend(fList)

    img_pathList.sort()
    
    for line in img_pathList:
        ofile.write(line+'\n')
    ofile.close()

    return True


def makeProjection_stack(dateStamp,pathDir,stackType='tflds'):
    if stackType == 'tflds':  ofname = os.path.join(pathDir,'targetImages_list.tflds.txt')
    else: ofname = os.path.join(pathDir,'targetImages_list.trace.txt')

    ofile = open(ofname,'w')
    img_pathList = list()
    
    for subDir in next(os.walk(pathDir))[1]:
        if stackType in subDir:
            subDirPath = os.path.join(pathDir,subDir)
            if 'both' in subDir or 'Both' in subDir:
                for ch in [1,2]:
                    fList = [os.path.join(subDirPath,f) for f in os.listdir(subDirPath) if
                         f.endswith('c'+str(ch)+'.tif')]
                    imgStack = ImageStack(fList)
                    maxImg, offsetDict = imgStack.stack(projection='maximum')
                    fImgdest,fOffsetDest = imgStack.saveProjected(maxImg,
                                                                  offsetDict,type='maximum',level=1)
                    img_pathList.append(fImgdest)
            else:
                    fList = [os.path.join(subDirPath,f) for f in os.listdir(subDirPath)]
                    imgStack = ImageStack(fList)
                    maxImg, offsetDict = imgStack.stack(projection='maximum')
                    fImgdest,fOffsetDest = imgStack.saveProjected(maxImg,
                                                              offsetDict,type='maximum',level=1,fname=subDir)
                    img_pathList.append(fImgdest)

    img_pathList.sort()
    for line in img_pathList:
        ofile.write(line+'\n')
    ofile.close()

    return True    
        

if __name__ == '__main__':
    monthIs =  {'01':'Jan','05':'May','06':'June','07':'July','08':'Aug','09':'Sept','10':'Oct','11':'Nov','12':'Dec'}

    [ARG,dateStamp] = sys.argv[1:]
    sourceDir = '/home/jaggu/marcotte_project/boulgakov/microscope'
    month = monthIs[dateStamp.split('-')[1]]
    year = dateStamp.split('-')[0]
    pathDir = os.path.join(sourceDir,year+'-'+month,dateStamp)

    t0 = time.clock()
    if ARG == 'TFLDS':
        print "Creating Maximum Projection images and Offset Dictionary ..."
        makeProjection_stack(dateStamp,pathDir,'tflds')
        print "Making Image list for peak fitting"
    elif ARG == 'TRACE':
        print "Creating Maximum Projection images for the time trace ..."
        makeProjection_stack(dateStamp,pathDir,'trace')
    
    elif ARG == 'FLDS':
        print "Finding only images for FLDS and not TFLDS"
        makeImageList_flds(dateStamp,pathDir)
    elif ARG == 'TEST':
        print "Testing...."
    else:
        raise SystemExit("Incorrect argument")

    t1 = time.clock()
    print ("Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0,time.strftime("%d %b %Y %H:%M:%S",time.localtime())))

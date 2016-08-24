#! /home/jaggu/anaconda/bin/python2.7

import os
import sys
import time
import collections
import numpy as np
import cPickle as pickle
import re

sys.path.append('/project/marcotte/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis')
from commonFunctions import locate,makeDir


"""CHANNEL,FIELD,H,W,CATEGORY,FRAME 0,FRAME 1,FRAME 2,FRAME 3,FRAME 4,FRAME 5,FRAME 6,FRAME 7
ch1,0,438.1,165.05,"(False, True, True, False, True, False, False,True)",4186.5,22855.0,20894.0,13362.5,17673.0,-1219.5,2956.5,19083.0
"""

def savePkl(db,f):
    pklPath = os.path.join(os.path.split(outputDir)[0],'pklFiles')
    if not os.path.exists(pklPath): makeDir(pklPath)
    fpath = os.path.join(pklPath,f)
    ofile = open(fpath,'w')
    pickle.dump(db,ofile)
    ofile.close()
    return True


def loadPkl(f):
    pklPath = os.path.join(os.path.split(outputDir)[0],'pklFiles')   
    fpath = os.path.join(pklPath,f)
    ifile = open(fpath,'r')
    db = pickle.load(ifile)
    ifile.close()
    return db


def makeIntensityCategory(trackFile,ch=1):
    ifile = open(trackFile,'r')
    lines = ifile.readlines()
    ifile.close()
    for line in lines[1:]:
        [catStr] = re.findall(r"\(.+\)",line)
        c
        category_tup = eval(catStr)
        intensityList = map(float,line.strip().split(',')[-8:])
        category_trackDetails_dict[category_tup].append(intensityList) 
        
    return category_trackDetails_dict


def divideIntensityList(channel=1):
    pklF = 'category_trackDetails.ch'+str(channel)+'.dict.pkl'
    db = loadPkl(pklF)
    
    
    desired_category = [[True]*(i+1) + [False]*(7-i) for i in range(8)]
    i = 0
    category = tuple(desired_category[i])
    intensityList = db[category]
    allIntensity = list()
    for l in intensityList:
        allIntensity.append(np.mean(l[0:i+1]))
    print len(allIntensity)
    print np.mean(allIntensity)
    print np.median(allIntensity)
    print np.min(allIntensity)




def makeCategoryDict():
    channel = 1 
    trackPhotometryFile = 'track_photometries_nx8t5v.csv'
    trackFile = os.path.join(outputDir,trackPhotometryFile)
    category_trackDetails_dict = makeIntensityCategory( trackFile,ch=channel)
    pklF = 'category_trackDetails.ch'+str(channel)+'.dict.pkl'
    savePkl(category_trackDetails_dict,pklF)
    print "Pickled : ",pklF

outputDir = '/project/boulgakov/microscope/2015-Oct/2015-10-27_hcomplete/output_2015-11-02'
category_trackDetails_dict = collections.defaultdict(list)

makeCategoryDict()
#divideIntensityList(channel=1)

    



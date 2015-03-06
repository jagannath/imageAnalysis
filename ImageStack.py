#! /home/jaggu/anaconda/bin/python2.7
from __future__ import division
import sys
import os
import time
import cv2
import numpy as np
import pickle
sys.path.append('/project/marcotte/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis')
import commonFunctions as cfunc


"""
This programme/class manipulates an image stack. It can generate - (a) maximum
intensity projection, (b) median intensity projection 
I will try to enhance it to have a stack of images from a trace and do some kind
of trace file analysis.
Outputs a tif file with the appropriate tags for further analysis. 
Also if not asked - no offset will be calculated. But an offset dictionary
pickle is saved that contains (x,y) = (0,0) offset against every file path in
the image stack. 
"""

class ImageStack:


    """
    Stack of images. Currently set up to do simple things like maximum intensity
    projection. In future will do a lot more single molecule tracking
    """
    def __init__(self,imagesList,IDX=0):
        assert isinstance(imagesList,list)
        self.imagesList = imagesList
        self.imagesList.sort()
        self.IDX = IDX

    def stack(self,projection='maximum',offset=None,IDX=0):
        if projection is 'maximum':
            projection_img, offset_dict = self._maxIntensity(offset,IDX)
        elif projection is 'median':
            projection_img, offset_dict = self._medIntensity(offset,IDX)
        else:
            raise SystemExit('Invalid projection on Image List')
        return projection_img,offset_dict

    def calculateOffset(template,img):
        # Need to edit it appropriately. 
        return ( 0,0)
    
    def _maxIntensity(self,offset=None,IDX=0):
        template = cv2.imread(self.imagesList[IDX],-1) 
        maxi = np.zeros([ 512,512],dtype=np.uint16)  
        offsetDict = dict()

        for f in self.imagesList:
            img = cv2.imread(f,-1)
            if offset: 
                delta = self.calculateOffset(template,img)
                print f, delta
            else:
                delta = ( 0,0)                               
                print f, delta
            (x,y) = delta
            (h,w) = (y,x)
            offsetDict[f] = delta
            M = np.float32([[1,0,w],[0,1,h]])
            rimg = cv2.warpAffine(img,M,template.shape)
            ix = rimg > maxi #Index where val[ar1]>val[ar2]
            maxi[ix] = rimg[ix]
        
        return maxi,offsetDict
   
    def _medIntensity(self,offset=None,IDX=0):
        template = cv2.imread(self.imagesList[IDX],-1)
        offsetDict = dict()
        allImgs = list()
        
        for f in self.imagesList:
            img = cv2.imread(f,-1)
            if offset: 
                delta = self.calculateOffset(template,img)
                print f, delta
            else:
                delta = ( 0,0)                               
                print f, delta
            (x,y) = delta
            (h,w) = (y,x)
            offsetDict[f] = delta
            M = np.float32([[1,0,w],[0,1,h]])
            rimg = cv2.warpAffine(img,M,template.shape)
            
            allImgs.append(rimg)
        
        final_stack = np.dstack(allImgs)
        medi = np.median(final_stack, axis=2)
        return medi,offsetDict


    def saveProjected(self,img,offsetDict,type='max',dest_dir=None,fname=None,level=1):
        if not dest_dir:
            dest_dir = os.path.split(self.imagesList[self.IDX])[0]
        if not fname:
            fname = '.'.join([self.imagesList[self.IDX][:-4],type+'_projection'])
        
        target_dir = os.path.join(dest_dir,'projectionImages')
        cfunc.makeDir(target_dir)
        f_img_name = os.path.join(target_dir,fname+'.level_'+str(level)+'.image.tif')
        f_offset_name = os.path.join(target_dir,fname+'.level_'+str(level)+'.offset_dict.pkl')
    
        cv2.imwrite(f_img_name,img)
        with open(f_offset_name,'wb') as ofile:
            pickle.dump(offsetDict,ofile)
        
        return f_img_name, f_offset_name













def test():
    pathDir = "/project2/marcotte/boulgakov/microscope/2015-Jan/2015-01-01"
    subDir = "150101_DiAS2_PR011_2nM_Mock1_TFA1h_A_Both_tflds002"

    subDirPath = os.path.join(pathDir,subDir)
    for ch in [1,2]:
        fList = [os.path.join(subDirPath,f) for f in os.listdir(subDirPath) if
                 f.endswith('c'+str(ch)+'.tif')]
        imgStack = ImageStack(fList)
        pmedImg, offsetDict = imgStack.stack(projection='median')
        fImgdest,fOffsetDest = imgStack.saveProjected(pmedImg,offsetDict,type='median')

    return True




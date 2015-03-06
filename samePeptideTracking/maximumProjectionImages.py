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
#from rigidRegistration_xy import RigidRegistrationXY_guess
from pylab import *#Remember this stupid 'where' 

class MaximumProjection:
    """
    This class is meant to obtain the maximum value of the pixel from a stack
    (list) of images. Currently there is no averaging or background subtraction,
    but it is a possibility. Finally there will be functions to save the tiff
    image or return the numpy array. Before the images are compared, it is
    aligned using the simple RigidRegistration function
    A pkl file is generated that contains the offset information between
    the image and the template image (typically the first one).
    Although offset is not always calculated (especially for the trace files),
    the function is set up. Offset calculation is time consuming and does not
    always give the best results. 
    If imgs are t0,t1,...; Pkl file will contain dictionary of the offset for
    the image when compared to the target image (maxProjection.tif)
    """

    def __init__(self,fList,startIDX= 0):
        self.fList = fList
        self.f0 = fList[startIDX]
        self.startIDX = 0
        self.orig,self.grayimg,self.cimg = self._openImages(self.fList[startIDX])

    def _openImages(self,fname,flag=-1):
        orig = cv2.imread(fname,flag)
        corig = cv2.imread(fname,0)#Always to 8bit color
        cimg = cv2.cvtColor(corig,cv2.COLOR_GRAY2BGR)
        return (orig, corig,cimg)

    def MeasureErr(self,img1,img2):
        diff = (img1-img2) # This is MSE
        # Need to investigate this later
        return sum(diff**2)

    def RigidRegistrationXY_guess(self,img,ximg):
        '''
        @author: Christian Rossmann, PhD
        @license:  Public Domain
        @blog: http://scientificcomputingco.blogspot.com/
        '''
    
        # Perform initial guess rotation & translation 
        v_range =  np.array(xrange(-10,10))
    
        err = np.array([self.MeasureErr(img,ndimage.shift(ximg,(v,0))) for v in v_range])
        x = v_range[np.where(err==err.min())[0]] 
        # Need to check what is this where from pylab. It gives the correct
        # result different from np.where    
        err = np.array([self.MeasureErr(img,ndimage.shift(ximg,(0,v))) for v in v_range])
        y = v_range[np.where(err==err.min())[0]]

        # List contains displacement in x and y and rotation
        param = [x,y]

        # Scipy optimize 
        def _ErrFunc(param,img=img,ximg=ximg):                              
        
            # Perform rotational and translational transformation
            _img = ximg.copy()
            _img = ndimage.shift(_img,param)
            return self.MeasureErr(img,_img)

        param = optimize.fmin(_ErrFunc,param,maxiter=100)
        return param,err.min()

    def calculateOffset(self,img,ximg):
        #rimg,param,guess = RigidRegistrationXY(tif1.orig,tif2.orig)           
        offset = ( 0, 0)
        guess1,errmin1 = self.RigidRegistrationXY_guess(img,ximg)
        guess2,errmin2 = self.RigidRegistrationXY_guess(ximg,img)
        if errmin1>errmin2:
            offset = [i*-1 for i in guess2]
        else: 
            offset = guess1
        return offset

    def obtainMaxProjection(self,offset=True):
        template = self.orig.copy()
        maxi = self.orig.copy()
        offsetDict = dict()
        for f in self.fList:
            img,grayimg, cimg = self._openImages(f)
            if offset is True: 
                delta = self.calculateOffset(template,img)
                print f, delta
            else:
                delta = ( 0,0)                               
                print f, delta
            (h,w) = delta
            offsetDict[f] = delta
            M = np.float32([[1,0,w],[0,1,h]])
            rimg = cv2.warpAffine(img,M,template.shape)
            ix = rimg > maxi #Index where val[ar1]>val[ar2]
            maxi[ix] = rimg[ix]
        return maxi,offsetDict

    def applyOffset(self,offsetDict,exptFnames,ch=1):
        f1Template = self.f0[:-5]+str(ch)+'.tif'
        orig1,gray1,cimg1 = self._openImages(f1Template)
        templateImg = orig1.copy()
        maxi = orig1.copy()
        for expt,(f1,f2) in exptFnames.items():
            offset = offsetDict[f2]
            (h,w) = offset#It returns x,y or w,h
            print f1,f2, h,w            
            img,gray,cimg = self._openImages(f1)
            M = np.float32([[1,0,w],[0,1,h]])
            rimg = cv2.warpAffine(img,M,templateImg.shape)
            #Place where can be buggy; Offset is x,y; cv2 wants h,w
            ix = rimg > maxi #Index where val[ar1]>val[ar2]
            maxi[ix] = rimg[ix]
        return maxi

def traceFiles(pathDir):
    """
    This function goes through all the trace folders in the directory. 
    """   
    type = 'trace'
    offset = False

    fwrite = os.path.join(pathDir,'targetImages_list.txt')
    ofile = open(fwrite,'w')
    for subDir in next(os.walk(pathDir))[1]:
        if type in subDir:
            traceDir = os.path.join(pathDir,subDir)
            fList = [os.path.join(traceDir,f) for f in os.listdir(traceDir) 
                    if not (f.endswith('maxProjImage.tif') or
                            f.endswith('.pkl'))]
            fList = [f for f in fList if os.path.isfile(f)]
            fList.sort() #First file is template by default
            t = MaximumProjection(fList)
            maxImg, offsetDict = t.obtainMaxProjection(offset=offset)
            destDir = os.path.join(pathDir,subDir,'maxProjection')
            if not os.path.exists(destDir): os.makedirs(destDir)
            f = os.path.join(destDir,subDir+'.maxProjImage.tif')
            cv2.imwrite( f,maxImg)
            pklFile = f + '.offsetDict.pkl'
            pickle.dump(offsetDict,open(pklFile,'w'))
            ofile.write(f+"\n")
    ofile.close()
    print "List created in %s"%(fwrite)
    

def getFpaths(frame):
    allFpaths = list() #tuple of (561,647)..Can extend
    pathList = list()

    pathDir1 = '/project2/marcotte/boulgakov/microscope/2014-Dec/test/A/mock1'
    subDir1 = '141213_DiAS1_PR011_2pM_Mock1_TFA1h_Both_A_flds004'

    pathDir2 = '/project2/marcotte/boulgakov/microscope/2014-Dec/test/A/mock2'
    subDir2 = '141213_DiAS1_PR011_2pM_Mock2_TFA1h_Both_A_flds007' 

    pathDir3 = '/project2/marcotte/boulgakov/microscope/2014-Dec/test/A/edman3'   
    subDir3 = '141213_DiAS1_PR011_2pM_Edman1_TFA1h_Both_A_flds008' 

    pathList = [(pathDir1,subDir1),(pathDir2,subDir2),(pathDir3,subDir3)]
    exptPaths = {'mock1':(pathDir1,subDir1),
                  'mock2':(pathDir2,subDir2),
                  'edman':(pathDir3,subDir3)}
    exptFnames_dict = dict()

    for expt,(pathDir,subDir) in exptPaths.items():
        # NOTE; PUT THIS IN THE ORDER YOU WANT
        fnamec1 = ''.join([subDir,frame,'c1','.tif'])
        fnamec2 = ''.join([subDir,frame,'c2','.tif'])
        fc1 = os.path.join(pathDir,fnamec1)
        fc2 = os.path.join(pathDir,fnamec2)
        allFpaths.append((fc1,fc2))
        exptFnames_dict[expt] = (fc1,fc2)
    return allFpaths,exptFnames_dict



def makeMaxProjection():
    # Same channel for reference; Same frame
    #ch = 'c2' #REF; THIS BECOMES agnostic
    
    frame = 'xy14'
    sameFld = 'A'
    exptOrder = ('mock1','mock2','edman')
    
    allFpaths,exptFnames_dict = getFpaths(frame)
    c2Fpaths = list()
    for expt in exptOrder:
        c2Fpaths.append(exptFnames_dict[expt][1])
    test = MaximumProjection(c2Fpaths)
    maxImg2,offsetDict = test.obtainMaxProjection(offset=True) 
    maxImg1 = test.applyOffset(offsetDict,exptFnames_dict)

    # Making image of the offsetdict
    destDir = '/project2/marcotte/boulgakov/microscope/2014-Dec/test/A/maxProjection'   
    f1 = 'testCycle0_1_2.c1.maxProjImage'
    f2 = 'testCycle0_1_2.c2.maxProjImage'

    f1 = os.path.join(destDir,f1+'.tif')
    f2 = os.path.join(destDir,f2+'.tif')
    cv2.imwrite(f1,maxImg1)
    cv2.imwrite(f2,maxImg2)
    
    ofile = open(f2+'.listImages.txt','w')
    ofile.write(f1+'\n'+f2)
    ofile.close()
    pklFile = f2+'offsetDict.pkl'
    pickle.dump(offsetDict,open(pklFile,'w'))


if __name__ == '__main__':
    monthIs =  {'05':'May','06':'June','07':'July','08':'Aug','09':'Sept','10':'Oct','11':'Nov','12':'Dec'}

    [ARG,dateStamp] = sys.argv[1:]
    sourceDir = '/project2/marcotte/boulgakov/microscope'
    month = monthIs[dateStamp.split('-')[1]]
    pathDir = os.path.join(sourceDir,"2014-"+month,dateStamp)

    t0 = time.clock()
    if ARG == 'MAXTRACE':
        print "Creating Maximum Projection images ..."
        traceFiles(pathDir)
    elif ARG == 'TEST':
        print "Testing...."
        test_case()
    else:
        raise SystemExit("Incorrect argument")

    t1 = time.clock()
    print ("Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0,time.strftime("%d %b %Y %H:%M:%S",time.localtime())))
      

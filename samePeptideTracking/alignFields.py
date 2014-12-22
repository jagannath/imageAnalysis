#! /home/jaggu/anaconda/bin/python2.7
from __future__ import division
import sys
import os
import time
import cv2
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
from scipy import ndimage,signal,misc
import pickle
sys.path.append('/project/marcotte/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis')
from commonFunctions import locate,simpleShow,heatMapImage,colorImage
#from rigidRegistration_xy import RigidRegistrationXY,RigidRegistrationXY_guess
from pylab import * 
#FOR THE LIFE OF ME; I CANT UNDERSTAND HOW THIS EFFECTS np.where

def bounds(x):
    if x < 0: x = 0
    if x > 511: x = 511
    return x


class TifImage:
    """
    This class handles a tif image. It opens them,retrieves the pkl and
    the associated peak information. It manipulates the image to create a black
    and white image (making the 3x3 peak as black. 
    
    PKL FILE
   {#KEY(h_0 rounded to nearest integer (i.e. to nearest pixel),
         w_0 rounded to nearest integer (i.e. to nearest pixel)):
    #VALUE
    #0---#1---#2-#3--#4-----#5-------#6------#7------#8------#9----#10--#11
    (h_0, w_0, H, A, sigma_h,sigma_w, theta, sub_img,fit_img, rmse, r_2, s_n),
     ...
     }
    
    NOTE: opencv uses BGR; hw is a better way to find pixel. To draw use xy
    """
    def __init__(self,tifF):
        self.tifF = tifF
        self.subDir,self.fname = os.path.split(tifF) 
        self.orig,self.grayimg, self.cimg = self._openImages(self.tifF)
        self.dim = self.orig.shape
        self.peakPkl = self._openPkl()

    def _openImages(self,fname,flag=-1):
        orig = cv2.imread(fname,flag)
        corig = cv2.imread(fname,0)#Always to 8bit color
        cimg = cv2.cvtColor(corig,cv2.COLOR_GRAY2BGR)
        return (orig, corig,cimg)

    def _showFig(self,img):
        simpleShow(img)
    
    def _openPkl(self):
        # Temporary and assumes only one pkl for that file
        [pklF] = [f for f in os.listdir(self.subDir) if f.startswith(self.fname) and f.endswith('.pkl')]
        pklFpath = os.path.join(self.subDir,pklF)
        return pickle.load(open(pklFpath))

    def _getRect(self,p,l):
        (h,w) = p
        posRect = [(i,j) for i in range(h-l+2,h+l-1) for j in
                   range(w-l+2,w+l-1)]
        pxLoc = [(bounds(i),bounds(j)) for (i,j) in posRect]
        return pxLoc

    def getAllPeakPositions(self):
        return self.peakPkl.keys()

    def colorPeaks(self,color=( 255,0,0)):
        horig,worig = self.dim
        peakDim = 3 #3x3 square
        cimg = np.zeros((horig,worig,3),np.uint8)
        cimg[:] = ( 0,0,0) # Making the background black
        bwimg = np.zeros((horig,worig),np.uint8)
        allPos = self.getAllPeakPositions()
        for pos in allPos:
            h0,w0 = map(int,pos)
            pxReg = self._getRect((h0,w0),peakDim)
            for (h,w) in pxReg: 
                cimg[h,w] = color
                bwimg[h,w] = 255
        return bwimg,cimg

    def colorMapImage(self,color='gist_heat',stretch=True):
        self.orig_cmap = colorImage(self.orig,color=color,stretch=stretch)
        return self.orig_cmap


class AlignImages:
    """
    The major class that should be imported to calculate offset. The most this
    should do is return h,w displacement of the first image with respect to the
    second. Currently it does not even call any optimization routine. It however
    calls a simple function currently called RigidRegistrationXY_guess from a file
    rigidRegistration_xy. This is a simple MSE minimiazation. The entire file
    moves 10x10 displacement on either side and calculates the difference. The
    lowest is used for alignment. The offset calculated is (h,w); 
    [NOTE- not x,y]
    Further functions are to call dest images or overlayed images in the form
    of peaks,grayscale,colorscale
    """

    def __init__(self,f1,f2):
        self.f1,self.f2 = f1,f2 #Full file name path
        self.tif1 = TifImage(f1)
        self.tif2 = TifImage(f2)
        assert self.tif1.orig.shape == self.tif2.orig.shape
        self.dim = self.tif1.orig.shape
        self.imgx = self.dim[0]

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
        print "Initial Guess : %d, %d"%(x,y)
        guess = [x,y]
        return guess

    def calculateOffset(self):
        #rimg,param,guess = RigidRegistrationXY(tif1.orig,tif2.orig)           
        guess = self.RigidRegistrationXY_guess(self.tif1.orig,self.tif2.orig)
        offset = guess
        h,w = np.around(offset)
        h,w = map(int,[h,w])
        return (h,w)



class FollowPeaks:
    """
    From the calculated offset and the two files that differ in the offset; i.e
    in the form (h,w); Also I will make a box around the peak and follow it in
    the second image after applying the offset and giving it a little wiggle
    room or 3 px radius. 
    # Someday I need to inherit the TifImage class
    """
    def __init__(self,f1,f2,offset):
        self.f1,self.f2 = f1,f2 
        #Full path of the file; Also assumes that the pklData is generated
        self.offset = offset; #(h,w) form
        self.tif1 = TifImage(f1)
        self.tif2 = TifImage(f2)
        assert self.tif1.orig.shape == self.tif2.orig.shape
        self.peakData1 = self.tif1.peakPkl
        self.peakData2 = self.tif2.peakPkl
        self.peakDataL = [self.peakData1,self.peakData2]

    def getPeakSummary(self):
        f1allPeakPos =  self.tif1.getAllPeakPositions()
        f2allPeakPos = self.tif2.getAllPeakPositions()
        nbrf1,nbrf2 = map(len, [f1allPeakPos, f2allPeakPos])
        nbrCommon = int()
        nbrNotMapped = int()

        for f1peak in f1allPeakPos:
            p1Info, p2Info = self.findSamePeak(f1peak)
            if p2Info[2]:
                nbrCommon+=1
            else:
                nbrNotMapped+=1
        return nbrf1,nbrf2,nbrNotMapped,nbrCommon

    
    def findNeighbouringPeak(self,p,rad=3):
        allPeaks = list()
        (h,w) = p
        posRect = [(i,j) for i in range(h-rad+2,h+rad-1) for j in
                   range(w-rad+2,w+rad-1)]
        pxLoc = [(bounds(i),bounds(j)) for (i,j) in posRect]
        for pPos in pxLoc:
            if pPos in self.peakData2.keys():
                allPeaks.append(pPos)
        if len(allPeaks)==1:return allPeaks[0]
        else: return None

    def _getExactPeakInfo(self,peakPos,displacement=None,frame=0):
        # Frame tells which peakData to search
        # Default is no posChange. If posChange=True; then have offset
        if displacement: 
            delH,delW = self.offset
        else:
            delH,delW = ( 0,0)
        dispH,dispW = map(int,[peakPos[0]-delH,peakPos[1]-delW])
        try:
            pInfo = self.peakDataL[frame][(dispH,dispW)]
            newH,newW = map(int,[dispH,dispW])
        except KeyError:
            nearPeakPos = self.findNeighbouringPeak((dispH,dispW))
            if nearPeakPos:
                newH,newW = nearPeakPos
                pInfo = self.peakDataL[frame][(newH,newW)]
            else:
                newH,newW = 0,0
                pInfo = None
#        print "Looked for peak at %d,%d"%(dispH,dispW)
#        print "Found peak at %d,%d"%(newH,newW)
        return [(dispH,dispW),(newH,newW),pInfo]

    def findSamePeak(self,peakPos):
        # Assumes that this peakPos was found somehow; It then uses the
        # peakData1 and looks at peakData2
        # pInfo = (dispH,dispW),(newH,newW),pInfo
        
        p1Info = None
        p2Info = None
        
        p1Info = self._getExactPeakInfo(peakPos)
        p2Info = self._getExactPeakInfo(peakPos,displacement=True,frame=1)
       
        pos1 = p1Info[1]
        pos2 = p2Info[1]

        return p1Info,p2Info

    def drawBoxImg(self,pos1,pos2,img):
        # This is the newPos; It should map to a peak
        if len(img.shape)==2: # This is grayscale;
            cimg = colorImage(img,color='gray')
        else: cimg = img

        y1,x1 = pos1
        y2,x2 = pos2
        
        cv2.rectangle(cimg,(x1-3,y1-3),(x1+3,y1+3),color=( 255,0,0),thickness=1)
        cv2.rectangle(cimg,(x2-3,y2-3),(x2+3,y2+3),color=( 0,255,0),thickness=1)
        return cimg

def overlayImages(f1,f2,offset):
    # f1,f2; Full path; Offset is (dh,dw)
    
    tif1 = TifImage(f1)
    tif2 = TifImage(f2)
    dim = tif1.orig.shape

    def _addMoveImg(offset):
        h,w = np.around(offset)
        M = np.array([[1,0,w],[0,1,h]],dtype=np.float32)

        newf2col = cv2.warpAffine(colimg2,M,dim)
        dest_peaks = cv2.add(colimg1,newf2col)

        newf2_gray = cv2.warpAffine(tif2.orig,M,dim)
        dest_gray = cv2.add(tif1.orig,newf2_gray)
        
        newf2_cmap = cv2.warpAffine(f2_cmap,M,dim)
        dest_cmap = cv2.addWeighted(f1_cmap,0.5,newf2_cmap,0.7,0)
        return dest_peaks,dest_gray,dest_cmap

    def _overlayImg(offset):
        overlay_peaks = cv2.add(colimg1,colimg2)
        overlay_gray = cv2.add(tif1.orig,tif2.orig)
        overlay_cmap = cv2.addWeighted(f1_cmap,0.5,f2_cmap,0.7,0)
        return overlay_peaks,overlay_gray,overlay_cmap

        
    bwimg1,colimg1 = tif1.colorPeaks()
    f1_cmap = tif1.colorMapImage(color='black_red',stretch='linear')
    bwimg2,colimg2 = tif2.colorPeaks(color=( 0,255,0))
    f2_cmap = tif2.colorMapImage(color='black_green',stretch='linear')

    aligned_peaks,aligned_gray,aligned_cmap = _addMoveImg(offset)
    overlay_peaks,overlay_gray,overlay_cmap = _overlayImg(offset)
        
    return ([aligned_peaks,aligned_gray,aligned_cmap],
            [overlay_peaks,overlay_gray,overlay_cmap])

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

    for pathDir,subDir in pathList:
        # NOTE; PUT THIS IN THE ORDER YOU WANT
        fnamec1 = ''.join([subDir,frame,'c1','.tif'])
        fnamec2 = ''.join([subDir,frame,'c2','.tif'])
        fc1 = os.path.join(pathDir,fnamec1)
        fc2 = os.path.join(pathDir,fnamec2)
        allFpaths.append((fc1,fc2))

    return allFpaths

def getAllPeaks(fpaths,pepCh):
    allPepFpaths = [f[pepCh] for f in fpaths]
    allPeakPos = list()
    for fc1 in allPepFpaths:
        tif = TifImage(fc1)
        allPeakPos.extend(tif.getAllPeakPositions())

    return allPeakPos



def test_case():
    # Same channel for reference; Same frame
    #ch = 'c2' #REF; THIS BECOMES agnostic
    frame = 'xy14'
    
    allFpaths = getFpaths(frame)

    ref = 1
    pepCh = 0
    compareCycles = [( 0,1),( 0,2),( 1,2),( 1,0),( 2,0),( 2,1)]
    totalPeaksPos = getAllPeaks(allFpaths,pepCh)
    print "Total peaks identified in total : %d"%len(totalPeaksPos)

    for cycle1,cycle2 in compareCycles:
        print "xxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
        print "Mapping Cycle : %d on Cycle : %d"%(cycle1,cycle2)
        f1c1,f2c1 = allFpaths[cycle1][pepCh],allFpaths[cycle2][pepCh]
        f1c2,f2c2 = allFpaths[cycle1][ref],allFpaths[cycle2][ref]
        twoCycles = AlignImages(f1c2,f2c2)
        offset = twoCycles.calculateOffset()
        aligned,overlay = overlayImages(f1c2,f2c2,offset)
        simpleShow(aligned[1])
        sys.exit(1)
        fp = FollowPeaks(f1c1,f2c1,offset)
        nbrf1c1,nbrf2c1,notMappedPeaks,mappedPeaks = fp.getPeakSummary()
        print "Summary in channel : %d"%(pepCh)
        print "Total Peaks in f1 : %d"%(nbrf1c1)
        print "Total Peaks in f2 : %d"%(nbrf2c1)
        print "Peaks in f1 not mapped to f2 : %d" %(notMappedPeaks)
        print "Peaks in f1 mapped to f2 : %d"%(mappedPeaks)
        print "#############################"

if __name__ == '__main__':
    monthIs = {'05':'May','06':'June','07':'July','08':'Aug','09':'Sept','10':'Oct','12':'Dec'}
    
    [ARG,dateStamp] = sys.argv[1:]

    t0 = time.clock()
    if ARG == 'TEST':
        print "Testing...."
        test_case()
    else:
        raise SystemExit("Incorrect argument")
    t1 = time.clock()
    print ("Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
          )

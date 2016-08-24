#! /home/jaggu/anaconda

import sys
import cPickle as pickle
import os
import time
import cv2
import numpy as np

# General functions
def savePkl(db,fpath):
    ofile = open(fpath,'w')
    pickle.dump(db,ofile)
    ofile.close()
    return True

def loadPkl(fpath):
    ifile = open(fpath,'r')
    db = pickle.load(ifile)
    ifile.close()
    return db

def openImage(fname,flag=-1):
    orig = cv2.imread(fname,flag)
    corig = cv2.imread(fname,0)#Always to 8bit color
    cimg = cv2.cvtColor(corig,cv2.COLOR_GRAY2BGR)
    return (orig, corig, cimg)

def saveImage(fname,img):
    cv2.imwrite(fname,img)
    return fname

def simpleShow(img):
    # Shows img when it is a matrix
    plt.imshow(img,cmap='gray')
    plt.show()
    return True

def makeDir(d):
    if not os.path.exists(d): os.makedirs(d)
    return d


"""
##
# PEAK 
         (h_0 rounded to nearest integer (i.e. to nearest pixel),
          w_0 rounded to nearest integer (i.e. to nearest pixel))
# VALUE
         (h_0, w_0, H, A, sigma_h, sigma_w, theta, sub_img, fit_img,
          rmse, r_2, s_n)
          
s_n = (max(sub_img) - mean(img_edge)) / stddev(img_edge)
##

Going to make a pickle path and analyzing the peaks through the traceDir..
Going to make an output directory in project2. A directory for each traceDir_[imgProjection] is created
pkl_files
    peak_intensity.pkl;
    peak_bgIntensity.pkl;
    peak_SNR.pkl;
    peak_Status.pkl;
"""

class Peak:
    def __init__(self,pos,subImage):
        self.pos = pos
        self.subImage = subImage
    def intensity(self):
        return np.sum(self.subImage[1:4,1:4])
    def bg_median(self):
        peakIndexList = [(r,c) for r in [1,2,3] for c in [1,2,3]]
        img_edge = [v for rc,v in np.ndenumerate(self.subImage) if not rc in peakIndexList]
        return (np.median(img_edge)*9)
    def bg_mean(self):
        peakIndexList = [(r,c) for r in [1,2,3] for c in [1,2,3]]
        img_edge = [v for rc,v in np.ndenumerate(self.subImage) if not rc in peakIndexList]
        return np.mean(img_edge)*9
    def snr(self):
        peakIndexList = [(r,c) for r in [1,2,3] for c in [1,2,3]]
        img_edge = [v for rc,v in np.ndenumerate(self.subImage) if not rc in peakIndexList]
        snr = ( np.max(self.subImage) - np.mean(img_edge) )/ np.std(img_edge)
        return snr

class PeakTrace:
    """
    The class handles an identified peak (after peak fitting) for the time trace files. 
    The main function is to calculate the subImage for every image in the traceDirectory. It calls the Peak class 
    and gets the statistics involved for the subImage identified. 
    """
    def __init__(self,pos,traceDir):
        self.pos = map(int, pos)
        try: 
            assert type(self.pos[0]) is int and type(self.pos[1]) is int
        except AssertionError:
            print self.pos
            sys.exit(1)
        self.traceFiles = sorted([os.path.join(traceDir,f) for f in os.listdir(traceDir) if f.endswith('.tif')])
    
    def bounds(self,x):
        if x < 0: x = 0
        if x > 511: x = 511
        return x
    
    def getRect(self,p,l=4):
        # 4 is for a 5x5 box
        (h,w) = p
        posRect = [(i,j) for i in range(h-l+2,h+l-1) for j in
                   range(w-l+2,w+l-1)]
        pxLoc = [(self.bounds(i),self.bounds(j)) for (i,j) in posRect]
        return pxLoc
    
    def getSubImage(self,pxL,f,dim=5):
        subimgL = list() 
        orig,gray8bit,cimg = openImage(f)
        for h,w in pxL:
            val = orig[(h,w)]
            subimgL = np.append(subimgL,val)
        sub_image = subimgL.reshape( dim,dim)
        return sub_image
    
    def runThroughFrames(self):
        pxLoc = self.getRect(self.pos)
        intensity_list = list()
        bg_median_list = list()
        snr_list = list()
        status_list = list()
        
        for f in self.traceFiles:
            sub_image = self.getSubImage(pxLoc,f)
            p = Peak(self.pos,sub_image)
            intensity_list.append(p.intensity())
            bg_median_list.append(p.bg_median())
            snr_list.append(p.snr())
            if p.snr()>4: status = True
            else: status = False
            status_list.append(status)
        
        return intensity_list, bg_median_list, snr_list, status_list
            


def analyzeFrames(sourceDir,outputDir,peakFilePath):
    # peakFile is the '.pkl' filepath which contains the peak information after image processing. 
    # an output directory structure is created based on the file name
    peak_dict = loadPkl(peakFilePath)
   
    peakFile = os.path.split(peakFilePath)[1]
    traceDir_name = peakFile.split('.')[0][:-4]# This is a hack you may have to do

    traceDir = os.path.join(sourceDir,traceDir_name)
    print "Analyzing directory ",traceDir_name,"...",len(peak_dict)

    peakOutput = peakFile.split('.')[0]+'_'+peakFile.split('.')[1]+'_output'
    peak_outputDir = makeDir(os.path.join(outputDir,peakOutput))
    
    pklPath = os.path.join(peak_outputDir,'pklFiles')
    if not os.path.exists(pklPath): makeDir(pklPath)
    
    peak_intensity_dict = dict()
    peak_snr_dict = dict()
    peak_bg_median_dict = dict()
    peak_status_dict = dict()

    for pos,val in peak_dict.items():
        # For each peak pos, it runs through all the frames in the time trace directory
        ptrace = PeakTrace(pos,traceDir)
        intensity_list, bg_median_list, snr_list, status_list = ptrace.runThroughFrames()
        peak_intensity_dict[pos] = intensity_list
        peak_bg_median_dict[pos] = bg_median_list
        peak_snr_dict[pos] = snr_list
        peak_status_dict[pos] = status_list

    savePkl(peak_intensity_dict,os.path.join(pklPath,'peak_intensity.dict.pkl'))
    savePkl(peak_bg_median_dict,os.path.join(pklPath,'peak_bg_median.dict.pkl'))
    savePkl(peak_snr_dict,os.path.join(pklPath,'peak_snr.dict.pkl'))
    savePkl(peak_status_dict,os.path.join(pklPath,'peak_status.dict.pkl'))
    
    return traceDir, pklPath










# In[ ]:




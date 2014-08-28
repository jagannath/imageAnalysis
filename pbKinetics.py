#! /home/jaggu/anaconda/bin/python2.7

"""
This programme does the image analysis of time trace files only. First it uses
Alexander's peak identifying algorithm on the first frame to identifying peaks.
Then (a) I measure the intensity of the 5x5px area of the peaks in the other
frames. I count whether the peak is still on based on the background intensity
of that frame. 
"""

import scipy.misc
import numpy as np
import scipy as sp
import os
import sys
from subprocess import call
import pflibMP
import cv2
import fnmatch
import re
import cPickle as pickle
import collections

sourceDir = '/project2/marcotte/boulgakov/microscope/2014-July/2014-07-27/AS2_Atto647NSE_2aM_16hWash_647_OSS2_50Perc_trace002_flds001'
imTIFF = 'AS2_Atto647NSE_2aM_16hWash_647_OSS2_50Perc_trace002_flds001t0001.tif'


class StarterFrame(object):
    """
    This is the peak finder algorithm (Alexander's) applied to the first frame.
    The list of peaks at the coordinates and the corresponding 5x5 pixel area
    for each of the peaks (h,w) is returned eventually
    """
    def __init__(self, inputTiff):
        self.inputTiff = inputTiff
        self. correlation_matrix=np.array([[-5935,-5935,-5935,-5935,-5935],
                             [-5935,8027,8027,8027,-5935],
                             [-5935,8027,30742,8027,-5935],
                             [-5935,8027,8027,8027,-5935],
                             [-5935,-5935,-5935,-5935,-5935]])

    def convertPNG(self):
        outputPNG = self.inputTiff + '.png'
        cmd = "convert " + self.inputTiff + " " + outputPNG
        call(cmd.split(),shell=False)
        return outputPNG

    def find_starterPeaks(self):
        self.outputPNG = self.convertPNG()
        img = scipy.misc.imread(self.outputPNG)
        results = pflibMP.fit_peaks([img],correlation_matrix=self.correlation_matrix)
        self.allPeakResults = results[0][1]
        return self.allPeakResults

    def getPeakXY(self,peakResult):
        h,w = peakResult[0],peakResult[1]
        return (h,w)
    
    def picklePeakCoords(self,allPeakCoords):
        ofname = self.inputTiff+'.pkl'
        ofile = open(ofname,'w')
        pickle.dump(allPeakCoords,ofile)
        ofile.close()
        return ofname

  
class TiffStack(object):
    """
    This gets an array of peaks with its 5x5 pixel number area. It runs through
    the tiff stack and returns back the intensity under that area and the
    background for each Tiff. Various outputs like number of peaks for each
    frame can be obtained. 
    @ arg: fname - This is the filename for the first frame Tiff file
    @ arg: allPeakInfo - This is the list of information about all the peaks
    """
    
    def __init__(self,fname,allPeakCoords):
        self.fname = fname #Name of 1st frame file
        self.allPeakCoords = allPeakCoords
        self.traceDir = os.path.dirname(self.fname)
        self.allTiffs = [os.path.join(self.traceDir,f) for f in
                         os.listdir(self.traceDir) if f.endswith('.tif')]
        self.allPeaks_dict = collections.defaultdict(list)

    def isItPeak(self,avg_intensity,background_intensity,stdev_intensity):
        snr = (avg_intensity - background_intensity)/stdev_intensity
        if snr > 2: return 1
        else: return 0

    def allPeakTraceIntensity(self,f,frameNbr):
        img = cv2.imread(f,-1)
        def _getPeakIntensity(peakAreaXY):
            peaklist = []
            for (i,j) in peakAreaXY: peaklist.append(img[i,j])
            return sum(peaklist)/len(peakAreaXY)            
        
        def _getBackgroundIntensity(bgAreaXY):
            bglist =  []
            for (i,j) in bgAreaXY:
                bglist.append(img[i,j])
            arr = np.array(bglist)
            return arr.mean(), arr.std()
        

        for (h,w) in self.allPeakCoords:

            sum_intensity = 0
            x,y = int(np.around(h)), int(np.around(w))
            if not (x in range(2,509)  and y in range(2,509)): continue; # Takes care of borders

            all_peakAreaXY = [(i,j) for i in xrange(x-2,x+3) for j in xrange(y-2,y+3)]
            peakAreaXY = [(i,j) for i in xrange(x-1,x+2) for j in xrange(y-1,y+2)]
            bgAreaXY = [(i,j) for (i,j) in all_peakAreaXY if not (i,j) in peakAreaXY]
            avg_peak_intensity = _getPeakIntensity(peakAreaXY)
            avg_background_intensity,background_stdev = _getBackgroundIntensity(bgAreaXY)
            
            peakStatus = self.isItPeak(avg_peak_intensity,avg_background_intensity,background_stdev)
            #1 if snr>2 else 0
            
            #Peak ID - (x,y); Full File path, Fname, FrameNbr, peakStatus (1,0),
            # avgPeakIntensity (avg of 3x3 pxl), avgBgIntensity (16 pxl around 3x3), stdev of the
            # Background
            self.allPeaks_dict[(x,y)].append([f,f.split('/')[-1],frameNbr,peakStatus, avg_peak_intensity,avg_background_intensity,background_stdev])
        return True 

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and
    below supplied root directory.'''
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            allFiles.append(os.path.join(path,filename))
    return allFiles

def pickleIt(fname,obj):
    ofile = open(fname,'w')
    pickle.dump(obj,ofile)
    ofile.close()


def getPeakCoords(fname):
    allPeakInfo = list()
    frame = StarterFrame(fname)
    results = frame.find_starterPeaks()
    allPeakCoords = [frame.getPeakXY(peak) for peak in results]
    pklFname = frame.picklePeakCoords(allPeakCoords)
    return allPeakCoords, pklFname

def getPeakTraceInfo(allPeakCoords,fname):
    traceFiles = TiffStack(fname,allPeakCoords)
    allTiffs = traceFiles.allTiffs
    allTiffs.sort()
    frameNbr = 0
    for frameTiff in allTiffs:
        frameNbr+=1
        traceFiles.allPeakTraceIntensity(frameTiff,frameNbr)
    ofname = fname+'.allPeaks_trace.dict.pkl'
    ofile = open(ofname,'w')
    pickle.dump(traceFiles.allPeaks_dict,ofile)
    ofile.close()
    return traceFiles.allPeaks_dict

def getPeakInttraces(pklFname):
    pklFname= (
        '/project2/marcotte/boulgakov/microscope/2014-July/2014-07-27/AS2_Atto647NSE_2aM_16hWash_TFA1_647_Wat_50Perc_trace002_flds002/AS2_Atto647NSE_2aM_16hWash_TFA1_647_Wat_50Perc_trace002_flds002t0001.tif.allPeaks_trace.dict.pkl')    
    peakInfo = pickle.load(open(pklFname))
    """
#Peak ID - (x,y); Full File path, Fname, FrameNbr, peakStatus (1,0),
# avgPeakIntensity (avg of 3x3 pxl), avgBgIntensity (16 pxl around 3x3), stdev
# of the background                                                          
    """
    ofname = pklFname+'.peakIntensityTraces.csv'
    ofile = open(ofname,'w')
    for k,v in peakInfo.items():
        frameNbr, peakStatus,peakInt,bgInt = v[2],v[3],v[4],v[5]





def analysePeaks():
    import shelve

    pklFname= (
        '/project2/marcotte/boulgakov/microscope/2014-July/2014-07-27/AS2_Atto647NSE_2aM_16hWash_TFA1_647_Wat_50Perc_trace002_flds002/AS2_Atto647NSE_2aM_16hWash_TFA1_647_Wat_50Perc_trace002_flds002t0001.tif.allPeaks_trace.dict.pkl')    
    peakInfo = pickle.load(open(pklFname))
    """
#Peak ID - (x,y); Full File path, Fname, FrameNbr, peakStatus (1,0),
# avgPeakIntensity (avg of 3x3 pxl), avgBgIntensity (16 pxl around 3x3), stdev
# of the background
    """
    traceDict = collections.defaultdict()
    # PBleaching rates
    totalPeaks = list() #Total active peaks/frame acrorss stack
    avgPeakIntensity = list() # Avg peak Intensity/frame across stack
    avgBGIntensity = list() # Avg background intensity/frame across stack
    for k, val in peakInfo.items():
        for v in val:
            fname,peakStatus,peakInt,bgInt = v[1],v[3],v[4],v[5]
            print fname, peakStatus, peakInt, bgInt

    sys.exit(1)
       


        





def main():
    sourceDir = '/project2/marcotte/boulgakov/microscope/2014-July/2014-07-27'
    allStarterTiffs = locate('*01.tif',sourceDir)
    firstFrameImages = [f for f in allStarterTiffs if
                        re.search(r't[0]*01.tif',f)] #all time traces have a t###01.tif
    
    for fname in firstFrameImages:
        print fname
        allPeakCoords,pklFname = getPeakCoords(fname)
        # After the Peak Coords are generated 
        allPeaks_dict = getPeakTraceInfo(allPeakCoords,fname)





def test2():
    fname = '/project2/marcotte/boulgakov/microscope/2014-July/2014-07-27/AS2_Atto647NSE_2aM_16hWash_647_OSS2_50Perc_trace002_flds002/AS2_Atto647NSE_2aM_16hWash_647_OSS2_50Perc_trace002_flds002t0001.tif.pkl'
    allPeakCoords = pickle.load(open(fname))
    traceFiles = TiffStack(fname,allPeakCoords)
    allTiffs = traceFiles.allTiffs
    allTiffs.sort()
    frameNbr = 0
    for frameTiff in allTiffs:
        frameNbr+=1
        traceFiles.allPeakIntensity(frameTiff,frameNbr)
    ofname = fname+'.dict.pkl'
    
def test3():
    fname = (
    '/project2/marcotte/boulgakov/microscope/2014-July/2014-07-27/AS2_Atto647NSE_2aM_16hWash_647_OSS2_50Perc_trace002_flds003/AS2_Atto647NSE_2aM_16hWash_647_OSS2_50Perc_trace002_flds003t0001.tif.pkl'
    )

    print fname
    
    allPeakCoords = pickle.load(open(fname))
    allPeaks_dict = getPeakTraceInfo(allPeakCoords,fname) 




def test():
    inputTIFF = sourceDir + '/' + imTIFF
    outputPNG = inputTIFF + '.png'
    cmd = "convert " + inputTIFF + " " + outputPNG
    call(cmd.split(),shell=False)
    img = scipy.misc.imread(outputPNG) 
    correlation_matrix=np.array([[-5935,-5935,-5935,-5935,-5935],
                             [-5935,8027,8027,8027,-5935],
                             [-5935,8027,30742,8027,-5935],
                             [-5935,8027,8027,8027,-5935],
                             [-5935,-5935,-5935,-5935,-5935]])

    print correlation_matrix
    results = pflibMP.fit_peaks([img],correlation_matrix = correlation_matrix)
    print results[0][1][1]

#allPeaks_dict = collections.defaultdict(list)
#test()
#main()
analysePeaks()
#test3()

#itest2()

#!/home/boulgakov/localpython2/bin/python
#! /home/jaggu/anaconda/bin/python2.7
#!/home/boulgakov/localpython2/bin/python


"""
This programme does the image analysis of time trace files only. First it uses
Alexander's peak identifying algorithm on the first frame to identifying peaks.
Then (a) I measure the intensity of the 5x5px area of the peaks in the other
frames. I count whether the peak is still on based on the background intensity
of that frame. 
"""
import sys
sys.path.append('/project/marcotte/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis')
import scipy.misc
import numpy as np
import os
from subprocess import call
import pflibMP
import cv2
import fnmatch
import re
import cPickle as pickle
import collections
import shelve
import time

class Peaks(object):
	
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

	def find_Peaks(self):
		self.outputPNG = self.convertPNG()
	        img = scipy.misc.imread(self.outputPNG)
	        print "Finding Peaks in %s ...."%(self.outputPNG)
            #peak fitting result tuple
            #(0  1  2--  3------  4-------  5--------  6---  7--, 8--)
            #(h, w, fit, fit_err, subimage, fit_image, rmse, r_2, s_n)
            #agpy.gaussfit parameters --
            # 0-----  1--------  2  3  4------  5------  6---
            #(height, amplitude, x, y, width_x, width_y, rota)
        	results = pflibMP.fit_peaks([img],correlation_matrix=self.correlation_matrix)
	        allPeakResults = results[0][1]
        	return allPeakResults

	def getPeakXY(self,peakResult):
        	h,w = peakResult[0],peakResult[1]
	return (h,w)
    
	def picklePeakCoords(self,allPeakCoords):
       		ofname = self.inputTiff+'.pkl'
	        ofile = open(ofname,'w')
	       	pickle.dump(allPeakCoords,ofile)
	        ofile.close()
        return ofname

class ImageDir(object):
	    """
	    This is a class calling the individual tiff image. The idea is to be able to
	    get the statistics of the single molecules
	    """
	def __init__(self,sourceDir,subDir):
		self.sourceDir = sourceDir #Full path
	        self.subDir = subDir #Dir name
        	self.fullPathSubDir = os.path.join(self.sourceDir,self.subDir)
	        self.files = [os.path.join(self.fullPathSubDir,f) for f in
        	              os.listdir(self.fullPathSubDir) if f.endswith('.tif')]
    
	def getPeaks(self,fTif):
        	imgTif = Peaks(fTif)
	        allPeakResults = Peaks.find_Peaks(imgTif)
	        return allPeakResults


	
	

	





sourceDir = "/project2/marcotte/boulgakov/microscope/2014-Oct/2014-10-02"
subDir = "ThiolSilane_RhodBNHS_2fM_flds004"
fname = "ThiolSilane_RhodBNHS_2fM_flds004xy02.tif"

img = ImageDir(sourceDir,subDir)
fTiff = img.files[4]
peakResults = img.getPeaks(fTiff)
fpkl = fTiff+'.allPeaks.pkl'
pickle.dump(peakResults,open(fpkl,'w'))




print "Yes"

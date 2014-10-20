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
import numpy as np
import os
from subprocess import call
import cPickle as pickle
import collections

class Peaks(object):
    """
    This is the peak finder algorithm (Alexander's) applied to the first frame.
    The list of peaks at the coordinates and the corresponding 5x5 pixel area
    for each of the peaks (h,w) is returned eventually
    This is the peak finder algorithm (Alexander's) applied to the first frame.    """
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
        sys.path.append('/project/marcotte/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis')
        import scipy.misc
        import pflibMP

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

class ImagePeaks(object):
    """
    This class will do all the image statistics for peaks in this image. Passed
    is the processed list of all the peaks in the image. This may be by opening
    the pkl list with the tiff stack.
    @Parameter: fname (tiff name), list(peak result tuple)
            #peak fitting result tuple
            #(0  1  2--  3------  4-------  5--------  6---  7--, 8--)
            #(h, w, fit, fit_err, subimage, fit_image, rmse, r_2, s_n)
    """
    def __init__(self,fname,peakList):
        self.fname = fname
        self.peakList = peakList

    def nbrPeaks(self):
        return len(self.peakList)

    def calcPeakStats(self):
        allPeakMean, allBgMean,allSNR = list(), list(),list()
        peakIndexList = [(r,c) for r in [1,2,3] for c in [1,2,3]]
        for peak in self.peakList:
            subImage,snr = peak[4],peak[8]
            allSNR.append(snr)
            peakImage = subImage[[[1],[2],[3]],[1,2,3]]
            allPeakMean.append(np.mean(peakImage))
            bgVals = [v for rc,v in np.ndenumerate(subImage) if not rc in peakIndexList]
            allBgMean.append(np.mean(bgVals))
        return map(np.mean,[allPeakMean,allBgMean,allSNR]),map(np.std,[allPeakMean,allBgMean,allSNR])


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

####################### BEGINNING OF CLASS #####################################

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
        print "Processing STARTER"
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
    Creates a shelve file
    """
    
    def __init__(self,fname,allPeakCoords,nbrFrames=600):
        self.fname = fname #Name of 1st frame file
        self.allPeakCoords = allPeakCoords
        self.nbrFrames = nbrFrames #Total number of frames to consider
        self.traceDir = os.path.dirname(self.fname)
        self.allTiffs = [os.path.join(self.traceDir,f) for f in
                         os.listdir(self.traceDir) if f.endswith('.tif')]
        self.allTiffs.sort() #sort all the TIFF files
        # Default dict
        self.peak_frameInfo_dict = collections.defaultdict(list)
        self.frame_peakInfo_dict = collections.defaultdict(list)

    def makeShelveFiles(self):
        shelveDir = os.path.join(self.traceDir,'shelve')

        if not os.path.exists(shelveDir): os.makedirs(shelveDir)
        makeShelveFile = lambda f: shelve.open(os.path.join(shelveDir,f))
        frame_peakInfo_shelve, peak_frameInfo_shelve = map(makeShelveFile,
                                                       ['frame_peakInfo_dict.shelve',
                                                        'peak_frameInfo_dict.shelve'])

        for k, v in self.peak_frameInfo_dict.items(): peak_frameInfo_shelve[k] = v
        for k, v in self.frame_peakInfo_dict.items(): frame_peakInfo_shelve[k] = v
        return shelveDir, frame_peakInfo_shelve, peak_frameInfo_shelve

    def runThruFrames(self):
        """ The function that routes the call of all functions to generate the
        relevant two dictionaries.
        """
        frameNbr = 0
        for tiffFile in self.allTiffs[0:self.nbrFrames]:
        #for tiffFile in self.allTiffs[0:10]:
            # tiffFile = Full Path of the tiff File
            print "Processing : %s"%(tiffFile)
            frameNbr += 1
            allPeakInfo_frame = self.allPeakTraceIntensity(tiffFile,frameNbr)
            self.frame_peakInfo_dict[tiffFile] = self.frameInfo(allPeakInfo_frame)

    def frameInfo(self,allPeakInfo):
        """ This function runs through the information of all peaks in that
        particular frame and obtains statistics about the Frame like number of
        active peaks, inactive peaks, avg peak intensity etc..
        """
            #allPeakInfo == [(x,y), FrameNbr, Peak Status (1,0), Full File path, Fname
            # avgPeakIntensity (avg of 3x3 pxl), peakIntensity stdev, 
            # avgBgIntensity (16 pxl around 3x3), stdev of the Background
            # snr],[..]...
        frameNbr, fname = allPeakInfo[0][1], allPeakInfo[0][4]
        activePeaks = list()
        inactivePeaks = list()

        def _meanCalc(peakInt_list):
            meanVal = lambda x: np.mean([item[x] for item in peakInt_list])
            return meanVal(0),meanVal(1),meanVal(2)

        for peak in allPeakInfo:
            peakStatus = peak[2]
            if peakStatus:
                activePeaks.append([peak[5],peak[7],peak[9]])
            else:
                inactivePeaks.append([peak[5],peak[7],peak[9]])

        nbrActivePeaks = len(activePeaks)
        nbrInactivePeaks = len(inactivePeaks)
            # frameInfo = [frameNbr, number of Active Peaks, number of Inactive Peaks,
            # (Active : mean Peak Intensity, mean Background Intensity, mean snr),
            # (Inactive : ....)
        return [frameNbr, nbrActivePeaks,nbrInactivePeaks,_meanCalc(activePeaks),_meanCalc(inactivePeaks)]

    def allPeakTraceIntensity(self,f,frameNbr):
        """ This function runs through all the peaks and then obtains the
        details for each peak in the assigned frame...
        """
        allPeakInfo_frame = list()
        img = cv2.imread(f,-1)
        
        def _isItPeak(avg_intensity,background_intensity,stdev_intensity):
            snr = (avg_intensity - background_intensity)/stdev_intensity
            if snr > 2: return ( 1, snr)
            else: return ( 0, snr)

        def _getPeakIntensity(peakAreaXY):
            peaklist = []
            for (i,j) in peakAreaXY: 
                peaklist.append(img[i,j])
            return np.mean(peaklist), np.std(peaklist)
        
        def _getBackgroundIntensity(bgAreaXY):
            bglist =  []
            for (i,j) in bgAreaXY:
                bglist.append(img[i,j])
            arr = np.array(bglist)
            return arr.mean(), arr.std()
      

        for (h,w) in self.allPeakCoords:
            sum_intensity = 0
            x,y = int(np.around(h)), int(np.around(w))
            peakID = 'ID_'+str(x)+'_'+str(y)
            if not (x in range(2,509)  and y in range(2,509)): continue; # Takes care of borders

            all_peakAreaXY = [(i,j) for i in xrange(x-2,x+3) for j in xrange(y-2,y+3)]
            peakAreaXY = [(i,j) for i in xrange(x-1,x+2) for j in xrange(y-1,y+2)]
            bgAreaXY = [(i,j) for (i,j) in all_peakAreaXY if not (i,j) in peakAreaXY]
            avg_peak_intensity,peak_stdev = _getPeakIntensity(peakAreaXY)
            avg_bg_intensity,bg_stdev = _getBackgroundIntensity(bgAreaXY)
            
            (peakStatus,snr) = _isItPeak(avg_peak_intensity,avg_bg_intensity,bg_stdev)
            #Peak ID - ID_x_y: [(x,y), FrameNbr, Peak Status (1,0), Full File path, Fname
            # avgPeakIntensity (avg of 3x3 pxl), peakIntensity stdev, 
            # avgBgIntensity (16 pxl around 3x3), stdev of the Background
            # snr
            peakInfo = (
                        [(x,y),frameNbr, peakStatus, f,f.split('/')[-1], 
                         avg_peak_intensity,peak_stdev,
                         avg_bg_intensity, bg_stdev, snr])
            allPeakInfo_frame.append(peakInfo)
            self.peak_frameInfo_dict[peakID].append(peakInfo)

        return allPeakInfo_frame

################## END OF CLASS #########################################                                      

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and
    below supplied root directory.'''
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            allFiles.append(os.path.join(path,filename))
    return allFiles

def test_findLastPeaks(dateStamp):
    """ FINDING LAST PEAKS ON THE FRAME """
    pathDir = os.path.join(sourceDir,dateStamp)
    allFinalTiffs = locate('*599.tif',pathDir)
    print allFinalTiffs




def findInitialPeaks(dateStamp):
    """ FINDING INITIAL PEAKS FUNCTION """

    def _getPeakCoords(fname):
        allPeakInfo = list()
        frame = StarterFrame(fname)
        results = frame.find_starterPeaks()
        allPeakCoords = [frame.getPeakXY(peak) for peak in results]
        pklFname = frame.picklePeakCoords(allPeakCoords)
        return allPeakCoords, pklFname

    pathDir = os.path.join(sourceDir,dateStamp)
    allStarterTiffs = locate('*01.tif',pathDir)
    firstFrameImages = [f for f in allStarterTiffs if re.search(r't[0]*01.tif',f)] #all time traces have a t###01.tif
    for fname in firstFrameImages:
        print "Processing Trace File : %s ...."%(fname)
        touchFile = os.path.join(os.path.split(fname)[0],'NOAUTO.txt')
        cmd = "touch " + touchFile
        call(cmd.split(),shell=False)
        allPeakCoords,pklFname = _getPeakCoords(fname)
        print len(allPeakCoords)

def runThroughFrames(fname):
    """ USES THE PKL FILE OF THE FIRST STARTER FRAME. THEN RUNS THROUGH THE REST
    OF THE TIFF FILE TO OBTAIN ALL THE PEAK AND FRAME INFORMATION """
    allPeakCoords = pickle.load(open(fname))
    traceFiles = TiffStack(fname,allPeakCoords)
    traceFiles.runThruFrames()
    shelveDir, frame_peakInfo_shelve, peak_frameInfo_shelve =  traceFiles.makeShelveFiles()
    print "Number of Frames : %d"%(len(frame_peakInfo_shelve))
    print "Number of Peaks Identified and Processed : %d"%(len(peak_frameInfo_shelve))

def makeBashScript(dateStamp):
    """ MAKES A TXT FILE THAT CONTAINS COMMAND LINE TO RUN THIS PYTHON PROGRAMME
    WITH THE 'ANALYSE' PARAMETER. IT PASSES THE IDENTIFIIED PKL FILE
    """
    def _chuck(lst,n):
        return [lst[i::n] for i in xrange(n)]
    pattern = "*.tif.pkl"
    pathDir = os.path.join(sourceDir,dateStamp)
    allPklFiles = locate(pattern,pathDir)
    bashDir = os.path.join(pathDir,'bashDir')
    if not os.path.exists(bashDir): os.makedirs(bashDir)
    chuckedList = _chuck(allPklFiles,4)
    for i,item in enumerate(chuckedList):
        ofname = os.path.join(bashDir,'ANALYSE_PYTHON_SCRIPTS.'+str(i)+'.sh')
        print ofname
        ofile = open(ofname,'w')
        for pklFile in item:
            cmd = ("python " 
                    "/project/marcotte/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis/pbKinetics.py "
                    "ANALYSE "+
                    pklFile + "\n"
                  )
            ofile.write(cmd)   


if __name__ == '__main__':
    [ARG, FILE] = sys.argv[1:]
    t0 = time.clock()
    sourceDir = '/project2/marcotte/boulgakov/microscope/2014-Oct'
    if ARG == 'STARTER':
        dateStamp = FILE 
        findInitialPeaks(dateStamp)
    elif ARG == 'ANALYSE': 
        fname = FILE
        print "In Analyse: PKL file"
        runThroughFrames(fname)
    elif ARG == 'TEST':
        print "Test"
        dateStamp = FILE
        test_findLastPeaks(dateStamp)
        print FILE
    elif ARG == 'GENERATE':
        dateStamp = FILE
        print "Generating a bash script to Analyse the pkl files"
        makeBashScript(dateStamp)
    else: SystemExit("Invalid Argument") 

    t1 = time.clock()
    print ("Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
          )

"""
##############################################
TEST CASES

def test_shelve(fname):
    db = shelve.open(fname)
    for k in db.keys()[0:10]:
        print k, len(db[k])
    db.close()

def runThroughFrames_test(fname):
    allPeakCoords = pickle.load(open(fname))
    traceFiles = TiffStack(fname,allPeakCoords)
    traceFiles.runThruFrames()
    shelveDir, f1, f2 =  traceFiles.makeShelveFiles()


def getPeakCoords(fname):
    allPeakInfo = list()
    frame = StarterFrame(fname)
    results = frame.find_starterPeaks()
    allPeakCoords = [frame.getPeakXY(peak) for peak in results]
    pklFname = frame.picklePeakCoords(allPeakCoords)
    return allPeakCoords, pklFname


def firstFrame_test():
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


##############################################
"""

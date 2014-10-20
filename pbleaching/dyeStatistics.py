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
import pflibMP
import numpy as np
import os
from subprocess import call
#import pflibMP
import cPickle as pickle
import collections
import time
from matplotlib import pyplot as plt



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
        self.tiffFiles = self.getFiles('.tif')
    
    def getFiles(self,pattern):
        fList = [os.path.join(self.fullPathSubDir,f) for f in
                 os.listdir(self.fullPathSubDir) if f.endswith(pattern)]
        return fList

    def getPeaks(self,fTif):
        imgTif = Peaks(fTif)
        allPeakResults = Peaks.find_Peaks(imgTif)
        return allPeakResults

    def analyseAllPeaks(self,pklFname):
        allPeakResults = pickle.load(open(pklFname))
        return allPeakResults

    def getAllIntensity(self,peakData):
        allIntensity = list()
        allSNR = list()
        for peak in peakData:
            subImage,snr = peak[4],peak[8]
            peakImage = subImage[[[1],[2],[3]],[1,2,3]]
            intensity = np.mean(peakImage)
            allIntensity.append(intensity)
            allSNR.append(snr)
        return allIntensity,allSNR

    def getPeakIntensity(self):
        # Gets the peak Intensity of all the peaks in all the files in the
        # subDir
        allIntensityList,allSNRList = list(),list()
        pklFiles = self.getFiles('allPeaks.pkl')
        for f in pklFiles:
            peakData = pickle.load(open(f))
            intensityList,snrList = self.getAllIntensity(peakData)
            allIntensityList.extend(intensityList)
            allSNRList.extend(snrList)
        return allIntensityList, allSNRList



def drawHistogram(data_dict,combineExpt=True,type='Intensity'):   
    fig = plt.Figure(figsize=( 5,5),dpi=300)
    if type is 'Intensity':  bins = range( 4000,30000,100)
    else: bins = np.arange(0,30,0.2)    
#    colorList = ['k','b','g','r','c','m','y','w']
    colorList = ['k','#552100','#370353','#050708','#AA7200','#FF8500','#901400']
    i=0
    for dye,dataDir in sorted(data_dict.items()):
        allVals = list()
        i+=1
        for subDir,val in dataDir:
            print dye,subDir,colorList[i]
            if not combineExpt: plt.hist(val,bins,label=dye,facecolor=colorList[i])
            else: allVals.extend(val)
        if combineExpt:
            if type is 'Intensity': plt.hist(allVals,bins,label=dye,color=colorList[i],histtype='stepfilled',alpha=0.7,linewidth=2,edgecolor='black')
            elif type is 'SNR': 
                plt.hist(allVals,bins,label=dye,facecolor=colorList[i],histtype='step',alpha=0.7,linewidth=2,edgecolor=colorList[i])
    name = 'Rhod_AF555_Atto647N_2_Water_1s.'+type+'.png'
    if type is 'Intensity':
        plt.xlabel('Intensity (Counts)')
        plt.ylabel('Frequency of Single Dye Molecules')
        plt.ylim( 0,750)
        plt.xlim( 4000,20000)
        plt.legend(loc='upper right')
        f = os.path.join(destDir,name)
        plt.savefig(f)
        plt.close()
    if type is 'SNR':
        plt.xlabel('SNR (illumina method)')
        plt.ylabel('Frequncy')
        plt.legend(loc='upper right')
        f = os.path.join(destDir,name)
        plt.savefig(f)
        plt.close()
    return True   

def processPeakFit(img):


    for fTif in img.tiffFiles:
        print "Processing %s ...."%(fTif)
        peakResults = img.getPeaks(fTif)
        fpkl = fTif + '.allPeaks.pkl'
        ofile = open(fpkl,'w')
        pickle.dump(peakResults,ofile)
        ofile.close()
    return img

def processDyeHistogram(dye_dict):
    allDyesINT_dict = collections.defaultdict(list)
    allDyesSNR_dict = collections.defaultdict(list)
    for dye,pathDirs in dye_dict.items():
        for path in pathDirs:
            sourceDir,subDir = os.path.split(path)
            dyeImg = ImageDir(sourceDir,subDir)
            print "Processing %s"%(subDir)
            intensityList, snrList = dyeImg.getPeakIntensity()
            allDyesINT_dict[dye].append([subDir,intensityList])
            allDyesSNR_dict[dye].append([subDir,snrList])

    drawHistogram(allDyesINT_dict,type='Intensity')
    drawHistogram(allDyesSNR_dict,type='SNR')




sourceDir = "/project2/marcotte/boulgakov/microscope/2014-July/2014-07-12"
subDir = "Glass_Alexa555SE_200nM_Overnite_561_40Perc_flds001"
destDir = "/project/marcotte/jagannath/projectfiles/compilationData/dyeStatistics/rhodamineComparisons"



dye_dict ={'Blank':
           ["/project2/marcotte/boulgakov/microscope/2014-Oct/2014-10-02/ThiolSilane2_Blank_flds005"],
          'RhodamineB':
           ["/project2/marcotte/boulgakov/microscope/2014-Oct/2014-10-02/ThiolSilane_RhodBNHS_2fM_flds004"],
          'Rhodamine':
           ["/project2/marcotte/boulgakov/microscope/2014-Oct/2014-10-02/ThiolSilane2_Blank_Wash_2pM_flds007"],
          'carboxyRhodamine':
           ["/project2/marcotte/boulgakov/microscope/2014-Oct/2014-10-09/ThiolSilane5_carbRhodB_2pM_flds004",
           "/project2/marcotte/boulgakov/microscope/2014-Oct/2014-10-09/ThiolSilane4_carbRhodB_20pM_flds002"],
           'Atto647N':
           ["/project2/marcotte/boulgakov/microscope/2014-Aug/2014-08-28/AS_Atto647N_2zM_DMFWash_647_flds"],
           'Alexa555':
           ["/project2/marcotte/boulgakov/microscope/2014-July/2014-07-12/Glass_Alexa555SE_200nM_Overnite_561_40Perc_flds001"]
          }

img = ImageDir(sourceDir,subDir)
#img = processPeakFit(img)
processDyeHistogram(dye_dict)

print "Completed"





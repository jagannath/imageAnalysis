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
#from Peaks import Peaks
import scipy.misc
import numpy as np
import os
from subprocess import call
import cPickle as pickle
import collections
import time
from matplotlib import pyplot as plt


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
        for fname,allPeaks in peakData.items():
            print "Retrieving peakData for %s"%(fname)
            for peak in allPeaks:
                subImage,snr = peak[4],peak[8]
                peakImage = subImage[[[1],[2],[3]],[1,2,3]]
                intensity = np.mean(peakImage)
                allIntensity.append(intensity)
                allSNR.append(snr)
        print len(allIntensity)
        return allIntensity,allSNR

    def getPeakIntensity(self):
        # Gets the peak Intensity of all the peaks in all the files in the
        # subDir
        allIntensityList,allSNRList = list(),list()
        pklFiles = self.getFiles('AllPeaks.pkl')
        for f in pklFiles:
            peakData = pickle.load(open(f))
            intensityList,snrList = self.getAllIntensity(peakData)
            allIntensityList.extend(intensityList)
            allSNRList.extend(snrList)
        return allIntensityList, allSNRList


def drawBoxPlots(data_dict,combineExpt=True,type='Intensity'):
    fig = plt.Figure(figsize=( 5,5),dpi=300)
    name = "Mean_fluorescence_Atto647N_"+type
    fpng = os.path.join(destDir,name+'.png')
    fsvg = os.path.join(destDir,name+'.svg')
    i = 0
    print type

    allData = list()
    allLabels = list()
    for dye,dataDir in sorted(data_dict.items()):
        allVals = list()
        i+=1
        for subDir,val in dataDir:
            if not combineExpt: print "check"
            else: allVals.extend(val)
        allData.append(allVals)
        allLabels.append(dye)
   
    print np.amax(allData[0])
    plt.boxplot(allData,notch=True,whis=2.0,patch_artist=True)
    plt.xticks(range(1,len(allData)+1),(allLabels),fontsize=8)
    plt.xlabel("Imaging Solvents")
    if type is 'Intensity':
        plt.ylim( 5000,25000)
        plt.ylabel("Mean Intensity of the peak")
    else:
        plt.ylim( 0,30)
        plt.ylabel("SNR of the peak")
    plt.savefig(fpng)
    plt.savefig(fsvg)
    plt.close()









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

def processDyeList(dye_dict):
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
    #drawHistogram(allDyesINT_dict,type='Intensity')
    #drawHistogram(allDyesSNR_dict,type='SNR')
    drawBoxPlots(allDyesINT_dict,type='Intensity')
    drawBoxPlots(allDyesSNR_dict,type='SNR')




sourceDir = "/project2/marcotte/boulgakov/microscope/2014-Oct/2014-10-26"
subDir = "Glass_Alexa555SE_200nM_Overnite_561_40Perc_flds001"
#destDir = "/project/marcotte/jagannath/projectfiles/compilationData/dyeStatistics/rhodamineComparisons"
destDir = "/project2/marcotte/jaggu/dataAnalysis/microscope1/2014-Oct/2014-10-26/graphs"
if not os.path.exists(destDir): os.makedirs(destDir)



"""dye_dict ={'Blank':
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
"""

dye_dict = {'Blank':
            [os.path.join(sourceDir,"GlassAS3_Blank_flds001")],
            'degasMeOH':
            [os.path.join(sourceDir,"GlassAS3_Atto647_20betaM_deGasMeOH_flds006")],
            'water':
            [os.path.join(sourceDir,"GlassAS3_Atto647_20betaM_DMFWash_Water_flds004")],
            'MeOH_\n10mMTrolox':
            [os.path.join(sourceDir,"GlassAS3_Atto647_20betaM_MeOH_10mMTrolox_flds")],
            'degasMeOH_\n2mMTrolox':
            [os.path.join(sourceDir,"GlassAS3_Atto647_20betaM_deGasMeOH_2mMTrolox_flds001")],
            'degasMeOH_\n10mMTrolox':
            [os.path.join(sourceDir,"GlassAS3_Atto647_20betaM_deGasMeOH_10mMTrolox_flds002")],
            'MeOH_\n2mMTrolox':
            [os.path.join(sourceDir,"GlassAS3_Atto647_20betaM_MeOH_2mMTrolox_flds007")],
            'MeOH':
            [os.path.join(sourceDir,"GlassAS3_Atto647_20betaM_MeOH_flds005")]
           }


processDyeList(dye_dict)

print "Completed"





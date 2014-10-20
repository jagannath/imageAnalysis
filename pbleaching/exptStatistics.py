#! /home/jaggu/anaconda/bin/python2.7

import os
import sys
from Peaks import ImagePeaks
import cPickle as pickle
import numpy as np
import time


def processSubDir(sourceDir,subDir):
    pklF = os.path.join(sourceDir,subDir,subDir+'.AllPeaks.pkl')
    subAnalysisDir = os.path.join(sourceDir,subDir,"subAnalysis")
    if not os.path.exists(subAnalysisDir): os.makedirs(subAnalysisDir)
    ofname = os.path.join(subAnalysisDir,subDir+'.ALL_FILES_SUMMARY.csv')

    ofile = open(ofname,'w')
    header = ["SUBDIR","TIFF FNAME","NBR PEAKS","MEAN PEAKS","MEAN PEAK BG","MEAN SNR","\n"]
    ofile.write("\t".join(header))
    with open(pklF) as ifile:
        allPeaks_dict = pickle.load(ifile)
    for fname, peakList in allPeaks_dict.items():
        print "Processing Tiff file ... %s"%(fname)
        img = ImagePeaks(fname,peakList)
        nbrPeaks = img.nbrPeaks()
        if nbrPeaks>1:
            [meanPeaks,meanBG,meanSNR],[stdevPeaks,stdevBg,stdevSNR] = img.calcPeakStats()
            ofile.write("\t".join([subDir,fname,str(nbrPeaks),str(meanPeaks),str(meanBG),str(meanSNR),"\n"]))
    ofile.close()


def checkSubDir(subDir):
    flag = True
    skipList = ['trace','stitch','zStep']
    for skip in skipList: 
        if skip in subDir: flag = False
    return flag

def summarize(sourceDir,analysisDir):
    monthStamp = analysisDir.split('/')[-1]
    if not os.path.exists(analysisDir): os.makedirs(analysisDir)
    ofLong = open(os.path.join(analysisDir,monthStamp+'_LONG_SUMMARY.csv'),'w')
    ofshortname = os.path.join(analysisDir,monthStamp+'_SHORT_SUMMARY.csv')
    ofShort = open(ofshortname,'w')
    header = (
        ["SUBDIR","MEAN NBR PEAKS","STDEV NBR PEAKS","MEANs PEAKS","STDEVs PEAKS",
         "MEANs PEAK BG","STDEVs PEAK BG","MEANs SNR","STDEVs SNR", "\n"]
        )
    ofShort.write("\t".join(header))
    print "Summarizing data for %s ..."%(sourceDir)
    for subDir in next(os.walk(sourceDir))[1]:
        print subDir
        sys.exit(1)
        if checkSubDir(subDir):
            nbrPeaksList,meansPeakList,meansBGList,meansSNRList = list(),list(),list(),list()
            fname = os.path.join(sourceDir,subDir,"subAnalysis",subDir+'.ALL_FILES_SUMMARY.csv')
            with open(fname) as ifile: readData = ifile.read()
            ofLong.write(readData)
            ofLong.write("\n")
        
            # Short summary
            data = readData.split('\n')[1:-1]# Header and last line is empty line
            for line in data:
                peakData = line.split('\t')[:-1]
                subDir,nbrPeaks,peakIntensity,BG,SNR = peakData[0],int(peakData[2]),float(peakData[3]),float(peakData[4]),float(peakData[5])
                nbrPeaksList.append(int(peakData[2]))
                meansPeakList.append(float(peakData[3]))
                meansBGList.append(float(peakData[4]))
                meansSNRList.append(float(peakData[5]))
        
            meanNbrPeaks, meansPeak, meansBG, meansSNR = map(np.mean,[nbrPeaksList,meansPeakList,meansBGList,meansSNRList])
            stdevNbrPeaks,stdevPeak, stdevBG, stdevSNR = map(np.std,[nbrPeaksList,meansPeakList,meansBGList,meansSNRList])
            ofShort.write("\t".join([subDir,str(meanNbrPeaks),str(stdevNbrPeaks),str(meansPeak),str(stdevPeak),str(meansBG),
                      str(stdevBG),str(meansSNR),str(stdevSNR),"\n"]))
        

    ofLong.close()
    ofShort.close()


def main(sourceDir,analysisDir):
    for subDir in next(os.walk(sourceDir))[1]:
        if checkSubDir(subDir):
            print "Processing %s ..."%(subDir)
            processSubDir(sourceDir,subDir)

    summarize(sourceDir,analysisDir)


if __name__ == '__main__':
    month = {'10':'Oct'}
    [ARG, dateStamp] = sys.argv[1:]
    pathDir = "/project2/marcotte/boulgakov/microscope"

    #dateStamp = "2014-10-15"
    monthStamp = "2014-"+month[dateStamp.split('-')[1]]
    sourceDir = os.path.join(pathDir,monthStamp,dateStamp)

    destDir = os.path.join("/project2/marcotte/jaggu/dataAnalysis/microscope1",monthStamp,dateStamp)

    t0 = time.clock()
    main (sourceDir,destDir)
    
    
    
    t1 = time.clock()
    print "Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0,
                                                       time.strftime("%d %b %Y  %H:%M:%S",time.localtime()))


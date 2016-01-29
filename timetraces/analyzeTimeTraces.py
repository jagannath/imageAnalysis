#! /home/jaggu/anaconda

"""
Analyzing time traces from a source Dir

"""
import sys
import cPickle as pickle
import os
import time
sys.path.append(os.path.join('/home/jaggu/marcotte_project/current/project/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis/timetraces'))
import picklePeakInfo 

def projectImages(sourceDir):
    print "maximum intensity projection done"

def pickle_peakInfo_allTraces(sourceDir,outputDir,peakID_files):
    for peakFilePath in peakID_files:
        t0 = time.ctime()
        traceDir, pklPath = picklePeakInfo.analyzeFrames(sourceDir,outputDir,peakFilePath)        
        print t0, time.ctime()


        
        
        
        
def main():
    sourceDir = '/home/jaggu/marcotte_project/boulgakov/microscope/2016-Jan/2016-01-20'
    peakSourceDir = '/home/jaggu/marcotte_project/boulgakov/microscope/2016-Jan/2016-01-20/projectedImages/first'
    outputDir = '/home/jaggu/marcotte_project/current/project2/jaggu/dataAnalysis/microscope1/2016-Jan/2016-01-20/output'
    peakID_files = [os.path.join(peakSourceDir,f) for f in os.listdir(peakSourceDir) if f.endswith('pkl')]
    pickle_peakInfo_allTraces(sourceDir, outputDir, peakID_files)


main()

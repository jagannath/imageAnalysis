#! /home/jaggu/anaconda/bin/python2.7
"""
This is a script that uses Alexander's method to identify the peaks. This script
runs it in a serial fashion and writes a dictionary pickle with 
filename: [peakResults];  
"""


from Peaks import Peaks
import os
import sys
import cPickle as pickle


def checkSubDir(subDir):
    flag = True
    skipList = ['trace','stitch','zStep']
    for skip in skipList: 
        if skip in subDir: flag = False
    return flag





def main(sourceDir):
    for subDir in next(os.walk(sourceDir))[1]:
        if checkSubDir(subDir):
            allPeaks_dict = dict()
            print "Processing %s ..."%(subDir)
            for f in next(os.walk(os.path.join(sourceDir,subDir)))[2]:
                if f.endswith('tif'): 
                    tifFname = os.path.join(sourceDir,subDir,f)
                    print "Finding peaks in %s "%(tifFname)
                    p = Peaks(tifFname)
                    res = p.find_Peaks()
                    allPeaks_dict[tifFname] = res
            pklFname = os.path.join(sourceDir,subDir,subDir+'.AllPeaks.pkl')
            print "Making Dictionary in pickle file : %s "%(pklFname)
            ofile = open(pklFname,'w')
            pickle.dump(allPeaks_dict,ofile)
            ofile.close()
            



if __name__ == '__main__':
    month = {'10':'Oct'}
    [ARG, dateStamp] = sys.argv[1:]
    pathDir = "/project2/marcotte/boulgakov/microscope"
    #dateStamp = "2014-10-15"
    monthStamp = "2014-"+month[dateStamp.split('-')[1]]
    sourceDir = os.path.join(pathDir,monthStamp,dateStamp)
    main (sourceDir)



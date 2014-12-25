#!/home/anaconda/python2.7

"""
The script aims to understand the intensity dynamics across the channels. I hope
to be able to reuse the script for understanding FRET dynamics etc. 
v1: Look for overlap in coordinates across the two channels after pairing up the
images. Then have a scatterplot of green intensity vs red intensity for the
identified peaks. 
"""

import os
import sys
import time
import fnmatch
import collections

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below supplied root directory.'''
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            allFiles.append(os.path.join(path, filename))
    return allFiles

class Date(object):
    def __init__(self,dateStamp):
        self.dateStamp = dateStamp
        self.month = self.month(dateStamp)
        self.DIR = "/project2/marcotte/boulgakov/microscope/2013-"+self.month+"/"+self.dateStamp
        self.analyseDIR = "/project2/marcotte/jaggu/dataAnalysis/2013-"+self.month+"/"+self.dateStamp
        self.checkDir(inputD = self.DIR,outputD = [self.analyseDIR]) #There can be many more analyseDir
    def month(self,dateStamp):
        monthName = \
        {'01':'Jan','02':'Feb','03':'March','04':'April','05':'May','06':'June','07':'July','08':'Aug','09':'Sept','10':'Oct','11':'Nov','12':'Dec'}
        return monthName[dateStamp.split('-')[1]]
    def checkDir(self,**kwargs):
        if not os.path.exists(kwargs['inputD']): raise SystemExit("Input Directory doesn't exist")
        [os.makedirs(d) for d in kwargs['outputD'] if not os.path.exists(d)]


class MultipleChannels(Date):
    """
    >> f = MultipleChannels(fname) #fname = *FRAME*CHANNEL COMPARISON.csv
    >> f.counts
    (#,#,#..) corresponds to channel 0, 1, 2.. 
    """
    def __init__(self,fname):
        self.fname = fname
        with open(fname) as f:
            self.header = f.readlines()[0]
            self.lines = f.readlines()[1:]
        self.counts = tuple()
    def getCounts(self): # Parses lines to extract channel information
        print header
        print len(self.lines)
        sys.exit(1)


def main():
    DIR = "/project2/marcotte/boulgakov/microscope/2013-"+month+dateStamp

    allFiles = 
    




if __name__ == '__main__': 
    
    t0 = time.clock()
    main()
    t1 = time.clock()
    print "Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))

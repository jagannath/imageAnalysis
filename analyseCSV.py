#! /home/jaggu/anaconda/bin/python2.7

# AUTHOR : Jagannath S; 
# Date Created : 2013-06-21
# Date Last Modified : 2013-06-21

"""
The aim of the method is to analyse the output of Alexander's peak finding
algorithm for the images. The goal is to tidy up the analysis of the csv files.
Batch: for batch analysis of the csv files in the directory 
TimeTrace: Getting the counts of the peptides/dyes at different times and making
a directory

"""
import numpy as np
import os.path
import os
import fnmatch
import csv
import collections
import time
import re
import shutil
import sys
from random import sample, choice
from subprocess import call

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below supplied root directory.'''
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            allFiles.append(os.path.join(path, filename))
    return allFiles

class Date(object):
    def __init__(self,dateStamp,default=True):
        self.dateStamp = dateStamp
        self.month = self.month(dateStamp)
        self.year = self.dateStamp[0:4]
        if default: 
            self.DIR = "/project2/marcotte/boulgakov/microscope/"+self.year+"-"+self.month+"/"+self.dateStamp
            self.analyseDIR = "/project2/marcotte/jaggu/dataAnalysis/microscope1/"+self.year+"-"+self.month+"/"+self.dateStamp
        else:
            self.DIR = "/project2/marcotte/boulgakov/microscope/PACJ_walker/2013-"+self.month+"/"+self.dateStamp
            self.analyseDIR = "/project2/marcotte/jaggu/dataAnalysis/PACJ_walker/2013-"+self.month+"/"+self.dateStamp
        self.checkDir(inputD = self.DIR,outputD = [self.analyseDIR]) #There can be many more analyseDir
    def month(self,dateStamp):
        monthName = \
                {'01':'Jan','02':'Feb','03':'Mar','04':'April','05':'May','06':'June','07':'July','08':'Aug','09':'Sept','10':'Oct','11':'Nov','12':'Dec','04':'Apr'}
        return monthName[dateStamp.split('-')[1]]
    def checkDir(self,**kwargs):
        if not os.path.exists(kwargs['inputD']): raise SystemExit("Input Directory doesn't exist")
        [os.makedirs(d) for d in kwargs['outputD'] if not os.path.exists(d)]

class Batch(Date):
    '''Handles the data within the dateStamp folder. It has an ignore  '''
    def __init__(self,dateStamp):
        Date.__init__(self,dateStamp,default=option)
        self.output = str()
    def analyse(self):
        allCSVs = locate('*SUMMARY.csv',self.DIR)
        ignoreList = ["expt","Trace","Time","trace","time","aaron","stitch"]
        def ignored(f):                     
            for ig in ignoreList:            
                if re.findall(ig,f): return True 
            return False                       
        for fname in allCSVs:
            if not ignored(fname):
                with open(fname,'r') as f: self.output += f.read()
        return self.output
    def writeCSV(self):
        fname = self.analyseDIR+"/"+self.dateStamp+".Counts.csv"
        with open(fname,'w') as f: f.write(self.output)
        return fname
    def generateImages(self):
        image = Image(self.dateStamp)
        print self.DIR
        print "in images"
        image.convertPNG(dir=self.DIR,expt='BATCH')

class TimeTrace(Date):
    '''Handles the Time trace data when images were gathered at random positions
    at a given time interval. Needs a dateStamp as argument. Functions are
    analyse, writeCSV and generateImages '''
    def __init__(self,dateStamp):
        Date.__init__(self,dateStamp) 
    def findFolders(self,pattern):
        timeTraceFolders = [self.DIR+'/'+d for d in os.listdir(self.DIR) if re.findall(pattern,d)]
        return timeTraceFolders
    def analyse(self,d):#Analyses only one timeTrace folder                    
        self.traceDIR = d
        allCSVs = locate('peak*.csv',self.traceDIR)
        self.timeTraceDICT = collections.defaultdict(list)
        def getRow(fname):
            with open(fname) as f:
                return f.readlines()[1]
        for f in allCSVs:
            row = getRow(f).split('\t')
            name = row[1]
            exptNbr = int(name.split('/')[-3]) #This is name of the directory above the place where the nd2 file was saved or two levels above tif file
            replicate = int(name.split('/')[-1].split('_')[-1].split('.')[0])
            count = int(row[2])
            self.timeTraceDICT[exptNbr].append([replicate,count])
        return self.timeTraceDICT
    def writeCSV(self):
        DICT = self.timeTraceDICT
        replNbr = len(DICT.values()[0])
        fname = self.analyseDIR+'/'+self.dateStamp+'.TimeTraces.csv'
        with open(fname,'w') as f:
            writer = csv.writer(f,delimiter='\t',quoting=csv.QUOTE_NONE)
            header = ['EXPT NBR']
            header.extend(['COUNT '+str(i) for i in range(replNbr)])
            writer.writerow(header)
            for exptNbr,countList in DICT.items():
                row = [exptNbr]
                row.extend([count[1] for count in countList])
                writer.writerow(row)
        return fname
    def generateImages(self):
        image = Image(self.dateStamp)
        image.convertPNG(dir=self.traceDIR,expt='TIMETRACE')

class Trace(Date):
    def __init__(self,dateStamp):
        Date.__init__(self,dateStamp)
        self.output = str()
    def getTrace(self):
        allFiles = locate("*SUMMARY.csv",self.DIR)
        allPbleachFiles = [f for f in allFiles if 'trace' in f]
        for fname in allPbleachFiles:
            with open(fname,'r') as f: self.output += f.read() 
    def writeCSV(self):
        fname = self.analyseDIR+"/"+self.dateStamp+".Counts.trace.csv"
        with open(fname,'w') as f: f.write(self.output)
        return fname


class Channels(Date):
    def __init__(self,dateStamp):
        Date.__init__(self,dateStamp)
    def getCounts(self):
        allFrameFiles = locate("*FRAME*CHANNEL COMPARISON.csv",self.DIR)
        ofname = self.analyseDIR+"/"+self.dateStamp+".Channels.Counts.csv"
        ofile = open(ofname,'a')
        ofile.write('FileName\tChannel-1 Counts\tChannel-2 Counts\tBothChannels Counts\n')

        def _getLines(fname):
            with open(fname) as f:
                lines = f.readlines()[1:]
            return lines

        for fname in allFrameFiles:
            lines = _getLines(fname)
            ch1Count, ch2Count, chBothCount = 0, 0, 0
            for line in lines:
                ch1Status, ch2Status = line.split('\t')[2],line.split('\t')[16]
                if ch1Status == "Not Available": ch2Count += 1
                elif ch2Status == "Not Available": ch1Count +=1
                else: chBothCount +=1
            ofile.write(fname+'\t'+str(ch1Count)+'\t'+str(ch2Count)+'\t'+str(chBothCount)+'\n')

        ofile.close()
    


class Image(Date):
    '''Handles Images. Currently the main function is to convert the tif file to
    a png of an appropriate format.'''
    def __init__(self,dateStamp):
        Date.__init__(self,dateStamp,default=option)
    def __call__(self):
        convertPNG()
    def convertPNG(self,dir,expt='BATCH'):
        ''' This converts one random tif file in the directory to a png
        according to specifications '''
        DIRS = collections.defaultdict(list)
        allTIFs = locate('*.tif',dir)
        [DIRS[f.split('/')[-2]].append(f) for f in allTIFs]
        convertList = [choice(flist) for flist in DIRS.values()]
        if expt=='TIMETRACE':convertList = sample(convertList,5)
        if expt=='BATCH': convertList = [im for im in convertList if not ('TimeTraceFiles' in im or 'stitch' in im)]
        for f in convertList: 
            destFile = self.analyseDIR+'/'+f.split('/')[-1]+'.png'
            cmd = "convert -contrast-stretch 0.15x0.05% "+f+" "+destFile
            print destFile
            call(cmd.split(),shell=False)
        
def main(dateStamp,ARG):
    if ARG == 'BATCH':
        batch = Batch(dateStamp)
        print batch.DIR
        batch.analyse()
        batch.writeCSV()
        batch.generateImages()

    if ARG == 'TIMETRACE':
        trace = TimeTrace(dateStamp)
        [dir] = trace.findFolders('TimeTraceFiles')
        trace.analyse(dir)           
        trace.writeCSV()
        trace.generateImages()
    if ARG == 'TRACE':
        pbleach = Trace(dateStamp)
        pbleach.getTrace()
        pbleach.writeCSV()
    if ARG == 'IMAGES':
        image = Image(dateStamp)
    if ARG == 'CHANNELS':
        ch = Channels(dateStamp)
        ch.getCounts()
        #for expt in ['BATCH','TIMETRACE']: image.convertPNG('expt'=expt)

def getArguments(argv):
    def invalidArg(): raise SystemExit("Incorrect Argument : python analyseCSV.py [dateStamp in format yyyy-mm-dd] [BATCH,TIMETRACE,IMAGES or TEST]")
    try: 
        [dateStamp,ARG] = sys.argv[1:]
        assert re.match(r"\d{4}-\d{2}-\d{2}$",dateStamp) and ARG in ['BATCH','TIMETRACE','IMAGES','TEST','CHANNELS','TRACE'] 
        return (dateStamp,ARG)
    except (ValueError,NameError, AssertionError): 
        print invalidArg()
    else:
        raise
    
if __name__ == '__main__': 
    dateStamp,ARG = getArguments(sys.argv[1:])
    #option = False
    option = True
    
    t0 = time.clock()
    main(dateStamp,ARG)
    t1 = time.clock()
    print "Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
   

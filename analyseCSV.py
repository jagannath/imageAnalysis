#!/usr/bin/python

# AUTHOR : Jagannath S; 8th May 2013

"""
The aim of this script is to use Alexander's algorithm to walk through every Tiff file in the directory - including subdirectory. As each file is being processed, a png file and a png file with the peaks identified is output in the same directory. 
An output file (csv file) is generated DIR.csv. This will have subdirectory as the row and the columns will be counts of the files in the subdirectory. The output function will be made modular to encompass other ways to output. Right now (all files within a directory must be averaged).  
I have extended the script to now include counting peptides in files within a directory (including subdirectory) and generate a txt file with the file names and peptide counts.
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
from random import choice
from subprocess import call

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
    supplied root directory.'''
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
	for filename in fnmatch.filter(files, pattern):
	    allFiles.append(os.path.join(path, filename))
    return allFiles

''' TIME TRACE EXPTS  '''
def getCountRow(f):
    ifile = open(f,'r')
    row = ifile.readlines()[1] # This is typically the important line; The first is the header
    ifile.close()
    return row.split('\t') # In earlier algorithms it was ',' as a delimiter

def writeTimeCSV(fname,timeDICT):
    '''The Dict has Expt number, replicates, counts and SNR. Outputing a csv file '''
    ofile = open(fname,'wb')
    writer = csv.writer(ofile,delimiter='\t',quoting=csv.QUOTE_NONE)
    header = ["EXPERIMENT NUMBER","FIELD1","COUNTS1","SNR1","FIELD2","COUNTS2","SNR2"]
    writer.writerow(header)
    for k, v in timeDICT.items():
        row = [k]
        exptNbr = k
        for flds in sorted(v):
            fld,[count,snr] = flds
	    row.extend([fld,count,snr])
        writer.writerow(row)
    ofile.close()

def timeTraceExpts(dateStamp,DIR,fID):
    '''The Time trace directory is indicated. It runs through all the CSV computed earlier beginning as peak*.csv (1 per directory). It then retrieves the counts and the SNR and produces a final csv file        with the details of all the time trace expts. '''
    allCounts = list()
    allCSVs = locate('peak*.csv',DIR+fID)
    timeDICT = collections.defaultdict(list)
    outfile = DIR+fID+dateStamp +'_'+'_WashBuffer_AUTO.csv'
    for f in allCSVs:
        row = getCountRow(f)
        # Columns - Fname NbrPeaks MeanPeaks StdPeaks MeanSNR StdSNR
        name = row[0].replace('/project2/marcotte/jaggu/mnicroscopeFiles','') 
        exptNbr = int(name.split('/')[-3]) #This is name of the directory above the place where the nd2 file was saved or two levels above where the .tif is.
        replicate = int(name.split('/')[-1].split('_')[-1].split('.')[0])
        count,snr  = int(row[1]),float(row[4])
        print "Expt : %d %d Count = %d SNR = %f"%(exptNbr, replicate,count,snr)
        allCounts.append((exptNbr,name,count))
	timeDICT[exptNbr].append([replicate,[count,snr]])
    writeTimeCSV(outfile,timeDICT)        
    dest1 = '/project2/marcotte/jaggu/dataAnalysis/2013-May'
    dest2 = DIR + 'analysis'
    map(lambda dest: shutil.copy(outfile,dest),[dest1,dest2])

''' BATCH EXPERIMENTS '''
         
def didIgnore(f,ignoreList):
    flag = True
    for ig in ignoreList:
	if re.findall(ig,f):flag = False
    return flag

def convertPNG(dateStamp,mthName,DIR):
    ''' This converts one random tif file in the directory to a png and copies this to the analyse Directory '''
    allDIRs = collections.defaultdict(list)
    analyseDir = '/project2/marcotte/jaggu/dataAnalysis/2013-'+mthName + '/' + dateStamp 
    if not os.path.exists(analyseDir): os.makedirs(analyseDir)
    allTIFs = locate('*.tif',DIR)
    convertList = list()
    for f in allTIFs:
	dirs = f.split('/')[-2]
        allDIRs[dirs].append(f)
    
    convertList = [choice(flist) for dirs,flist in  allDIRs.items()]
    
    for f in convertList:
	destFile = analyseDir + '/' + f.split('/')[-1]+'.png' #Just the file name in the destination dir
	print destFile
	command = "convert -contrast-stretch 0.15x0.05% "+f+" "+destFile
        call(command.split(),shell=False)

def batchAnalysis(dateStamp,mthName,DIR):
    ''' The Batch analysis is done per directory given. It has an ignore list - time traces, expts (which is for timeTrace one) and then outputs a single csv file     
    '''
    allCSVs = locate('peak*.csv',DIR)
    allSubDir = collections.defaultdict(list)
    analyseDir = '/project2/marcotte/jaggu/dataAnalysis/2013-'+mthName + '/' + dateStamp 
    analyseDir2 ="project/marcotte/jagannath/projectfiles/Single Molecule Microscopy/dataAnalysis/2013-"+mthName+'/'+dateStamp
    map(lambda DIR: DIR if os.path.exists(DIR) else os.makedirs(DIR),[analyseDir,analyseDir2]) #Make two analysis Directory. 
    ignoreList = ["expt","Trace","Time","trace","time","aaron"]
    outfile = DIR+dateStamp+'_peptideCounts_AUTO.txt'
    ofile = open(outfile,'w')
    ofile.write('FILE \t COUNTS \n')
    for f in allCSVs:
	if didIgnore(f,ignoreList):
            ofile.write(open(f,'r').read())

    ofile.close()
    dest1 = analyseDir 
    dest2 = DIR + 'analysis'                                    
    dest3 = analyseDir2
    map(lambda dest: shutil.copy(outfile,dest),[dest1,dest2,dest3])

def test():
    count = countPeptides.simplecount(THRESHOLD,2,fname)
    assert count == 385

def main(dateStamp,ARG):
    monthName = {'01':'Jan','02':'Feb','03':'March','04':'April','05':'May','06':'June'}
    month = dateStamp.split('-')[1]
    if ARG == 'TIMETRACES':
    	DIR = "/project2/marcotte/boulgakov/microscope/2013-"+monthName[month]+"/"+dateStamp+"/"
        fID = findFolder('TimeTraceFiles',DIR)+"/"
        timeTraceExpts(dateStamp,DIR,fID)
    elif ARG == 'BATCH':
        DIR = "/project2/marcotte/boulgakov/microscope/2013-"+monthName[month]+"/" + dateStamp + '/'
        batchAnalysis(dateStamp,monthName[month],DIR)
        convertPNG(dateStamp,monthName[month],DIR)
    elif ARG == 'TEST':
	fname = 'photobleaching005t716.tif'
        test()
    else: print invalidArg()
   
def findFolder(pattern,DIR):
    '''This expects only one folder '''
    allDir = os.listdir(DIR) 
    for d in allDir:
	if re.findall(pattern,d):
	    return d

def invalidArg():
    return "Incorrect Argument : python analyseCSV.py [dateStamp in format yyyy-mm-dd] [BATCH,TIMETRACE or TEST]"
    
if __name__ == '__main__': 
    #ARG = 'BATCH'
    #ARG = 'TIMETRACES'
    #ARG = 'TEST'
    try: [dateStamp,ARG] = sys.argv[1:]
    except NameError or ValueError: print invalidArg()
    try: assert re.match(r"\d{4}-\d{2}-\d{2}$",dateStamp) and ARG in ['BATCH','TIMETRACES','TEST']
    except AssertionError: print invalidArg()

    t0 = time.clock()
    main(dateStamp,ARG)
    t1 = time.clock()
    print "Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
   

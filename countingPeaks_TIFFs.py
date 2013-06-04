#!/usr/bin/python

# AUTHOR : Jagannath S; 8th May 2013

"""
The aim of this script is to use Alexander's algorithm to walk through every Tiff file in the directory - including subdirectory. As each file is being processed, a png file and a png file with the peaks identified is output in the same directory. 
An output file (csv file) is generated DIR.csv. This will have subdirectory as the row and the columns will be counts of the files in the subdirectory. The output function will be made modular to encompass other ways to output. Right now (all files within a directory must be averaged).  
I have extended the script to now include counting peptides in files within a directory (including subdirectory) and generate a txt file with the file names and peptide counts.
"""


import numpy as np
import os.path
import countPeptides
#import countPeptides_astd as countPeptides
import fnmatch
import csv
import collections
import time
import re
import shutil
import sys

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
    supplied root directory.'''
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
	for filename in fnmatch.filter(files, pattern):
	    allFiles.append(os.path.join(path, filename))
    return allFiles

def timeTraceExpts():
    allCounts = list()
    allTIFFs = locate('*.tif',DIR)
    outfile = DIR+'/'+dateStamp +'_' + fID+'_WashBuffer_'+str(THRESHOLD)+'.csv'
    for f in allTIFFs:
	count = countPeptides.simplecount(THRESHOLD,f)
	name = f.replace('/project2/marcotte/jaggu/microscopeFiles','')
	if not re.findall('\#',name):
	    exptNbr = int(name.split('/')[-3]) #This is name of the directory above the place where the nd2 file was saved or two levels above where the .tif is.
	    replicate = int(name.split('/')[-1].split('_')[-1][:-4])
	    print "Expt : %d %d Count = %d"%(exptNbr, replicate,count)
	    allCounts.append((exptNbr,name,count))
	    arrayExpts[exptNbr,replicate] = int(count)

    np.savetxt(outfile,arrayExpts,fmt='%d',delimiter='\t')
    dest = '/project2/marcotte/jaggu/dataAnalysis/2013-May'
    shutil.copy(outfile,dest)
         
def batchAnalysis():
    # This function looks through all the TIFF ured_UTJJ1200pM_TimeTraceFiles
    allTIFFs = locate('*.tif',DIR)
    outfile = DIR+dateStamp+'_PeptideCounts.txt'
    ofile = open(outfile,'w')
    ofile.write('FILE \t COUNTS \n')
    for f in allTIFFs:
	if not (re.findall(r'expt',f) or re.findall(r'Trace',f) or re.findall(r'Time',f) or re.findall(r'trace',f)):
	    count = countPeptides.simplecount(THRESHOLD,f)
	    print "Processing %s: Counts = %d"%(f,count)        
	    ofile.write(f+'\t'+str(count)+'\n')
            
    ofile.close()
    dest = '/project2/marcotte/jaggu/dataAnalysis/2013-May'
    shutil.copy(outfile,dest)
    
    
def test():
    count = countPeptides.simplecount(THRESHOLD,2,fname)
    assert count == 385
    

def main():
    if ARG == 'BATCH': batchAnalysis()
    elif ARG == 'TIMETRACES': timeTraceExpts()
    elif ARG == 'TEST':test()
    else:
	print "Incorrect Argument- TEST, BATCH or TIMETRACES with appropriate directory and file indicated"
    
    
if __name__ == '__main__':
    #ARG = 'BATCH'
    #ARG = 'TIMETRACES'
    [ARG] = sys.argv[1:]
    #ARG = 'TEST'
    if ARG == 'TIMETRACES':
	dateStamp = '2013-05-20'
        fID = "AS-N2CuredCy5Dye_TimeTraceFiles"
	DIR = "/project2/marcotte/boulgakov/microscope/2013-May/"+dateStamp+"/"+fID
    elif ARG == 'BATCH':
	dateStamp = '2013-05-23'
	DIR = "/project2/marcotte/boulgakov/microscope/2013-May/" + dateStamp + '/'
    elif ARG == 'TEST':
	fname = 'photobleaching005t716.tif'
    
    THRESHOLD = 4
    arrayExpts = np.zeros((180,6))

    t0 = time.clock()
    main()
    t1 = time.clock()
    print "Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
   

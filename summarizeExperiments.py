#! /home/jaggu/anaconda/bin/python2.7

import os
import sys
import time
import collections
import numpy as np
import cPickle as pickle

sys.path.append('/project/marcotte/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis')
from commonFunctions import locate,makeDir

"""
Aimed at summarizing the processed image directory (in boulgakov/microscope/..)
. It now seems to contain multiple types of files, - trace; - flds; - tflds;
Appropriate types of summary will need to be done. 
"""

def summarize_fld(dateStamp,pathDir,destDir):
    """
    Summary of flds; This typically would be the count density/field;
    """
    def _skip(f):
        skip = ['offset','maxProj','proj','dict','tflds']
        for word in skip:
            if word in f: return False
        return True

    def _checkChannel(f):
        ch = None
        if '561' in f: ch = '561'
        if '647' in f: ch = '647'
        if 'c1.tif' in f: ch = '561' #Must take precautions
        if 'c2.tif' in f: ch = '647'
        return ch

    def _summary(ch,listToWrite,ofile):
        ofile.write('\n CHANNEL \t '+ch+'\n')
        for l in listToWrite:
            ofile.write('\t'.join(l)+'\n')
        return True


    out_long =  os.path.join(destDir,dateStamp+'_flds.long.summary.csv')
    out_short = os.path.join(destDir,dateStamp+'_flds.short.summary.csv')

    ofile_long = open(out_long,'w')
    header = '\t'.join(['FULL PATH','SUBDIR','CH','PEAK COUNT'])
    ofile_long.write(header+'\n')

    ofile_short = open(out_short,'w')
    header = '\t'.join(['FULL PATH','SUBDIR','CH','AVERAGE PEAK COUNT', 
                       'STDEV PEAK COUNT'])
    ofile_short.write(header)

    list_561 = list()
    list_647 = list()
    dict_561 = collections.defaultdict(list)
    dict_647 = collections.defaultdict(list)
    

    pattern = "*flds*.tif.png*.pkl"
    allfPkl = locate(pattern, pathDir)
    fList = [f for f in allfPkl if not _skip(f)]
    for f in fList:
        print "Loading file : %s"%(os.path.split(f)[1])
        with open(f) as ifile:
            peakInfo = pickle.load(ifile)
        nbrPeaks = len(peakInfo) #Just this for now. Simplest
        exptPath,fname = os.path.split(f)
        expt = os.path.split(exptPath)[1]
        ch = _checkChannel(fname)
        if ch == '561': 
            list_561.append([f,expt,ch,str(nbrPeaks)])
            dict_561[expt].append([exptPath,expt,nbrPeaks])
        elif ch == '647':
            list_647.append([f,expt,ch,str(nbrPeaks)])
            dict_647[expt].append([exptPath,expt,nbrPeaks])
        else: raise SystemExit("Something wrong with Channels")

    _summary('561',list_561,ofile_long)
    _summary('647',list_647,ofile_long)
    ofile_long.close()

    list_561 = list()
    list_647 = list()
    list_ch = [list(),list()]
    for idx,d in enumerate([dict_561,dict_647]):
        for expt,val in d.items():
            exptPath = val[0][0]
            nbrPeaksL = zip(*val)[2]
            avgPeaks,stdPeaks = np.mean(nbrPeaksL),np.std(nbrPeaksL)
            list_ch[idx].append([exptPath,expt,str(idx+1),str(avgPeaks),str(stdPeaks)])

    _summary('561',list_ch[0],ofile_short)
    _summary('647',list_ch[1],ofile_short)
    ofile_short.close()

    return True


if __name__ == '__main__':
    month = {'01':'Jan','02':'Feb','10':'Oct','11':'Nov','12':'Dec'}
    [ARG, dateStamp] = sys.argv[1:]
    sourceDir = "/project2/marcotte/boulgakov/microscope"
    
    yearStamp = dateStamp.split('-')[0]
    monthStamp = yearStamp+"-"+month[dateStamp.split('-')[1]]
    pathDir = os.path.join(sourceDir,monthStamp,dateStamp)
    
    destDir = os.path.join("/project2/marcotte/jaggu/dataAnalysis/microscope1",monthStamp,dateStamp)
    makeDir(destDir)
    t0 = time.clock()   

#    allDir = [subDir for subDir in next(os.walk(pathDir))[1] if 'trace' in subDir]

    if ARG == 'FLDS':
        summarize_fld(dateStamp,pathDir,destDir)
    
    t1 = time.clock()
    print "Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0,
                                                       time.strftime("%d %b %Y  %H:%M:%S",time.localtime()))

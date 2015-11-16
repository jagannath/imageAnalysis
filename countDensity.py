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

def testImport():
    return True

def getMeanIntensity(db):
    allPEAKIntensity = list()
    for k,v in db.items():
        subImage = v[7]
        peak = subImage[1:4,1:4]
        peakInt = np.sum(peak)
        allPEAKIntensity.append(peakInt)

    meanIntensity = np.mean(allPEAKIntensity)
    return meanIntensity

def summarize_fld(dateStamp,pathDir,destDir,typeVal='count'):
    """
    Summary of flds; This typically would be the count density/field;
    """
    def _skip(f):
        skip = ['offset','maxProj','proj','dict','tflds']
        for word in skip:
            if word in f: return False
        return True

    def _checkChannel(f):
        ch = 'default'
        if '561' in f: ch = '561'
        if '647' in f: ch = '647'
        if '561_647' in f: ch = '561_647'
        if 'c1.tif' in f: ch = '561' #Must take precautions
        if 'c2.tif' in f: ch = '647'
        if 'c3.tif' in f: ch = '561_647'
        return ch

    def _summary(ch,listToWrite,ofile):
        ofile.write('\n CHANNEL \t '+ch+'\n')
        for l in listToWrite:
            ofile.write('\t'.join(l)+'\n')
        return True
    
    if typeVal == 'count':
        out_long =  os.path.join(destDir,dateStamp+'_flds.long.count.summary.csv')
        out_short = os.path.join(destDir,dateStamp+'_flds.short.count.summary.csv')
    if typeVal == 'intensity':
        out_long =  os.path.join(destDir,dateStamp+'_flds.intensity.long.summary.csv')
        out_short = os.path.join(destDir,dateStamp+'_flds.intensity.short.summary.csv')

    ofile_long = open(out_long,'w')
    if typeVal == 'count':  header = '\t'.join(['FULL PATH','SUBDIR','CH','PEAK COUNT'])
    if typeVal == 'intensity': header = '\t'.join(['FULL PATH', 'SUBDIR','CH','MEAN PEAK INTENSITY'])
    ofile_long.write(header+'\n')

    ofile_short = open(out_short,'w')
    if typeVal == 'count':  header = '\t'.join(['FULL PATH','SUBDIR','CH','AVERAGE PEAK COUNT', 'STDEV PEAK COUNT'])
    if typeVal == 'intensity': header = '\t'.join(['FULL PATH','SUBDIR','CH', 'AVERAGE of MEAN PEAK INTENSITY','STDEV MEAN INTENSITY'])
    ofile_short.write(header)

    list_561 = list()
    list_647 = list()
    list_561_647 = list()
    list_default = list()
    dict_561 = collections.defaultdict(list)
    dict_647 = collections.defaultdict(list)
    dict_561_647 = collections.defaultdict(list)
    dict_default = collections.defaultdict(list)
    
    pattern = "*fld*.tif.png*.pkl"
    allfPkl = locate(pattern, pathDir)
    print pattern, pathDir
    fList = [f for f in allfPkl if not _skip(f)]
	
    for f in fList:
        print "Loading file : %s"%(os.path.split(f)[1])
        with open(f) as ifile:
            peakInfo = pickle.load(ifile)
        if typeVal == 'count':  colVal = len(peakInfo) #Just this for now. Simples
        if typeVal == 'intensity': colVal = getMeanIntensity(peakInfo)
        exptPath,fname = os.path.split(f)
        expt = os.path.split(exptPath)[1]
        ch = _checkChannel(fname)
        if ch == '561': 
            list_561.append([f,expt,ch,str(colVal)])
            dict_561[expt].append([exptPath,expt,colVal])
        elif ch == '647':
            list_647.append([f,expt,ch,str(colVal)])
            dict_647[expt].append([exptPath,expt,colVal])
        elif ch == '561_647':
            list_561_647.append([f,expt,ch,str(colVal)])
            dict_561_647[expt].append([exptPath,expt,colVal])
        elif ch == 'default':
            list_default = list()
            dict_default[expt].append([exptPath,expt,colVal])
        else: raise SystemExit("Something wrong with Channels")

    _summary('561',list_561,ofile_long)
    _summary('647',list_647,ofile_long)
    _summary('561_647',list_561_647,ofile_long)
    _summary('default',list_default,ofile_long)
    ofile_long.close()

    list_561 = list()
    list_647 = list()
    list_561_647 = list()
    list_ch = [list(),list(),list()]
    list_default = list()
    
    for idx,d in enumerate([dict_561,dict_647,dict_561_647]):
        for expt,val in d.items():
            exptPath = val[0][0]
            colValL = zip(*val)[2]
            avgPeaks,stdPeaks = np.mean(colValL),np.std(colValL)
            list_ch[idx].append([exptPath,expt,str(idx+1),str(avgPeaks),str(stdPeaks)])

    _summary('561',list_ch[0],ofile_short)
    _summary('647',list_ch[1],ofile_short)
    _summary('561_647',list_ch[2],ofile_short)
    ofile_short.close()

    return True


if __name__ == '__main__':
    month = {'01':'Jan','02':'Feb','03':'Mar','04':'Apr','06':'June','07':'July','10':'Oct','11':'Nov','12':'Dec'}
    [ARG, dateStamp] = sys.argv[1:]

    yearStamp = dateStamp.split('-')[0]
    monthStamp = yearStamp+"-"+month[dateStamp.split('-')[1]]

    microscope = 1 
    if microscope is 1:                                                                                                                     
        sourceDir = "/project/boulgakov/microscope"                                                                                         
        pathDir = os.path.join(sourceDir,monthStamp,dateStamp)                                                                              
        destDir = os.path.join("/project/current/project2/jaggu/dataAnalysis/microscope1",monthStamp,dateStamp)     
    elif microscope is 2:                                                                                                                   
        sourceDir = "/project/boulgakov/microscope2/jagannath/rawFiles"                                                                     
        pathDir = os.path.join(sourceDir,monthStamp,dateStamp)                                                                              
        destDir = os.path.join("/project/current/project2/jaggu/dataAnalysis/microscope2",monthStamp,dateStamp)     
    
    elif microscope is 3:
        sourceDir = "/project/boulgakov/microscope3/rawFiles/jagannath"
        pathDir = os.path.join(sourceDir,monthStamp,dateStamp)
        destDir = os.path.join("/project/current/project2/jaggu/dataAnalysis/microscope3",monthStamp,dateStamp)
    
    else:                                                                                                                                   
        raise SystemExit("Incorrect microscope number. Quiting ..")                                                                         
    
    makeDir(destDir)
    t0 = time.clock()   

#    allDir = [subDir for subDir in next(os.walk(pathDir))[1] if 'trace' in subDir]

    val = 'count'
    #val = 'intensity'

    if ARG == 'FLDS':
        summarize_fld(dateStamp,pathDir,destDir,typeVal=val)
    
    t1 = time.clock()
    print "Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0,
                                                       time.strftime("%d %b %Y  %H:%M:%S",time.localtime()))

#! /home/jaggu/anaconda/bin/python2.7

import os
import sys
import time
import collections
import numpy as np

sys.path.append('/project/marcotte/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis')
from commonFunctions import locate,makeDir






def summarize_photobleaching(dateStamp,pathDir,destDir):
    def _summary(ofile, listLines):
        norm_trueList_561,norm_trueList_647,norm_falseList_561,norm_false_647 = listLines

        ofile.write('\n NORMALIZATION OF PHOTOBLEACHING CURVE \n')
        ofile.write('\n 561 CHANNEL \n')
        for line in norm_trueList_561: ofile.write(line+'\n')
        ofile.write('\n 647 CHANNEL \n')
        for line in norm_trueList_647: ofile.write(line+'\n')

        ofile.write('\n NO NORMALIZATION OF PHOTOBLEACHING \n')
        ofile.write('\n 561 CHANNEL \n')
        for line in norm_falseList_561: ofile.write(line+'\n')
        ofile.write('\n 647 CHANNEL \n')
        for line in norm_falseList_647: ofile.write(line+'\n')

        ofile.close()
        print "Summary of Photobleaching written "
        return True

    def _calcMeans(expt,val):
        h1L,h2L,r2L,rmseL = zip(*val)
        h1_mean = str(np.mean(h1L))
        h2_mean = str(np.mean(h2L))
        h1_std = str(np.std(h1L))
        h2_std = str(np.std(h2L))
        r2L_mean = str(np.mean(r2L))
        r2L_std = str(np.std(r2L))
        rmseL_mean = str(np.mean(rmseL))
        rmseL_std = str(np.std(rmseL))
        line = \
            '\t'.join([expt,h1_mean,h2_mean,h1_std,h2_std,
                       r2L_mean,r2L_std,rmseL_mean,rmseL_std])
        return line

    def _addTrueDict(line):
        _expt,_h1,_h2,_r2,_rmse = [line.split('\t')[i] for i in[1,2,3,5,6]]
        expt = '_'.join(_expt.split('_')[:-1])
        expt_true_dict[expt].append([float(_h1),float(_h2),float(_r2),float(_rmse)])
    
    def _addFalseDict(line):
        print line,"yes"

        _expt,_h1,_h2,_r2,_rmse = [line.split('\t')[i] for i in[1,2,3,5,6]]
        expt = '_'.join(_expt.split('_')[:-1])
        expt_false_dict[expt].append([float(_h1),float(_h2),float(_r2),float(_rmse)])




# Destination directory is same as the pathdir for photobleaching. 
    out_long = os.path.join(destDir,dateStamp+'_photobleaching.long.summary.csv')
    out_short = os.path.join(destDir,dateStamp+'_photobleaching.short.summary.csv')
    
    ofile_short = open(out_short,'w')
    
    ofile_long = open(out_long,'w')
    header_long = '\t'.join(['FULL PATH','SUBDIR',
                         'HALF-LIFE-1','HALF-LIFE-2',
                         'PARAMS','r_2','RMSE','NORMALIZATION','\n'])
    ofile_long.write(header_long)

    norm_trueList_561 = list()
    norm_trueList_647 = list()
    norm_falseList_561 = list()
    norm_falseList_647 = list()
    expt_true_dict = collections.defaultdict(list)
    expt_false_dict = collections.defaultdict(list)

    pattern = '*.pbleaching.csv'
    all_pbleaching_files = locate(pattern, destDir)
    """
    path = \
        '/project2/marcotte/jaggu/dataAnalysis/microscope1/2014-Dec/2014-12-22/images/141222_DiAS1_PR011_200fM_100uMTXMeOH_561_500ms_40C_trace001/_frameAnalysis/pbleaching'
    f = \
        os.path.join(path,'141222_DiAS1_PR011_200fM_100uMTXMeOH_561_500ms_40C_trace001.two_exp.pbleaching.csv')
    all_pbleaching_files = [f]
    """
    for f in all_pbleaching_files:
        ifile = open(f,'r')
        lines = ifile.readlines()
        for line in lines:
            line = line.rstrip()
            if (len(line.split('\t')))==8:
                print line
                if line.split('\t')[-1] == 'True':
                   if '561' in line:
                       norm_trueList_561.append(line)
                       _addTrueDict(line)
                   if '647' in line: 
                       _addTrueDict(line)                       
                       norm_trueList_647.append(line)
                elif line.split('\t')[-1] == 'False':
                   if '561' in line:
                       print line
                       _addFalseDict(line)                                              
                       norm_falseList_561.append(line)
                   if '647' in line: 
                       _addFalseDict(line)                                                                     
                       norm_falseList_647.append(line)
                else: pass
    
    listLines = [norm_trueList_561,norm_trueList_647,norm_falseList_561,norm_falseList_647]
    _summary(ofile_long,listLines)

    
    header_short = '\t'.join(['SUBDIR',
                              'HALF-LIFE-1','HALF-LIFE-2','STD-1','STD-2',
                              'r_2','STD-r_2','RMSE','STD_RMSE','\n'])
    ofile_short.write(header_short)

    norm_trueList_561 = list()
    norm_trueList_647 = list()
    norm_falseList_561 = list()
    norm_falseList_647 = list()

    for expt,val in expt_true_dict.items():
        line = _calcMeans(expt,val)
        if '561' in expt:
            norm_trueList_561.append(line)
        if '647' in expt:
            norm_trueList_647.append(line)

    for expt,val in expt_false_dict.items():
        line = _calcMeans(expt,val)
        if '561' in expt:
            norm_falseList_561.append(line)
        if '647' in expt:
            norm_falseList_647.append(line)

    listLines = [norm_trueList_561,norm_trueList_647,norm_falseList_561,norm_falseList_647]
    _summary(ofile_short,listLines)


    return True



if __name__ == '__main__':
    month = {'10':'Oct','11':'Nov','12':'Dec'}
    [ARG, dateStamp] = sys.argv[1:]
    sourceDir = "/project2/marcotte/boulgakov/microscope"

    monthStamp = "2014-"+month[dateStamp.split('-')[1]]
    pathDir = os.path.join(sourceDir,monthStamp,dateStamp)
    
    destDir = os.path.join("/project2/marcotte/jaggu/dataAnalysis/microscope1",monthStamp,dateStamp)
    t0 = time.clock()   

#    allDir = [subDir for subDir in next(os.walk(pathDir))[1] if 'trace' in subDir]

    if ARG == 'PHOTOBLEACH':
        summarize_photobleaching(dateStamp,pathDir,destDir)
    
    t1 = time.clock()
    print "Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0,
                                                       time.strftime("%d %b %Y  %H:%M:%S",time.localtime()))

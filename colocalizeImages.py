#! /home/jaggu/anaconda/bin/python2.7

import os
import sys
import time
import collections
import numpy as np
import cPickle as pickle
import itertools

sys.path.append('/project/marcotte/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis')
from commonFunctions import locate,makeDir


class Colocalize:
    def __init__(self,f1,f2):
        f1Peaks_dict = pickle.load(open(f1))
        f2Peaks_dict = pickle.load(open(f2))
        self.f1Pos_list,self.f2Pos_list = map(self.pos_list,[f1Peaks_dict,f2Peaks_dict])
        self.f1NbrPeaks, self.f2NbrPeaks = map(len, [f1Peaks_dict,
                                                     f2Peaks_dict])

    def pos_list(self,d):
        return d.keys()

    def wigglePos(self,pos):
        bound = 2
        nearbyPos = [(pos[0]+i,pos[1]+j) for i in range(-1*bound,bound+1,1) for j in range(-1*bound,bound+1,1)]
        return nearbyPos
    
    def checkPosOverlap(self):
        overlapPeaks = list() 
        pL1 = self.f1Pos_list
        pL2 = self.f2Pos_list
        for pos in pL1:
            for npos in self.wigglePos(pos):
                if npos in pL2:
                    overlapPeaks.append([pos, npos])
        return overlapPeaks

def colocalizeImages(dateStamp,pathDir,destDir):
    print dateStamp
    print pathDir
    print destDir
    processed_dir = os.path.join(destDir,"colocalizePKL")
    makeDir(processed_dir)

    pklFname = os.path.join(processed_dir,dateStamp + '_colocalize.pkl')
    ofile = open(pklFname,'w')
    subDir_colocal_dict = collections.defaultdict(list)

    allSubDir = [x[0] for x in os.walk(pathDir) if not x[0].endswith('PKL')]
    for subDir in allSubDir:
        print "Processing directory %s ..."%(subDir)
        f = str()
        allFile_pairs = list() 
        f = os.path.split(subDir)[1]
        allPkl = [os.path.join(subDir,f) for f in os.listdir(subDir) if
                  f.endswith('.pkl') and 'c1.tif.png' in f]
        for pklF in allPkl:
            file_list = list()
            suf = ".".join(pklF.split('.')[1:])
            pre = pklF.split('.')[0]
            if pre.endswith('c1'):
                for nbr in range( 1,5):
                    f = pre[:-2]+'c'+str(nbr)+'.'+suf
                    if os.path.isfile(f): file_list.append(f)
                    else: 
                        file_couple = list(itertools.combinations(file_list,2))
                        allFile_pairs.append(file_couple)
                        break
            else: raise SystemExit("Need atleast two files to compare")
            
        for subDir_list in allFile_pairs:
            for f1,f2 in subDir_list:
                colocal = Colocalize(f1,f2)
                overlapPeaks = colocal.checkPosOverlap()
                subDir_colocal_dict[subDir].append([f1,f2,colocal.f1NbrPeaks,colocal.f2NbrPeaks,len(overlapPeaks)])

    pickle.dump(subDir_colocal_dict,ofile)

    return subDir_colocal_dict


def writeToFile(colocal_dict,dateStamp,destDir):
    ofname1 = os.path.join(destDir,dateStamp+"_colocalize.long.csv")
    ofile1 = open(ofname1,'w')
    ofile1.write("\t".join(["SUBDIR","\n"]))
    ofile1.write("\t".join(["FILE 1","FILE 2","Nbr Peaks in File 1",
                           "Nbr Peaks in File 2", "Nbr Peaks Colocalized","\n"]))
    ofname2 = os.path.join(destDir,dateStamp+"_colocalize.short.csv")
    ofile2 = open(ofname2,'w')
    ofile2.write("\t".join(["FULLPATH","SUBDIR","CHANNEL 1--2","CHANNEL 1--3","CHANNEL 2--3",
                            "STDEV 1--2","STDEV 1--3", "CHANNEL 2--3","\n"]))
    
    od_colocal_dict = collections.OrderedDict(sorted(colocal_dict.items()))

    for subDir,val in od_colocal_dict.items():
        d = collections.defaultdict(list)
        ofile1.write(subDir+"\n")
        for col in val:
            ofile1.write("\t".join(map(str,col)))
            ofile1.write("\n")

            ch1 = col[0].split('.')[0][-2:]
            ch2 = col[1].split('.')[0][-2:]
            k = ch1 + '-'+ ch2
            d[k].append(col[-1])

        ofile1.write("\n")

        odDict = collections.OrderedDict(sorted(d.items()))
        avg_list = list()
        std_list = list()
        for chOverlap,v_list in odDict.items():
            avg_list.append(np.mean(v_list))
            std_list.append(np.std(v_list))
        
        newList = [subDir] + [os.path.split(subDir)[1]] + avg_list + std_list
        writeShortLine = "\t".join(map(str,newList))
        ofile2.write(writeShortLine)
        ofile2.write("\n")
   
    ofile1.close()
    ofile2.close()


def test():
    addThis = 'psfs_nr0wyr'
    fname = os.path.join(pathDir,subDir)
    xy = '10'
    f1 = fname+'xy'+str(xy)+'c1'+'.tif.png_'+addThis+'.pkl'
    f2 = fname+'xy'+str(xy)+'c2'+'.tif.png_'+addThis+'.pkl'
    f3 = fname+'xy'+str(xy)+'c3'+'.tif.png_'+addThis+'.pkl'
    
    for f1,f2 in [[f1,f2],[f2,f3],[f1,f3]]:
        colocal = Colocalize(f1,f2)
        overlapPeaks = colocal.checkPosOverlap()
        print len(overlapPeaks)



if __name__ == '__main__':
    month = {'01':'Jan','02':'Feb','03':'Mar','04':'Apr','06':'June','07':'July','10':'Oct','11':'Nov','12':'Dec'}
    [ARG, dateStamp] = sys.argv[1:]

    yearStamp = dateStamp.split('-')[0]
    monthStamp = yearStamp+"-"+month[dateStamp.split('-')[1]]

#    sourceDir = "/project/boulgakov/microscope3/rawFiles/jagannath/2015-July/2015-07-02"  
#    subDir = "150702_A2S_CycleMock4Edman4_OSS_BOTH_FRET_flds004"
#    pathDir = os.path.join(sourceDir,subDir)

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

    if ARG == 'COLOCALIZE':
        colocal_dict = colocalizeImages(dateStamp,pathDir,destDir)
        writeToFile(colocal_dict,dateStamp,destDir)

    
    t1 = time.clock()
    print "Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0,
                                                       time.strftime("%d %b %Y  %H:%M:%S",time.localtime()))

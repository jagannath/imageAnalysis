#! /home/jaggu/anaconda/bin/python2.7
from __future__ import division
import sys
import os
import time
import cv2
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
from scipy import ndimage,signal,misc,optimize
import cPickle as pickle
import collections
from matplotlib_venn import venn3, venn3_circles
sys.path.append('/project/marcotte/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis')
from commonFunctions import locate,simpleShow,bounds,simpleOpen


from matplotlib import rcParams
rcParams['axes.labelsize'] = 9
rcParams['xtick.labelsize'] = 9
rcParams['ytick.labelsize'] = 9
rcParams['legend.fontsize'] = 9
rcParams['font.size'] = 12 


""" This should be very close to the tracking the trace files. It obtains the
peak info of the maximumprojection image and the image,offset values associated
(with the c2(ref) version of file). It then extracts the subimage info for each
of the peaks from the file. This is then the information for tracking the
peptides. This file may be similar in a lot of way to the classes built for the
tracking file. But I need to do a lot more book keeping for the eventual venn
diagram and so must not be copied"""

class TrackPeptides:
    """
    Restricted to obtaining peak information of the maximum intensity projection
    file. Then it obtains the offset dict and looks for the corresponding image.
    Retrieves the subimage for each peak after the offset application.
    """
    def __init__(self,pathDir,subDir,fname,ch):
        self.pathDir = pathDir
        self.subDir = subDir
        self.fname = fname
        self.ch = ch
        self.maxProjDir = os.path.join(pathDir,subDir,'maxProjection')
        self.analysisDir = os.path.join(pathDir,subDir,'trackingAnalysis')
        if not os.path.exists(self.analysisDir): os.makedirs(self.analysisDir)
        self.offsetDict = self._getOffsetDict()

    def _getOffsetDict(self):
        [offsetDictFname] = [f for f in os.listdir(self.maxProjDir) if 
                           f.startswith(self.fname) 
                           and f.endswith('.offsetDict.pkl')]
        offsetDict = \
                pickle.load(open(os.path.join(self.maxProjDir,offsetDictFname)))
        return offsetDict

    def _openPeakPkl(self,ch):
        [peakF] = [f for f in os.listdir(self.maxProjDir) if
                   f.endswith('.pkl') and
                   ch in f and # This will crash later when the fname is c2
                   not f.endswith('.offsetDict.pkl')
                ]
        print peakF
        assert len(peakF)>0
        self.peakInfo = pickle.load(open(os.path.join(self.maxProjDir,peakF)))
        return self.peakInfo

    def _savePkl(self,db,small=True):
        if small:
            pklFname = \
                os.path.join(self.analysisDir,self.fname+'.c'+str(self.ch)+'.allPeaksTrace.small.pkl')
        else:
            pklFname = \
                os.path.join(self.analysisDir,self.fname+'.c'+str(self.ch)+'.allPeaksTrace.pkl')
        pickle.dump(db,open(pklFname,'w'),-1)
        return pklFname

    def loadPeakInfo(self):
        pklFname = \
            os.path.join(self.analysisDir,self.fname+'.allPeaksTrace.pkl')
        print "Loading file :%s"%(pklFname)
        self.peakInfo = pickle.load(open(pklFname,'rb'))
        assert len(self.peakInfo) > 0
        return self.peakInfo


class Peak:
    """
    !!$$ NOTE _ ALMOST COPIED FROM TRACING PEAKS FILE
    This is the information of the identified peak position. It performs
    functions like - (a) get subimage of peak given the position (b) find the
    offset (c) calculate SNR (d) is it a peak etc
    """
    def __init__(self,pos):
        self.pos = pos

    def applyOffset(self,pos,offset):
        delx, dely = offset #CAPABLE OF BUG; OFFSET RETURNS X,Y or W,H
        delh,delw = dely,delx
        new_pos = map(int,np.around([pos[0]-delh,pos[1]-delw]))
        return new_pos

    def getRect(self,p,l=4):
        # 4 is for a 5x5 box
        (h,w) = p
        posRect = [(i,j) for i in range(h-l+2,h+l-1) for j in
                   range(w-l+2,w+l-1)]
        pxLoc = [(bounds(i),bounds(j)) for (i,j) in posRect]
        return pxLoc
    
    def getSubImage(self,pxL,f):
        subimgL = list() 
        orig,gray8bit,cimg = simpleOpen(f)
        for h,w in pxL:
            val = orig[(h,w)]
            subimgL = np.append(subimgL,val)
        subimage = subimgL.reshape(5,5)
        return subimage
    
    def getPeakInfo(self,subImage):
        status = False
        peak = subImage[1:4,1:4]
        edge = \
            np.concatenate((subImage[0],subImage[4],subImage[:,0][1:-1],subImage[:,4][1:-1]))
        peakIntensity = np.sum(peak)
        bgIntensity = np.sum(edge)
        snr = (np.max(subImage) - np.mean(edge))/np.std(edge)
        if snr>2: status = True
        return status,peakIntensity,bgIntensity, snr



class AnalyseAllExperiments:
    """
    The idea is that this class will analyse all the peaks associated with each
    of the experiment and make some kind of truth table
    """
    def __init__(self,pathDir,subDir,fname,frame,ch):
        self.pathDir = pathDir
        self.subDir = subDir
        self.fname = fname
        self.frame = frame
        self.ch = ch
        self.maxProjDir = os.path.join(pathDir,subDir,'maxProjection')
        self.analysisDir = os.path.join(pathDir,subDir,'trackingAnalysis')
        self.destDir = os.path.join(self.analysisDir,frame)
        if not os.path.exists(self.destDir): os.makedirs(self.destDir)
    
    def _loadPeakTrack(self,ch):
        # Peak tracking across experiments
        [peakF] = [f for f in os.listdir(self.analysisDir) if
                   f.endswith('.pkl') and not
                   'rong' in f and 
                   ch in f  # This will crash later when the fname is c2
                ]
        print peakF

        assert len(peakF)>0
        self.peakInfo = pickle.load(open(os.path.join(self.analysisDir,peakF)))
        return self.peakInfo

    def makeTruthTable(self,nbrSets=3):
        # Someday based on nbrSets, I will recurse
        truth = ['map','not']
        truthTable = [(t1,t2,t3) for t1 in truth 
                      for t2 in truth 
                      for t3 in truth]
        return truthTable

    def getPeakStatus(self,val,snr=2):
        subImage = val[6]
        peaksum = val[3]
        bgsum = val[4]
        floor = bgsum*( 9/16)
        ht_overBg = peaksum - floor
        edge = \
            np.concatenate((subImage[0],subImage[4],subImage[:,0][1:-1],subImage[:,4][1:-1]))
        peak = subImage[1:4,1:4]
        calcsnr = (np.max(subImage) - np.mean(edge))/np.std(edge)
        if calcsnr>snr:
            status = True
        else:
            status = False
        return status,calcsnr,ht_overBg,peaksum,bgsum

    def populateTruthTable(self,peakInfo,truthTable,snr=2):
        truthDict = collections.defaultdict(int)
        peakTruthDict = collections.defaultdict(list)
        for peakPos,allInfo in peakInfo.items():
            truthCond = list()
            peakInfo = list()
            for val in allInfo:
                status,calcSNR,ht_overBg,peaksum,bgsum = self.getPeakStatus(val,snr)#can do more funny thinsg 
                peakInfo.append([calcSNR,ht_overBg,peaksum,bgsum])
                if calcSNR>snr:
                    truthCond.append('map')
                else:
                    truthCond.append('not')
            truthDict[tuple(truthCond)]+=1
            peakTruthDict[tuple(truthCond)].append(peakInfo)
            print peakTruthDict
            sys.exit(1)
        return truthDict,peakTruthDict

    def makeVenn3(self,truthDict,snr,addon=''):
        type1 = '.svg'
        type2 = '.png'
        f = \
            os.path.join(self.destDir,self.fname+'.c'+str(self.ch)+'.snr_'+str(snr)+addon)


        reqdVennOrder = [('map','not','not'),
                         ('not','map','not'),
                         ('map','map','not'),
                         ('not','not','map'),
                         ('map','not','map'),
                         ('not','map','map'),
                         ('map','map','map')]
        valSet = list()
        for vSet in reqdVennOrder:
            val = truthDict[vSet]
            valSet.append(val)
        
        unmapped = truthDict[('not','not','not')]
        # Making venn diagram
        plt.figure(figsize=( 5,5))
        v = \
            venn3(subsets=valSet,set_labels=('Cycle1:Mock','Cycle2:Mock','Cycle3:Edman'))
        c = venn3_circles(subsets=valSet,ls='solid')
        txt = 'unmapped='+str(unmapped)+'\n SNR='+str(snr)
        plt.title('Peak Mapping :'+self.fname + '.'+self.frame+'\n channel:'+str(self.ch))
        plt.figtext( 0.7,0.1,txt)
        plt.savefig(f+type1,dpi=300)        
        plt.savefig(f+type2,dpi=300)
        plt.close()

def getFpaths(frame):
    allFpaths = list() #tuple of (561,647)..Can extend
    pathList = list()

    pathDir1 = '/project2/marcotte/boulgakov/microscope/2014-Dec/test/A/mock1'
    subDir1 = '141213_DiAS1_PR011_2pM_Mock1_TFA1h_Both_A_flds004'

    pathDir2 = '/project2/marcotte/boulgakov/microscope/2014-Dec/test/A/mock2'
    subDir2 = '141213_DiAS1_PR011_2pM_Mock2_TFA1h_Both_A_flds007' 

    pathDir3 = '/project2/marcotte/boulgakov/microscope/2014-Dec/test/A/edman3'   
    subDir3 = '141213_DiAS1_PR011_2pM_Edman1_TFA1h_Both_A_flds008' 

    pathList = [(pathDir1,subDir1),(pathDir2,subDir2),(pathDir3,subDir3)]
    exptPaths = {'mock1':(pathDir1,subDir1),
                  'mock2':(pathDir2,subDir2),
                  'edman':(pathDir3,subDir3)}
    exptFnames_dict = dict()

    for expt,(pathDir,subDir) in exptPaths.items():
        # NOTE; PUT THIS IN THE ORDER YOU WANT
        fnamec1 = ''.join([subDir,frame,'c1','.tif'])
        fnamec2 = ''.join([subDir,frame,'c2','.tif'])
        fc1 = os.path.join(pathDir,fnamec1)
        fc2 = os.path.join(pathDir,fnamec2)
        allFpaths.append((fc1,fc2))
        exptFnames_dict[expt] = (fc1,fc2)
    return allFpaths,exptFnames_dict


def analysePeaks_truthTable(pathDir):
    subDir = 'A'
    fname = 'testCycle0_1_2'
    frame = 'xy14'
    fname = fname
    exptOrder = ('mock1','mock2','edman')
    ref = 2
    ch = 2 

    SNR = 2                  
    i = AnalyseAllExperiments(pathDir,subDir,fname,frame,ch)
    peakInfo = i._loadPeakTrack('c'+str(ch))
    truthTable = i.makeTruthTable()
    truthDict,peakTruthDict = i.populateTruthTable(peakInfo,truthTable,snr=SNR)










def findTruthTable(pathDir):
    print pathDir
    subDir = 'A'
    fname = 'testCycle0_1_2'
    frame = 'xy14'
    fname = fname
    exptOrder = ('mock1','mock2','edman')
    ref = 2
    ch = 2 

    SNR = 2                  
    addon = '.calcSNR'
    e = AnalyseAllExperiments(pathDir,subDir,fname,frame,ch)
    peakInfo = e._loadPeakTrack('c'+str(ch))
    # The order of experiments must match with the truth table 
    truthTable = e.makeTruthTable(nbrSets = len(exptOrder))
    truthDict,peakTruthDict = e.populateTruthTable(peakInfo,truthTable,snr=SNR)
    e.makeVenn3(truthDict,SNR,addon)


def track_thru_expts(pathDir):
    print pathDir
    subDir = 'A'
    fname = 'testCycle0_1_2'
    frame = 'xy14'
    exptOrder = ('mock1','mock2','edman')
    ref = 2
    ch = 1 

    c = TrackPeptides(pathDir,subDir,fname,ch)
    peakDict = c._openPeakPkl('ch'+str(ch))#willchange later
    offsetDict = c.offsetDict
    allFpaths,exptFnames_dict = getFpaths(frame)
    allPeakTrace_dict = dict()
    for i, pos in enumerate(peakDict.keys()):
        print "Processing Peak # %d at position %s"%(i,str(pos))
        peakInfo = list()

        p = Peak(pos)
        for order in exptOrder:
            fpath = exptFnames_dict[order][ch-1] #Ch is 2-1;
            fref = fpath[:-5]+str(ref)+'.tif'
            offset = offsetDict[fref]           

            # Need to get the peakPos in subpixel
            peakPos =  peakDict[pos][0],peakDict[pos][1]
            r2 = peakDict[pos][10]
            newPos = p.applyOffset(peakPos,offset)
            pxL = p.getRect(newPos)
            subImage = p.getSubImage(pxL,fpath)
            status,peakInt, bgInt, snr = p.getPeakInfo(subImage)
            peakInfo.append([fpath,offset,status,peakInt,bgInt,snr,subImage,r2])
        allPeakTrace_dict[pos] = peakInfo        
    pklF = c._savePkl(allPeakTrace_dict,small=False)
    print "Generated Pickle dictionary of the peak across experiments for \
            channel %s in %s"%('ch'+str(ch),pklF)

    SNR = 2                  
    addon = '.rectSNR'
    e = AnalyseAllExperiments(pathDir,subDir,fname,frame,ch)
    peakInfo = e._loadPeakTrack('c'+str(ch))
    # The order of experiments must match with the truth table 
    truthTable = e.makeTruthTable(nbrSets = len(exptOrder))
    truthDict = e.populateTruthTable(peakInfo,truthTable,snr=SNR)
    e.makeVenn3(truthDict,SNR,addon)


if __name__ == '__main__':
    monthIs =  {'05':'May','06':'June','07':'July','08':'Aug','09':'Sept','10':'Oct','11':'Nov','12':'Dec'}

    [ARG,dateStamp] = sys.argv[1:]
    sourceDir = '/project2/marcotte/boulgakov/microscope'
    month = monthIs[dateStamp.split('-')[1]]
    pathDir = os.path.join(sourceDir,"2014-"+month,dateStamp)
    
    pathDir = '/project2/marcotte/boulgakov/microscope/2014-Dec/test'

    t0 = time.clock()
    if ARG == 'TEST':
        analysePeaks_truthTable(pathDir)
    elif ARG == 'TRACK':
        track_thru_expts(pathDir)
        findTruthTable(pathDir)
    else:
        raise SystemExit("Incorrect argument")

    t1 = time.clock()
    print ("Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0,time.strftime("%d %b %Y %H:%M:%S",time.localtime())))

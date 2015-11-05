#! /home/jaggu/anaconda/bin/python2.7
from __future__ import division
import os
import sys
from Peaks import ImagePeaks
import cPickle as pickle
import numpy as np
import time
import collections
import cv2
import shelve
import random
import multiprocessing
from  matplotlib import pyplot as plt
from matplotlib import gridspec as gridspec
from matplotlib.patches import Rectangle

sys.path.append('/project/current/project/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis')
from commonFunctions import locate,simpleShow,simpleOpen,bounds,makeDir
import my_math 


# MATPLOTLIB DEFAULT PARAMS
from matplotlib import rcParams
rcParams['axes.labelsize'] = 8 
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['legend.fontsize'] = 8
rcParams['font.size'] = 8

# TICKS
from matplotlib.ticker import MaxNLocator
my_xlocator = MaxNLocator(5)
my_ylocator = MaxNLocator(3)



class TimeTracePeaks:
    """
    I will try to restrict this to only obtaining peak information by looking
    through the tiff stacks (sorting it based). It will determine if it is a
    peak if the SNR calculated is >threshold
    """
    def __init__(self,pathDir,fname,snr_cutoff=4):
        self.pathDir = pathDir
        self.fname = fname
        self.subDir = os.path.join(pathDir,self.fname)
        self.snr_cutoff = 4
        self.maxProjDir = os.path.join(self.subDir,'maxProjection')
        self.offsetDict,self.peakDict = self._openPkl()
        self.analysisDir = os.path.join(self.subDir,'traceAnalysis')
        if not os.path.exists(self.analysisDir): os.makedirs(self.analysisDir)
        
    
    def makeShelve(self,writeback=False):
        shelveFname = \
                os.path.join(self.analysisDir,self.fname+'.allPeaksTrace.shelve')
        db = shelve.open(shelveFname,writeback=False)
        return db

    def _openPkl(self):
        # Temporary and assumes only two pkl for that file
        pklF = [f for f in os.listdir(self.maxProjDir) if
                f.startswith(self.fname) and f.endswith('.pkl')]
        assert len(pklF) == 2
        for f in pklF:
            if 'png' in f:
                fpath = os.path.join(self.maxProjDir,f)
                peakDict = pickle.load(open(fpath,'r'))
            else: 
                fpath = os.path.join(self.maxProjDir,f)
                offsetDict = pickle.load(open(fpath,'r'))
                offsetDict_ordered = \
                    collections.OrderedDict(sorted(offsetDict.items()))
        return offsetDict_ordered, peakDict

    def _savePkl(self,db,small=True):
        if small:
            pklFname = \
                os.path.join(self.analysisDir,self.fname+'.allPeaksTrace.small.pkl')
        else:
            pklFname = \
                os.path.join(self.analysisDir,self.fname+'.allPeaksTrace.pkl')
        pickle.dump(db,open(pklFname,'w'),-1)
        return pklFname

    def loadPeakInfo(self):
        pklFname = \
            os.path.join(self.analysisDir,self.fname+'.allPeaksTrace.pkl')
        print "Loading file :%s"%(pklFname)
        self.peakInfo = pickle.load(open(pklFname,'rb'))
        return self.peakInfo

    def frameKinetics(self,peakInfo_dict,nbrFrames=300):
        
        countDict = collections.defaultdict(int)
        peak_intensity = collections.defaultdict(list)
        bg_intensity = collections.defaultdict(list)
        snr_dict = collections.defaultdict(list)
        ht_overFloor = collections.defaultdict(list)
        frame_peaksFound = collections.defaultdict(list)
        frame_peaksNotFound = collections.defaultdict(list)
        nbr_onOff_dict = collections.defaultdict(list)
        binSize = 20

        for pos,peakInfo in peakInfo_dict.items():
            print "Processing peak %s"%(str(pos))
            nbr_on,nbr_off = 0,0 
            binNbr = 1 
            allonRate_freq = list()
            
            for i in xrange(0,nbrFrames):
                frame = peakInfo[i]
                peak_sum,bg_sum,snr_val = frame[3],frame[4],frame[5]
                if snr_val > self.snr_cutoff:
                    frame_peaksFound[i].append(pos)
                    countDict[i]+=1
                    ht = (peak_sum - (bg_sum*9/16))
                    ht_overFloor[i].append(ht)
                    peak_intensity[i].append(peak_sum)
                    bg_intensity[i].append(bg_sum)
                    snr_dict[i].append(snr_val)
                    if i == binSize*binNbr:
                        allonRate_freq.append([binNbr,nbr_on,nbr_off])
                        nbr_on = 0
                        nbr_off = 0
                        binNbr +=1
                
                    nbr_on+=1

                else:
                    frame_peaksNotFound[i].append(pos)
                    if i == binSize*binNbr:
                        allonRate_freq.append([binNbr,nbr_on,nbr_off])
                        nbr_on = 0
                        nbr_off=0
                        binNbr +=1
                    nbr_off+=1
            nbr_onOff_dict[pos].extend(allonRate_freq)
        
        frameDict = {'counts':countDict,
                     'peak_intensity':peak_intensity,
                     'bg_intensity':bg_intensity,
                     'snr_list':snr_dict,
                     'ht_overFloor':ht_overFloor,
                     'frame_peaksFound':frame_peaksFound,
                     'frame_peaksNotFound':frame_peaksNotFound,
                     'nbr_onOff':nbr_onOff_dict
                   }
        f = os.path.join(self.analysisDir,self.fname+\
                         '.snrcutoff'+str(self.snr_cutoff)+'.frameInfo.dict.pkl')
        pickle.dump(frameDict,open(f,'wb'),-1)
        return frameDict

    def getNewPeaks_frame(self,frame_peaksDict,nbrFrames=20):
        #First frame is reference
        all_peaks_found = set()
        all_peaks_notfound = set()
        newPeaks_frame = list()
        missingPeaks_frame = list()
        for i,allPeaks in frame_peaksDict.items()[0:nbrFrames]:
            pos_set = set(allPeaks)
            if i == 0: 
                all_peaks_found = pos_set
                all_peaks_notfound = pos_set
            unique_peaks = pos_set - all_peaks_found
            missing_peaks = all_peaks_found - pos_set
            all_peaks_found = pos_set.union(all_peaks_found)
            newPeaks_frame.append(len(unique_peaks))
            missingPeaks_frame.append(len(missing_peaks))
        return newPeaks_frame,missingPeaks_frame

    def graphDelta_peaks(self,newPeaks,missingPeaks):
        type1 = '.png'
        type2 = '.svg'
        self.destDir = \
            os.path.join('/project2/marcotte/jaggu/dataAnalysis/microscope1/2014-Dec',dateStamp,'images')

        imageDestDir = os.path.join(self.destDir,self.fname,'_frameAnalysis')
        makeDir(imageDestDir) 
        
        fig = plt.figure(figsize=( 5,5),dpi=300)
        gs = gridspec.GridSpec(nrows=2,ncols=1) 
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
        xL = xrange(0,len(newPeaks))
        yL1 = newPeaks
        yL2 = [-1*i for i in missingPeaks]
        ax1.bar(xL,yL1,color='#0E3D59',width=1)
        ax2.bar(xL,yL2,color='#F25C05',width=1)
        
        totalPeaks = len(self.peakDict.keys())
        s = 'Total peaks : '+str(totalPeaks)
        ax1.text( 0.7,0.8,s,fontsize=6,                                   
            horizontalalignment='left',transform=ax1.transAxes)        
        
        ax1.set_ylim([0,np.amax(yL1)*1.2])
        ax1.set_ylabel('Number of new peaks ')
        ax2.set_xlabel('Number of frames /500 ms')
        ax2.set_ylabel('Number of missing peaks ')
        ax1.set_xticks([])
        ax2.set_ylim([np.amin(yL2)*1.2,0])
        ax2.xaxis.set_major_locator(MaxNLocator(5))
        ax2.yaxis.set_major_locator(MaxNLocator(3))

        f = os.path.join(imageDestDir,'delta_changePeaks')       
        plt.tight_layout()
        fig.savefig(f+type1)
        fig.savefig(f+type2)
        plt.close()
        return True

    def plotPbleachingCurve(self,ydata):
        # Doing one fit exponential and two fit exponentials.

        def _drawPlot(typeDecay):
            type1 = '.svg'
            type2 = '.png'

            # Drawing
            fig = plt.figure(figsize=( 5,5),dpi=300)
            ax = fig.add_subplot(111)
            xL = range(0,len(y_data))

            pform = ["{0:.3f}".format(n) for n in param_opt]
            r2_format = "{0:.3f}".format(r_2)
            rmse_format = "{0:.3f}".format(rmse)
            t_half_1_format = "{0:.3f}".format(t_half_1)
            t_half_2_format = "{0:.3f}".format(t_half_2)
                                               
            s1 = 'A * exp(-b*t) + A*exp(-d*t) + e '
            s2 = 'A:'+pform[0] + '; b:'+pform[1]+'; d:'+pform[2]+'; e:'+pform[3]
            s3 = 'RMSE: '+rmse_format + '\nR2: '+r2_format
            s4 = 'T(1/2)_b: '+t_half_1_format+'\nT(1/2)_d: '+t_half_2_format
            s = s1+'\n'+s2+'\n'+s3+'\n'+s4
            
            ax.plot(xL,y_data,color='#6E1E90',marker='o',markersize=2,ls='None')
            ax.plot(xL,y_fit,color='#88A61B',lw=2)
            ax.xaxis.set_major_locator(my_xlocator)
            ax.yaxis.set_major_locator(MaxNLocator(5))
            ax.set_xlabel('Number of frames /500 ms')
            if normalization: 
                ax.set_ylim([0,1.2])
                ax.set_ylabel('Normalized number of peaks identified')
            else:
                ax.set_ylim([0,np.amax(y_data)*1.2])
                ax.set_ylabel('Number of peaks identified')
            ax.text( 0.5,0.6,s,fontsize=6,
                    horizontalalignment='left',transform=ax.transAxes)
            
            fig.savefig(os.path.join(imageDestDir,f+'.'+str(normalization)+type1))
            fig.savefig(os.path.join(imageDestDir,f+'.'+str(normalization)+type2))
            plt.close()
            return True

        self.destDir = \
            os.path.join('/project2/marcotte/jaggu/dataAnalysis/microscope1/2014-Dec',dateStamp,'images')
        imageDestDir = \
                os.path.join(self.destDir,self.fname,'_frameAnalysis','pbleaching')            
        makeDir(imageDestDir)                                                 

        typeDecay = 'two_exp'
        f = self.fname + '.'+typeDecay+'.pbleaching'
        ofile = open(os.path.join(imageDestDir,f+'.csv'),'w')
        ofile.write('\t'.join(['FULL PATH','SUBDIR',
                                'HALF-LIFE-1','HALF-LIFE-2',
                                'PARAMS','r_2','RMSE','NORMALIZATION','\n']))

        for normalization in [True,False]:
            if normalization: y_data = np.array(ydata)/np.amax(ydata)
            else: y_data = ydata
            param_opt, rmse, r_2, y_fit = my_math.two_exponential(y_data,norm=normalization)
            a,b,d,e = param_opt
            t_half_1 = np.log(2)/b
            t_half_2 = np.log(2)/d
            _drawPlot(typeDecay)
            ofile.write('\t'.join([self.subDir,self.fname,str(t_half_1),str(t_half_2),str(param_opt),str(r_2),str(rmse),str(normalization),'\n']))
        ofile.close()
        return True

    def plotMeanHT_SNRCurves(self,yData,ySTD,stdFill=False):
        type1 = '.png'
        type2 = '.svg'
        self.destDir = \
            os.path.join('/project2/marcotte/jaggu/dataAnalysis/microscope1/2014-Dec',dateStamp,'images')

        imageDestDir = os.path.join(self.destDir,self.fname,'_frameAnalysis')
        makeDir(imageDestDir) 
        
        fig = plt.figure(figsize=( 5,5),dpi=300)
        ax1 = fig.add_subplot(111)
        
        xL = xrange(0,len(yData[0]))
        yHT,ySNR = yData
        yHT_std,ySNR_std = ySTD

        ax1.xaxis.set_major_locator(MaxNLocator(4))
        ax1.set_xlabel('Number of frames /500 ms')
        ax1.plot(xL,yHT,color='#6E1E90',lw=2)

        ax1.yaxis.set_major_locator(my_xlocator)
        ax1.set_ylabel('Mean Intensity of Peaks (over floor) \n \
                       (Sum of counts in 3x3 pixels)',color='#6E1E90',fontsize=8)
        ax1.set_ylim([0,np.amax(yHT)*1.2])
        ax1.set_xlim([0,300])

        ax2 = ax1.twinx()
        ax2.plot(xL,ySNR,color='#F25C05',lw=2)
        ax2.yaxis.set_major_locator(MaxNLocator(5))
        ax2.set_ylabel('Mean Signal to Noise (SNR)',color='#F25C05',fontsize=8)
        ax2.set_ylim([0,np.amax(ySNR)*1.2])
        ax2.set_xlim([0,300]) 

        if stdFill:
            yHT_plus = np.array(yHT)+np.array(yHT_std) 
            yHT_minus = np.array(yHT)-np.array(yHT_std)
            ySNR_plus = np.array(ySNR)+np.array(ySNR_std)
            ySNR_minus = np.array(ySNR)-np.array(ySNR_std)

            ax1.fill_between(xL,yHT_plus,yHT_minus,color='#6E1E90',alpha=0.1)
            ax2.fill_between(xL,ySNR_plus,ySNR_minus,color='#F25C05',alpha=0.1)       
            ax1.set_ylim([0,np.amax(yHT)*2.0])
            ax2.set_ylim([0,np.amax(yHT)*2.0])
            f = 'HeightOverFloor_SNR.withSTD'

        f = 'HeightOverFloor'
        plt.tight_layout()

        fig.savefig(os.path.join(imageDestDir,f+type1))
        fig.savefig(os.path.join(imageDestDir,f+type2))
        plt.close()
        return True


    def analyseOnRates(self,onOff_freqList):
        type1 = '.png'                                                                                             
        type2 = '.svg'
        self.destDir = \
            os.path.join('/project2/marcotte/jaggu/dataAnalysis/microscope1/2014-Dec',dateStamp,'images')

        imageDestDir = \
            os.path.join(self.destDir,self.fname,'_frameAnalysis','dwellStates')
        makeDir(imageDestDir)
        fig = plt.figure(figsize=( 5,5),dpi=300)
        gs = gridspec.GridSpec(nrows=2,ncols=2)           

        def _plotOnRateHist(nbrs,ax,nbrFrames,normalized=False):
            
            if normalized:
                norm_weights = np.ones_like(nbrs)/len(nbrs)#Interesting
                n,bins,patches = \
                    ax.hist(nbrs,bins=50,color='#0E3D59',weights=norm_weights)
                ax.set_ylabel('Normalized Frequency')
            else:
                n,bins,patches = \
                    ax.hist(nbrs,bins=50,color='#0E3D59')
                ax.set_ylabel('Frequency')                   
            
            
            for i,rate in enumerate(bins):
                if rate <0.5: patches[(i-1)].set_fc('#D92525')
            ax.set_xlabel('Fraction of ON states')
            ax.xaxis.set_major_locator(MaxNLocator(3))
            ax.set_xticks([0,0.5,1])
            ax.yaxis.set_major_locator(MaxNLocator(4))

            s = 'First '+str(nbrFrames) + ' Frames'
            ax.text( 0.5,0.7,s,fontsize=6,                              
                    horizontalalignment='left',transform=ax.transAxes)
            return ax

        bin_0 = list()
        bin_100 = list()
        bin_200 = list()
        bin_300 = list()

        for freqL in onOff_freqList:
            forBin_0 = freqL[0][1]/(freqL[0][1]+freqL[0][2])
            bin_0.append(forBin_0)
            forBin_100 = np.mean([on/(on+off) for (b,on,off) in freqL[0:5]])
            bin_100.append(forBin_100)
            forBin_200 = np.mean([on/(on+off) for (b,on,off) in freqL[0:10]])
            bin_200.append(forBin_200)
            forBin_300 = np.mean([on/(on+off) for (b,on,off) in freqL])
            bin_300.append(forBin_300)

        for norm in [True,False]:
            ax0 = plt.subplot(gs[0])
            _plotOnRateHist(bin_0,ax0,20,normalized=norm)
            ax1 = plt.subplot(gs[1])
            _plotOnRateHist(bin_100,ax1,100,normalized=norm)
            ax2 = plt.subplot(gs[2])
            _plotOnRateHist(bin_200,ax2,200,normalized=norm)
            ax3 = plt.subplot(gs[3])
            _plotOnRateHist(bin_300,ax3,300,normalized=norm)
            plt.tight_layout()
            
            f = os.path.join(imageDestDir,'onRate.hist.'+str(norm))
            fig.savefig(f+type1)
            fig.savefig(f+type2)
        return True

    def analyseFrameKinetics(self):

        f = os.path.join(self.analysisDir,self.fname+\
                         '.snrcutoff'+str(self.snr_cutoff)+'.frameInfo.dict.pkl')
        frameDict = pickle.load(open(f,'rb'))
        print "Frame Dictionary : %s  loading "%f 
        
        onOff_freqL = frameDict['nbr_onOff'].values()
        self.analyseOnRates(onOff_freqL)
        
        new_peaks_frame,missingPeaks_frame = self.getNewPeaks_frame(frameDict['frame_peaksFound'])
        self.graphDelta_peaks(new_peaks_frame,missingPeaks_frame)
        print "Delta Peaks Graph plotted "
        
        xL = frameDict['counts'].keys()
        nbrCounts = frameDict['counts'].values()
        self.plotPbleachingCurve(nbrCounts)
        print "Photobleaching Curve plotted "
        htList = frameDict['ht_overFloor'].values()
        htL = [np.mean(h) for h in htList]
        htL_std = [np.std(h) for h in htList]
        snrList = frameDict['snr_list'].values()       
        snrL = [np.mean(s) for s in snrList]
        snrL_std = map(np.std,snrList)
        self.plotMeanHT_SNRCurves([htL,snrL],[htL_std,snrL_std])
        print "Height over floor and SNR Curves plottted"

        return True


class Peak:
    """
    This is the information of the identified peak position. It performs
    functions like - (a) get subimage of peak given the position (b) find the
    o ffset (c) calculate SNR (d) is it a peak etc
    """
    def __init__(self,pos):
        if isinstance(pos,str):
            self.pos = eval(pos)
        else: 
            self.pos = pos
        self.destDir = \
            os.path.join('/project2/marcotte/jaggu/dataAnalysis/microscope1/2014-Dec',dateStamp,'images')
        if not os.path.exists(self.destDir):os.makedirs(self.destDir)

    def applyOffset(self,pos,offset):
        delh, delw = offset
        new_pos = map(int,[pos-delh,pos-delw])
        return new_pos
    
    def getRect(self,p,l=4):
        # 4 is for a 5x5 box
        (h,w) = p
        posRect = [(i,j) for i in range(h-l+2,h+l-1) for j in
                   range(w-l+2,w+l-1)]
        pxLoc = [(bounds(i),bounds(j)) for (i,j) in posRect]
        return pxLoc
    
    def getSubImage(self,pxL,f,dim=5):
        subimgL = list() 
        orig,gray8bit,cimg = simpleOpen(f)
        for h,w in pxL:
            val = orig[(h,w)]
            subimgL = np.append(subimgL,val)
        subimage = subimgL.reshape( dim,dim)
        return subimage
    
    def getPeakInfo(self,subImage):
        status = False
        peak = subImage[1:4,1:4]
        edge = \
            np.concatenate((subImage[0],subImage[4],subImage[:,0][1:-1],subImage[:,4][1:-1]))
        peakIntensity = np.sum(peak)
        bgIntensity = np.sum(edge)
        snr = (np.max(peak) - np.mean(edge))/np.std(edge)
        if snr>4: status = True
        return status,peakIntensity, bgIntensity, snr

    def intensityTrace(self,peakTrace,nbrFrames=300):
        #Peak Info is a list
        # --0----------1------2--------3------------------4-------------
        #[Filename, offset, status,  peak Intensity, Background Intensity
        # ---5--------6-----
        # SNR        subImage
        peakIntensity_traceL = list()
        floorIntensity_traceL = list()#3x3 pixel
        snr_traceL = list()
        ht_floorL = list()
        for i in xrange(0,nbrFrames):
            peakIntensity,bgIntensity,snr =\
                peakTrace[i][3],peakTrace[i][4],peakTrace[i][5]
            ht_floor = (peakIntensity) - (bgIntensity*9/16)
            peakIntensity_traceL.append(peakIntensity)
            floorIntensity_traceL.append(bgIntensity*9/16)
            snr_traceL.append(snr)
            ht_floorL.append(ht_floor)
       
        return peakIntensity_traceL,floorIntensity_traceL,snr_traceL,ht_floorL


    def linePlotHistogram(self,intensityList,traceDir,label,color='#6E1E90'):
        type1 = '.png'
        type2 = '.svg'
        posName = 'PeakPos_'+str(self.pos)
        imageDestDir = os.path.join(self.destDir,traceDir,posName)
        makeDir(imageDestDir) 
        
        fig = plt.figure(figsize=( 5,5),dpi=300)
        gs = gridspec.GridSpec(nrows=1,ncols=2,width_ratios=[ 2,1],wspace=0.02)
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
        # Line Plot
        xL1 = range(0,len(intensityList))
        ax1.plot(xL1,intensityList,lw=2,label=label)
        ax1.axhline(y=0,color='#000000')
        ax1.set_xlabel('Number of frames/500 ms')
        ax1.set_ylabel('Intensity sum (3x3pixel area) (Counts)')
        ax1.yaxis.set_major_locator(my_ylocator)
        ax1.xaxis.set_major_locator(MaxNLocator(3))
        #ax1.xaxis.set_ticklabels(xTicks)
        s = 'Position (h,w)'+str(self.pos)
        ax1.text( 150,80000,s,fontsize=8)
        
        # Dwell time Histogram
        nbrBins = 50
        max_intensity = np.amax(intensityList)
        n,bins,patches =  ax2.hist(intensityList,bins=nbrBins,
                                   orientation='horizontal',fc='#F25C05')
        ax2.xaxis.set_major_locator(MaxNLocator(4))
        ax2.yaxis.set_ticks([])
        
        f = os.path.join(imageDestDir,str(self.pos)+'_'+label+'.DwellTimesHist')
        fig.savefig(f+type1)
        fig.savefig(f+type2)
        plt.close()
        return True

    def drawOnePlot(self,intensityList,traceDir,label,color='#6E1E90'):
        type1 = '.png'
        type2 = '.svg'
        posName = 'PeakPos_'+str(self.pos)
        imageDestDir = os.path.join(self.destDir,traceDir,posName)
        makeDir(imageDestDir) 

        fig = plt.figure(figsize=( 5,5),dpi=300)
        ax = fig.add_subplot(111)
        if len(intensityList)==2:
            xL = range(0,len(intensityList[0]))
            peakLabel,floorLabel = label.split('_')
            peakL = intensityList[0]
            floorL = intensityList[1]
            ax.plot(xL,peakL,color='#0E3D59',lw=2,label=peakLabel)
            ax.plot(xL,floorL,color='#F25C05',lw=2,label=floorLabel)
        if len(intensityList)>2:
            xL = range(0,len(intensityList))
            ax.plot(xL,intensityList,color='#6E1E90',lw=2,label=label)
            if label is 'SNR': ax.axhline(y=2,color='#D92525')
            else: ax.axhline(y=0,color='#000000')
        
        s = 'Position (h,w)'+ str(self.pos)
        ax.text( 150,80000,s,fontsize=10)
        ax.legend()
        ax.set_xlabel('Number of frames /500 ms')
        ax.set_ylabel('Intensity sum (3x3 area) (Counts)')
        ax.xaxis.set_major_locator(my_xlocator)
        ax.yaxis.set_major_locator(my_ylocator)
        plt.tight_layout()

        f = os.path.join(imageDestDir,str(self.pos)+'_'+label)
        fig.savefig(f+type1)
        fig.savefig(f+type2)
        plt.close()

    def makePeakImageStack(self,pathDir,traceDir,nbrFrames=300,autoscale=False):
        imgsLoc = os.path.join(pathDir,traceDir)
        allImages = [os.path.join(imgsLoc,f) for f in os.listdir(imgsLoc) if
                     os.path.isfile(os.path.join(imgsLoc,f)) 
                     if not 'maxProj' in f]
        allImages.sort() 
        #Hoping this is the only way I get the right trace like before
        posName = 'PeakPos_'+str(self.pos)
        if autoscale:
            imageDestStackDir = \
                os.path.join(self.destDir,traceDir,posName,'imageStack_scaled')
        else:
            imageDestStackDir = os.path.join(self.destDir,traceDir,posName,'imageStack')
        makeDir(imageDestStackDir)

        for i,f in enumerate(allImages[0:nbrFrames]):
            fig = plt.figure(figsize=( 3,3),dpi=300)
            ax = fig.add_subplot(111)
            frameNbr = f[-7:-4]
            print "Processing Frame : %s"%frameNbr
            orig,gray,cimg = simpleOpen(f)
            h,w = map(int,self.pos)
            roi = orig[h-10:h+11,w-10:w+11]
            if i == 0:
                self.v_min = int(np.amin(roi))
                self.v_max = int(np.amax(roi))
                self.v_med = int(np.median(roi))
                self.v_75 = int( 0.75 * self.v_max)
                    
                self.ctickRange = [self.v_min,self.v_med,self.v_75,self.v_max]
                self.ctickLabels = [self.v_min,self.v_med,self.v_75,self.v_max]
            xticklabels = range(w-11,w+11,5)
            yticklabels = range(h-11,h+11,5)#It seems only 4 ticks set
            ax.xaxis.set_ticklabels(xticklabels)
            ax.yaxis.set_ticklabels(yticklabels)
            t = "Frame Number : %s"%(frameNbr)
            ax.set_title(t,fontsize=8)
            if autoscale:cax = ax.imshow(roi,cmap='gist_gray')
            else:cax = ax.imshow(roi,cmap='gist_gray',vmin=self.v_min,vmax=self.v_max)
            cbar = fig.colorbar(cax) #Color map jumps a lot around
            cbar.set_ticks(self.ctickRange)
            cbar.set_ticklabels(self.ctickLabels)
             
            ax.add_patch(Rectangle(( 8,8),5,5,fill=False,color='#F29F05',lw=1))
            fname = 'PeakStack'+'_'+str(self.pos)+frameNbr+'.tif'
            fig.savefig(os.path.join(imageDestStackDir,fname))
            plt.close()
        
        return True

######## END OF CLASS ############################

def plot_frameKinetics(traceDir):

    SNR_CUTOFF = 4 
    print "Processing %s"%traceDir 

    t = TimeTracePeaks(pathDir,traceDir,snr_cutoff=SNR_CUTOFF)
    #allPos = t.peakDict.keys()
    
    #peakTraceInfo = t.loadPeakInfo()
    #frameDict = t.frameKinetics(peakTraceInfo)
    t.analyseFrameKinetics()

    return "Done" 


def make_plots_individualPeakTrace(traceDir):
#    allDir = ['141222_DiAS1_PR011_200fM_AfTFA_MeOH_561_10x500ms_40C_trace001',
#              '141222_DiAS1_PR011_200fM_deGasMeOH_561_500ms_40C_trace001',
#              '141222_DiAS1_PR011_200fM_10mMTXMeOH_561_500ms_40C_trace002']
    
#    traceDir = allDir[2]

    def _inside(x):
        return True if (x<450 and x>50) else False


    SNR_CUTOFF = 4
    #def chunks(l,n):
    #    for i in xrange(0,len(l),n):
    #        yield l[i:i+n]

    #splitDir = [subDir for subDir in next(os.walk(pathDir))[1] if 'trace' in
    #            subDir]
    #allDir = list(chunks(splitDir,10))

    t = TimeTracePeaks(pathDir,traceDir,snr_cutoff=SNR_CUTOFF)                    
    #Opens it from the calculated peaks; Should match the peak pos from the
    print "Processing Trace Dir: %s"%(traceDir)

    t = TimeTracePeaks(pathDir,traceDir,snr_cutoff=SNR_CUTOFF)
    #Opens it from the calculated peaks; Should match the peak pos from the
    allPos = t.peakDict.keys()
    peakTraceInfo = t.loadPeakInfo()

    random_pos_list = random.sample(allPos,5)

    allIntensityDict = collections.defaultdict(list)
    nbrFrames = 300
    for randomPos in random_pos_list:
        h,w = randomPos
        if not _inside(h) or not _inside(w):continue;

        print "Peak %s"%str(randomPos)
        p = Peak(randomPos)
        try:
            pTrace = peakTraceInfo[randomPos]
        except KeyError: continue
        traceInfoL = p.intensityTrace(pTrace)
        peak,floor,snr,ht_overFloor = traceInfoL
        p.linePlotHistogram(ht_overFloor,traceDir,
                            'Peak Intensity over Floor',randomPos)
        p.drawOnePlot(snr,traceDir,'SNR',randomPos)
        p.drawOnePlot([peak,floor],traceDir,'PeakIntensity_floor',randomPos)
        p.drawOnePlot(ht_overFloor,traceDir,'Peak Intensity over Floor',randomPos)

        p.makePeakImageStack(pathDir,traceDir,autoscale=True)
        p.makePeakImageStack(pathDir,traceDir)
        print "Completed plots for Trace Dir : %s"%(traceDir)
    return "Done" 

def make_peakInfo_dict(pathDir):
#    traceDir = '141222_DiAS1_PR011_200fM_10mMTXMeOH_561_500ms_40C_trace002'

#    allDir = ['141222_DiAS1_PR011_200fM_AfTFA_MeOH_561_10x500ms_40C_trace001',
#              '141222_DiAS1_PR011_200fM_deGasMeOH_561_500ms_40C_trace001']
    SNR_CUTOFF = 4

    def chunks(l,n):
        for i in xrange(0,len(l),n):
            yield l[i:i+n]

    splitDir = [subDir for subDir in next(os.walk(pathDir))[1] if 'trace' in
                subDir]
    allDir = list(chunks(splitDir,8))
    
    for traceDir in allDir[6]:
        print "Processing folder %s"%traceDir
        if not 'trace' in traceDir:
            continue;
        t = TimeTracePeaks(pathDir,traceDir,snr_cutoff=SNR_CUTOFF)
        allPeakTrace_dict = dict() 
        allPeakTrace_small_dict = dict()
        #allPeakTrace_shelve = t.makeShelve()

        for i,pos in enumerate(t.peakDict.keys()):
            pos = tuple(map(int,pos)) #Ensures it is always an integer
            print "Tracing peak %d at position %s"%(i, str(pos))
            peakInfo = list()
            peakInfo_small = list()
            p = Peak(pos)
            pxL = p.getRect(pos) #The position of this should be changed if offset
            for f,offset in t.offsetDict.items():
                subImage = p.getSubImage(pxL, f)
                status, peakInt, bgInt, snr = p.getPeakInfo(subImage)
                peakInfo_small.append([f,offset,status,peakInt,bgInt,snr])
                peakInfo.append([f,offset,status,peakInt,bgInt,snr,subImage])
            allPeakTrace_small_dict[pos] = peakInfo_small
            allPeakTrace_dict[pos] = peakInfo
            #allPeakTrace_shelve[str(pos)] = peakInfo #for shelve it is str

        pklFname = t._savePkl(allPeakTrace_small_dict,small=True)
        pklFname = t._savePkl(allPeakTrace_dict,small=False)
#        allPeakTrace_shelve.close()
        return "Done" 

def start_process():
    print "Starting ",multiprocessing.current_process().name

if __name__ == '__main__':
    month = {'10':'Oct','11':'Nov','12':'Dec'}
    [ARG, dateStamp] = sys.argv[1:]
    sourceDir = "/project2/marcotte/boulgakov/microscope"

    monthStamp = "2014-"+month[dateStamp.split('-')[1]]
    pathDir = os.path.join(sourceDir,monthStamp,dateStamp)
    
    destDir = os.path.join("/project2/marcotte/jaggu/dataAnalysis/microscope1",monthStamp,dateStamp)
    t0 = time.clock()   

    allDir = [subDir for subDir in next(os.walk(pathDir))[1] if 'trace' in subDir]


    if ARG == 'TEST':
        traceDir = '141222_DiAS1_PR011_200fM_AfTFA_MeOH_647_10x500ms_40C_trace001'
        plot_frameKinetics(traceDir)
        sys.exit(1)

    if ARG == 'FRAMEANALYSE':
        #pathDir = "/project2/marcotte/boulgakov/microscope/2014-Dec/traceTest"
        #allDir = ['141222_DiAS1_PR011_200fM_AfTFA_MeOH_561_10x500ms_40C_trace001',
        #          '141222_DiAS1_PR011_200fM_deGasMeOH_561_500ms_40C_trace001',
        #          '141222_DiAS1_PR011_200fM_10mMTXMeOH_561_500ms_40C_trace002']
        for traceDir in allDir:
            plot_frameKinetics(traceDir)


        #traceDir = allDir
        #allDir = allDir[12:]
        #pool_size = 6 
        #pool = multiprocessing.Pool(processes=pool_size, initializer=start_process)        
        #pool_outputs = pool.map(plot_frameKinetics,allDir)
        #pool.close()
        #pool.join()
        #print 'Pool :',pool_outputs

    if ARG == 'MAKEINFO':
        pathDir = os.path.join(sourceDir,monthStamp,dateStamp)
        make_peakInfo_dict(pathDir)
    if ARG == 'MAKEPLOTS':

        make_plots_individualPeakTrace(pathDir)    
        #pool_size = 6 
        #pool = multiprocessing.Pool(processes=pool_size, initializer=start_process)        
        #pool_outputs = pool.map(make_plots_individualPeakTrace,allDir)
        #pool.close()
        #pool.join()

        #print 'Pool :',pool_outputs
    
    t1 = time.clock()
    print "Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0,
                                                       time.strftime("%d %b %Y  %H:%M:%S",time.localtime()))

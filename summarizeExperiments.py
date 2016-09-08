#!/home/jaggu/anaconda/bin/python

"""
Combing all the summarizing scripts to make something more consistent and huge file
Three major classes or functions are  - 
(a) Converting the images - autoscaling them - random images as well as converting a field or all fields
(b) Tab delimited file for summarizing the count density
(c) Tab delimited file for summarizing the track counts
"""

import os
import sys
import time
import collections
import numpy as np
import cPickle as pickle
import re
import argparse
import socket

hostname = socket.gethostname()
if hostname == 'canopus': headDir = "/home/jaggu/marcotte_project"
else: headDir = "/project"

sys.path.append(os.path.join('marcotte/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis'))
from commonFunctions import savePkl, loadPkl
import convertingImages as images
import peptideTracks as peptracks 
import countDensity as cdensity
import intensityTracks as iTracks 

if images.testImport(): print "Import Successful"
else: raise SystemExit("Import error of convertingImages. Quiting ..")

if peptracks.testImport(): print "Import Successful"
else: raise SystemExit("Import error of peptideTracks. Quiting ..")

if cdensity.testImport(): print "Import Successful"
else: raise SystemExit("Import error of countDensity. Quiting ..")

def getDirectory(hostname,dateStamp,microscopeNbr):
    month = {'01':'Jan','02':'Feb','03':'Mar','04':'Apr', '05':'May', '06':'June','07':'July','08':'Aug','09':'Sept', '10':'Oct','11':'Nov','12':'Dec'}
    #[ARG, dateStamp] = sys.argv[1:]
    yearStamp = dateStamp.split('-')[0]
    monthStamp = yearStamp+"-"+month[dateStamp.split('-')[1]]
    
    if microscopeNbr is 1:
        sourceDir = os.path.join(headDir,"boulgakov/microscope")
        pathDir = os.path.join(sourceDir,monthStamp,dateStamp)
        destDir = os.path.join(headDir,"jagannath/projectfiles/singleMoleculeMicroscopy/dataAnalysis/microscope1",monthStamp,dateStamp)
    elif microscopeNbr is 2:
        sourceDir = os.path.join(headDir,"boulgakov/microscope2/jagannath/rawFiles")
        pathDir = os.path.join(sourceDir,monthStamp,dateStamp)
        destDir = os.path.join(headDir,"jagannath/projectfiles/singleMoleculeMicroscopy/dataAnalysis/microscope2",monthStamp,dateStamp)
    elif microscopeNbr is 3:                                                                                                                 
        sourceDir = os.path.join(headDir,"boulgakov/microscope3/rawFiles/jagannath")                                                                   
        pathDir = os.path.join(sourceDir,monthStamp,dateStamp)                                                                             
        destDir = os.path.join(headDir,"jagannath/projectfiles/singleMoleculeMicroscopy/dataAnalysis/microscope3",monthStamp,dateStamp)    
    else:
        raise SystemExit("Incorrect microscope number. Quiting ..")

    print "Source folder      : ",pathDir
    print "Destination folder : ",destDir
    if not os.path.exists(destDir): os.makedirs(destDir)
    return sourceDir,pathDir,destDir

class ConvertImages:
    """
    Class that imports the convertingimages.py and does the appropriate images.
    """
    def __init__(self,pathDir, destDir):
        self.img_destDir = os.path.join(destDir,"images")
        if not os.path.exists(self.img_destDir): os.makedirs(self.img_destDir)
        self.pathDir = pathDir
        
    def convertPattern(self,pattern):
        """
        Walks through all the tifs in the pathDir and converts it to png, with autoscaling;This is generic and can cause the entire directories to be converted
        """
        images.convertSameFldImages(pattern,self.pathDir,self.img_destDir)
        return True
   
    def convertFrames(self,frames,totalFrames=100):
        #frames - this is a list; It can contain one or many frames. 
        #Need to change totalFrames if the frames are 1000 and increments of 10s. 
        for frame in frames:
            framePattern = (len(str(totalFrames))-(len(str(frame)))) * "0" + str(frame) # To account for 001 in frame number
            pattern = "*Cycles*xy"+framePattern+"*.tif"
            print "Processing pattern ...",pattern
            images.convertPatternImages(pattern,self.pathDir,self.img_destDir)
        return True

    def randomFields(self):
        # Chooses a random field from each step of the cycle and prints a random image (c1,c2 and c3 if so). 
        for d in os.walk(self.pathDir).next()[1]:
            print "Processing dir ...",d
            srcImg_dir = os.path.join(self.pathDir,d)
            if not (d.startswith('output') or d.startswith('pkl')):
                images.convertOneRandomImage(srcImg_dir,self.img_destDir) 
        return True 

class PeptideTracks:
    """
    Class that summarizes results from peptide density
    """
    def __init__(self, pathDir, destDir, trackFile):
        self.trackFile = trackFile
        self.outputDirpath = os.path.split(trackFile)[0]
        outputDirname = os.path.split(self.outputDirpath)[1]
        self.peptracks_destDir = os.path.join(destDir,outputDirname)
        if not os.path.exists(self.peptracks_destDir): os.makedirs(self.peptracks_destDir)       
        self.pathDir = pathDir

    def pickle_categoryDict(self,category_dict,fname,ch):
        pklF = fname + str(ch)+'.dict.pkl'
        savePkl(category_dict,pklF,self.pathDir)
        print "Pickled : ",pklF
        return pklF  
    
    def countTracks(self, pickle_exists=False,nbrframes=7,nbrchannels=2):
        trackFile = self.trackFile
        if pickle_exists:
            category_dict_ch1 = loadPkl('category_trackDetails.ch1.dict.pkl',self.pathDir)
            category_dict_ch2 = loadPkl('category_trackDetails.ch2.dict.pkl',self.pathDir)
            print "Categories Dictionaries loaded "
        else:
            category_dict_ch1,category_peakDetails_ch1 = peptracks.makeIntensityCategory(trackFile,ch=1)
            category_dict_ch2,category_peakDetails_ch2 = peptracks.makeIntensityCategory(trackFile,ch=2)
            self.pickle_categoryDict(category_dict_ch1,'category_trackDetails.ch',1)
            self.pickle_categoryDict(category_dict_ch2,'category_trackDetails.ch',2)
            self.pickle_categoryDict(category_peakDetails_ch1,'category_peakDetails.ch',1)
            self.pickle_categoryDict(category_peakDetails_ch2,'category_peakDetails.ch',2)
            print "Categories Dictionaries pickled"
        return category_dict_ch1, category_dict_ch2

    def writeTracks(self,category_dict_ch1,category_dict_ch2,nbrframes=7,nbrchannels=2):
        if len(category_dict_ch2) <= 1: 
            assert nbrchannels == 1
        desired_category = peptracks.getDesiredCategory('staggered',nbrframes)
        if nbrchannels == 2:
            self.desired_categoryList_ch1 = peptracks.getIntensityList_category(category_dict_ch1,desired_category)
            self.desired_categoryList_ch2 = peptracks.getIntensityList_category(category_dict_ch2,desired_category)
            ofname = peptracks.writeCountsFile([self.desired_categoryList_ch1,self.desired_categoryList_ch2],self.peptracks_destDir,nbrchannels)
        else:
            self.desired_categoryList_ch1 = peptracks.getIntensityList_category(category_dict_ch1,desired_category)           
            self.desired_categoryList_ch2 = None
            ofname = peptracks.writeCountsFile([self.desired_categoryList_ch1],self.peptracks_destDir,nbrchannels)       
        print "Peptide track counts and intensity file created ",ofname
        return self.desired_categoryList_ch1, self.desired_categoryList_ch2,ofname

    def edmanEfficiency(self,ofname,edmanframe_ch1=3,edmanframe_ch2=4,nbrchannels=2,calcType='aboveDecay'):
        if nbrchannels == 2:
            counts_ch1 = [item[1] for item in self.desired_categoryList_ch1]
            counts_ch2 = [item[1] for item in self.desired_categoryList_ch2]
            eff1_ch1, eff2_ch1,totalPeptides_ch1 = peptracks.calculateEdmanEfficiency(counts_ch1,edmanframe_ch1,calcType)
            eff1_ch2, eff2_ch2,totalPeptides_ch2 = peptracks.calculateEdmanEfficiency(counts_ch2,edmanframe_ch2,calcType)
            ofname = peptracks.writeEdmanEfficiency(ofname,eff1_ch1,eff2_ch1,totalPeptides_ch1,edmanframe_ch1,'ch1',calcType)
            ofname = peptracks.writeEdmanEfficiency(ofname,eff1_ch2,eff2_ch2,totalPeptides_ch2,edmanframe_ch2,'ch2',calcType)
        else:
            counts_ch1 = [item[1] for item in self.desired_categoryList_ch1]
            eff1_ch1, eff2_ch1,totalPeptides_ch1 = peptracks.calculateEdmanEfficiency(counts_ch1,edmanframe_ch1,calcType)
            ofname = peptracks.writeEdmanEfficiency(ofname,eff1_ch1,eff2_ch1,totalPeptides_ch1,edmanframe_ch1,'ch1',calcType)
        return ofname

    def sequencingEfficiency(self,ofname,edmanframe_ch1=3,edmanframe_ch2=4,nbrchannels=2,calcType='conservative'):
        if nbrchannels == 2:
            counts_ch1 = [item[1] for item in self.desired_categoryList_ch1]
            counts_ch2 = [item[1] for item in self.desired_categoryList_ch2]
            eff_ch1, cleaved_ch1,totalPeptides_ch1 = peptracks.calculateSequencingEfficiency(counts_ch1,edmanframe_ch1,calcType)
            eff_ch2, cleaved_ch2,totalPeptides_ch2 = peptracks.calculateSequencingEfficiency(counts_ch2,edmanframe_ch2,calcType)
            ofname = peptracks.writeSequencingEfficiency(ofname,eff_ch1,cleaved_ch1,totalPeptides_ch1,'ch1',calcType)
            ofname = peptracks.writeSequencingEfficiency(ofname,eff_ch2,cleaved_ch2,totalPeptides_ch2,'ch2',calcType)
        else:
            counts_ch1 = [item[1] for item in self.desired_categoryList_ch1]
            eff_ch1, cleaved_ch1,totalPeptides_ch1 = peptracks.calculateSequencingEfficiency(counts_ch1,edmanframe_ch1,calcType)
            ofname = peptracks.writeSequencingEfficiency(ofname,eff_ch1,cleaved_ch1,totalPeptides_ch1,'ch1',calcType)
        return ofname

    def writeTracks_quartilesCutoff(self,ofname,nbrframes=7,nbrchannels=2,cycle=3):
        cycle = 2  
        desired_category = peptracks.getDesiredCategory('staggered',nbrframes)
        ofname = ofname[:-4]+'.quartileCutoff.tab'
        ofile = open(ofname,'w')
        ofile.write('#Category \t Counts \n')
        ofile.write('\n Cycle used for cutoff - \t '+str(cycle)+'\n')
        if nbrchannels == 2:
            print len(self.desired_categoryList_ch1)
            print len(self.desired_categoryList_ch2)
        else:
            all_trackCounts = list()
            for q_range in ([0,25],[25,50],[50,75],[75,100]):
                
                trackCounts = peptracks.applyCutoff(self.desired_categoryList_ch1,desired_category,q_range, cycle)
                ofile.write('\n \n Quartile '+str(q_range)+ '\n') 
                for cat, counts_qTracks,q_intensityList in trackCounts:
                    ofile.write(str(cat)+'\t'+str(counts_qTracks)+'\n')
            ofile.close()
        print ofname
        return True

    def draw_boxPlotPeptideIntensity(self,category_intensityList,exptCycle,nbrframes=7,nbrchannels=2):
        """
        Makes a series of boxplot where the intensity of the peptides are applied.  
        """
        edmanCycle = exptCycle - 2 #(nbrMock=4;)
        desired_category = peptracks.getDesiredCategory('staggered',nbrframes)       
        category = desired_category[exptCycle]
        all_peptrack_intensity = category_intensityList[tuple(category)]
        title = 'intensity_edman_'+str(edmanCycle)
        iTracks.plot_boxPlots(self.pathDir,self.peptracks_destDir,all_peptrack_intensity,title)
        return True

class CountDensity:
    """
    Class that summarizes results from count density. This is probably the old way of computing the stats. So may be deprecated in future
    """
    def __init__(self,pathDir,destDir,dateStamp):
        self.pathDir = pathDir
        self.destDir = destDir
        self.dateStamp = dateStamp
        
    def fldDensity(self,typeVal='count'):
        cdensity.summarize_fld(self.dateStamp,self.pathDir,self.destDir,typeVal='count')
        return True

## TEST
def test_images():
    c = ConvertImages(pathDir,destDir)
    pattern = "*Cycles*.tif"
    c.convertPattern(pattern)
    c.convertFrames([10])
    c.randomFields()
    
def test_tracks():
    p = PeptideTracks(pathDir,destDir)
    db_ch1,db_ch2, ofname =p.countTracks(pickle_exists=True)
    p.edmanEfficiency(ofname)

def test_density():
    dateStamp = '2015-10-10'
    sourceDir,pathDir,destDir = getDirectory(hostname,dateStamp,1)
    d = CountDensity(pathDir,destDir)
    d.fldDensity()
## TEST

def runTrackExperiment(trackFile, pathDir,destDir,edmanframes,nbrframes=7,nbrchannels=2):
    # Need to make this an argument for the future
    edmanPos = 5
    nbrMock = 4
    exptCycle = (nbrMock-2)+edmanPos
    #    
    if len(edmanframes)==2:
        edmanframe_ch1,edmanframe_ch2 = map(int,edmanframes)
    else:
        edmanframe_ch1 = int(edmanframes[0])
        edmanframe_ch2 = 4
    p = PeptideTracks(pathDir,destDir,trackFile)
    db_ch1,db_ch2 = p.countTracks(pickle_exists=False,nbrframes=nbrframes,nbrchannels=nbrchannels)
    desired_categoryList_ch1,desired_categoryList_ch2, ofname = p.writeTracks(db_ch1,db_ch2,nbrframes=nbrframes,nbrchannels=nbrchannels)
    #p.draw_boxPlotPeptideIntensity(db_ch1,exptCycle,nbrframes=nbrframes,nbrchannels=nbrchannels)
    #p.writeTracks_quartilesCutoff(ofname, nbrframes=nbrframes,nbrchannels=nbrchannels,cycle=3)
    p.edmanEfficiency(ofname,edmanframe_ch1,edmanframe_ch2,nbrchannels=nbrchannels)
    p.sequencingEfficiency(ofname,edmanframe_ch1,edmanframe_ch2,nbrchannels=nbrchannels)

def runImageExperiment(pathDir,destDir,imgType,imgframes,imgpattern):
    i = ConvertImages(pathDir,destDir)
    if imgType == 'random':
        i.randomFields()
    elif imgType == 'frames':
        i.convertFrames(imgframes)
    elif imgType is 'pattern':
        i.convertPattern(imgpattern)
    else:
        raise SystemExit("Choose random,frames or pattern as arguments for images")

def runPeakDensityExperiment(pathDir,destDir,dateStamp,peakdensity_type):
    c = CountDensity(pathDir,destDir,dateStamp)
    c.fldDensity(peakdensity_type)

def runExperimentType(dateStamp,microscopeNbr,trackFile,ARGLIST=[]):
    edmanframes,nbrframes,nbrchannels,imgTypeList,peakdensity_type= ARGLIST
    sourceDir,pathDir,destDir = getDirectory(hostname,dateStamp,microscopeNbr)   
    if edmanframes:
        runTrackExperiment(trackFile, pathDir,destDir,edmanframes,nbrframes,nbrchannels)
    if imgTypeList[0]: #None is the default
        imgType,imgframes,imgpattern = imgTypeList
        print "Processing Images :%s", imgType,imgframes,imgpattern
        runImageExperiment(pathDir,destDir,imgType,imgframes,imgpattern)
    if peakdensity_type:
        print "Processing peak density : ",peakdensity_type
        runPeakDensityExperiment(pathDir,destDir,dateStamp,peakdensity_type) 

def test():
    test_tracks()
    test_images()
    test_density()


parser = argparse.ArgumentParser(description="""This is the script meant to summarize the results after computing the basic_experiment.py and/or basic_imagescript.py. There is an output destination folder\
        created (defaults to \project\jagannath\projectfiles\singleMoleculeMicroscopy\dataAnalysis)\
                                 ck_10ImageEdman_OSS_Both_Gold_Cycles*
        There are currently three components - (a) Converting Images: random images from each sub directory and entering a frame number \
                                               (b) Peptide tracks: Outputing the counts for each of the desired tracks. Computes the edman efficiency in the two channels \ 
                                               Can be extended to change the desired category in later versions \
                                               (c) Count density: The peak counts for each field (typically random) is computed and output txt file created 
        For example: python -u summarizeExperiments.py --edmanframes 3 4 --nbrcycles 11 --nbrchannels 2 --images frames --image_frames 18 49 66  
                                              """ )
parser.add_argument('trackPhotometry_file',action="store",type=str)
parser.add_argument('-nbr','--nbrcycles',action="store",dest='nbrcycles',type=int,help="Number of Cycles performed (Edman + Mock); Just count #cycles")
parser.add_argument('--edmanframes',action="store", dest='edmanframes',nargs='+',default=[3, 4],help="Edman frames for Channel-1 and/or Channel-2. eg --edmanframes 3 4")
parser.add_argument('--nbrchannels',action="store", dest='nbrchannels',default=2,type=int,help="Number of channels (1 or 2) e.g --nbrchannels 1")
parser.add_argument('--images',action="store", dest='imgType',default='frames',help="Converts the images into [random,all,frames]")
parser.add_argument('--peakdensity',action="store",dest='peakdensity',nargs=1,default=None,help="Adds the option of [counts or intensity] for each field in the experiment")
parser.add_argument('--image_frames',action="store",dest='imgframes',nargs='+',default=[13,31,79],help="Images whose Frames are converted to pngs: eg [10,20]")
parser.add_argument('--image_pattern',action="store",dest='imgpattern',default='*Cycle*',help="Converts images of the given pattern")
args = parser.parse_args()

trackFile = os.path.abspath(args.trackPhotometry_file)
dateStamp = trackFile.split('/')[-3]
scopeName = trackFile.split('/')[-5]

if scopeName == 'microscope': microscopeNbr = 1
else: "Microscope Number not 1"; sys.exit(1)

nbrframes = args.nbrcycles - 1

ARGLIST = [args.edmanframes,nbrframes,args.nbrchannels,[args.imgType,args.imgframes,args.imgpattern],args.peakdensity]
runExperimentType(dateStamp,microscopeNbr,trackFile,ARGLIST)

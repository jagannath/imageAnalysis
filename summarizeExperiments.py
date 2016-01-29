#! /home/jaggu/anaconda

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

if images.testImport(): print "Import Successful"
else: raise SystemExit("Import error of convertingImages. Quiting ..")

if peptracks.testImport(): print "Import Successful"
else: raise SystemExit("Import error of peptideTracks. Quiting ..")

if cdensity.testImport(): print "Import Successful"
else: raise SystemExit("Import error of countDensity. Quiting ..")

def getDirectory(hostname,dateStamp,microscopeNbr):
    month = {'01':'Jan','02':'Feb','03':'Mar','04':'Apr','06':'June','07':'July','10':'Oct','11':'Nov','12':'Dec'}
    #[ARG, dateStamp] = sys.argv[1:]
    yearStamp = dateStamp.split('-')[0]
    monthStamp = yearStamp+"-"+month[dateStamp.split('-')[1]]
    
    if microscopeNbr is 1:
        sourceDir = os.path.join(headDir,"boulgakov/microscope")
        pathDir = os.path.join(sourceDir,monthStamp,dateStamp)
        destDir = os.path.join(headDir,"current/project2/jaggu/dataAnalysis/microscope1",monthStamp,dateStamp)
    elif microscopeNbr is 2:
        sourceDir = os.path.join(headDir,"boulgakov/microscope2/jagannath/rawFiles")
        pathDir = os.path.join(sourceDir,monthStamp,dateStamp)
        destDir = os.path.join(headDir,"current/project2/jaggu/dataAnalysis/microscope2",monthStamp,dateStamp)
    elif microscopeNbr is 3:                                                                                                                 
        sourceDir = os.path.join(headDir,"boulgakov/microscope3/rawFiles/jagannath")                                                                   
        pathDir = os.path.join(sourceDir,monthStamp,dateStamp)                                                                             
        destDir = os.path.join(headDir,"current/project2/jaggu/dataAnalysis/microscope3",monthStamp,dateStamp)    
    else:
        raise SystemExit("Incorrect microscope number. Quiting ..")
    
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
    def __init__(self, pathDir, destDir):
        outputDirs = [item for item in os.walk(pathDir).next()[1] if item.startswith('output')]
        try: assert len(outputDirs)==1
        except AssertionError: raise SystemExit("There is no output directory for peptide tracks. Quiting ...")
        [outputDir] = outputDirs
        self.outputDirpath = os.path.join(pathDir,outputDir)
        self.peptracks_destDir = destDir
        self.pathDir = pathDir

    def pickle_categoryDict(self,category_dict,ch):
        pklF = 'category_trackDetails.ch'+str(ch)+'.dict.pkl'
        savePkl(category_dict,pklF,self.pathDir)
        print "Pickled : ",pklF
        return pklF  
    
    def countTracks(self,pickle_exists=False,nbrframes=7):
        trackFileList = [f for f in os.walk(self.outputDirpath).next()[2] if f.startswith('track_photometries')]
        assert len(trackFileList)==1
        [trackPhotometryFile] = trackFileList
        trackFile = os.path.join(self.outputDirpath,trackPhotometryFile)
        if pickle_exists:
            category_dict_ch1 = loadPkl('category_trackDetails.ch1.dict.pkl',self.pathDir)
            category_dict_ch2 = loadPkl('category_trackDetails.ch2.dict.pkl',self.pathDir)
            print "Categories Dictionaries loaded "
        else:
            category_dict_ch1 = peptracks.makeIntensityCategory(trackFile,ch=1)
            category_dict_ch2 = peptracks.makeIntensityCategory(trackFile,ch=2)
            self.pickle_categoryDict(category_dict_ch1,1)
            self.pickle_categoryDict(category_dict_ch2,2)
            print "Categories Dictionaries pickled"
    
        desired_category = peptracks.getDesiredCategory('staggered',nbrframes)
        self.desired_categoryList_ch1 = peptracks.getIntensityList_category(category_dict_ch1,desired_category)
        self.desired_categoryList_ch2 = peptracks.getIntensityList_category(category_dict_ch2,desired_category)
        
        ofname = peptracks.writeCountsFile([self.desired_categoryList_ch1,self.desired_categoryList_ch2],self.peptracks_destDir)
        print "Peptide track counts and intensity file created ",ofname
        return self.desired_categoryList_ch1, self.desired_categoryList_ch2,ofname

    def edmanEfficiency(self,ofname,edmanframe_ch1=3,edmanframe_ch2=4,calcType='aboveDecay'):
        counts_ch1 = [item[1] for item in self.desired_categoryList_ch1]
        counts_ch2 = [item[1] for item in self.desired_categoryList_ch2]
        eff1_ch1, eff2_ch1,totalPeptides_ch1 = peptracks.calculateEdmanEfficiency(counts_ch1,edmanframe_ch1,calcType)
        eff1_ch2, eff2_ch2,totalPeptides_ch2 = peptracks.calculateEdmanEfficiency(counts_ch2,edmanframe_ch2,calcType)
        ofname = peptracks.writeEdmanEfficiency(ofname,eff1_ch1,eff2_ch1,totalPeptides_ch1,edmanframe_ch1,'ch1',calcType)
        ofname = peptracks.writeEdmanEfficiency(ofname,eff1_ch2,eff2_ch2,totalPeptides_ch2,edmanframe_ch2,'ch2',calcType)
        return ofname

    def sequencingEfficiency(self,ofname,edmanframe_ch1=3,edmanframe_ch2=4,calcType='conservative'):
        counts_ch1 = [item[1] for item in self.desired_categoryList_ch1]
        counts_ch2 = [item[1] for item in self.desired_categoryList_ch2]
        eff_ch1, cleaved_ch1,totalPeptides_ch1 = peptracks.calculateSequencingEfficiency(counts_ch1,edmanframe_ch1,calcType)
        eff_ch2, cleaved_ch2,totalPeptides_ch2 = peptracks.calculateSequencingEfficiency(counts_ch2,edmanframe_ch2,calcType)
        ofname = peptracks.writeSequencingEfficiency(ofname,eff_ch1,cleaved_ch1,totalPeptides_ch1,'ch1',calcType)
        ofname = peptracks.writeSequencingEfficiency(ofname,eff_ch2,cleaved_ch2,totalPeptides_ch2,'ch2',calcType)
        return ofname


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


def runTrackExperiment(pathDir,destDir,edmanframes,nbrframes=7):
    if len(edmanframes)==2:
        edmanframe_ch1,edmanframe_ch2 = map(int,edmanframes)
    p = PeptideTracks(pathDir,destDir)
    db_ch1,db_ch2,ofname = p.countTracks(pickle_exists=False,nbrframes=nbrframes)
    p.edmanEfficiency(ofname,edmanframe_ch1,edmanframe_ch2)
    p.sequencingEfficiency(ofname,edmanframe_ch1,edmanframe_ch2)

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


def runExperimentType(dateStamp,microscopeNbr,ARGLIST=[]):
    edmanframes,nbrframes,imgTypeList,peakdensity_type= ARGLIST
    sourceDir,pathDir,destDir = getDirectory(hostname,dateStamp,microscopeNbr)   
    if edmanframes:
        runTrackExperiment(pathDir,destDir,edmanframes,nbrframes)
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
        created (defaults to \project2\jaggu\dataAnalysis)\
        There are currently three components - (a) Converting Images: random images from each sub directory and entering a frame number \
                                               (b) Peptide tracks: Outputing the counts for each of the desired tracks. Computes the edman efficiency in the two channels \ 
                                               Can be extended to change the desired category in later versions \
                                               (c) Count density: The peak counts for each field (typically random) is computed and output txt file created                                              
                                              """ )

parser.add_argument('--date', action="store", dest='dateStamp')
parser.add_argument('--scope', action="store", dest='microscopeNbr',default=1,type=int)
parser.add_argument('--edmanframes',action="store", dest='edmanframes',nargs='+',default=[],help="Edman frames for Channel-1 and/or Channel-2. eg --edmanframes 3 4")
parser.add_argument('--nbrframes',action="store", dest='nbrframes',default=7,type=int,help="Edman frames for Channel-1 and/or Channel-2. eg --edmanframes 3 4")
parser.add_argument('--images',action="store", dest='imgType',default=None,help="Converts the images into [random,all,frames]")
parser.add_argument('--peakdensity',action="store",dest='peakdensity',nargs=1,default=None,help="Adds the option of [counts or intensity] for each field in the experiment")
parser.add_argument('--image_frames',action="store",dest='imgframes',nargs='+',default=[5],help="Images whose Frames are converted to pngs: eg [10,20]")
parser.add_argument('--image_pattern',action="store",dest='imgpattern',default='*Cycle*',help="Converts images of the given pattern")
args = parser.parse_args()
ARGLIST = [args.edmanframes,args.nbrframes,[args.imgType,args.imgframes,args.imgpattern],args.peakdensity]
runExperimentType(args.dateStamp,args.microscopeNbr,ARGLIST)

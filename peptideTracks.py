#! /home/jaggu/anaconda

"""
This script is to consolidate the many graphs (counts etc) generated from running basic_experiment.py and getting the number of tracks that changed.
"""
from __future__ import division
import os
import sys
import collections
import re
import socket
import numpy as np




hostname = socket.gethostname()
if hostname == 'canopus': headDir = "/home/jaggu/marcotte_project"
else: headDir = "/project"

sys.path.append(os.path.join(headDir,'/current/project/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis'))
from commonFunctions import loadPkl, savePkl
import my_math


def makeIntensityCategory(trackFile,ch=1):
    """
    track_photometries_[dateCode].csv is trackFile. The header is 
    [0]CHANNEL,[1]FIELD,[2]H,[3]W,[4]CATEGORY,[5]FRAME 0,FRAME 1,FRAME 2,FRAME 3,FRAME 4,FRAME 5,FRAME 6,FRAME 7
    example of CATEGORY: "(True, False, False, False, True, False, False, False)"
    @function : It parses this file and makes a dictionary for the peptide tracks for each peak
    @return   : A dictionary in the form - CATEGORY: (Intensity List1),... 
    """
    query_channel = 'ch'+str(ch)
    print "Processing Channel : ",query_channel
    category_trackDetails_dict = collections.defaultdict(list)
    category_peakDetails_dict = collections.defaultdict(list) 
    ifile = open(trackFile,'r')
    lines = ifile.readlines()
    ifile.close()
    for line in lines[1:]:
        [catStr] = re.findall(r"\(.+\)",line)
        channel = line.split(',')[0]
        if channel == query_channel:
	    xy = int(line.split(',')[1])
	    hStr,wStr = line.split(',')[2], line.split(',')[3]
	    if hStr == 'None': hStr = 0
	    if wStr == 'None': wStr = 0
	    h, w = map(int, [float(hStr),float(wStr)])		
	    category_tup = eval(catStr)
	    intensity_split = re.findall('\)",[-|\d]\d*.*',line.strip())[0]
	    intensityList = map(float,intensity_split.split(',')[1:])
	    category_trackDetails_dict[category_tup].append(intensityList) 
	    category_peakDetails_dict[category_tup].append((xy, (h,w),intensityList))
    return category_trackDetails_dict,category_peakDetails_dict    


def getIntensityList_category(db,desired_categoryList):
    """
    getIntensityList_category: The dictionary containing all the category is parsed to output only the desired categories. 
    For example - the category such as [True, False,....],[True,True,....]
    The number of peaks assigned to that category and the mean of the intensity for each of the peaks when "live" is stored as a dictionary
    @ return : A list in the form - (category,number of peptide tracks, list of mean intensity of all the tracks),....
    """
    categories_intensity_list = list()
    # [Categories, counts, meanIntensityList]
    for idx, category in enumerate(desired_categoryList):
        category = tuple(category)
        intensityList = db[category]
	allMeanIntensity = list()
        for l in intensityList:
            allMeanIntensity.append(np.mean(l[0:idx+1]))
            # Calculates the mean of all the peaks which were true till that frame 
            # (like True True False): Only true true intensity values
        categories_intensity_list.append([category,len(intensityList),allMeanIntensity])
    return categories_intensity_list

def getDesiredCategory(category_type='staggered',nbrFrames=7):
    if category_type is 'staggered': 
	    desired_category = [[True]*(i+1) + [False]*(nbrFrames-i) for i in range(nbrFrames+1)]
    return desired_category

def applyCutoff(desired_categoryList,desired_category,q_range,cycle=3):
	category = tuple(desired_category[cycle])
	cutoff_intensityList = np.array(desired_categoryList[cycle][2])
	q1_intensity = np.percentile(cutoff_intensityList,q_range[0])
	q2_intensity = np.percentile(cutoff_intensityList,q_range[1])
	trackCounts = list()
	for cat,counts,intensityList in desired_categoryList:
		intensity_cutoffList = [item for item in intensityList if item > q1_intensity and item < q2_intensity]
		trackCounts.append([cat,len(intensity_cutoffList),intensity_cutoffList])
	return trackCounts

def writeCountsFile(desired_category_details_chList,destDir,nbrchannels):
    """
    @ function : Outputs the list of counts for each of the category (desired) and the mean of the mean intensity of the peaks
    @ return   : outputfile name
    """
    if nbrchannels == 2:
	    category_details_list_ch1 = desired_category_details_chList[0]
	    category_details_list_ch2 = desired_category_details_chList[1]
	    ofname = os.path.join(destDir,'category_counts.traditional.tab')
	    ofile = open(ofname,'w')
	    ofile.write('Category \t Counts \t Mean Intensity \n')
	    ofile.write('CHANNEL 1 \n')
	    print "Processing ch1"
	    for category,counts,intensityList in category_details_list_ch1:
        	ofile.write(str(category)+'\t'+str(counts)+ '\t' + str(np.mean(intensityList))+'\n')
	    ofile.write('CHANNEL 2 \n')
	    print "Processing ch2"
	    for category,counts,intensityList in category_details_list_ch2:
        	ofile.write(str(category)+'\t'+str(counts)+'\t' + str(np.mean(intensityList))+'\n')
	    ofile.close()

    else:
	    category_details_list_ch1 = desired_category_details_chList[0]
	    ofname = os.path.join(destDir,'category_counts.traditional.tab')
	    ofile = open(ofname,'w')
	    ofile.write('Category \t Counts \t Mean Intensity \n')
	    ofile.write('CHANNEL 1 \n')
	    print "Processing ch1"
	    for category,counts,intensityList in category_details_list_ch1:
        	ofile.write(str(category)+'\t'+str(counts)+ '\t' + str(np.mean(intensityList))+'\n')
	    ofile.close()
    return ofname

def calculateEdmanEfficiency(countList,edmanframe,calcType='aboveDecay'):
    """
    Calculates the edman efficiency. It requires the track counts for the truth table (works for the staggered type)
    For now I have two types - 
    
    (a) conservative: The total peptides is the counts of all counts except the first
    efficiency 1 = (lossCounts[edman frame]) / total peptides
    efficiency 2 = (lossCounts[edman frame + 1]) / total peptides

    (b) aboveDecay: The counts are modeled as an exponential decay curve. The edman frame and the next frame skipped in the modeling. 
    The number of loss at that step is the proxy if Edman did not happen. So the efficiency is -
    total peptides = absolute value of difference between the true counts - expected counts
    expected loss1 = expFit(edman frame); expected loss2 = expFit(edman frame + 1)
    efficiency  1 = (lossCounts[edman frame] - expected loss1)/total peptides
    efficicency 2 = (lossCounts[edman frame + 1]  - expected loss2)/total peptides
    returns (efficiency 1,efficiency 2)
    """
   
    def model_func(t, A, K, C):
        return A * np.exp(K * t) + C

    def fit_exp_linear(t,y,C=0):
	y = np.array(y,dtype=float)
	t = np.array(t,dtype=int)
	y = y - C
	y = np.log(y)
	K, A_log = np.polyfit(t,y,1)
	A = np.exp(A_log)
	y_fit = model_func(t,A,K,C)
	return A,K,C, y_fit

    if calcType is 'conservative':
	    edmanCleaved = countList[edmanframe]
	    totalPeptides = sum(countList[1:])
	    eff1 = countList[edmanframe]/totalPeptides
	    eff2 = countList[edmanframe+1]/totalPeptides
    if calcType is 'aboveDecay':
	    index_toRemove = [0,len(countList)-1,edmanframe,edmanframe+1]
	    ydata = [(countList[idx],idx) for idx in range(len(countList)) if not idx in index_toRemove]
	    y, t = zip(*ydata)
	    A,K,C,y_fit = fit_exp_linear(t,y)
	    y_edmanframe = countList[edmanframe]
	    restFrames = range(edmanframe,len(countList))
	    totalPeptides_edmanframe = sum([abs(countList[i]-model_func(i,A,K,C)) for i in restFrames])
	    eff1 = y_edmanframe/totalPeptides_edmanframe
	    y_nextEdmanframe = countList[edmanframe+1]
	    restFrames = range(edmanframe+1,len(countList))
	    totalPeptides_nextEdmanframe = sum([abs(countList[i] - model_func(i,A,K,C)) for i in restFrames])
	    eff2 = y_nextEdmanframe/totalPeptides_nextEdmanframe
	    #totalPeptides = sum([abs(countList[i]-model_func(i,A,K,C)) for i in range(len(countList))])
	    #eff1 = y_edmanframe/totalPeptides
	    #print eff1, eff2
	    #sys.exit(1)
	    #eff2 = y_nextEdmanframe/totalPeptides
    return (eff1,eff2,totalPeptides_edmanframe)

def writeEdmanEfficiency(fname,eff1,eff2,totalPeptides, edmanframe,chName,calcType='aboveDecay'):
    """
    Writes the Edman efficiency on the output file of the counts; It takes the output file as an appendable form
    return : file name
    """
    ofile = open(fname,'a')
    ofile.write("\n \n \n")
    ofile.write("Edman efficiency \n ")
    ofile.write("Calculation Type : "+ calcType + "\n")
    ofile.write("Channel " + chName + "\n")
    ofile.write("Edman frame  : " + "\t" + str(edmanframe) + "\n")
    ofile.write("Total peptides : "+ "\t" + str(totalPeptides) + "\n")
    ofile.write("Efficiency 1 : " + "\t" + str(eff1) + "\n")
    ofile.write("Efficiency 2 : " + "\t" + str(eff2) + "\n")
    ofile.close()
    return fname

def calculateSequencingEfficiency(countList,edmanframe,calcType='conservative'):
    # This is the similar as the conservative estimate of edman efficiency; I calculate the cleaved peptide (edman frame and next frame) / total peptides
    if calcType is 'conservative':
	    cleaved = countList[edmanframe] + countList[edmanframe+1]
	    totalPeptides = sum(countList[1:])
	    eff = cleaved/totalPeptides
	
    return (eff, cleaved, totalPeptides)

def writeSequencingEfficiency(fname,eff,cleaved,totalPeptides,chName,calcType='conservative'):
    """
    Writes the Sequencing efficiency on the output file of the counts; 
    """
    ofile = open(fname,'a')
    ofile.write("\n \n \n")
    ofile.write("Sequencing efficiency \n ")
    ofile.write("Calculation Type : "+calcType + "\n")
    ofile.write("Channel " + chName + "\n")
    ofile.write("Peptides cleaved " + str(cleaved) + "\n")
    ofile.write("Total peptides " + str(totalPeptides) + "\n" )
    ofile.write("Efficiency " + str(eff) + "\n" )
    ofile.close()
    return fname 

def testImport():
	return True

def test():
    pathDir = '/home/jaggu/marcotte_project/boulgakov/microscope/2015-Oct/2015-10-27_forHopper'
    outputDirname = "output_2015-11-02"                                                                         	
                                                                                                            	
    outputDir = os.path.join(pathDir,outputDirname)                                                             	
    destDir = '/home/jaggu/marcotte_project/current/project2/jaggu/dataAnalysis/microscope1/2015-Oct/2015-10-27'	
	
    trackPhotometryFile = 'track_photometries_nx8t5v.csv'
    trackFile = os.path.join(outputDir,trackPhotometryFile)
    #category_trackDetails_ch1 = makeIntensityCategory(trackFile,ch=1)
    #category_trackDetails_ch2 = makeIntensityCategory(trackFile,ch=2)
    #print "Categories Dictionaries pickled"

    category_dict_ch1 = loadPkl('category_trackDetails.ch1.dict.pkl',pathDir)
    category_dict_ch2 = loadPkl('category_trackDetails.ch2.dict.pkl',pathDir)
    print "Categories Dictionaries loaded "
    
    desired_category = getDesiredCategory('staggered')
    desired_categoryList_ch1 = getIntensityList_category(category_dict_ch1,desired_category)
    desired_categoryList_ch2 = getIntensityList_category(category_dict_ch2,desired_category)

    writeCountsFile([desired_categoryList_ch1,desired_categoryList_ch2])
    print "Peptide track counts and intensity file created"

    edmanframe_ch1 = 3
    edmanframe_ch2 = 4
    allCounts_ch1 = [item[1] for item in desired_categoryList_ch1]
    allCounts_ch2 = [item[1] for item in desired_categoryList_ch2]
    
    eff1_ch1,eff2_ch1 = calculateEdmanEfficiency(allCounts_ch1,edmanframe_ch1)
    eff1_ch2,eff2_ch2 = calculateEdmanEfficiency(allCounts_ch2,edmanframe_ch2)
    print "Overall Edman efficiency - \n Ch1 : ",eff1_ch1,eff2_ch1," \n Ch2 : ",eff1_ch2,eff2_ch2







#test()

#! /home/jaggu/anaconda

"""
This script is aimed at plotting boxplots of intensity of peptides for a given category; For e.g. category is (True,True,True,False); 
This applies to an experimental cycle; The category is chosen in the summarizeExperiment file and what is passed is the intensityList of each peptideTracks
"""
from __future__ import division
import os
import sys
import collections
import re
import socket
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
plt.style.use('ggplot')

hostname = socket.gethostname()
if hostname == 'canopus': headDir = "/home/jaggu/marcotte_project"
else: headDir = "/project"

sys.path.append(os.path.join(headDir,'/current/project/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis'))
from commonFunctions import loadPkl, savePkl

def applyIntensityCutoff(quartile,cycle,all_intensity_perPeptide):
	all_intensity_perPeptide = np.array(all_intensity_perPeptide)
	cutoff_intensity_perPeptide = list()
	q1,q2 = quartile
	all_intensity_perFrame = all_intensity_perPeptide.transpose()
	intensity_atFrame = all_intensity_perFrame[cycle]
	#print np.median(intensity_atFrame)
	cutoff1 = np.percentile(intensity_atFrame,q1)
	cutoff2 = np.percentile(intensity_atFrame,q2)
	for intensity_list in all_intensity_perPeptide:
		if intensity_list[cycle] > cutoff1 and intensity_list[cycle] < cutoff2:
			cutoff_intensity_perPeptide.append(intensity_list)
	return np.array(cutoff_intensity_perPeptide)

def plot_boxPlots(sourceDir,dest_imgDir,all_pepTrack_intensity,title):
	print sourceDir, dest_imgDir, len(all_pepTrack_intensity), title
	cycle = 2 #Making the 4th Mock as the quartile for cutoff
	# Need to change how I enter labels (later)
	labels = ['Mock']*3 +['Edman']* 9
	for quartile in [(0,100),(0,25),(25,50),(50,75),(75,100)]:
		title = title + '_' + str(quartile)
		intensity_perPeptide = applyIntensityCutoff(quartile,cycle,all_pepTrack_intensity)
		bp = plt.boxplot(intensity_perPeptide,patch_artist=True, showfliers=False)
		for box in bp['boxes']:
			# change outline color
			#box.set( color='#7570b3', linewidth=2)
			# change fill color
			box.set( facecolor = '#0072b2')
		plt.xticks(range(1,len(labels)+1),labels,rotation='vertical')
		plt.title(title)
		plt.ylabel('Intensity of peptides')
		#plt.ylim([0,90000])
		plt.tight_layout()
	        new_imgDir = os.path.join(dest_imgDir,'allIntensity')
	        if not os.path.exists(new_imgDir): os.makedirs(new_imgDir)
		fname = os.path.join(new_imgDir,title)
		plt.savefig(fname + '.png')
		plt.savefig(fname + '.svg')
		plt.close()
	return True

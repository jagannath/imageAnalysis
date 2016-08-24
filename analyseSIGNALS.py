#! /home/jaggu/anaconda

"""
The idea of the script is to perform analysis on track**.SIGNALS.pkl
I will need to pass this file as argument. As second argument (optional) will output the results from the relevant function.
As a starting point, it will make a new directory in the SIGNALS directory and output the csv or txt.
"""
import sys
import os
import argparse
import cPickle as pickle
import collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def freq_stepDrops(signal_filepath):
    freqDrops_dict = collections.defaultdict(int)
    signal_dict = pickle.load(open(signal_filepath,'r'))
    headDir,signal_fname = os.path.split(signal_filepath)
    destDir = os.path.join(headDir,'signals_analysis')
    if not os.path.exists(destDir): os.makedirs(destDir)
    for drops, counts in signal_dict.items():
        nbrDrops = len(drops[0])
        status = drops[1] #True means drop is to background intensity
        # Counting the remainders only.Drop is not to background. So - ((A,0),False),#)
        if nbrDrops is 1: # Calculating 0 drop case
            dropPos = drops[0][0][1]
            if dropPos == 0:
                nbrDrops = 0
        
        freqDrops_dict[nbrDrops]+=counts
    
    ofname = os.path.join(destDir,signal_fname+'.freqDrops.tab')
    ofile = open(ofname,'w')
    ofile.write('Drop Number \t Frequency \n')
    for nbrDrops, freq in freqDrops_dict.items():
        ofile.write(str(nbrDrops) + '\t' + str(freq) + '\n')

    ofile.close()
    return ofname

def get_stepCounts(signal_dict):
    """Filters the signals to a list of one and two Step drop counts """
    twoStepCounts_list = list()
    oneStepCounts_list = list()
    for drops, counts in signal_dict.items():
        nbrDrops = len(drops[0])
        status = drops[1]
	if nbrDrops is 1 and status is True:
	    oneStepCounts_list.append([drops,counts])
	if nbrDrops is 2 and status is True:
            twoStepCounts_list.append([drops,counts])
    return oneStepCounts_list, twoStepCounts_list

def heatMap_twoStep(twoStepCounts_list):
    """ Creates a heat map in the signal_analysis directory"""
    # Loading labels
    nbrEdman = args.nbrEdman
    nbrMock = args.nbrMock
    nbrOmit = args.nbrOmit
    sequence = args.sequence
    channel = args.channel

    nrows = ncols = (nbrMock - nbrOmit ) + nbrEdman + 1
    twoDrops_mat = np.zeros([nrows,ncols])
    for drops,counts in twoStepCounts_list:
        first_drop = drops[0][0][1]
        second_drop = drops[0][1][1]
        twoDrops_mat[first_drop,second_drop] = counts
  
    if channel is 1:
	cmap_color = 'YlOrRd'
	channelColor = "561"

    if channel is 2:
        cmap_color = 'YlGnBu'
	channelColor = "647"
    
    fig, ax = plt.subplots()
    heatmap_array = twoDrops_mat[1:,1:]
    labels = ['E'+str(item+1) for item in range(nrows-1)]
    #labels.extend('R')
    text_limit = np.amax(heatmap_array)
    for (i,j),counts in np.ndenumerate(heatmap_array):
        if counts > (text_limit)*0.75: textColor = 'white'
        else: textColor = 'black'
        ax.text(j,i, '{:0.0f}'.format(counts),ha='center',va='center',color=textColor)
    

    # Color map changing accordingly
    count_min = np.amin(heatmap_array)
    count_max = np.amax(heatmap_array)

    cax = ax.matshow(heatmap_array,cmap=cmap_color,vmin=count_min,vmax=count_max)
    fig.colorbar(cax)
    
    ax.set_xticks(range(nrows))
    ax.set_yticks(range(nrows))
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)
    title = 'Double drop' + '('+channelColor+'):'+ sequence
    ax.set_title(title)

    fname = os.path.join(destDir,title+'_'+channelColor+'_twoStepDrops')
    plt.savefig(fname+'.png')
    plt.savefig(fname+'.svg')
    plt.close()

    return fname




parser = argparse.ArgumentParser(description="""This is the script meant to summarize the results from "track_photometries_NO_NONES*csv_ch2_SIGNALS.pkl""")
#parser.add_argument('--signal_file','-f', action="store", dest='signal_file',type=str)
parser.add_argument('signal_file',action="store",type=str)
parser.add_argument('--mock','-m',action="store",dest="nbrMock",type=int,default=4)
parser.add_argument('--omit','-o',action="store",dest="nbrOmit",type=int,default=4)
parser.add_argument('--edman','-e',action="store",dest="nbrEdman",type=int,default=4)
parser.add_argument('--channel','-c',action="store",dest="channel",type=int,default=2)
parser.add_argument('--sequence','-s',action="store",dest="sequence",type=str,default="Peptide")
args = parser.parse_args()

signal_filepath = os.path.abspath(args.signal_file)

signal_dict = pickle.load(open(signal_filepath,'r'))
headDir,signal_fname = os.path.split(signal_filepath)
destDir = os.path.join(headDir,'signals_analysis')
if not os.path.exists(destDir): os.makedirs(destDir)

#ofname = freq_stepDrops(os.path.abspath(args.signal_file))
sequence = 'GC*[A647N]AGC*[A647N]AGAG'
sequence = args.sequence

oneStepCounts_list, twoStepCounts_list = get_stepCounts(signal_dict) 

#heatMap_oneStep(oneStepCounts_list)
heatMap_twoStep(twoStepCounts_list)






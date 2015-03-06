#! /home/jaggu/anaconda/bin/python2.7 

from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
import os
import sys
sys.path.append('/project/marcotte/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis')
from commonFunctions import locate,simpleShow

import cPickle as pickle
import collections
import time

# MATPLOTLIB DEFAULT PARAMS
from matplotlib import rcParams
rcParams['axes.labelsize'] = 9
rcParams['xtick.labelsize'] = 9
rcParams['ytick.labelsize'] = 9
rcParams['legend.fontsize'] = 9
rcParams['font.size'] = 10


# TICKS
from matplotlib.ticker import MaxNLocator
my_xlocator = MaxNLocator(5)
my_ylocator = MaxNLocator(3)





class RadialProfile:
    """
    This is the major class where the radial profile of the intensity is
    plotted. It ensures that the right size, dpi of the figure is maintained. 
    """
    def __init__(self,pathDir):
        plt.close()
        self.fig = plt.figure(figsize = ( 5,5),dpi=300)
        self.pathDir = pathDir
        colors = (('purple','#ae1e6f'),
                  ('red','#d92525'),
                  ('orange','#f25c05'),
                  ('yellow','#f29f05'),
                  ('green','#88a61b'),
                  ('blue','#0e3d59'))
        self.selcolor = dict(colors)

    def writeNbrCircles(self,nbrCircles,ax):
        s = 'n = '+str(nbrCircles)+' beads'
        ax.text( 0.7,0.9,s,
                horizontalalignment='center',
                verticalalignment='center',
                transform = ax.transAxes
               )
        return ax
    
    def drawProfile(self,xyDict,fname,nbrCircles=int(),col='blue'):
        type1 = '.png'
        type2 = '.svg'
        ax = self.fig.add_subplot(111)
        xList = xyDict.keys()
        normXList = [x/50 for x in xList]
        yList = xyDict.values()
        ax.plot(normXList,yList,color=self.selcolor[col],lw='4')
        if nbrCircles: ax = self.writeNbrCircles(nbrCircles,ax)
        ax.set_ylim([0,60000])
        ax.xaxis.set_major_locator(my_xlocator)
        ax.yaxis.set_major_locator(my_ylocator)
        ax.set_xlabel('Normalized Radial Distance')
        ax.set_ylabel('Mean Intensity (Counts/pixel)')
        ax.set_adjustable('box')
        print "Saving File %s"%(fname+type1)
        self.fig.savefig(fname+type1) 
        self.fig.savefig(fname+type2)
        plt.close()

    def drawProfileWstdev(self,xList,yList,yStd_plus,yStd_minus,fname,nbrCircles=int(),ylim=[0,60000],col='red'):
        type1 = '.png'
        type2 = '.svg'
        ax = self.fig.add_subplot(111)
        ax.plot(xList,yList,color=self.selcolor[col],lw='4')
        ax.fill_between(xList,yStd_plus,yStd_minus,color=self.selcolor[col],alpha=0.5)
        if nbrCircles: ax = self.writeNbrCircles(nbrCircles,ax)        
        ax.set_ylim(ylim)
        ax.xaxis.set_major_locator(my_xlocator)
        ax.yaxis.set_major_locator(my_ylocator)
        ax.set_xlabel('Normalized Radial Distance')
        ax.set_ylabel('Normalized Intensity ')
        ax.set_adjustable('box')
        print "Saving File %s"%(fname+type1)
        self.fig.savefig(fname+type1) 
        self.fig.savefig(fname+type2)
        plt.close()


    def combinePlots(self,packList,fname):
        type1 = '.png'
        type2 = '.svg'
        ax = self.fig.add_subplot(111)
        for unpack in packList:
            pklF,ls,col = unpack
            nBeads,xyDict = pickle.load(open(pklF))
            label = pklF.split('_')[4]+'_'+pklF.split('_')[5]
            xList = xyDict.keys()
            yList = xyDict.values()
            normXList = [x/50 for x in xList]
            print max(yList)
            ax.plot(normXList,yList,linestyle=ls,lw='4',color=self.selcolor[col],label=label)

        ax.set_ylim([0,60000])
        ax.xaxis.set_major_locator(my_xlocator)
        ax.yaxis.set_major_locator(my_ylocator)
        ax.set_xlabel('Normalized Radial Distance')
        ax.set_ylabel('Mean Intensity (Counts/pixel)')
        ax.set_adjustable('box')
        handles,labels = ax.get_legend_handles_labels()
        ax.legend(handles,labels)
        print "Saving File %s"%(fname+type1)
        self.fig.savefig(fname+type1) 
        self.fig.savefig(fname+type2)
        plt.close()


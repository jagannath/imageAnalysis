#! /home/jaggu/anaconda

from __future__ import division
import sys
import os
import matplotlib.pyplot as plt
import numpy as np


def normalizeMatrix(mat):
    newmat = np.array(mat)
    normVal = 152
    normMat = np.zeros((8,8),dtype=int)
    maxVal = np.max(newmat)
    for (i,j),counts in np.ndenumerate(newmat):
        if not counts == 0:
            newCounts = (counts/maxVal)*normVal
            normMat[(i,j)]=newCounts
    return normMat


mat = np.zeros((8,8))
mat = [[0,161,69,50,83,30,16,12],
        [0,0,120,94,213,76,45,35],
        [0,0,0,40,62,32,16,11],
        [0,0,0,0,48,17,10,8],
        [0,0,0,0,0,14,8,6],
        [0,0,0,0,0,0,1,1],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0]]

normMat = normalizeMatrix(mat)

heatmap_array = np.array(normMat)
print heatmap_array.size

fig, ax = plt.subplots()
labels = ['E'+str(item) for item in range(9)]
text_limit = np.amax(heatmap_array)
for (i,j),counts in np.ndenumerate(heatmap_array):
    if counts > (text_limit)*0.75: textColor = 'white'
    else:
        textColor = 'black'
    ax.text(j,i, '{:0.0f}'.format(counts),ha='center',va='center',color=textColor)

cmap_color = 'YlGnBu'

# Color map changing accordingly
count_min = np.amin(heatmap_array)
count_max = np.amax(heatmap_array)
count_min = 0
count_max = 150

cax = ax.matshow(heatmap_array,cmap=cmap_color,vmin=count_min,vmax=count_max)
fig.colorbar(cax)
ax.set_xticklabels(labels)
ax.set_yticklabels(labels)
title = 'Double drop' + '(ch2):'
ax.set_title(title)

destDir ='/home/jaggu/marcotte_project/jagannath/projectfiles/alexanders_results' 
fname = os.path.join(destDir,title+'_simulatedTwoStepsDrops')
plt.savefig(fname+'.svg')

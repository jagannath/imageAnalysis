#! /home/jaggu/anaconda/python2.7

"""
Parses the results of the peak csv files for the count of peptides. Then makes a
heat map on the count vs xy
"""

import os
import sys
import fnmatch
import numpy as np
import matplotlib.pyplot as plt

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below supplied root directory.'''
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            allFiles.append(os.path.join(path, filename))
    return allFiles

def getCount(out):
    count = out.split('\n')[1].split('\t')[2]
    return count

def getXY(xy):
    # This is complicated because the delimiter was '-'. The negative signs gets
    # split as well. So this elaborate logic
    xysplit = xy.split('-')
    if len(xysplit)>2:
        if 'x=' in xysplit:
            xidx = xysplit.index('x=')
            xx='-'+xysplit[xidx+1]
        else: xx = xysplit[0].split('=')[-1]
        if 'y=' in xysplit:
            yidx = xysplit.index('y=')
            yy = '-'+xysplit[yidx+1]
        else: yy = xysplit[1].split('=')[-1]
    else:
        xx,yy = map(lambda foo:foo.split('=')[-1],xysplit)
    return (xx,yy)


def getPos(out):
    fname = out.split('\n')[1].split('\t')[1].split('/')[-1]
    (i,j) = map(lambda x:x[2:],fname.split('.')[1].split('-'))
    xy = fname.split('.')[0].split('_')[1]
    (x,y) = getXY(xy)
    return [(x,y),(i,j)]


def sumFiles(DIR):
    allCSVS = locate('peak*.csv',DIR)
    ofname = analyseDIR + '/' + 'heatMap.txt'
    allOutput = str()
    ijCountDICT = dict()
    for csv in allCSVS:
        with open(csv) as f:
            out = f.read()
            allOutput+=out
            count = getCount(out)
            [(x,y),(i,j)] = getPos(out)
            ijCountDICT[(i,j)] = [(x,y),count]
    with open(ofname,'w') as f:
        f.write(allOutput)
    return ijCountDICT

def createMatrix(ijCountDICT):
    mat = np.zeros(shape=( 50, 20)) #Making a matrix of zeros
    for k,v in ijCountDICT.items():
        (i,j) = map(int,[k[0],k[1]])
        count = v[1]
        mat[i,j] = count
    tmat = np.transpose(mat)
    for index, count in np.ndenumerate(tmat):
        if count == 0:
            count = np.mean(mat)
        tmat[index] = count
    return tmat

def heatMap(mat,fname):
    fig1 = plt.figure(figsize=( 4,4))
    ax = fig1.add_subplot(111)
    im = ax.imshow(mat)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    im.set_cmap('YlGnBu_r')
    cb = plt.colorbar(im,orientation='horizontal',ticks=[0,500,1000,1500,2000])
    cb.set_label('Counts/field')
    plt.setp(cb.ax.get_xticklabels(),visible=True)
    plt.savefig(fname,dpi=1000)
    return plt




analyseDIR = '/project/marcotte/jagannath/projectfiles/SingleMoleculeMicroscopy/dataAnalysis/2013-July/2013-07-04'
inputDIR = '/project2/marcotte/boulgakov/microscope/2013-July/2013-07-04/B2111_20pMCy5Dye/'
ijCountDICT = sumFiles(inputDIR)
countMat = createMatrix(ijCountDICT)
print "Total Counts : %s "%str(np.sum(countMat))
print "Average Counts/Field : %s "%str(np.mean(countMat))
print "Standard Deviation Counts : %s"%str(np.std(countMat))
fname = analyseDIR+'/heatmap.png'
heatMap(countMat,fname)




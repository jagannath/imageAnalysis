#! /home/jaggu/anaconda/bin/python2.7

"""
This file contains common custom functions that were made and used repeatedly
"""
from __future__ import division
import os
import fnmatch 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and
    below supplied root directory.'''
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            allFiles.append(os.path.join(path,filename))
    return allFiles

def simpleShow(img):
    # Shows img when it is a matrix
    plt.imshow(img,cmap='gray')
    plt.show()
    return True

def heatMapImage(img,cmap='Reds',transpose=False):
    # Given an image; It is a grayscale; Makes a heatmap; It doesnt transform
    # the image; The numpy array reads as (h,w); 
    plt.imshow(img,cmap)
    plt.show()
    return True

def colorImage(img,color='Reds',stretch='linear'):
    # Making my own color map
    if color is 'black_green':
        cdict = {'red': ((0.0,  0.0, 0.0),
                        (0.5,  0.0, 0.0),
                        (1.0,  0.0, 0.0)),
          
              'green': ((0.0,  0.0, 0.0),
                        (0.25, 0.0, 0.0),
                        (0.75, 1.0, 1.0),
                        (1.0,  1.0, 1.0)),
          
              'blue':  ((0.0,  0.0, 0.0),
                        (0.5,  0.0, 0.0),
                        (1.0,  0.0, 0.0))}
        colormap = LinearSegmentedColormap('black_green',cdict)
        
    elif color is 'black_red':
        cdict = {'red': ((0.0,  0.0, 0.0),
                        (0.25, 0.0, 0.0),
                        (0.75, 1.0, 1.0),
                        (1.0,  1.0, 1.0)),
                    
              'green': ((0.0,  0.0, 0.0),
                        (0.5,  0.0, 0.0),
                        (1.0,  0.0, 0.0)),
      
              'blue':  ((0.0,  0.0, 0.0),
                        (0.5,  0.0, 0.0),
                        (1.0,  0.0, 0.0))}
        colormap = LinearSegmentedColormap('black_red',cdict)

    else: colormap = color

    hdim,wdim = img.shape
    if stretch is 'linear': img_map = (img-img.min())/(img.max()-img.min())
    elif stretch is 'sigma':
        alpha= (img.max() - img.min())
        beta = alpha/8
        # Inew = (newMax-newMin)*(1/(1+exp(-(I-beta)/alpha))) + newMin
        img_map = 1/(1 + np.exp(-1*(img-beta)/alpha))
    else: img_map = img/img.max()


    # Need to normalize the img; 0-1 for the cmap to work
    cmap = plt.get_cmap(colormap)
    rgba_img = cmap(img_map,bytes=True)
    rgb_img = np.delete(rgba_img,3,2)
    return rgb_img

def makeCmap():
    print "making cmap"
    cdict = {'green': ((0.0,  0.0, 0.0),
                    (0.5,  0.0, 0.0),
                    (1.0,  0.0, 0.0)),
          
          'red': ((0.0,  0.0, 0.0),
                    (0.25, 0.0, 0.0),
                    (0.75, 1.0, 1.0),
                    (1.0,  1.0, 1.0)),
          
          'blue':  ((0.0,  0.0, 0.0),
                    (0.5,  0.0, 0.0),
                    (1.0,  0.0, 0.0))}


    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.cm as cm
    import sys
    newcdict = cm.get_cmap('gist_heat')._segmentdata
    #print newcdict['red'][0]

    x = np.arange(0, np.pi, 0.1)
    y = np.arange(0, 2*np.pi, 0.1)
    X, Y = np.meshgrid(x,y)
    Z = np.cos(X) * np.sin(Y) * 10

    my_cmap = LinearSegmentedColormap('new_cmap',cdict)
    plt.imshow(Z,cmap=my_cmap)
    plt.colorbar()
    plt.show()

makeCmap()



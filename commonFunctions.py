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
import cv2

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and
    below supplied root directory.'''
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            allFiles.append(os.path.join(path,filename))
    return allFiles

def simpleOpen(fname,flag=-1):
    orig = cv2.imread(fname,flag)
    corig = cv2.imread(fname,0)
    cimg = cv2.cvtColor(corig,cv2.COLOR_GRAY2BGR)
    return (orig, corig, cimg)

def makeDir(dirname):
    if not os.path.exists(dirname):os.makedirs(dirname)
    return True

def bounds(x):
    if x < 0: x = 0
    if x > 511: x = 511
    return x


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
        cdict = {'red': ((0.0,  0.0, 1.0),
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

    elif color is 'black_orange':
        cdict = {'green': ((0.0,  0.0, 0.0),       
                        (0.15, 0.0, 0.0),
                        (1.0,  0.57, 0.57)),
          
              'red':   ((0.0,  0.0, 0.0),
                       (0.15, 0.0, 0.0),
                       (0.75, 1.0, 1.0),
                       (1.0,  1.0, 1.0)),
          
              'blue':  ((0.0,  0.0, 0.0),
                       (0.5,  0.0, 0.0),
                       (1.0,  0.0, 0.0))}
        colormap = LinearSegmentedColormap('black_orange',cdict)
    elif color is 'black_blue':
        cdict = {'green': ((0.0,  0.0, 0.0),
                           (0.15, 0.0, 0.0),
                           (0.75, 0.637, 0.637),
                           (1.0,  0.637, 0.630)),
            
              'red':   ((0.0,  0.0, 0.0),
                        (0.15, 0.0, 0.0),
                        (0.75, 0.05, 0.05),
                        (1.0,  0.055, 0.055)),
          
              'blue':  ((0.0,  0.0, 0.0),
                        (0.15, 0.0, 0.0),
                        (0.75, 0.95, 0.95),
                       (1.0,  0.975, 0.975))}
        colormap = LinearSegmentedColormap('black_blue',cdict)

    elif color is 'white_blue':
        cdict = {'green': ((0.0,  0.0, 0.0),       
                           (0.05, 0.49, 0.49),
                            (0.75, 0.0, 0.0),
                            (1.0,  0.0, 0.0)),
          
                  'red':   ((0.0,  0.0, 0.0),
                           (0.05, 0.29, 0.29),
                           (0.75, 0.0, 0.0),
                           (1.0,  0.0, 0.0)),
          
                  'blue':  ((0.0,  0.0, 0.0),
                            (0.05, 0.98, 0.98),
                            (0.75, 1.0, 1.0),
                           (1.0,  1.0, 1.0))}
        colormap = LinearSegmentedColormap('white_blue',cdict)       
    elif color is 'white_orange':
        cdict = {'green': ((0.0, 0.0, 0.0),              
                        (0.05, 0.92, 0.92),
                        (0.75, 0.49,0.49),
                        (1.0,  0.49, 0.49)),
          
              'red':   ((0.0, 0.0, 0.0),
                       (0.05, 1.0, 1.0),
                       (0.75, 1.0, 1.0),
                       (1.0,  1.0, 1.0)),
          
              'blue':  ((0.0,  0.0, 0.0),
                       (0.05,  0.4, 0.4), 
                       (1.0, 0.0, 0.0))}                       
        colormap = LinearSegmentedColormap('white_orange',cdict)          
    elif color is 'black_yellow':
        cdict = {   'red':      ((0.0,  0.0, 0.0),       
                                (0.65, 1.0, 1.0),
                                (0.95, 1.0, 1.0),
                                (1.0,  1.0, 1.0)),
          
                    'green':    ((0.0,  0.0, 0.0),
                                (0.65, 0.655, 0.655),#FFAA00
                                (0.95, 0.749, 0.749),#FFBF00
                                (1.0,  1.0, 1.0)),#FFFF00
          
                    'blue':     ((0.0,  0.0, 0.0),
                                (0.65, 0.0, 0.0),
                                (0.95, 0.0, 0.0),
                                (1.0,  0.0, 0.0))}
        colormap = LinearSegmentedColormap('black_yellow',cdict)          
    elif color is 'black_purple':
        cdict = {   'red':      ((0.0,  0.0, 0.0),       
                                (0.75, 0.0, 0.0),
                                (0.95, 0.082, 0.082),
                                (1.0,  0.51, 0.51)),
          
                    'green':    ((0.0,  0.0, 0.0),
                                (0.75, 0.0, 0.0),
                                (0.95, 0.14, 0.14),
                                (1.0,  0.0, 0.0)),
          
                    'blue':     ((0.0,  0.0, 0.0),
                                (0.75, 1.0, 1.0),#0000F
                                (0.95, 0.584, 0.584),
                                (1.0,  0.29, 0.29))} 
        colormap = LinearSegmentedColormap('black_purple',cdict)          
    
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
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.cm as cm
    import sys

    print "making cmap"
    print "making black-yellow"
    cdict = {   'red':      ((0.0,  0.0, 0.0),       
                            (0.25, 1.0, 1.0),
                            (0.75, 1.0, 1.0),
                            (1.0,  1.0, 1.0)),
          
                'green':    ((0.0,  0.0, 0.0),
                            (0.25, 0.655, 0.655),#FFAA00
                            (0.75, 0.749, 0.749),#FFBF00
                            (1.0,  1.0, 1.0)),#FFFF00
          
                'blue':     ((0.0,  0.0, 0.0),
                            (0.25, 0.0, 0.0),
                            (0.75, 0.0, 0.0),
                            (1.0,  0.0, 0.0))}

    x = np.arange(0, np.pi, 0.1)
    y = np.arange(0, 2*np.pi, 0.1)
    X, Y = np.meshgrid(x,y)
    Z = np.cos(X) * np.sin(Y) * 10

    my_cmap = LinearSegmentedColormap('new_cmap',cdict)
    plt.imshow(Z,cmap=my_cmap)
    plt.colorbar()
    plt.show()

    print "making black-purple"
    cdict = {   'red':      ((0.0,  0.0, 0.0),       
                            (0.25, 0.0, 0.0),
                            (0.75, 0.082, 0.082),
                            (1.0,  0.51, 0.51)),
          
                'green':    ((0.0,  0.0, 0.0),
                            (0.25, 0.0, 0.0),
                            (0.75, 0.14, 0.14),
                            (1.0,  0.0, 0.0)),
          
                'blue':     ((0.0,  0.0, 0.0),
                            (0.25, 1.0, 1.0),#0000F
                            (0.75, 0.584, 0.584),
                            (1.0,  0.29, 0.29))} 
    plt.close()
    my_cmap = LinearSegmentedColormap('new_cmap',cdict)
    plt.imshow(Z,cmap=my_cmap)
    plt.colorbar()
    plt.show()








#makeCmap()



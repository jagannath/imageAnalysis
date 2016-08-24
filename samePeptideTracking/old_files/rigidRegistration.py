#! /home/jaggu/anaconda/bin/python2.7
'''
@author: Christian Rossmann, PhD
@license:  Public Domain
@blog: http://scientificcomputingco.blogspot.com/
'''

import numpy as np
import Image
from scipy import ndimage
from scipy import optimize
from scipy import misc

from pylab import *

def MeasureErr(img1,img2):
    diff = (img1-img2)
    return sum(diff**2)

def RigidRegistration(img,ximg):
    
    # Perform initial guess rotation & translation 
    v_range =  np.array(xrange(-10,10))
    
    err = np.array([MeasureErr(img,ndimage.shift(ximg,(v,0))) for v in v_range])
    x = v_range[where(err==err.min())[0]]
    
    err = np.array([MeasureErr(img,ndimage.shift(ximg,(0,v))) for v in v_range])
    y = v_range[where(err==err.min())[0]]

    err = np.array([MeasureErr(img,ndimage.rotate(ximg,v,reshape=0)) for v in v_range])
    r = v_range[where(err==err.min())[0]]

    # List contains displacement in x and y and rotation
    
    param = [x,y,r]
    
    def ErrFunc(param,img=img,ximg=ximg):
        
        # Perform rotational and translational transformation
        
        _img = ximg.copy()
        _img = ndimage.rotate(_img,param[2],reshape=0)
        _img = ndimage.shift(_img,param[:2])
        
        return MeasureErr(img,_img)

    
    param = optimize.fmin(ErrFunc,param)
    
    #Final transformation
    _img = ximg.copy()
    _img = ndimage.rotate(_img,param[2],reshape=0)
    _img = ndimage.shift(_img,param[:2])
    
    return (_img,param)


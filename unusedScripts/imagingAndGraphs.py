#! /home/jaggu/anaconda/bin/python2.7

"""
This script performs advanced image drawing and representation and draws graphs
as need be. 
"""

import scipy.misc
import numpy as np
import os
import sys
from subprocess import call
#import pflibMP
import cv2
import fnmatch
import re
import cPickle as pickle
import collections
import shelve
import time
from PIL import Image, ImageDraw, ImageFilter

dateStamp = "2014-08-28"
sourceDir = "/project2/marcotte/boulgakov/microscope/2014-Aug/2014-08-28/AS_Atto647N_2zM_48hlater_647_OSS_COT40_40Perc_flds_trace001"
#name = 'AS_Atto647N_2zM_48hlater_647_OSS_COT40_40Perc_flds_trace001t001.tif_scaled.png'
name = 'AS_Atto647N_2zM_48hlater_647_OSS_COT40_40Perc_flds_trace001t001.tif'
inputF = os.path.join(sourceDir,name)

destDir = os.path.join("/project2/marcotte/jaggu/dataAnalysis/microscope1/2014-Aug/",dateStamp,"peakResults","images",name)
if not os.path.exists(destDir): os.makedirs(destDir)

x,y = ( 250, 256)
def test_draw():
    img = Image.open(inputF)
    print img.format, img.size, img.mode
    img2 = img.convert('I')
    print img2.format, img2.size, img2.mode
    tup = (x-10,y-20,x+10,y+20)
    img3 = img2.crop(tup).resize(( 200,200))
    print img3.format, img3.size, img3.mode
    draw = ImageDraw.Draw(img3)
    a,b = ( 100,100)
    draw.rectangle([a-5,b-5,a+5,b+5],fill = "white", outline="blue")
#    draw.rectangle([x-5,y-5,x+5,y+5],fill="white", outline="blue")

    fname = os.path.join(destDir,name)
    #img.save(fname,"PNG") 
    img3.save(fname+'RECT_CROP.PNG',"PNG")

test_draw()


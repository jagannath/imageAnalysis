#!/home/jaggu/anaconda/bin/python

"""
The function of the script is to extract the raw images from the different
cycles that corresponds to the defined peak (position). 
:w:The main input parameters required are 
1. Channel - ch1, ch2,.. ;Name defined in the track photometry csv files
2. Field (frame) = 0,1...99..; This is the number defined by Alex and does not
typically correspond to the naming convention used in the file. In my case, the
file name starts from ..xy001.. and thus field will correspond to xy+1 image
file
3. Peak position (h,w); This is defined by Alex in his output files and
corresponds to the peak position at the first frame. 
4. Name pattern: This is the file name pattern that will be used to retrieve all
the images corresponding to the Edman cycle [such as 160621_4Mock4Edman_cycles*.tif]
5. offset dictionary: This is one of the output from the basic_experiment_script
and corresponds to the offset values that needs to be applied for the
field after successive cycle. 

A number of other ancillary inputs are defined such as
1. boxsize - Dimension of the cropped image (in terms of pixels); Default of 50
2. peakbox - A golden colored square is drawn around each peak; Default of 10x10px
3. output_directory - If not defined, it makes a directory in the current
directory as "extractedPeakImages"
4. 
"""
from __future__ import division
import sys
import cPickle as pickle
import os
import fnmatch
import numpy as np
import matplotlib.pyplot as plt
import cv2
import argparse
from matplotlib.patches import Rectangle
import math

# Header functions #
def savePkl(db,fpath):
    ofile = open(fpath,'w')
    pickle.dump(db,ofile)
    ofile.close()
    return True

def loadPkl(fpath):
    ifile = open(fpath,'r')
    db = pickle.load(ifile)
    ifile.close()
    return db

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and
    below supplied root directory.'''
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            allFiles.append(os.path.join(path,filename))
    return allFiles

def openImage(fname,flag=-1):
    orig = cv2.imread(fname,flag)
    corig = cv2.imread(fname,0)#Always to 8bit color
    cimg = cv2.cvtColor(corig,cv2.COLOR_GRAY2BGR)
    return (orig, corig, cimg)

def makeDir(d):
    if not os.path.exists(d): os.makedirs(d)
    return d

# End of header functions #

class ExtractPeakImages:
    def __init__(self,channel,field,(h0,w0),fnamePattern,offset_fpath,sourceDir,destDir):
        self.ch = channel
        self.fld = field #Field defined by the analysis (starts from 0)
        self.h0,self.w0 = h0,w0
        self.fnamePattern = fnamePattern
        self.offset_file = offset_fpath
        self.sourceDir = sourceDir
        self.destDir = destDir

    def get_allImagePath(self):
        """ Finds all images across cycles and creates a dictionary pair
        cycle number : imagePath
        """
        def _findAllFiles(xy_name):
            imageName_pattern = namePattern + 'xy'+xy_name + 'c'+str(ch)+'.tif'
            allCycles = locate(imageName_pattern, sourceDir)
            return allCycles

        allCycles = list()
        cycle_imagePath_dict = dict()

        namePattern = self.fnamePattern
        fld = self.fld
        ch = self.ch
        sourceDir = self.sourceDir
        
        # Retrieving files based on differences in the xy format 01 or 001
        xy_name = '{0:02}'.format(fld+1)
        allCycles = _findAllFiles(xy_name)
        if not allCycles:
            xy_name = '{0:03}'.format(fld+1)
            allCycles = _findAllFiles(xy_name)
        sorted_allCycle_images = sorted(allCycles)
        # Populating dictionary
        for cycle, fname in enumerate(sorted_allCycle_images):
            cycle_imagePath_dict[cycle] = fname

        return cycle_imagePath_dict

    def get_allOffsets(self):
        """ Uses the offset_dictionary to create a dictionary pair
        cycle number: offset (from the first image)
        """
        cycle_offset_dict = dict()
        offsetDict_file = self.offset_file

        offset_dict = loadPkl(offsetDict_file)
        xy = self.fld 
        ch = 'ch' + str(self.ch) #The analysed file reads as ch1 etc
        offset = ( 0, 0)

        for cycle in sorted(offset_dict.keys()):
            xy_dict = offset_dict[cycle]
            new_offset = xy_dict[xy][ch]
            offset = (offset[0]+new_offset[0], offset[1]+new_offset[1])
            cycle_offset_dict[cycle] = offset

        return cycle_offset_dict


    def getExactPosition(self,h0,w0,imagePath):
        """ Wiggles around to find the exact position of the peak """

        def _getClosest((h0,w0),allPeakPos):
            distances_list = [((math.sqrt((float(h)-float(h0))**2+(float(w)-float(w0))**2)),(h,w))
                             for h,w in allPeakPos]
            distance, nearest_hw = sorted(distances_list)[0]
            if distance > 2: nearest_hw = (h0,w0)
            return nearest_hw

        # Uses the pickle file to find all the peak positions and determine the
        # nearest
        cycleDir,imageFile = os.path.split(imagePath)
        psfs_file = [os.path.join(cycleDir,f) for f in os.listdir(cycleDir) 
                    if f.startswith(imageFile) and f.endswith('.pkl')]
        try:
            assert len(psfs_file) == 1
        except AssertionError: 
            systemExit("There are more psfs files to determine exact position")

        [psfs_fpath] = psfs_file
        allPeaks_db = loadPkl(psfs_fpath)
        allPeakPos = allPeaks_db.keys()
        new_h, new_w = _getClosest((h0,w0),allPeakPos)
        return new_h, new_w

    def getImageBox(self, imagePath,offset,title,cycle, intensityScaling, cropsize):
        """ Crops the image to the appropriate boxsize and rescales to a defined
        value. In further experiments,I will tinker with the intensity scaling
        to match the first of the images """

        destDir = self.destDir
        ch = self.ch
        h, w = map(float, [self.h0, self.w0])
        v_min, v_max = intensityScaling

        orig, corig, cimg = openImage(imagePath)
        
        new_h = h - offset[0]
        new_w = w - offset[1]

        h_exact, w_exact = self.getExactPosition(new_h, new_w, imagePath)

        crop_image = orig[h_exact-(cropsize/2):h_exact+1+(cropsize/2),w_exact-(cropsize/2):w_exact+1+(cropsize/2)]
        lowerleft = (-5 + cropsize/2,-5 + cropsize/2)

        # Drawing figure
        fig = plt.figure(figsize=(3,3),dpi=300)
        ax = fig.add_subplot(111)
        cax = ax.imshow(crop_image,cmap='gist_gray',interpolation='none',vmin=v_min,vmax=v_max)
        #cax = ax.imshow(crop_image,cmap='gist_gray',interpolation='none')
        ax.add_patch(Rectangle(lowerleft,10,10,fill=False,color='#FDBE46',lw=2))
        ax.axis('off')
        ax.set_title(title,fontsize=8)
        cbar = fig.colorbar(cax)
        
        # Saving figure
        f = os.path.join(destDir,title)
        plt.savefig(f+'.svg')
        plt.savefig(f+'.png')
        print "Saving file ...",f 
        plt.close()


def main(ch,fld,h,w,fpattern,offset_fpath,sourceDir,destDir,intensityScaling,cropsize):
    p = ExtractPeakImages(ch,fld,(h,w),fpattern,offset_fpath,sourceDir,destDir)
    cycle_imagePath_dict = p.get_allImagePath()
    cycle_offset_dict = p.get_allOffsets()
    
    for cycle, imagePath in cycle_imagePath_dict.items():
        title = 'cycle'+str(cycle)+'_ch'+str(ch)+'_fld'+str(fld)+'_hw'+str((h,w))
        imagePath, offset = cycle_imagePath_dict[cycle], cycle_offset_dict[cycle]
        p.getImageBox(imagePath, offset, title, cycle,intensityScaling,cropsize)
        




parser = argparse.ArgumentParser(description="""This script extracts the images surrounding the peaks. \n Saves the images in png and svg format \n
Example : /project/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis/extractPeakImages.py --file_pattern "2016_08_20_JSP045025_0*" --offset_file out_a02_hat_5_4/offsets_dict_ocbbvj.pkl \n --output_directory tempOutput \n --intensity_scaling 300 1000 --crop_size 24 -c 1 -fld 73 -hw 231 242 
                                 """)

parser.add_argument('-c','--channel',action="store",dest="ch",default=2,type=int,\
                    help="Channel of choice; Default=2")
parser.add_argument('-fld','--field',action="store",dest="fld",default=50,type=int,\
                    help="field; Note - the number may not correspond to the \
                    actual filename's xy")
parser.add_argument('-hw','--position',action="store",nargs=2,dest="hw",
                    help="peak position (h,w) as determined after the analysis; example: --position 400 192")
parser.add_argument('--file_pattern',nargs='+',action="store",dest="fnamePattern",type=str,
                    help="File name pattern (the directory structure where every cycle images are stored) ")
parser.add_argument('--offset_file',action="store",dest="offset_file",type=str)
parser.add_argument('--header_directory',action="store",dest="sourceDir",default=os.getcwd(),\
                    help="This is the source directory to start searching the pattern from ")
parser.add_argument('--output_directory',action="store",dest="destDir",default="extractedPeakImages",type=str,
                    help="Destination directory where the images are saved. ")
parser.add_argument('--intensity_scaling',action="store",dest="intensityScaling",nargs=2,default=[200,1000],
                    help="Intensity range to be used for scaling the final image")
parser.add_argument('--crop_size',action="store",dest="cropsize",default=12,type=int,
                    help="Size (pixels) to be cropped around the peak")

args = parser.parse_args()

try:
    assert len(args.fnamePattern) == 1
except AssertionError:
    raise SystemExit("****** Assertion Error *******\n. Please enter file pattern within quotes")

ch = args.ch
fld = args.fld
h,w = args.hw
[fnamePattern] = args.fnamePattern
sourceDir = os.path.abspath(args.sourceDir)
offset_fpath = os.path.abspath(args.offset_file)
destDir = os.path.abspath(args.destDir)
peakInfo = 'ch-'+str(ch)+'_fld-'+str(fld)+'_hw-'+str(h)+'_'+str(w)
destDir = os.path.join(destDir,peakInfo)
makeDir(destDir)
intensityScaling = map(int,args.intensityScaling)
cropsize = args.cropsize


main(ch,fld,h,w,fnamePattern,offset_fpath,sourceDir,destDir,intensityScaling,cropsize)



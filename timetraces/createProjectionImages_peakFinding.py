#! /home/jaggu/anaconda

"""
This is a file to cleanup the photobleaching timetrace files. The idea is to move through the sourceDir and retrieve the first
(and perhaps later maximum projection etc) image. Then I will create a new directory and copy the first tif file of the tflds.
"""

import os
import sys
import fnmatch
import shutil 


def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and
    below supplied root directory.'''
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            allFiles.append(os.path.join(path,filename))
    return allFiles

class ProjectedImage:
    """
    Creates a projected image given the traceDir
    """
    def __init__(self,traceDir,pattern='trace'):
        self.traceDir = traceDir
        self.pattern = '*'+pattern+'*.tif'

    def firstFrame(self):
        allTiff_files = locate(self.pattern, self.traceDir)
        firstFrame_file = sorted(allTiff_files)[0]
        return firstFrame_file




def firstFrame(traceDir):
    allTiff_files = os.path.listdir(traceDir)
    print allTiff_files


def createImages(sourceDir,proj_type='first',pattern='trace'):
    sourceDir = '/home/jaggu/marcotte_project/boulgakov/microscope/2014-Dec/2014-12-22t1'
    peakSourceDir = os.path.join(sourceDir,'projectedImages',proj_type)
    if not os.path.exists(peakSourceDir): os.makedirs(peakSourceDir)

    allSubDirPath = [os.path.join(sourceDir, subDir) for subDir in os.walk(sourceDir).next()[1] if pattern in subDir]
    for subDirPath in allSubDirPath:
        p = ProjectedImage(subDirPath,pattern)
        if proj_type is 'first':
            proj_file = p.firstFrame()
        shutil.copy(proj_file, peakSourceDir)
            
    


#createImages()
    


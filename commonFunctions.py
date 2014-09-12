#! /home/jaggu/anaconda/bin/python2.7

"""
This file contains common custom functions that were made and used repeatedly
"""

import os
import fnmatch 



def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and
    below supplied root directory.'''
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            allFiles.append(os.path.join(path,filename))
    return allFiles

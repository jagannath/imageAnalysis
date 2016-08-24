#! /home/jaggu/anaconda/python2.7

"""
I want to be able to convert and/or stitch images for the large number of images
gathered by tiling across the microscope slide. 
"""

from subprocess import call
import sys
import os
import collections
import fnmatch
from operator import itemgetter

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below supplied root directory.'''
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            allFiles.append(os.path.join(path, filename))
    return allFiles

class Images(object):
    """
    Converting the images in the directory to appropriate format for stitching
    it using fiji (for now).
    """
    def __init__(self,DIR):
        self.DIR = DIR
        self.fijiDIR = os.path.dirname(DIR) + '/' + 'forFIJI'
        self.checkDIR(self.DIR,self.fijiDIR)
        self.posFile = self.fijiDIR + '/' + 'TileConfiguration.txt'
        self.allPNGS = locate('*.png',self.DIR)
        self.allTIFS = locate('*.tif',self.DIR)
    def checkDIR(self,inputDIR,outputDIR):
        if not os.path.exists(inputDIR):raise SystemExit("Input Director doesn't exist")
        if not os.path.exists(outputDIR):os.makedirs(outputDIR)

    def path(self,(XTILES,YTILES)):
        xylist = [(x,y) for x in range(XTILES) for y in range(YTILES)]
        path = list()
        for j in range(YTILES):
            row = [xy for xy in xylist if xy[1]==j]
            if j%2 == 1: row = sorted(row,key=itemgetter(0),reverse=True)
            path.extend(row)
        return path

    def getPositions(self,fname):
        def getXY(fname):
            # This is temporary and complicated because the delimiter was '-'. The negative signs gets
            # split as well. So this elaborate logic
            xy = fname.split('.')[0].split('_')[1]
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
   
        fname = fname.split('/')[-1]
        (ii,jj) = map(lambda x:x[2:],fname.split('.')[1].split('-'))
        (xx,yy) = getXY(fname)
        return [(ii,jj),(xx,yy)]
                    
    def snakeByRows(self,order=('left','up')):        
        """For now left up is the default. may need to rethink for other options
            (50,20) ... <-
            |->....       |
           (50,1)...  <--(1,1)Pos=D
        This is snaking the way up from the bottommost position D; 
        """
        DIM = ( 50,20)
        path = self.path(DIM)#Interesting. Always from (0,0) -- (50,20) and snakes
        indexDICT = dict()
        for fname in self.allPNGS:
            [(ii,jj),(xx,yy)] = self.getPositions(fname)
            idx =  path.index((int(ii),int(jj)))
            newfname = self.fijiDIR + '/tile_%s.tiiiif'%(idx)
            cmd = ['convert',fname,newfname]
            call(cmd,shell=False)
                                   

def main():
    dir = '/project/marcotte/jagannath/projectfiles/SingleMoleculeMicroscopy/dataAnalysis/2013-July/2013-07-04/B2111_20pMCy5Dye/png' 
    im = Images(dir)
    im.snakeByRows()


def temp():
    DIR = '/project/marcotte/jagannath/projectfiles/SingleMoleculeMicroscopy/dataAnalysis/2013-July/2013-07-04/B2111_20pMCy5Dye/png'
    indexDIR = \
    '/project/marcotte/jagannath/projectfiles/SingleMoleculeMicroscopy/dataAnalysis/2013-July/2013-07-04/B2111_20pMCy5Dye/newfiji'
    allPNGS = locate('*.png',DIR)
    posfile = indexDIR + '/TileConfiguration.txt'
    with open(posfile,'w') as f:
        f.write('# Define the number of dimensions we are working on \n dim = 2 \n')
        f.write('# Define the image coordinates')
        counter = 0
        for fname in allPNGS:
            counter+=1
            i,j = fname.split('.')[-2].split('-')
            ii,jj = map(lambda x:x[-2:],[i,j])
            
            xy = fname.split('/')[-1].split('.')[0].split('_')[-1]
            xysplit = xy.split('-')
            if len(xysplit)>2:
                if 'x=' in xysplit: 
                    xidx = xysplit.index('x=')
                    xx='-'+xysplit[xidx+1] 
                if 'y=' in xysplit: 
                    yidx = xysplit.index('y=')
                    yy = '-'+xysplit[yidx+1]
            else:
                xx,yy = map(lambda foo:foo.split('=')[-1],xysplit)

            f.write("image%s.tif;;(%s,%s)\n"%(xx,yy))
#            f.write('_%s_%s;;(%s,%s,0.0)\n'%(str(counter),ii,jj))
            newfname = indexDIR + '/' + 'image%s.tif'%(counter)
            cmd = ['cp',fname,newfname]
            call(cmd,shell=False)
            print newfname


main()        
                





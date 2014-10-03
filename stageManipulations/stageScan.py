'''#! /home/jaggu/anaconda/bin/python2.7'''

# AUTHOR : Jagannath S; 
# Date Created : 2013-07-03
# Date Last Modified : 2013-07-03

"""
The script is designed to generate a list of xy coordinates when given the left
uppermost xy coordinate of the microscope slide.Essentially it is to raster scan
the entire slide. Not only will it provide a heat map but can also tell how much
of the peptide density (Slide is
A(4005.0,1355.0),B(-285.0,3355.0),C(-5005.0,1355.0),D(-5005.0,-2645.0),E(-285.0,-3645.0),F(4005,-2645))
The X,Y System however is completely reversed. A(+,+);B(-+),C(-,+) and so on.
The whole idea would be to flip over so that C can be the reference (0,0) and A
would be (+x,0) and so on. The diagram below is the real nature. 

WALL
___________________________
//////////////////////////

       B
    A--|--C
    |     |
    F--|--D
       E

    0
   -|- (ME)
   / \
(0,+y)  (+x,+y)
C-------A ^
|       | |+y
D-------F
(0,0)   (+x,0)
->+X
"""


check = 10

import numpy as np
import os.path
import os
import fnmatch
import csv
import collections
import time
import re
import shutil
import sys
from random import sample, choice
from subprocess import call
from operator import itemgetter



class StageScan(object):
    '''
    This class generates a grid given the position,origin,%overlap and the
    dimensions of the grid. 
    '''
    def __init__(self):
        self.ORIGIN = 10
        sys.exit(1) 
        try:
            self.ORIGIN = kwargs['ORIGIN']
            self.POS = kwargs['POS']
            self.OVERLAP = kwargs['OVERLAP']
            self.DIM = kwargs['DIM']
        except KeyError:
            raise SystemExit("Please define all: ORIGIN=(0,0),POS='D',OVERLAP=10,DIM=(200,20)")
        self.fovSize = 200 #Distance as 200 nm
        self.signX,self.signY = self.getSign()
        self.stagePath = list()
        self.indexPathDICT = dict()

    def getSign(self): return {'D':(+1,+1),'F':(-1,+1),'A':(-1,-1),'C':(+1,-1)}[self.POS]
    def stepsize(self): return self.fovSize - ((self.fovSize*self.OVERLAP)/100)
    def getPath(self,(XTILES,YTILES)):
        xylist = [(x,y) for x in range(XTILES) for y in range(YTILES)]
        path = list()
        for j in range(YTILES):
            row = [xy for xy in xylist if xy[1]==j]
            if j%2 == 1: row = sorted(row,key=itemgetter(0),reverse=True)
            path.extend(row)
        return path
    def getStagePath(self):
        '''
        By default, I provide the x and y position of the right lowermost corner-D.
        (see figure above in docstring for reference) and the percentage overlap
        needed. In later iteration maybe I will work out to include smaller regions
        and defined edges. In essence, the grid is a 45x20tiles on the most
        conservative estimate to ensure the PFS doesn't run off.
        (xi,yj) = (x0+i*dx,y0+j*dy) [i (0,45) and j (0,20)]
        If other positions, then just change the sign for dx and dy
        appropriately. 
        '''
        (x0,y0) = self.ORIGIN
        dx = dy = self.stepsize()
        print dx,dy
        path = self.getPath(self.DIM)
        self.stagePath = list()
        for xy in path:
            xVal = x0 + self.signX*xy[0]*dx
            yVal = y0 + self.signY*xy[1]*dy
            self.stagePath.append((xVal,yVal))
            self.indexPathDICT[(xVal,yVal)]=(xy)
        return self.stagePath
    def getVertices(self):
        '''
        D = (0,0); F = (200-1,0); A = (200-1,20-1); C = (0,20-1) [The first tile
        takes up a position.
        returns dict(A:aval,C:cval,D:dval,F:fval) # To confuse order
        '''
        xlim,ylim = self.DIM
        path = self.stagePath
        dpos = path[0]
        fpos = path[xlim-1]
        if ylim%2 == 0:
            cpos = path[xlim*(ylim-1)]
            apos = path[(xlim*ylim)-1]
        else:
            cpos = path[(xlim*ylim)-1]
            apos = path[xlim*(ylim-1)]
        return {'A':apos,'C':cpos,'D':dpos,'F':fpos}

    def getRandomPos(nbr=10): return sample(self.stagePath,nbr)




#def main():
#    ORIGIN = (-5005,-2645)
#    POS = 'D'
#    OVERLAP = 0
#    DIM = (200,20)
#
#    s = StageScan(ORIGIN=ORIGIN,POS=POS,OVERLAP=OVERLAP,DIM=DIM)
#    print len(s.getStagePath())
#    print s.getVertices()
#    print s.indexPathDICT


#if __name__ == '__main__': 
#    
#    t0 = time.clock()
#    main()
#    t1 = time.clock()
#    print "Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))

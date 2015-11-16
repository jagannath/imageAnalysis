"""
I have 8000; However after 4000 images, it got out of focus. It is 1000 (cols) x
40 (rows)
There is 40% overlap between two images. I am going to use FIJI's image
stitching software to join 
"""
import os
import re
from subprocess import call

sourceImgDir = '/project/boulgakov/microscope/2015-Oct/8000Images/151027_P004-005_20nM_Both_Gold_flds'
destDir = '/project/current/project2/jaggu/dataAnalysis/microscope1/2015-Oct/2015-10-27/stitchImages/4000Images'
destDir2 = '/project/current/project2/jaggu/dataAnalysis/microscope1/2015-Oct/2015-10-27/stitchImages/4000Images/scaledDown_25'

fnames = os.walk(sourceImgDir).next()[2]
for f in sorted(fnames):
    ch = f[-5]
    [nbr] = re.findall("xy\d+",f)
    fnbr = int(nbr.split('xy')[1])
     
    inputTif = os.path.join(sourceImgDir,f)
    if fnbr<=4001:
        
        outputTif1 = os.path.join(destDir,'ch'+str(ch),'ch'+str(ch)+'_tile_'+str(nbr)+'.tif')
        outputTif2 = os.path.join(destDir2,'ch'+str(ch),'scaled_ch'+str(ch)+'_tile_'+str(nbr)+'.tif')
        print "Renaming and shrinking to ..",outputTif1
        cmd1 = "convert " + inputTif + " " + outputTif1
        cmd2 = "convert -scale 25% " + inputTif + " " + outputTif2
        call(cmd1.split(),shell=False)
        call(cmd2.split(),shell=False)

print "Renaming Completed"

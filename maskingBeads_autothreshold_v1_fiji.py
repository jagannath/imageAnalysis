#! /python/anaconda

import os
import fnmatch
import sys, traceback
from ij import IJ, ImagePlus
from ij.process import ImageStatistics as IS 
from ij.plugin import ImageCalculator, Thresholder, Commands
from fiji.threshold import Auto_Threshold
from ij.plugin.filter import GaussianBlur, RankFilters, ParticleAnalyzer, Filters
from ij.gui import GenericDialog


class BeadStack(object):
	def __init__(self, fname, DIC_index=1):
		self.fname = fname # Note this has to be a TIFF stack
		self.DIC_index = DIC_index # Need to know which is the DIC_index
		self.img, self.stack, self.title = self._getInfo(fname)
		
	def checkStack(self):
		if self.img.getNSlices()>1: return True
		else: return None
			
	def _getInfo(self,fname):
		self.sourceDir = os.path.dirname(fname)
		img = IJ.openImage(fname)
		title = img.getTitle()
		stack = img.getImageStack()
		return img, stack, title

	def _dialog(self,comment):
		gd = GenericDialog(comment)
		gd.enableYesNoCancel()
		gd.showDialog()
		if gd.wasCanceled():
			print "user cancelled dialog"
			sys.exit(1)
			
		option = gd.wasOKed()
		return option

	def makeMask(self):
		"""
		This function makes the mask. The steps are (1) Minimum Filter - makes a darker boundary around beads (2) Autothresholding using the Huang algorithm - has some fuzzy logic and seems to work (3) Analyze particles with a size between 500-50000 and 
		circularity between 0.4 to 1.0; The mask generated is sometimes black on beads and white around. Then I need to invert the LUTs
		"""
	
		ipOriginal = self.stack.getProcessor(self.DIC_index)
		ip = ipOriginal.duplicate()
		imgUpdate = ImagePlus("New",ip)
		imgUpdate.setProcessor("Mask",ip)
		
		img0 = ImagePlus("Before",ipOriginal)
		img0.show()
		
		# Minimum filter
		RankFilters().rank(ip,2,RankFilters.MIN)
		img1 = ImagePlus("Filter",ip)
		# Autothreshold - Huang algorithm
		hist = ip.getHistogram()
		lowTH = Auto_Threshold.Huang(hist)
		ip.setThreshold(0,lowTH, ImageProcessor.BLACK_AND_WHITE_LUT)
		img3 = ImagePlus("Thresholded",ip)
		img3.show()

		# Making a binary mask
		IJ.run(img3,"Convert to Mask","")
		
		if self._dialog("Invert Mask ??") is True: IJ.run("Invert LUT")
		img3.updateAndDraw()
		
		# The true mask after Particle Analysis; Creates a mask image around the particles
		IJ.run(img3,"Analyze Particles...", "size=500-50000 circularity=0.40-1.00 show=Masks")
		img1.close()
		#img3.close()

		# Editing the masks (filling holes and dilating it twice)
		imgActive = IJ.getImage()
		IJ.run(imgActive,"Convert to Mask","")
		IJ.run(imgActive,"Fill Holes","")
		for i in range(8): IJ.run(imgActive,"Dilate","")
		ipActive = imgActive.getProcessor().convertToFloat()
		
		# Saving the mask
		maskFname = self.sourceDir + "\\" + self.title + '_mask'
		IJ.saveAs(imgActive, "PNG", maskFname)
				
		# Confirming that the image is masked and the histogram is correct
		#IJ.run(imgActive, "Histogram", "")
		
		#stats = ipActive.getStatistics()
		pixels = ipActive.getPixels()
		self.maskPixels = [pix/255 for pix in pixels]
		self.areaMask = self.maskPixels.count(1)

		# Checking if the image is fine. If not, returns option to skip
		ImageCalculator().calculate("zero create", img0, imgActive)

		skip = False
		if self._dialog("Skip Image ??") is True: skip = True

		IJ.run("Close All")
		return self.maskPixels, skip
		
	def _getStatsMask(self,pixels):
		pixUnderMask = list()
		for i in xrange(0,len(pixels)):
			pixUnderMask.append(pixels[i]*self.maskPixels[i])
		return sum(pixUnderMask)/self.areaMask
	
	def maskOtherChannels(self):
		outList = list()
		for i in xrange(1, self.img.getNSlices()+1): #It counts the area under DIC as well
			cp = self.stack.getProcessor(i).convertToFloat()
			pixels = cp.getPixels()
			intPerArea = self._getStatsMask(pixels)
			outList.append((i,intPerArea))
		return outList


## IMPORTANT: The PATH FROM LINUX WILL NOT WORK. '\' instead of '/'. Also confirm the escape characters. 
#fname = sourceDir + '\\' + 'NegControl_NoDye_TentagelBeads_Before_multichannel_1fps001.tif'

def locate(pattern, root=os.curdir):
	""" Locate all files matching supplied filename pattern in and below supplied root directory. """
	allFiles = list()
	
	for path, dirs, files in os.walk(root):
		for filename in fnmatch.filter(files, pattern):
			allFiles.append(os.path.join(path, filename))
	return allFiles

def main():
	pattern = "*.tif"
#	sourceDir = "Z:\\jagannath\projectfiles\SingleMoleculeMicroscopy\dataAnalysis\\2014-May\\2014-05-14\\rawImages\Water_Hilyte647Mal_TentagelBeadSH_24h_multichannel_1fps018"
#	sourceDir = "Z:\\jagannath\projectfiles\SingleMoleculeMicroscopy\dataAnalysis\\2014-May\\2014-05-11"
	#sourceDir = "Z:\\jagannath\projectfiles\EpiMicroscopy\\2014-June\\2014-06-27\\rawImages\TentagelNH2_BODIPYNHS_100ms_Before_flds107"
	#sourceDir = "Z:\\jagannath\projectfiles\EpiMicroscopy\\2014-June\\2014-06-27"
	#fname = sourceDir + "\\" + "TentagelNH2_Alexa405NHS_1fps_TFA0030h_flds007.tif0_noScale.png"
	#sourceDir = "Z:\\jagannath\projectfiles\EpiMicroscopy\\2014-June\\2014-06-24\\rawImages\TentagelNH2_Alexa405NHS_1fps_TFA0030h_flds007"
	dateStamp = sourceDir.split('\\')[-1]
	
	ofname = sourceDir + '\\' + dateStamp+'_beadsMaskedIntensity.csv'
	ofile = open(ofname, 'w')
	ofile.write('\t'.join(['ImagePath','ImageName','DIC_channel','AreaUnderMask (in Pixels)','Channel Number','Intensity/Pixel','Channel Number','Intensity/Pixel','Channel Number','Intensity/Pixel','Channel Number','Intensity/Pixel','Channel Number','Intensity/Pixel','\n']))
	allTiffs = locate(pattern,sourceDir)
	tiffStackFiles = sorted(allTiffs)
	# I am making a special loop structure here. It actually redoes the option of inverting mask, so that I dont waste the image. Also it skips to next image if I skipped thrice 
	restart = True
	i = 0
	nbrSkips = 0
	while restart and i<len(tiffStackFiles):
		fname = tiffStackFiles[i]
		i+=1
		tiffStack = BeadStack(fname,DIC_index=1)
		if tiffStack.checkStack(): #Checks if the tif 
			maskPixels,skip = tiffStack.makeMask()
			print len(maskPixels), maskPixels.count(1), skip
			if skip is True: 
				if nbrSkips == 2:
					nbrSkips = 0
					continue
				else:
					nbrSkips +=1
					i-=1
					continue;
			outList = tiffStack.maskOtherChannels()
			outLine = fname + '\t' + fname.split('\\')[-1] + '\t' + str(tiffStack.DIC_index) + '\t' + str(tiffStack.areaMask) + '\t' + '\t'.join(str(outList[i][j]) for i in xrange(0,len(outList)) for j in [0,1]) + '\n'
			ofile.write(outLine)
			print outLine
	ofile.close()
	print "Done"	

def test(): 
	imp = IJ.openImage(fname)
	ip = imp.getProcessor()
	RankFilters().rank(ip,2,RankFilters.MIN)

	IJ.run("Analyze Particles...", "size=500-50000 circularity=0.40-1.00 show=Masks");
	imgActive = IJ.getImage()
	ipActive = imgActive.getProcessor().convertToFloat()

	imgActive.show()

main()

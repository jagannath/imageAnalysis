#! /python/anaconda

import os
import fnmatch
import sys, traceback
from ij import IJ, ImagePlus, process
from ij.process import ImageStatistics as IS 
from ij.plugin import ImageCalculator, Thresholder, Commands
from fiji.threshold import Auto_Threshold
from ij.plugin.filter import GaussianBlur, RankFilters, ParticleAnalyzer, Filters
from ij.gui import GenericDialog

sourceDir = "W:\\boulgakov\microscope2\jagannath\\rawFiles\\2014-07\\2014-07-15\\AS_patterned_2uMAlexa555SE_TRITC_50ms_16hWash_flds2"
fname = sourceDir + "\\"+ "AS_patterned_2uMAlexa555SE_TRITC_50ms_16hWash_flds2xy02.tif"
print fname

class PatternDirectory():
	"""
	This is the directory containing all the images for the Photolithographed patterned fields. It should not have "Stitched" or "TimeTraced Images"
	"""
	def __init__(self, directory):
		self.directory = directory

	def getImages(self):
		imgFiles = os.listdir(self.directory)
		return imgFiles
	


def test():
	newImg = ImagePlus("GrayScaled",imp)
	newip = newImg.getProcessor()

	hist = newip.getHistogram()
	lowTH = Auto_Threshold.IsoData(hist)
	newip.setThreshold(lowTH, max(hist),ImageProcessor.BLACK_AND_WHITE_LUT)


	rt = ResultsTable()
	pa = ParticleAnalyzer(ParticleAnalyzer.SHOW_RESULTS | ParticleAnalyzer.SHOW_OVERLAY_OUTLINES, Measurements.AREA |Measurements.MEAN |\
		Measurements.MEDIAN | Measurements.STD_DEV | Measurements.MIN_MAX | Measurements.RECT, rt,50, 200000, 0.5, 1  )
	pa.setResultsTable(rt)
	pa.analyze(newImg)
	rt.show("Results")

sourceDir = "W:\\boulgakov\microscope2\jagannath\\rawFiles\\2014-07\\2014-07-15\\AS_patterned_2uMAlexa555SE_TRITC_50ms_16hWash_flds2"
patDir = PatternDirectory(sourceDir)
print patDir.getImages()
#img = IJ.openImage(fname)
#imp = img.getProcessor(1)






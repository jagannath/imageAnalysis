#!/usr/bin/python

# AUTHOR : Alexander Boulgakov
# Modified By: Jagannath S; 2nd May 2013

"""
The aim of this script is to use Alexander's algorithm to walk through every Tiff file in the directory - including subdirectory. As each file is being processed, a png file and a png file with the peaks identified is output in the same directory. 
An output file (csv file) is generated DIR.csv. This will have subdirectory as the row and the columns will be counts of the files in the subdirectory. The output function will be made modular to encompass other ways to output. Right now (all files within a directory must be averaged).  
"""

import sys
import scipy.ndimage.filters as spf
import scipy.misc
import scipy.signal
import numpy as np
import Image
import ImageOps
import ImageDraw
import datetime
import math
import subprocess
import os.path

def simplecount(c_std, imageFile):
    image_list = [imageFile]
    c_std = float(c_std)
    if image_list[0][-3:] == 'tif':
	for f, filename in enumerate(image_list):
	    subprocess.call(["convert", filename, filename + '.png'])
	    image_list[f] = filename + '.png'
    output_filename_hash = (datetime.datetime.now().minute * 60 + datetime.datetime.now().second)
    images = [scipy.misc.imread(filename).astype(np.uint16)
	    for filename in image_list]
    #correlation_matrix = np.array([[-1, -1, -1, -1, -1],
    #                               [-1, -1, -1, -1, -1],
    #                               [-1, -1, 24, -1, -1],
    #                               [-1, -1, -1, -1, -1],
    #                               [-1, -1, -1, -1, -1]])
    correlation_matrix = np.array(
					    [[-5935, -5935, -5935, -5935, -5935],
					    [-5935,  8027,  8027,  8027, -5935],
					    [-5935,  8027, 30742,  8027, -5935],
					    [-5935,  8027,  8027,  8027, -5935],
					    [-5935, -5935, -5935, -5935, -5935]])
    median_diameter = 5
    local_max = 3
    #c_std = 1
    processed_images = [np.copy(image).astype(np.int64) for image in images]
    processed_images = \
	[np.subtract(image, np.minimum(spf.median_filter(image, median_diameter),
					image))
	for image in images]
    processed_images  = \
	[np.maximum(scipy.signal.correlate(image, correlation_matrix, mode='same'),
		    np.zeros_like(image)).astype(np.int64)
	for image in processed_images]
    thresholded_images = [image > np.mean(image) + c_std * np.std(image) and images[i] > a_std  # Two changes done to Alexanders script
			for i,image in enumerate(processed_images)]
    for i, mask in enumerate(thresholded_images):
	for (h, w), valid in np.ndenumerate(mask):
	    if valid:
		local_slice = \
		    np.copy(processed_images[i][h - local_max:h + local_max + 1,
						w - local_max:w + local_max + 1])
		if (h + local_max >= mask.shape[0] or
		    w + local_max >= mask.shape[1] or
		    h - local_max < 0 or w - local_max < 0):
		    mask[h, w] = False
		    continue
		local_slice[local_max, local_max] = 0
		if np.amax(local_slice) >= processed_images[i][h, w]:
		    mask[h, w] = False
    size = 4
    for i, image in enumerate(images):
	i_std = np.std(image)
	i_median = np.median(image) - i_std
	i_max = np.amax(image)
	for (h, w), value in np.ndenumerate(image):
	    if image[h, w] > i_median:
		image[h, w] = int(np.around(math.sqrt(image[h, w] - i_median) /
					    math.sqrt(i_max - i_median)
					    * float(2**8 - 1)))
	    else:
		image[h, w] = 0
    output_images = [ImageOps.colorize(
			Image.fromstring('L', image.shape, image.astype(np.uint8)),
					(0,0,0), (255,255,255))
			for image in images]
    for i, image in enumerate(images):
	for (h, w), value in np.ndenumerate(image):
	    if thresholded_images[i][h, w]:
		box = ((w - size, h - size),(w + size, h + size))
		draw = ImageDraw.Draw(output_images[i])
		draw.rectangle(box, fill=None, outline='blue')
    # Editing the output file name; Will have the same as the filename + png 
    #[image.save(str(i).zfill(int(np.ceil(math.log(len(output_images), 10)))) +'_' + imageFile[:-4] + str(output_filename_hash) + '.png')
	#for i, image in enumerate(output_images)]
    [image.save(imageFile[:-4]+str(output_filename_hash)+'_peaks.png') for image in output_images]
    for image in thresholded_images:
	pfCounts = np.sum(np.sum(image))
    
    return pfCounts

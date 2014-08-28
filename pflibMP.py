#! /home/boulgakov/localpython2/bin/python
#!/home/boulgakov/localpython/bin/python
#!/home/boulgakov/anaconda/bin/python

"""
Jag and Alex's automated pipeline. Thank you Zack for help!

pflibMP stands for "peak fitting library -- multiple processor".

This file contains subpixel image analysis routines cobbled together to analyse
images resulting from single molecule peptide sequencing on our inverted TIRF
assembly. The algorithms used are fairly general, however some of the default
parameter values are based on experience with our specific datasets.
"""

import sys
import numpy as np
import scipy.ndimage.filters as spf
import scipy.signal
import scipy.misc
import scipy.stats
import agpy
#import Image#ALEX:THIS IS FOR OLD PIL...PILLOW USES THE NEXT LINE INSTEAD
from PIL import Image
from PIL import ImageOps#SAME FOR THIS
from PIL import ImageDraw#AND THIS
import itertools
import math
import time
import datetime
import multiprocessing
import cProfile
import logging
import os
import subprocess

def fit_peaks(images, image_tuples=None, correlation_matrix=None,
              median_diameter=5, r_2_threshold=0.5, c_std=4,
              consolidation_radius_sq=4**2, child_number=0,
              child_priority_increment=0, silent=True):
    """
    The core function of the algorithm. In each image, finds all peaks and
    evaluates their resemblance to a Gaussian PSF model.
    """
    #renice
    os.nice(child_priority_increment)
    #generate empty image metadata tuples if empty
    if image_tuples is None:
        if not silent:
            logging.info('pflib.fit_peaks child number ' + str(child_number) +
                         ': no image_tuples metadata passed, using empty ' + 
                         'tuples')
        image_tuples = [() for t in images]
    elif len(image_tuples) != len(images):
        if not silent:
            logging.info('pflib.fit_peaks child number ' + str(child_number) +
            ': number of image metadata tuples does not match number of ' +
            'images, using empty tuples')
        image_tuples = [() for t in images]
    #perform squareroot just this once for the entire function call
    consolidation_box = np.ceil(math.sqrt(consolidation_radius_sq))
    #copy images to prevent tampering with originals by any function called
    images = [np.copy(image).astype(np.int64) for image in images]
    images_s3 = images
    if not silent:
        logging.info("pflib.fit_peaks child number " + str(child_number) +
                     ": copied images")
        sys.stdout.flush()
    #TODO:the following line does not check that images remain within bounds as
    #spalttiefpass and gaussian may increase positive values
    #subtract median filter
    images_s3 = \
      [np.subtract(image, np.minimum(spf.median_filter(image, median_diameter),
                                     image))
       for image in images_s3]
    if not silent:
        logging.info("pflib.fit_peaks child number " + str(child_number) +
                     ": applied median filter")
        sys.stdout.flush()
    #perform correlation
    if correlation_matrix is not None:
        images_s3 = [np.maximum(scipy.signal.correlate(image,
                                                       correlation_matrix,
                                                       mode='same'),
                                np.zeros_like(image)).astype(np.int64)
                     for image in images]
        if not silent:
            logging.info("pflib.fit_peaks child number " + str(child_number) +
                         ": applied correlation filter")
            sys.stdout.flush()
    consolidated_peak_stack = [[] for image in images]
    for i, image in enumerate(images_s3):
        start_time = time.clock()
        #number of fit attempts
        fit_counter = 0
        #no dictionary comprehensions in python 2.6 :(
        #pixel_bins = {(h,w): [] for h in range(image.shape[0])
        #                        for w in range(image.shape[1])}
        #pixel_bins = dict([((h,w), []) for h in range(image.shape[0])
        #                               for w in range(image.shape[1])])
        pixel_bins = [[[] for w in range(image.shape[1])]
                          for h in range(image.shape[0])]
        #threshold beneath which a pixel needs not to be checked
        candidate_threshold = np.mean(image) + c_std * np.std(image)
        for h, w in itertools.product(range(image.shape[0]),
                                      range(image.shape[1])):
            #keep track of progress within each image
            peak_time = time.clock()
            if (not 1 < h < image.shape[0] - 2 or
                not 1 < w < image.shape[1] - 2):
                continue
            if image[h, w] < candidate_threshold:
                continue
            fit_counter = fit_counter + 1
            subimage = images[i][h - 2:h + 3, w - 2:w + 3]
            #agpy.gaussfit parameters --
            #(height, amplitude, x, y, width_x, width_y, rota)
            ((fit, fit_err), fit_image) = \
                     agpy.gaussfit(subimage,
                                   params=(np.mean(subimage),
                                           np.amax(subimage),
                                           2.5, 2.5, 1, 1, 0),
                                   return_all=True,
                                   limitedmin = [True] * 7,
                                   limitedmax = [False, False,
                                                 True, True, True, True, True],
                                   minpars=np.array([0.00,
                                                     (np.amax(subimage) -
                                                      np.mean(subimage)) / 3.0,
                                                     2.00, 2.00, 0.75, 0.75,
                                                     0.00]),
                                   maxpars=np.array([0.00, 0.00, 3.00, 3.00,
                                                     2.00, 2.00, 360.00]),
                                   returnfitimage=True)
            #rmse
            rmse = math.sqrt(
                         sum([(subimage[x,y] - fit_image[x,y])**2
                              for x,y in itertools.product(range(5),range(5))])
                         / 9.0)
            #coefficient of determination r**2
            r_2 = (1.0 -
                        sum(np.reshape((subimage - fit_image)**2, -1)) /
                        sum((np.reshape(subimage, -1) - np.mean(subimage))**2))
            if r_2 < r_2_threshold:
                continue
            #Illumina S/N = 
            #([peak pixel] - [mean of 16 outer edge pixels]) /
            #    (std[16 outer pixels])
            op = [subimage[x, y]
                      for x in range(subimage.shape[0])
                          for y in range(subimage.shape[1])
                              if x == 0 or
                                 x == subimage.shape[0] - 1 or
                                 y == 0 or
                                 y == subimage.shape[1] - 1]
            s_n = (np.amax(subimage) - np.mean(op)) / np.std(op)
            #peak fitting result tuple
            #(0  1  2--  3------  4-------  5--------  6---  7--, 8--)
            #(h, w, fit, fit_err, subimage, fit_image, rmse, r_2, s_n)
            #agpy.gaussfit parameters --
            # 0-----  1--------  2  3  4------  5------  6---
            #(height, amplitude, x, y, width_x, width_y, rota)
            peak = (fit[2] + h - 2, fit[3] + w - 2, fit, fit_err, subimage,
                    fit_image, rmse, r_2, s_n)
            #for all fits, consolidate them into the best one within a radius;
            #check peaks in the consolidated peak stack
            #default is that this is a new peak
            if not (0 <= peak[0] <= image.shape[0] and
                    0 <= peak[1] <= image.shape[1]):
                continue
            (pixel_bins[int(np.around(peak[0]))]
                       [int(np.around(peak[1]))].append(peak))
            if (not silent
                and
                (h + w) % int(image.shape[0] * image.shape[1] / 1000) == 0):
                logging.info("pflib.fitpeaks child number " +
                             str(child_number) +
                             ": done fitting at h,w " + str((h,w)) + " in " +
                             str(time.clock() - peak_time) + " sec; layer " +
                             str(i) + ' of ' + str(len(images_s3)) +
                             "; cumulative time " +
                             str(time.clock() - start_time) +
                             "; fit_counter " + str(fit_counter))
                sys.stdout.flush()
        consolidation_time = time.clock()
        for h, w in itertools.product(range(image.shape[0]),
                                      range(image.shape[1])):
            if pixel_bins[h][w]:
                pixel_bins[h][w] = max(pixel_bins[h][w],
                                       key=lambda x: x[7])
            else:
                pixel_bins[h][w] = None
        for h, w in itertools.product(range(image.shape[0]),
                                      range(image.shape[1])):
            if not pixel_bins[h][w]:
                continue
            local_max = True
            for j, q in itertools.product(
                        range(int(max(h - consolidation_box - 2, 0)),
                              int(min(h + consolidation_box + 3,
                                  image.shape[0]))),
                        range(int(max(w - consolidation_box - 2, 0)),
                              int(min(w + consolidation_box + 3,
                                  image.shape[1])))):
                if not pixel_bins[j][q] or (j == h and q == w):
                    continue
                if ((pixel_bins[j][q][0] - pixel_bins[h][w][0])**2 +
                    (pixel_bins[j][q][1] - pixel_bins[h][w][1])**2 <
                    consolidation_radius_sq):
                    if pixel_bins[j][q][7] < pixel_bins[h][w][7]:
                        pixel_bins[j][q] = None
                    else:
                        local_max = False
            if not local_max:
                pixel_bins[h][w] = None
        consolidated_peak_stack[i] = \
                         (image_tuples[i],
                          [pixel_bins[h][w]
                           for h, w in itertools.product(range(image.shape[0]),
                                                         range(image.shape[1]))
                           if pixel_bins[h][w]])
        if not silent:
            logging.info("pflib.fit_peaks child number " + str(child_number) +
                         ": consolidation accomplished in " +
                         str(time.clock() - consolidation_time) + " sec; " +
                         str(len(consolidated_peak_stack[i][1])) +
                         " peaks remaining of " + str(fit_counter) +
                         " fit_counter")
        if not silent:
            logging.info("pflib.fit_peaks child number " + str(child_number) +
                         ": image " + str(i + 1) + " of " +
                         str(len(images_s3)) + " processed in " +
                         str(time.clock() - start_time) + " sec;" +
                         " fit_counter " + str(fit_counter))
    return consolidated_peak_stack
    
result_wrapper = []
def wrapper(child_images, child_tuples, correlation_matrix, r_2_threshold,
            c_std, c, child_priority_increment, silent):
    result_wrapper.append(fit_peaks(child_images, child_tuples,
                                    correlation_matrix=correlation_matrix,
                                    r_2_threshold=r_2_threshold, c_std=c_std,
                                    child_number=c,
                                    child_priority_increment=child_priority_increment,
                                    silent=silent))

def run_fit_peaks(image_tuples, crop=0,
                  correlation_matrix=np.array(
                                        [[-5935, -5935, -5935, -5935, -5935],
                                         [-5935,  8027,  8027,  8027, -5935],
                                         [-5935,  8027, 30742,  8027, -5935],
                                         [-5935,  8027,  8027,  8027, -5935],
                                         [-5935, -5935, -5935, -5935, -5935]]),
                  r_2_threshold=0.5, c_std=4,
                  output_filename_hash=datetime.datetime.now().minute,
                  drift=(0.0, 0.0), max_processors=0,
                  child_priority_increment=0, silent=False, profile=False):
    if not silent:
        logging.info("\npflib.run_fit_peaks: staring datetime " +
              str(datetime.datetime.now()) + "\n")
        logging.info("\npflib.run_fit_peaks: fitting peaks for images in\n")
        logging.info(image_tuples)
    #def image_multiplexer#####################################################
    def image_multiplexer(child_tuples, c, q):
        child_images = [
                    scipy.misc.imread(metadata[3]) if crop is 0
                    else scipy.misc.imread(metadata[3])[crop:-crop, crop:-crop]
                        for metadata in child_tuples]
        if not silent:
            logging.info("pflib.run_fit_peaks child number " + str(c) +
                  ": loaded " + str(len(child_images)) + " images")
        if profile:
            profile_file = ("profile run " + str(output_filename_hash) +
                            " child " + str(c) + ".p")
            cProfile.runctx('wrapper(child_images, child_tuples, ' +
                            'correlation_matrix, r_2_threshold, c_std, c, ' +
                            'child_priority_increment, silent)',
                            globals(),
                            {'child_images': child_images,
                             'child_tuples': child_tuples,
                             'correlation_matrix': correlation_matrix,
                             'r_2_threshold': r_2_threshold,
                             'c_std': c_std,
                             'c': c,
                             'child_priority_increment': child_priority_increment,
                             'silent': silent},
                            profile_file)
            q.put(result_wrapper[0])
        else:
            print('pflibMP.run_fit_peaks: nonprofiling branch not ' + 
                  'maintained; exiting')
            sys.exit()
    #end image_multiplexer#####################################################
    child_count = (min(max_processors, multiprocessing.cpu_count(),
                       len(image_tuples))
                   if max_processors > 0
                   else min(multiprocessing.cpu_count(), len(image_tuples)))
    result_queue = multiprocessing.Queue()
    child_process_list = []
    child_image_tuples = [[] for x in range(child_count)]
    cil_cycle = itertools.cycle(child_image_tuples)
    for i, image_tuple in enumerate(image_tuples):
        cil_cycle.next().append(image_tuple)
    if not silent:
        logging.info("\nnumber of images: " + str(len(image_tuples)) +
              "; cpu count: " + str(multiprocessing.cpu_count()) +
              "; max processors: " + str(max_processors) +
              "; spawning children with " +
              str([len(x) for x in child_image_tuples]) +
              " images each")
    for c, child_tuples in enumerate(child_image_tuples):
        p = multiprocessing.Process(target=image_multiplexer,
                                    args=(child_tuples, c, result_queue))
        p.start()
        child_process_list.append(p)
    #stored as [((image metadata tuple), layer)]
    consolidated_peak_stack = []
    for process in child_process_list:
        consolidated_peak_stack += result_queue.get()
    for process in child_process_list:
        process.join()
    #sort first by channel x[0][4], then by index number x[0][2]
    consolidated_peak_stack = sorted(sorted(consolidated_peak_stack,
                                            key=lambda x: x[0][4]),
                                     key=lambda x: x[0][2])
    if not silent:
        logging.info("\npflib.run_fit_peaks: end datetime " +
              str(datetime.datetime.now()) + "\n")
    return consolidated_peak_stack

def linked_analysis(consolidated_peak_stack, wiggle=3, silent=False):
    #(0  1  2--  3------  4-------  5--------  6---  7--)
    #(h, w, fit, fit_err, subimage, fit_image, rmse, r_2)
    #agpy.gaussfit parameters --
    # 0-----  1--------  2  3  4------  5------  6---
    #(height, amplitude, x, y, width_x, width_y, rota)
    #flash vs persistent;mean, std
    #height, amplitude, width_x, width_y, volume, rota, rmse, r_2
    #{(bin name, variable name): [ordered data])}
    #TODO:drift correction
    phantom_hashes = {}
    linked_consolidated_peak_stack = [[] for i in consolidated_peak_stack]
    for i, (m, layer) in enumerate(consolidated_peak_stack):
        for p, peak in enumerate(layer):
            linked_consolidated_peak_stack[i].append((peak, [], [], i,
                                                      peak[0], peak[1]))
            phantom_hashes.setdefault((linked_consolidated_peak_stack[i][p][4],
                                       linked_consolidated_peak_stack[i][p][5],
                                       i),
                                      (i, p))
    for i, layer in enumerate(linked_consolidated_peak_stack[:-1]):
        start_time = time.time()
        next_layer = linked_consolidated_peak_stack[i + 1]
        for p, peak in enumerate(layer):
            for n_peak in next_layer:
                if ((peak[4] - n_peak[4])**2 + (peak[5] - n_peak[5])**2 <=
                    wiggle**2):
                    n_peak[1].append((i, peak[4], peak[5]))
                    peak[2].append((i + 1, n_peak[4], n_peak[5]))
        if not silent:
            logging.info("pflibMP.linked_analysis: connected peaks in layer " +
                  str(i) + " of " +
                  str(len(linked_consolidated_peak_stack[:-1]) - 1) + " in " +
                  str(time.time() - start_time) + " sec")
    phantom_dictionary = {}
    phantom_inverse_dictionary = {}
    phantom_int = 0
    for i, layer in enumerate(linked_consolidated_peak_stack):
        if i == 0:
            for peak in layer:
                phantom_dictionary.setdefault((peak[4], peak[5], 0), phantom_int)
                phantom_inverse_dictionary.setdefault(phantom_int, [(peak[4], peak[5], 0)])
                phantom_int = phantom_int + 1
            continue
        start_time = time.time()
        last_layer = linked_consolidated_peak_stack[i - 1]
        for p, peak in enumerate(layer):
            b = [False for x in last_layer]
            for L, L_peak in enumerate(last_layer):
                if ((peak[4] - L_peak[4])**2 + (peak[5] - L_peak[5])**2 <=
                    wiggle**2):
                    phantom_dictionary.setdefault((peak[4], peak[5], i), phantom_dictionary[(L_peak[4], L_peak[5], i - 1)])
                    phantom_inverse_dictionary[phantom_dictionary[(L_peak[4], L_peak[5], i - 1)]].append((peak[4], peak[5], i))
                    b[L] = True
            if not any(b):
                phantom_int = phantom_int + 1
                phantom_dictionary.setdefault((peak[4], peak[5], i), phantom_int)
                phantom_inverse_dictionary.setdefault(phantom_int, [(peak[4], peak[5], i)])
                continue
            if len([x for x in b if x]) == 1:
                continue
            consolidation_group = list(set([phantom_dictionary[(last_layer[x][4], last_layer[x][5], i - 1)] for x, tv in enumerate(b) if tv]))
            leader = consolidation_group[0]
            for g in consolidation_group[1:]:
                for c in phantom_inverse_dictionary[g]:
                    phantom_dictionary[c] = leader
                    phantom_inverse_dictionary[leader].append(c)
                del phantom_inverse_dictionary[g]
            phantom_inverse_dictionary[leader] = list(set(phantom_inverse_dictionary[leader]))
        if not silent:
            logging.info("pflibMP.linked_analysis: int connected peaks in layer " +
                  str(i - 1) + " of " +
                  str(len(linked_consolidated_peak_stack[:-1])) + " in " +
                  str(time.time() - start_time) + " sec")
    e = []
    for phantom, members in phantom_inverse_dictionary.iteritems():
        e.append(len(members))
    return linked_consolidated_peak_stack, phantom_dictionary, phantom_inverse_dictionary, phantom_hashes

def save_linked_images(images, phantom_inverse_dictionary, phantom_hashes,
                       linked_consolidated_peak_stack, drift=(0.0,0.0),
                       output_filename_hash=("LINKED run " +
                                          str(datetime.datetime.now().minute)),
                       image_names=None, size=4, brightness_shift=4,
                       silent=False):
    #try using imagemagick instead of brightness_shift
    #brightness_shift = 0
    #TODO:implement drift, this variable is currently not used
    #drift = [float(d) for d in drift]
    start_time = time.time()
    output_images = [ImageOps.colorize(
                      Image.fromstring('L', image.shape,
                                       (np.copy(image) >> brightness_shift).astype(np.uint8)),
                                       (0,0,0), (255,255,255))
                     for image in images]
    if not silent:
        logging.info("pflib.save_linked_images: images converted into PIL in " +
              str(time.time() - start_time) + " sec")
    for phantom, peaks in phantom_inverse_dictionary.iteritems():
        if len(peaks) == 1:
            peak_h = peaks[0][0]
            peak_w = peaks[0][1]
            peak_i = peaks[0][2]
            h = peak_h
            w = peak_w
            box = ((w - size + np.around(i / drift[1])
                    if drift[1] != 0 else w - size,
                    h - size + np.around(i / drift[0])
                    if drift[0] != 0 else h - size),
                   (w + size + np.around(i / drift[1])
                    if drift[1] != 0 else w + size,
                    h + size + np.around(i / drift[0])
                    if drift[0] != 0 else h + size))
            draw = ImageDraw.Draw(output_images[peak_i])
            draw.rectangle(box, fill=None, outline='lightblue')
            continue
        degenerate = False
        for i in range(len(linked_consolidated_peak_stack)):
            count_i = 0
            for peak in peaks:
                if peak[2] == i:
                    count_i = count_i + 1
            if count_i > 1:
                degenerate = True
                break
        if degenerate:
            for peak in peaks:
                peak_h = peak[0]
                peak_w = peak[1]
                peak_i = peak[2]
                h = peak_h
                w = peak_w
                box = ((w - size + np.around(i / drift[1])
                        if drift[1] != 0 else w - size,
                        h - size + np.around(i / drift[0])
                        if drift[0] != 0 else h - size),
                       (w + size + np.around(i / drift[1])
                        if drift[1] != 0 else w + size,
                        h + size + np.around(i / drift[0])
                        if drift[0] != 0 else h + size))
                draw = ImageDraw.Draw(output_images[peak_i])
                draw.rectangle(box, fill=None, outline='red')
        else:
            for peak in peaks:
                peak_h = peak[0]
                peak_w = peak[1]
                peak_i = peak[2]
                h = peak_h
                w = peak_w
                box = ((w - size + np.around(i / drift[1])
                        if drift[1] != 0 else w - size,
                        h - size + np.around(i / drift[0])
                        if drift[0] != 0 else h - size),
                       (w + size + np.around(i / drift[1])
                        if drift[1] != 0 else w + size,
                        h + size + np.around(i / drift[0])
                        if drift[0] != 0 else h + size))
                draw = ImageDraw.Draw(output_images[peak_i])
                draw.rectangle(box, fill=None, outline='blue')
    if not silent:
        logging.info("\npflib.save_linked_images: saving images using " +
              output_filename_hash + "\n")
    for i, img in enumerate(output_images):
        img.save(str(i).zfill(int(np.ceil(math.log(len(output_images), 10)))) +
                 ' ' + output_filename_hash + ' -- ' + image_names[i] +
                 ' --.png')

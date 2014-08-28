#!/home/boulgakov/localpython/bin/python

import numpy as np
import sys
import pflibMP
import datetime
import cPickle
import glob
import IPython
import logging
import subprocess
#import os.path
import os
import csv

output_filename_hash = datetime.datetime.now().minute * 60 + datetime.datetime.now().second

print("output_filename_hash " + str(output_filename_hash))

logging.basicConfig(filename='run log ' + str(output_filename_hash) + '.log', level=logging.INFO, format='%(asctime)s %(message)s')
logging.info("output_filename_hash " + str(output_filename_hash))

drift = (0.0, 0.0)
#drift=(50.0,136.0)
c_std = 2
r_2_threshold = 0.7
child_priority_increment = 1
correlation_matrix = np.array(
                                        [[-5935, -5935, -5935, -5935, -5935],
                                         [-5935,  8027,  8027,  8027, -5935],
                                         [-5935,  8027, 30742,  8027, -5935],
                                         [-5935,  8027,  8027,  8027, -5935],
                                         [-5935, -5935, -5935, -5935, -5935]])

logging.info("drift = " + str(drift))
logging.info("r_2_threshold " + str(r_2_threshold))
logging.info("correlation_matrix " + str(correlation_matrix))
brightness_shift = 6

image_list = sorted(sys.argv[1:])
if image_list[0][-3:] == 'tif':
    for f, filename in enumerate(image_list):
        subprocess.call(["convert", filename, filename + '.png'])
        image_list[f] = filename + '.png'

consolidated_peak_stack, peak_cube, peak_data_cube, images = pflibMP.run_fit_peaks(image_list, crop=0, correlation_matrix=correlation_matrix, r_2_threshold=r_2_threshold, c_std=c_std, output_filename_hash=output_filename_hash, drift=drift, child_priority_increment=child_priority_increment, silent=False, profile=True)

#arrays for making histograms
#height, amplitude, width_x, width_y, rmse, r_2, s_n
#histograms = tuple([[peak[2][0] for layer in consolidated_peak_stack for peak in layer],
#                    [peak[2][1] for layer in consolidated_peak_stack for peak in layer], 
#                    [peak[2][4] for layer in consolidated_peak_stack for peak in layer], 
#                    [peak[2][5] for layer in consolidated_peak_stack for peak in layer], 
#                    [peak[6] for layer in consolidated_peak_stack for peak in layer],
#                    [peak[7] for layer in consolidated_peak_stack for peak in layer],
#                    [peak[8] for layer in consolidated_peak_stack for peak in layer]])
#cPickle.dump(histograms, open('histograms ' + str(output_filename_hash) + '.pkl', 'w'))

#coallate all pixels in each image into either background or peaks
background_pixels = [[] for x in consolidated_peak_stack]
peak_pixels = [[] for x in consolidated_peak_stack]
#partition masks for background vs peak
partition_masks = [[] for x in consolidated_peak_stack]
for L, layer in enumerate(consolidated_peak_stack):
    for peak in layer:
        partition_masks[L].append((np.around(peak[0]), np.around(peak[1])))
    partition_masks[L] = set(partition_masks[L])
    for h in range(images[L].shape[0]):
        for w in range(images[L].shape[1]):
            if (h, w) in partition_masks[L]:
                peak_pixels[L].append(images[L][h, w])
            else:
                background_pixels[L].append(images[L][h, w])

ss = csv.writer(open('peak stats ALGORITHM v1 ' + str(output_filename_hash) + '.csv', 'wb'), dialect='excel-tab')
ss.writerow(['filename'] + ['number of peaks'] + ['mean(peak fit quality r^2)'] + ['std dev(peak fit quality r^2)'] + ['mean(peak illumina snr)'] + ['std dev(peak illumina snr)'] +
            ['mean(background pixels)'] + ['std dev(background pixels)'] + ['mean(1x1 peak pixel value)'] + ['std dev(1x1 peak pixel value)'])
for L, layer in enumerate(consolidated_peak_stack):
    ss.writerow([image_list[L]] + [str(len(layer))] + [str(np.mean([peak[7] for peak in layer]))] + [str(np.std([peak[7] for peak in layer]))] + [str(np.mean([peak[8] for peak in layer]))] + [str(np.std([peak[8] for peak in layer]))] + 
                [str(np.mean(background_pixels[L]))] + [str(np.std(background_pixels[L]))] + [str(np.mean(peak_pixels[L]))] + [str(np.std(peak_pixels[L]))])


#IPython.embed()
#summary for all files
ss.writerow(['all files'] + [str(sum([len(layer) for layer in consolidated_peak_stack]))] + [str(np.mean([peak[7] for layer in consolidated_peak_stack for peak in layer]))] +
[str(np.std([peak[7] for layer in consolidated_peak_stack for peak in layer]))] + [str(np.mean([peak[8] for layer in consolidated_peak_stack for peak in layer]))] +
[str(np.std([peak[8] for layer in consolidated_peak_stack for peak in layer]))] +
[str(np.mean([x for layer in background_pixels for x in layer]))] +
[str(np.std([x for layer in background_pixels for x in layer]))] +
[str(np.mean([x for layer in peak_pixels for x in layer]))] +
[str(np.std([x for layer in peak_pixels for x in layer]))])

#(persistent_bin, flash_bin, unplotted_bin, persistent_bin2,
#            flash_bin2, last_starts, last_stops, cumulative_blinking,
#            persistent_boolean, restart_wiggle, continued_wiggle,
#            blinking_cube, start_stop_cache) = pflibMP.blink_analysis(peak_cube, peak_data_cube, r_2_threshold=r_2_threshold, wiggle=3,
#                   silent=False)

#pflibMP.save_labeled_images(images, blinking_cube, drift=drift,
#                          output_filename_hash=("run " +
#                                          str(output_filename_hash)),
#                          size=4, silent=False)


#data_nexus = pflibMP.correlation_shufti(persistent_bin, flash_bin, unplotted_bin,
#                       persistent_bin2, flash_bin2, cumulative_blinking,
#                       last_starts, last_stops, restart_wiggle,
#                       continued_wiggle, drift=drift, wiggle=3)

#cPickle.dump(data_nexus, open("data_nexus " + str(output_filename_hash) + ".pkl", "w"))

linked_consolidated_peak_stack, phantom_dictionary, phantom_inverse_dictionary, phantom_hashes = pflibMP.linked_analysis(consolidated_peak_stack, wiggle=3, silent=False)

#cPickle.dump((linked_consolidated_peak_stack, phantom_dictionary, phantom_inverse_dictionary, phantom_hashes), open("linked_analysis " + str(output_filename_hash) + ".pkl", "w"))

pflibMP.save_linked_images(images, phantom_inverse_dictionary, phantom_hashes, linked_consolidated_peak_stack, drift=drift,
                       output_filename_hash=("LINKED run " +
                                          str(output_filename_hash)), brightness_shift=brightness_shift)

#chain_length = 5
#culled_phantoms = [phantom for phantom in phantom_inverse_dictionary if len(phantom_inverse_dictionary[phantom]) >= chain_length]
#pk = [len([peak for peak in layer if phantom_dictionary[(peak[4], peak[5], L)] in culled_phantoms]) for L, layer in enumerate(linked_consolidated_peak_stack)]
#logging.info("number of peaks in each frame that last longer than " + str(chain_length) + " frames")
#for f, frame in enumerate(pk):
#    logging.info("frame number " + str(f) + " has " + str(pk[f]) + " peaks that last longer than " + str(chain_length) + " frames")
#chain_length = 1
#culled_phantoms = [phantom for phantom in phantom_inverse_dictionary if len(phantom_inverse_dictionary[phantom]) >= chain_length]
#pk = [len([peak for peak in layer if phantom_dictionary[(peak[4], peak[5], L)] in culled_phantoms]) for L, layer in enumerate(linked_consolidated_peak_stack)]
#logging.info("number of peaks in each frame that last longer than " + str(chain_length) + " frames")
#for f, frame in enumerate(pk):
#    logging.info("frame number " + str(f) + " has " + str(pk[f]) + " peaks that last longer than "  + str(chain_length) + " frames")

#IPython.embed()

for filename in image_list:
    os.remove(filename)

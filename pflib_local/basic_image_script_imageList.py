#!/home/boulgakov/anaconda2/bin/python

"""
Sample of the most basic image processing script.

Invoke with '--help' and see pflib.py for more documentation.
"""

import pflib
import logging
import os
import os.path
import datetime
import argparse
import time


#get timestamp that will be used throughout this run
timestamp_epoch = time.time()
timestamp_datetime = datetime.datetime.fromtimestamp(timestamp_epoch)

#parse arguments
script_description_string = ("Will traverse all target_directories and "
                             "process all found *.tif files through "
                             "pflib.parallel_image_batch. ")
parser = argparse.ArgumentParser(description=script_description_string)
num_processes_help_string = ("Number of processes to use. Default defined by "
                             "pflib.parallel_image_batch.")
parser.add_argument('-n', '--num_processes', type=int, nargs=1, default=[None],
                    help=num_processes_help_string)
default_log_directory = '/project2/marcotte/boulgakov/microscope/pf_log_files'
default_log_filename = 'basic_image_script_' + str(timestamp_datetime) + 'imageList.log'
log_path_help_string = \
                   ("Pathname for log. Default is " +
                    "/project2/marcotte/boulgakov/microscope/pf_log_files/"
                    "basic_image_script_ + str(datetime.datetime.now()) + .log"
                    ". If the log file already exists, further logging output "
                    "follows default Python behavior, which is currently to "
                    "append the file.")
default_log_path = os.path.join(default_log_directory, default_log_filename)
parser.add_argument('-L', '--log_path', nargs=1, default=[default_log_path],
                    help=log_path_help_string)
#target_directories_help_string = ("Directories to process. At least one must "
#                                  "be specified.")
#parser.add_argument('target_directories', nargs='+',
#                    help=target_directories_help_string)
target_imageList_help = ("A list of images to process. A full path must be"
                         "provided.")
parser.add_argument('target_list',nargs='+',
                    help=target_imageList_help)

args = parser.parse_args()

#normalize all target directories to absolute paths
#target_directories = [os.path.abspath(d) for d in args.target_directories]
target_images_lists = args.target_list

#setup logging for debug and error statements
logging.basicConfig(filename=args.log_path[0], level=logging.DEBUG)
logger = logging.getLogger()
logger.info("basic_image_scipt starting at " + str(timestamp_datetime))

#find all tif files in directories
target_images = []
for imageList_file in target_images_lists:
    with open(imageList_file) as ifile:
        lines = ifile.readlines()
    for f in lines:
        f = f.rstrip()
        if f.endswith('.tif'):
            target_images.append(f)

"""
Simple hack to just compile files with maxProjImage.tif'
for target_dir in target_directories:
    for root, subfolders, files in os.walk(target_dir):
        for f in files:
            if f.endswith('maxProjImage.tif'):
                target_images.append(os.path.join(root, f))
"""
                    
#confirm to log what files will be processed
logger.info("Scanned target directories\n"  + '\n'.join(target_images_lists))
logger.info("Will process target images\n" + '\n'.join(target_images))

#this does all the work!
processed_images = \
                pflib.parallel_image_batch(target_images,
                                           timestamp_epoch=timestamp_epoch,
                                           num_processes=args.num_processes[0])

#summarize current run
logger.info("Pathnames of images processed: " +
            str('\n'.join(processed_images.keys())))
logger.info("basic_image_scipt finished at " + str(datetime.datetime.now()))

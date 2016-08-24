#! /home/jaggu/anaconda
"""
# ### General image statistics for the time traces
# The statistics are - (a) Photobleaching curve (b) Lifetime histogram (i.e frequency of frames with peptides till the first OFF state) and (c) Number of frames ON. 
# (a) is generated as plots and a txt file. (b) and (c) are generated as list in the pickle file. 
"""
from __future__ import division
import collections
import sys
import cPickle as pickle
import os
import time
import numpy as np
import matplotlib.pyplot as plt
import cv2
import numpy as np
import subprocess

def savePkl(db, fpath):
    ofile = open(fpath, 'w')
    pickle.dump(db, ofile)
    ofile.close()
    return True


def loadPkl(pklPath, f):
    fpath = os.path.join(pklPath, f)
    ifile = open(fpath, 'r')
    db = pickle.load(ifile)
    ifile.close()
    return db


def simpleShow(img):
    plt.imshow(img, cmap='gray')
    plt.show()
    return True


def makeDir(d):
    if not os.path.exists(d):
        os.makedirs(d)
    return d


def model_func(t, A, K, C):
    return A * np.exp(K * t) + C


def fit_exp_linear(t, y, C = 0):
    y0 = np.array(y, dtype=float)
    y = np.array(y, dtype=float)
    t = np.array(t, dtype=int)
    y = y - C
    y = np.log(y)
    K, A_log = np.polyfit(t, y, 1)
    A = np.exp(A_log)
    y_fit = model_func(t, A, K, C)
    rmse = np.sqrt(((y_fit - y0) ** 2).mean())
    return (A, K, y_fit, rmse)


def computeFrameONs(pklPath):
    peak_status_fname = 'peak_status.dict.pkl'
    peak_status_dict = loadPkl(pklPath, peak_status_fname)
    frames_tillOff = list()
    nbrFrames = len(peak_status_dict.values()[0])
    sum_array = np.array([False] * nbrFrames, dtype=bool)
    for pos, status_list in peak_status_dict.items():
        try:
            frame_firstOff = next((i for i, status in enumerate(status_list) if status is False))
        except StopIteration:
            frame_firstOff = len(status_list)

        frames_tillOff.append(frame_firstOff)
        sum_array = np.vstack((sum_array, status_list))

    frameCounts = np.sum(sum_array, axis=0)
    peakONCounts = np.sum(sum_array, axis=1)
    return (frameCounts, peakONCounts, frames_tillOff)


def plot_pbleaching(frameCounts, trace_outputDir):

    def _plot(y, norm = False):
        t = range(len(y))
        A, K, y_fit, rmse = fit_exp_linear(t, y)
        plt.plot(t, y, ls='None', marker='o', markerfacecolor='#0072B2', markersize=6)
        plt.plot(t, y_fit, color='#E69F00', linewidth=2)
	plt.ylim(ymin=0)
        plt.xlabel('Frames')
        plt.ylabel('Number of peptides ON')
        if norm:
            fname1 = 'pbleaching.norm'
        else:
            fname1 = 'pbleaching'
        plt.savefig(os.path.join(pbleach_outputDir, fname1 + '.svg'))
        plt.savefig(os.path.join(pbleach_outputDir, fname1 + '.jpg'))
        plt.close()
        return (A, K, rmse)

    def _halfLife(K):
        t_half = np.log(2) / -K
        return t_half

    pbleach_outputDir = os.path.join(trace_outputDir, 'pbleaching')
    makeDir(pbleach_outputDir)
    param = _plot(frameCounts)
    t_half = _halfLife(param[1])
    norm_frameCounts = [ float(item / max(frameCounts)) for item in frameCounts ]
    param_norm = _plot(norm_frameCounts, norm=True)
    t_half_norm = _halfLife(param_norm[1])
    fname = os.path.split(os.path.split(pbleach_outputDir)[0])[1] + '.pbleach.tab'
    ofile = open(os.path.join(pbleach_outputDir, fname), 'w')
    ofile.write('#A \t #K \t #RMSE \t #HALF-LIFE (/FRAME) \n')
    ofile.write(str(param[0]) + '\t' + str(param[1]) + '\t' + str(param[2]) + '\t' + str(t_half) + '\n')
    ofile.write('#NORMALIZED \n')
    ofile.write(str(param_norm[0]) + '\t' + str(param_norm[1]) + '\t' + str(param_norm[2]) + '\t' + str(t_half_norm) + '\n')
    ofile.close()
    return True


def pickle_peakONs(peakONCounts, frames_tillOff, pklPath):
    pkl_f1 = os.path.join(pklPath, 'peakONCounts.list.pkl')
    pkl_f2 = os.path.join(pklPath, 'frames_tillOff.list.pkl')
    savePkl(peakONCounts, pkl_f1)
    savePkl(frames_tillOff, pkl_f2)
    return (pkl_f1, pkl_f2)


def analyze_pbleaching(traceFile, outputDir):
    # traceFile - the name of the subDir containing the directory to process 
    print "Analyzing photobleaching for ",traceFile, "..."
    #trace_outputDir = os.path.join(outputDir, traceFile + '_first_output')
    trace_outputDir = os.path.join(outputDir, traceFile + '_tif_output')
    pklPath = os.path.join(trace_outputDir, 'pklFiles')
    frameCounts, peakONCounts, frames_tillOff = computeFrameONs(pklPath)
    plot_pbleaching(frameCounts, trace_outputDir)
    pickle_peakONs(peakONCounts, frames_tillOff, pklPath)
    return True


def convert_timetraces(traceFile, sourceDir, outputDir, autoscale = False):
    scale = 'not_scaled'
    traceDir = os.path.join(sourceDir, traceFile)
    allInput = [ os.path.join(traceDir, f) for f in os.listdir(traceDir) if f.endswith('.tif') ]
    if autoscale:
        scale = 'scaled'
    trace_outputDir = os.path.join(outputDir, traceFile + '_first_output', scale + '_images')
    makeDir(trace_outputDir)
    for inputF in allInput:
        ofname = os.path.split(inputF)[1].split('.')[0] + '.' + scale + '.tif'
        outputF = os.path.join(trace_outputDir, ofname)
        if autoscale:
            cmd = 'convert -contrast-stretch 0.015x0.05% ' + inputF + ' ' + outputF
            subprocess.call(cmd.split(), shell=False)
        else:
            cmd = 'convert ' + inputF + ' ' + outputF
            FNULL = open(os.devnull, 'w')
            subprocess.call(cmd.split(), shell=False, stdout=FNULL, stderr=subprocess.STDOUT)

    return True

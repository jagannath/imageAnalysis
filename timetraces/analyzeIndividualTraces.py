#! /home/jaggu/anaconda

import os
import sys
import cPickle as pickle
import os
import time
import fnmatch
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
plt.style.use('ggplot')
sys.path.append('/project/current/project/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis/proteanseq/pflib')
import flexlibrary

def savePkl(db, fpath):
    ofile = open(fpath, 'w')
    pickle.dump(db, ofile)
    ofile.close()
    return True


def loadPkl(fpath):
    ifile = open(fpath, 'r')
    db = pickle.load(ifile)
    ifile.close()
    return db


def makeDir(d):
    if not os.path.exists(d):
        os.makedirs(d)
    return d


def locate(pattern, root = os.curdir):
    """Locate all files matching supplied filename pattern in and
    below supplied root directory."""
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            allFiles.append(os.path.join(path, filename))

    return allFiles


def plot_individualTraces(f, outputDir):

    def _yFit(plateau_trace):
        y_fit = list()
        for start, stop, ht in plateau_trace:
            intensity_fit = [ht] * (stop - start)
            y_fit.extend(intensity_fit)

        return y_fit

    def _figure(trace, y_fit, bg_trace, fig_outputDir):
        left, width = (0.1, 0.65)
        bottom, height = (0.1, 0.65)
        bottom_h = left_h = left + width + 0.02
        rect_trace = rect_scatter = [left,
         bottom,
         width,
         height]
        rect_histy = [left_h,
         bottom,
         0.2,
         height]
        plt.figure(1, figsize=(8, 8))
        ax_trace = plt.axes(rect_trace)
        ax_histy = plt.axes(rect_histy)
        nullfmt = NullFormatter()
        ax_histy.yaxis.set_major_formatter(nullfmt)
        ax_trace.plot(trace, color='#0072B2')
        ax_trace.plot(y_fit, color='#D55E00', lw=2)
        ax_trace.plot(bg_trace, color='black')
        ax_histy.hist(trace, 100, orientation='horizontal', facecolor='#E69F00')
        ax_trace.set_ylabel('Intensity (sum of 3x3 area)')
        ax_trace.set_xlabel('Frames (/s or /100 ms)')
        fname = 'trace_' + str(pos)
        f = os.path.join(fig_outputDir, fname + '.png')
        plt.savefig(f, dpi=300)
        plt.close()
        return True

    peak_intensity_dict = loadPkl(f)
    f_bg = os.path.join(os.path.split(f)[0], 'peak_bg_median.dict.pkl')
    peak_bg_median_dict = loadPkl(f_bg)
    traceFile = f.split('/')[-3]
    fig_outputDir = os.path.join(outputDir, traceFile, 'individualTraces')
    makeDir(fig_outputDir)
    pklDir = os.path.split(f)[0]
    peak_plateau_dict = dict()
    for pos, trace in peak_intensity_dict.items():
        bg_trace = peak_bg_median_dict[pos]
        h, w = pos
        t = flexlibrary.PhotometryTrace(trace, h, w)
        photometries, unmirrored_ck_filtered_photometries, unmirrored_plateaus, unmirrored_t_filtered_plateaus = t.stepfit_photometries(h, w)
        plateau_trace = unmirrored_t_filtered_plateaus.trace
        y_fit = _yFit(plateau_trace)
        diffTrace = [ trace[i] - bg_trace[i] for i in range(len(trace)) ]
        _figure(trace, y_fit, bg_trace, fig_outputDir)
        peak_plateau_dict[pos] = plateau_trace

    pkl_fname = os.path.join(pklDir, 'peak_plateau.dict.pkl')
    savePkl(peak_plateau_dict, pkl_fname)
    return True


def plotTraces(outputDir):
    pklFile_pattern = 'peak_intensity.dict.pkl'
    f_list = locate(pklFile_pattern, outputDir)
    for f in sorted(f_list):
        print 'Processing ', f.split('/')[-3], ' ...'
        plot_individualTraces(f, outputDir)

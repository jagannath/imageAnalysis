#!/home/jaggu/anaconda/bin/python2.7


"""
A python script that executes a number of analysis scripts. The basic image
script must have been completed (maybe TACC) and perhaps sync'ed back to the
appropriate directory and server (preferably hopper)
Each of the subsequent scripts can be run independantly. But it requires remembering the enter input arguments.
The script essentially calls the independent scripts by subprocess modules and
keeps tabs on the output file.
"""

import os
import sys
import subprocess
import cPickle as pickle
import argparse
import shlex 
import subprocess
import datetime

def getTimeNow():
    time_now = '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
    return time_now

def run_experimentScript(**kwargs):
    [fnamePattern] = kwargs['fnamePattern'] 
    nbrChannels = kwargs['nbrChannels']
    dateStamp = '{:%Y%m%d}'.format(datetime.datetime.now())
    output_directory = os.path.abspath("output_"+dateStamp)

    if nbrChannels is 3:
        alignment_file_pattern = os.path.join(fnamePattern,fnamePattern+'xy*c3.tif')
        peptide_file_pattern = os.path.join(fnamePattern,fnamePattern+'xy*c1.tif')
        second_channel_pattern = os.path.join(fnamePattern,fnamePattern+'xy*c2.tif')
    else:
        raise SystemExit("Number of channels is not 3.")
    cmdLine = "/project/boulgakov/microscope/basic_experiment_script.py \
            -L expt.log --output_directory "+output_directory + "\
            --alignment_files "+alignment_file_pattern + "\
            --peptide_files "+peptide_file_pattern + "\
            --second_channel "+second_channel_pattern +"\
            >expt.out 2>expt.err"
    
#    cmdLine = "/project/boulgakov/microscope/basic_experiment_script.py \
#            -h "
    args = shlex.split(cmdLine)
    file_out = open("expt3.out",'w')
    file_err = open("expt3.err",'w')
    subprocess.call(args,stdout=file_out,stderr=file_err)
    file_out.close()
    file_out.close()    file_err.close()
    file_out.close()    
    time_finished = getTimeNow()
    print "basic_experiment_script.py : completed on ",time_finished 

    return output_directory


def run_remainderScript(output_directory,time_hash=None):
    def _get_time_hash():
        category_stats_file = [f for f in os.walk(output_directory).next()[2] if f.startswith('category_stats_') and f.endswith('pkl')]
        try:
            assert len(category_stats_file) is 1
        except AssertionError:
            raise SystemExit("Ensure only one time_hash file in the output directory. Enter time_hash code")
        [fname] = category_stats_file
        time_hash = fname.split('_')[-1][:-4]
        return time_hash

    output_directory = os.path.join(os.path.split(output_directory)[0],'output160904')
    if not time_hash:
        time_hash = _get_time_hash()
    else:
        time_hash = time_hash
   
    track_photometry_file = os.path.join(output_directory,'track_photometries_NO_NONES_'+time_hash+'.csv')

    cmdLine = "/home/boulgakov/peptidesequencing/git/proteanseq/pflib/remainder_correction.py "\
            +track_photometry_file

    args = shlex.split(cmdLine)
    file_out = open("remainder.out",'w')
    file_err = open("remainder.err",'w')
    subprocess.call(args,stdout=file_out,stderr=file_err)
    file_out.close()
    file_err.close()

    time_finished = getTimeNow()
    print "remainder_correction.py : completed at ",time_finished

    remainder_adjusted_file = track_photometry_file + '_adjusted.csv'
    return remainder_adjusted_file

def run_logfitterScript(remainder_adjusted_file,**kwargs):
    def _runChannel(ch):
        if ch == '561':
            s = kwargs['sequence_561']
            c = 1
        elif ch == '647':
            s = kwargs['sequence_647']
            c = 2
        else:
            raise SystemExit("Something wrong with channels")

        cmdLine = "/home/boulgakov/peptidesequencing/git/proteanseq/pflib/lognormal_fitter.py " \
                + remainder_adjusted_file + \
                " -c "+str(c) + \
                " -m "+str(nbrMock) + \
                " -o "+str(nbrOmit) + \
                " -e "+str(nbrEdman) + \
                " -s "+s + \
                " --last_drop_method"
        args = shlex.split(cmdLine)
        file_out = open("lognormal."+str(ch)+".out",'w')
        file_err = open("lognormal."+str(ch)+".err",'w')
        subprocess.call(args,stdout=file_out,stderr=file_err)
        file_out.close()
        file_err.close()
        return True

    nbrMock = kwargs['nbrMock']
    nbrOmit = kwargs['nbrOmit']
    nbrEdman = kwargs['nbrEdman']
    nbrChannels = kwargs['nbrChannels'] 
    if not nbrChannels is 3:
        raise SystemExit("Number of channels is not 3.")
    
    _runChannel('561')
    _runChannel('647')
    time_finished = getTimeNow()
    print "lognormal_fitter.py : completed at ",time_finished
    return True


def run_summarizeExperiments(remainder_adjusted_file,**kwargs):
    nbrMock = kwargs['nbrMock']
    nbrEdman = kwargs['nbrEdman']
    nbrCycles = nbrMock + nbrEdman

    cmdLine = "/project/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis/summarizeExperiments.py " \
            + remainder_adjusted_file + \
            " --nbrCycles " + str(nbrCycles) \
    
    args = shlex.split(cmdLine)
    file_out = open("summary.out",'w')
    file_err = open("summary.err",'w')
    subprocess.call(args,stdout=file_out,stderr=file_err)
    file_out.close()
    file_err.close()
    
    time_finished = getTimeNow()
    print "summarizeExperiments.py : completed at ",time_finished
    return True



parser = argparse.ArgumentParser(description=""" Run all the analysis script; \n 
                                  [1] basic_experiment_script.py \n
                                  [2] remainder_correction.py \n
                                  [3] lognormal_fitter.py \n
                                  [4] summarizeExperiments.py 
                                  """)
parser.add_argument('--file_pattern',nargs='+',action="store",dest="fnamePattern",type=str,
                    help="File name pattern (the directory structure where every cycle images are stored) ")
parser.add_argument('--mock','-m',action="store",dest="nbrMock",type=int,default=4)
parser.add_argument('--omit','-o',action="store",dest="nbrOmit",type=int,default=4)
parser.add_argument('--edman','-e',action="store",dest="nbrEdman",type=int,default=4)
parser.add_argument('--nbr_channels','-nc',action="store",dest="nbrChannels",type=int,default=3)
parser.add_argument('--sequence1','-s1',action="store",dest="sequence_561",type=str,default="PEPTIDE_561")
parser.add_argument('--sequence2','-s2',action="store",dest="sequence_647",type=str,default="PEPTIDE_647")

args = parser.parse_args()

print "Starting at : ",getTimeNow()
output_directory = run_experimentScript(**vars(args))
remainder_adjusted_file = run_remainderScript(output_directory)
run_logfitterScript(remainder_adjusted_file,**vars(args))
run_summarizeExperiments(remainder_adjusted_file,**vars(args))


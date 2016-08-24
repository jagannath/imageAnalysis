#!/home/jaggu/anaconda

"""
Analyzing time traces from a source Dir

"""
import sys
import cPickle as pickle
import os
import time
import fnmatch
sys.path.append(os.path.join('/home/jaggu/marcotte_project/current/project/jagannath/projectfiles/proteomics_single_molecule/imageAnalysis/timetraces'))
import picklePeakInfo
import generalImageStatistics
import createProjectionImages_peakFinding as projectionImage

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and   
    below supplied root directory.'''
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            allFiles.append(os.path.join(path,filename))
    return allFiles


def projectImages(sourceDir,proj_type ='first',pattern='trace'):
	projectionImage.createImages(sourceDir,proj_type)
	print "first image projection done..Now run basic image script"

def pickle_peakInfo_allTraces(sourceDir,outputDir,peakID_files):
    for peakFilePath in peakID_files:
        t0 = time.ctime()
        traceDir, pklPath = picklePeakInfo.analyzeFrames(sourceDir,outputDir,peakFilePath)        
        print t0, time.ctime()

def computeGeneralStatistics(sourceDir, outputDir, peakID_files):
    for peakFilePath in peakID_files:
        print sourceDir, outputDir, peakFilePath
        traceFile = os.path.split(peakFilePath)[1].split('.')[0]
        generalImageStatistics.analyze_pbleaching(traceFile, outputDir)

def summarizePhotobleaching(outputDir):
    
    # Parses through the frames_tillOff.list.pkl and obtains the number of peaks that remainded after a given frame
    def _remainders(frames_tillOff_list, frameCutoff=8):
        remainders = 0
        for frames_tillOff in frames_tillOff_list:
            if frames_tillOff >= frameCutoff: remainders +=1 
        return remainders
    
    def summarize_framesTillOff(outputDir):
        outputF = os.path.join(outputDir,'SUMMARIZED_all_framesTillOff.tab')
        ofile = open(outputF, 'w')
        ofile.write('# Frame Cutoff = 8 \n')
        ofile.write('# Experiment Name \t Total Peaks \t Remainder Peaks \n')
        frames_tillOff_files = locate('frames_tillOff.list.pkl',outputDir)
        for frames_tillOff_f in frames_tillOff_files:
            exptName = frames_tillOff_f.split('/')[-3][:-13] #the file trails _first_output
            frames_tillOff_list = pickle.load(open(frames_tillOff_f))
            ofile.write('\t'.join([exptName,str(len(frames_tillOff_list)), str(_remainders(frames_tillOff_list)),'\n']))
        ofile.close()
        return outputF

    def _halfLife(f):
        ifile = open(f)
        lines = ifile.readlines()
        halfLife = lines[1].strip().split('\t')[3]
        ifile.close()
        return halfLife

    def summarize_halfLifeStats(outputDir):
        outputF = os.path.join(outputDir, 'SUMMARIZED_all_halfLife.tab')
        ofile = open(outputF, 'w')
        ofile.write('# Experiment Name \t Half Life \n')
        halfLife_files = locate('*pbleach.tab',outputDir)
        for halfLife_f in halfLife_files:
            exptName = os.path.split(halfLife_f)[1].split('.')[0][:-13]
            halfLife = _halfLife(halfLife_f)
            ofile.write(exptName + '\t' + str(halfLife) + '\n')
        ofile.close()
        return outputF
    
    outputF_framesTillOff = summarize_framesTillOff(outputDir)
    outputF_halfLifeStats = summarize_halfLifeStats(outputDir)
    return True  
        
        
def main():
    proj_type = 'first'
    pattern = 'trace'
    
    import socket
    hostname = socket.gethostname()                                    
    if hostname == 'canopus': headDir = "/home/jaggu/marcotte_project" 
    else: headDir = "/project"                                         
    
    sourceDir = os.path.join(headDir,'boulgakov/microscope/2014-Dec/2014-12-22t1')
    peakSourceDir = os.path.join(sourceDir,'projectedImages',proj_type)
    if not os.path.exists(peakSourceDir): os.makedirs(peakSourceDir)
    
    outputDir = os.path.join(headDir,'jagannath/projectfiles/singleMoleculeMicroscopy/dataAnalysis/microscope1/2014-Dec/2014-12-22t1/output')
    if not os.path.exists(outputDir): os.makedirs(outputDir)

    # Step 1
    #projectImages(sourceDir,proj_type='first',pattern='trace')
    # Step 2: Run basic_image_script.py
    # Step 3: Run through the identified peaks for each traceDir and pickle the information
    peakID_files = [os.path.join(peakSourceDir,f) for f in os.listdir(peakSourceDir) if f.endswith('pkl')]
    #pickle_peakInfo_allTraces(sourceDir, outputDir, peakID_files)
    # Step 4: From the pickled information, compute the general statistics;
    computeGeneralStatistics(sourceDir, outputDir, peakID_files)
    # Step 5: Summarize Photobleaching
    #summarizePhotobleaching(outputDir)

main()

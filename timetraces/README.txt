The idea of the timetrace folder is to do image processing of the time traces and summarize important results such as photobleaching etc.
0. When generating time traces, make it such that one directory has one channel only and the trace tifs. (exporting the nd2 file to directory)
1. The "analyzeTimeTraces.py" script has all the relevant modules that is loaded, but needs to be done in steps
2. Run ProjectImages function first; Currently it is set to only copy the first image of the traceDir (sorted). 
3. This creates a folder in the source Dir with projectedImages/[projection type]. Run basic_image_script.py on this directory
4. Perform analysis on the identified peaks and quantitate them through time. It uses picklePeakInfo module. One could make changes to the peak size etc in it. This function generates the pickled information of the peaks. 
5. The next function is computeGeneralStatistics where the general statistical information and photobleaching curves are generated. 
6. The photobleaching statistics are finally computed. 


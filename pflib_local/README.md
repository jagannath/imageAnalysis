Please do not modify pflib.py. If you wish to add to it, please branch first.
PFLIB has not been modified. This is a local copy where additional files are being made

basic_image_script.py is a sample script that uses pflib. Try it with the --help option and feel free to read the source.
basic_image_script_imageList.py is changes to the basic_image_script.py and has a basic_image_script_imageList.jag.py soft link in boulgakov/microscope/ folder. It accepts a file which contains a list of all image files that needs to be processed. Typically it is a set of tiff files with '.maxProjection.tif'

The last image pipeline was hard to develop and add features because it was not well modularized. I was becoming confused about my own code! If we can keep things modularized, the pipeline should hopefully grow naturally and easily.

The pipeline's modularity is based on independent peak fitting of every image. The first step of any analysis is to fit peaks for all images involved. In this step, the peaks (point spread functions, or PSFs) are found. The peak finding is agnostic to the image's underlying content or meaning. PSFs are stored independently for each image. PSF results can be stored in different ways. The three default formats are the Python pickle format, CSV, and PNG. Further formats are added as needed.

For example, given the file 'peptides.tif', the first step of the pipeline would be to generate (by default) the files:

1. peptides.tif.png._psfs_timestamp.pkl -- contains the PSFs found in peptides.tif as a pkl
2. peptides.tif.png._psfs_timestamp.csv -- contains the PSFs found in peptides.tif as a csv
3. peptides.tif.png._psfs_timestamp.png -- contains the PSFs found in peptides.tif as a png

All subsequent image analysis depends on combining information contained in these files through utility scripts. Ideally, these should be designed a la UNIX: each component does one thing and do it well, and they can be chained together.

For example, one of the utility scripts can be given as input a sequence of pkl files containing PSFs of fiduciary marker dyes in one field, and yields as output the list of translations/rotations aligning the image sequence.

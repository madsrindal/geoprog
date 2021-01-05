# TBA4251 - Programming in Geomatics

This is a description of how to run the code made for computations of geoidal heights and spherical harmonic coefficients. 

In order to create a geoid model based on EGM2008 from raw data, the data file needs to be downloaded from http://icgem.gfz-potsdam.de/tom_longtime. Select the EGM2008 model with GRS90 as the reference system, and save the file in the Datafiles/Part1 folder in the project. The file should be on format .gfc.

The code is also formatted to be able to run from pre-processed data. In order to do so, do the following:

1. Clone the locally and install the packages used in the project
2. Run the files located in the Plot folders in both part 1 and 2

This will create the same plots presented in the report.

There is also created two test scripts, one for each part, that presents most of the steps in the computation process.
These run with the use of raw data, but the script related to part 1 can be ran without the download of the EGM2008 data. In this case, change the boolean egm-value from 'True' to 'Flase' in the main function at the end. 

Both scripts are written with a much lower degree of n_max than used in the project paper in order to obtain faster results. The degree can also be changed if desired, but a higher degree may result in a sufficent increase in time complexity. 

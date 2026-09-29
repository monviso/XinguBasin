# Forest structure types in the southern Amazon Basin are influenced by anthropogenic edges and physiography.

This collection of R script allow to reproduce all the workflow of the "Forest structure types in the southern Amazon Basin are influenced by anthropogenic edges and physiography". Please cite the paper if you use this scripts applied to your research.

## 1 - Preliminary steps and disclaimers
- As a preliminary step download all the data from Zenodo.\
- Steps 01-09 where run locally on a Dell Precision 3680 , Processor:	Intel(R) Core(TM) i7-14700K, 3400 Mhz, 20 Core(s), 28 Logical Processor(s).\
- To run ALS steps (07, 08) ALS dataset that overlap spatially with GEDI observations in the Xingu basin must be downloaded locally. The file .... reports the download path of each file.
- Steps 13 -14 and 22 where implemented on the Northern Arizona University HPC Monsoon (https://in.nau.edu/arc/details/) due to the high computational demand. Step 13 and 14 can be run directly using the Xingu_GAM_prepared.rds file (Zenodo). We do not report here the predictor calculation scripts as they were extracted from a new Amazon-wide landscape physiography descriptor raster product that is currently under review. The .rds file already encapsulate all the required information to run the GAM fitting.
- Steps 15-17 and 23-24 can be run again locally.\
- Missing script numbers are ancillary scripts already integrated in the current script set version.
- All the scripts must be carefully edited to match your local machine/HPC directory structure.
- All intermidiate results and plots proidcued here can be direcly downlaoded from Zenodo.

## 2 - Steps 01-09
The scripts 01-09 cover the NMF and k-means classification, robustness tests, radiometric and ALS characterization. Again, to run ALS steps (07, 08) ALS dataset that overlap spatially with GEDI observations in the Xingu basin must be downloaded locally. The file .... reports the download path of each file. 

## 3 - Steps 

 

#### Small-scale fishing modelling using Bayesian networks (a.k.a. belief propagation networks)

## Author: Phillip P.A. Staniczenko  pstaniczenko@brooklyn.cuny.edu
## Please cite: 
## Social ties explain catch portfolios of small-scale fishers in the Caribbean
## Alexander, S.M., Staniczenko, P.P.A. & Bodin, O.
## Fish and Fisheries, vol. X, ppXX--XX (2019)

## Date: 7th October 2019

#### Notes
## Before running the R script, SSF_BNAnalysis_Staniczenko_distribute.R, you must first compile for your specific computer the C source file
## fNMLofDAG.c
## to generate the executable file
## fNMLGet
##
## To run the R script, ensure all required files are in your working directory, then in R type
## source("SSF_BNAnalysis_Staniczenko_distribute.R")

#### Outputs
## Total lengths for model-data combinations;
## prints results to screen and also saves results as a csv file in the working directory

#### Required program files
## ZipData.py
## fNMLGet executable (must be compiled for you computer from source file, fNMLofDAG.c)

#### Required data files
## manID.csv
## manTrait_metadata.csv (heirarchy is column 3)
## fishID.csv
## fishTrait_metadata.csv
## catchportfolio_manfishedgelist.csv
## socialnetwork_edgelist.csv
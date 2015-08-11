#
# 
# File for compiling all files into CM and DrEEMS
# 10aug2015
# Ashlee Jollymore's phd
#######################
## set working directory
rm(list = ls())
ls()

## necessary packages
library(reshape)
library(plyr)
library(gsubfn)

##########
all.EEMS <- prechlor.files <- "/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_all_corrected"

# CM directory
directoryCM <-"/Users/ashlee/Documents/MATLAB/CorrectedEEMS" 

project <- "DBPall"

##########################################################################
# Compile all files and make graph headings file for running CM on all EEMS
setwd(all.EEMS) 
filelist_EEMScor <- list.files(pattern = ".csv$")

ex.PARAFAC <- seq(240, 800, by = 2) #change if excitation wavlenegths are different

# call function
setwd("/Users/ashlee/SpecScripts") 
source("EEMCMtrimDBPall_function.R")
CMsave <- CMtrim(filedirectory = all.EEMS, filelist = filelist_EEMScor, project = project, exmin = "X240",
                 directoryCM = directoryCM, ex = ex.PARAFAC)

####################################################################################################################################
######### DOM Fluor
########
# Get files ready for DOMFLuor toolbox. 
# Need .csv file for ex, em and one csv file containing all of the fluorescence EEMS compiled
# Use Save Dr EEMS function
# Inputs include the filelist, the project, the vector containing sample names, and the excitation wavelength min you want to trim to
# File cuts EEMs from ex min that you want to

setwd("/Users/ashlee/SpecScripts") 
source("EEMSDrEEMsaveDBPall_function.R")

ex.DrEEMS = seq(240, 800, by = 2)
DrEEM.data = DrEEM(filelist = filelist_EEMScor, project = project, 
                   exmin = 'X240', filedirectory = all.EEMS, ex = ex.DrEEMS)

# check DrEEM.data. This is the compiled EEMS for DrEEM PARAFAC modelling
head(DrEEM.data)

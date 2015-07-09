#
# Script for creating delta DBP EEMS
# The aim of this script is to get the eems of the prechlorination minus the post chlorination
# DBP project, ashlee's phd
# 9July2015
# ashlee jollymore
#############

## set working directory
rm(list = ls())
ls()

## necessary packages

######
# load pre and post chlorination files, match according to sample ID, and then minus pre - post
# note that this is only for corrected files.

prechlor.files <- "/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_correctedEEMs"
postchlor.files <- "/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_correctedEEMs"

project <- "DBPdelta"

# directory where the delta eems will go
delta.eems <- "/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_delta"

# create file that matches the pre and post files together

# prechlor files
setwd(prechlor.files)
filelist_pre <- as.data.frame(list.files(pattern = "CorrectedNEW.csv$"))
#postchlor
setwd(postchlor.files)
filelist_post <- list.files(pattern = "CorrectedNEW.csv$")

#get sample ID from both files
n = length(filelist_pre)
for (i in 1:n){
  sampleID.pre.temp <- strapplyc(filelist_pre[i,], "(.*)_DBPPre", simplify = TRUE)
  filelist_pre$sampleID[i] <- sampleID.pre.temp 
}

sampleID.post <- 

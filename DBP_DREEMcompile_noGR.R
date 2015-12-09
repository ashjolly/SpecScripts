# SCript for compiling only streamwater samples for EEMs - removing greenroof and rainwater samples
# 8Dec2015
########

## clear workspace
rm(list = ls())
ls()
#######
# get file list for the blank, absorbance and EEMS file

library(reshape)
library(plyr)
library(gsubfn)
library(gplots)
library(R.matlab)
library(zoo)
library(pracma)
#


# take out the GR samples
GR <- c("DBP0056", "DBP0057", "DBP0058", "DBP0060", "DBP0061", "DBP0062", 
        "DBP0073", "DBP0074", "DBP0075", "DBP0091", "DBP0092", "DBP0093", "DBP0119", "DBP0013")
# chlorinated
GRc <- c("DBPChlor0056", "DBPChlor0057", "DBPChlor0058", "DBPChlor0060", "DBPChlor0061", "DBPChlor0062", 
         "DBPChlor0073", "DBPChlor0074", "DBPChlor0075", "DBPChlor0091", "DBPChlor0092", "DBPChlor0093", "DBPChlor0119", "DBPChlor0013")

#############
# Compile all of the prechlorinated EEMs
directoryCorrectedEEMS <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_correctedEEMS"
project = "DBPPre_noGR"
#
setwd(directoryCorrectedEEMS) 
filelist_EEMScor <- list.files(pattern = "_Corrected.csv$")
# change the sample name column to reflect the pre chlronation format - important for later merging
samples <- sapply(strsplit(as.character(filelist_EEMScor), split='_', fixed=TRUE), function(x) (x[1]))

filelist_EEMScor <- data.frame(cbind(filelist_EEMScor, samples))
# take out the GR samples from the prechlor
filelist_EEMScor <-  filelist_EEMScor[!(filelist_EEMScor$samples  %in% GR),] # from wq
filelist_EEMScor <- as.character(filelist_EEMScor[,1])

setwd("/Users/user/SpecScripts") 
source("EEMSDrEEMsave_function.R")

ex.DrEEMS = seq(240, 800, by = 2)
DrEEM.data = DrEEM(filelist = filelist_EEMScor, project = project, 
                   exmin = 'X240', filedirectory = directoryCorrectedEEMS, ex = ex.DrEEMS, 
                   DrEEMfolder = "/Users/user/Documents/MATLAB/toolbox/CorrEEMS")

# check DrEEM.data. This is the compiled EEMS for DrEEM PARAFAC modelling
head(DrEEM.data)

##########
# Compile all of the post chlorinated EEMs
directoryCorrectedEEMS <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_correctedEEMS"
project = "DBPPost_noGR"
#
setwd(directoryCorrectedEEMS) 
filelist_EEMScor <- list.files(pattern = "_Corrected.csv$")
# change the sample name column to reflect the pre chlorination format - important for later merging
samples <- sapply(strsplit(as.character(filelist_EEMScor), split='_', fixed=TRUE), function(x) (x[1]))

filelist_EEMScor <- data.frame(cbind(filelist_EEMScor, samples))
# take out the GR samples from the prechlor
filelist_EEMScor <-  filelist_EEMScor[!(filelist_EEMScor$samples  %in% GRc),] # from wq
filelist_EEMScor <- as.character(filelist_EEMScor[,1])

setwd("/Users/user/SpecScripts") 
source("EEMSDrEEMsave_function.R")

ex.DrEEMS = seq(240, 800, by = 2)
DrEEM.data = DrEEM(filelist = filelist_EEMScor, project = project, 
                   exmin = 'X240', filedirectory = directoryCorrectedEEMS, ex = ex.DrEEMS, 
                   DrEEMfolder = "/Users/user/Documents/MATLAB/toolbox/CorrEEMS")

# check DrEEM.data. This is the compiled EEMS for DrEEM PARAFAC modelling
head(DrEEM.data)

#################### DBP delta
# Get the corrected EEMS from the delta file
directoryCorrectedEEMS <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_delta"
project = "DBPdelta_noGR"

setwd(directoryCorrectedEEMS) 
filelist_EEMScor <- list.files(pattern = "_DBPdelta.csv$")
# change the sample name column to reflect the pre chlorination format - important for later merging
samples <- sapply(strsplit(as.character(filelist_EEMScor), split='_', fixed=TRUE), function(x) (x[1]))
filelist_EEMScor <- data.frame(cbind(filelist_EEMScor, samples))
# take out the GR samples from the prechlor
filelist_EEMScor <-  filelist_EEMScor[!(filelist_EEMScor$samples  %in% GR),] # from wq
filelist_EEMScor <- as.character(filelist_EEMScor[,1])

setwd("/Users/user/SpecScripts") 
source("EEMSDrEEMsave_function.R")

ex.DrEEMS = seq(240, 800, by = 2)
DrEEM.data = DrEEM(filelist = filelist_EEMScor, project = project, 
                   exmin = 'X240', filedirectory = directoryCorrectedEEMS, ex = ex.DrEEMS, 
                   DrEEMfolder = "/Users/user/Documents/MATLAB/toolbox/CorrEEMS")

# check DrEEM.data. This is the compiled EEMS for DrEEM PARAFAC modelling
head(DrEEM.data)
#################### DBP pre + post
directoryCorrectedEEMS_pre <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_correctedEEMS"
directoryCorrectedEEMS_post <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_correctedEEMS"
project = "DBPprepost_noGR"

setwd(directoryCorrectedEEMS_pre) 
filelist_EEMSpre <- list.files(pattern = "_Corrected.csv$")

setwd(directoryCorrectedEEMS_post) 
filelist_EEMSpost <- list.files(pattern = "_Corrected.csv$")

# change the sample name column to reflect the pre chlorination format - important for later merging
samples <- sapply(strsplit(as.character(filelist_EEMSpre), split='_', fixed=TRUE), function(x) (x[1]))
filelist_EEMSpre <- data.frame(cbind(filelist_EEMSpre, samples))
# take out the GR samples from the prechlor
filelist_EEMSpre <-  filelist_EEMSpre[!(filelist_EEMSpre$samples  %in% GR),] # from wq
filelist_EEMSpre <- as.character(filelist_EEMSpre[,1])

samples <- sapply(strsplit(as.character(filelist_EEMSpost), split='_', fixed=TRUE), function(x) (x[1]))
filelist_EEMSpost <- data.frame(cbind(filelist_EEMSpost, samples))
# take out the GR samples from the postchlor
filelist_EEMSpost <-  filelist_EEMSpost[!(filelist_EEMSpost$samples  %in% GRc),] # from wq
filelist_EEMSpost <- as.character(filelist_EEMSpost[,1])

filelist_EEMSall <- rbind(filelist_EEMSpre, filelist_EEMSpost)
setwd("/Users/user/SpecScripts") 
source("EEMSDrEEMsave_function.R")

ex.DrEEMS = seq(240, 800, by = 2)
DrEEM.data = DrEEM(filelist = filelist_EEMSall, project = project, 
                   exmin = 'X240', filedirectory = directoryCorrectedEEMS, ex = ex.DrEEMS, 
                   DrEEMfolder = "/Users/user/Documents/MATLAB/toolbox/CorrEEMS")

# check DrEEM.data. This is the compiled EEMS for DrEEM PARAFAC modelling
head(DrEEM.data)
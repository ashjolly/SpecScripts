#
# Script for creating delta DBP EEMS
# The aim of this script is to get the eems of the prechlorination minus the post chlorination
# Also, to do the same for the absorbance file and calculate the absorbance indicies from this
# DBP project, ashlee's phd
# 9July2015
# ashlee jollymore
#############

## set working directory
rm(list = ls())
ls()

## necessary packages
library(reshape)
library(plyr)
library(gsubfn)

######
# load pre and post chlorination files, match according to sample ID, and then minus pre - post
# note that this is only for corrected files.

prechlor.files <- "/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_correctedEEMs"
postchlor.files <- "/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_correctedEEMs"

# directory where the delta eems will go
delta.eems <- "/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_delta"

# CM directory
directoryCM <-"/Users/ashlee/Documents/MATLAB/CorrectedEEMS" 

project <- "DBPdelta"

##################
# create file that matches the pre and post files together
######## Prechlor files
setwd(prechlor.files)
filelist_pre <- list.files(pattern = "_Corrected.csv$")
filelist_preAbs <- list.files(pattern = "AbsCorrected.csv$")

# Get sample ID from pre chlor files - EEMS
sample.ID <- 0 #create sample ID variable
n = length(filelist_pre)
for (i in 1:n){
  sampleID.pre.temp <- strapplyc(filelist_pre[i], "DBP(.*)_DBPPre_Corrected.csv", simplify = TRUE)
  sample.ID[i] <- sampleID.pre.temp 
}

filelist_preEEMS <- as.data.frame(cbind(filelist_pre, sample.ID))

###
# Get sample ID from pre chlor files - Abs
sample.ID <- 0 #create sample ID variable
n = length(filelist_preAbs)
for (i in 1:n){
  sampleID.pre.temp <- strapplyc(filelist_preAbs[i], "DBP(.*)_DBPPre_AbsCorrected.csv", simplify = TRUE)
  sample.ID[i] <- sampleID.pre.temp 
}

filelist_preAbs <- as.data.frame(cbind(filelist_preAbs, sample.ID))

filelist_pre <- as.data.frame(merge(filelist_preEEMS, filelist_preAbs, by = "sample.ID"))

######################
# Get sample ID from post chlor files
#postchlor
setwd(postchlor.files)
filelist_post <- list.files(pattern = "_Corrected.csv$")
filelist_postAbs <- list.files(pattern = "AbsCorrected.csv$")

# get sample ID from EEMS
sample.ID <- 0 #reset sample ID variable
n = length(filelist_post)
for (i in 1:n){
  sampleID.post.temp <- strapplyc(filelist_post[i], "DBPChlor(.*)_DBPPost_Corrected.csv", simplify = TRUE)
  sample.ID[i] <- sampleID.post.temp 
}

filelist_postEEMS <- as.data.frame(cbind(filelist_post, sample.ID))

# get sample ID for post Abs
sample.ID <- 0 #reset sample ID variable

n = length(filelist_postAbs)
for (i in 1:n){
  sampleID.post.temp <- strapplyc(filelist_postAbs[i], "DBPChlor(.*)_DBPPost_AbsCorrected.csv", simplify = TRUE)
  sample.ID[i] <- sampleID.post.temp 
}

filelist_postAbs <- as.data.frame(cbind(filelist_postAbs, sample.ID))
filelist_post <- as.data.frame(merge(filelist_postEEMS, filelist_postAbs, by = "sample.ID"))

# merge files by sample ID
all <- merge(filelist_post, filelist_pre, by = "sample.ID", all = FALSE)

Spec.ind = data.frame(matrix(vector(), 0, 8))

#### Run a loop  that loads each file and do the subtraction

n = dim(all)[1]
for (i in 1:n){
  
  ##### Subtract EEM files
  # pre EEM files
  setwd(prechlor.files)
  pre.samplename <- toString(all[i,4])
  pre.temp <- as.data.frame(read.delim(pre.samplename, header= TRUE, sep = ",", stringsAsFactors=FALSE))
  
  #post EEM files
  setwd(postchlor.files)
  post.samplename <- toString(all[i,2])
  post.temp <- as.data.frame(read.delim(post.samplename, header= TRUE, sep = ",", stringsAsFactors=FALSE))
  
  # ensure that EEMS are the same size
  #trim so that exitation and emission goes from the same
  # call EEM trim function
  setwd("/Users/ashlee/SpecScripts") 
  source("EEMtrim_function.R")
  
  post.trim <- EEMtrim(eem = post.temp, minex = "X240")
  pre.trim <- EEMtrim(eem = pre.temp, minex = "X240")

  ## create delta eem variable
  delta.temp <- pre.trim - post.trim
  
  # get generalized sample ID - for example, DBP0001
  sample.ID <- strapplyc(pre.samplename, "(.*)_DBPPre_Corrected", simplify = TRUE)
  
  # save delta eem file
  deltapath <- file.path(delta.eems, paste(sample.ID,"_", project,".csv", sep = ""))
  write.table(delta.temp, file = deltapath, row.names = TRUE,col.names = TRUE, sep = ",")

  ##########
  # Plot delta eems as a contour plot
  # variables to change
  zmax = max(delta.temp,na.rm=TRUE) # put the max intensity of that you want to graph
  #EEMmax[i] <- zmax #to show the maximum fluorescence for all files
  xlimit <- range(300, 700, finite=TRUE)
  ylimit <- range(240, 450, finite = TRUE)
  
  numcont = 100 # number of contour levels you want: Change if you want
  
  ##### contour plotting function
  library(gplots)
  
  # call contour plot function
  setwd("/Users/ashlee/SpecScripts") 
  source("EEM_contour_v1.R")
  
  #Plot contours and save in correction file
  plotpath <- file.path(delta.eems, paste(sample.ID,"_", project,".jpeg", sep = ""))
  
  g <- length(delta.temp)
  EEMplot <- data.frame(delta.temp) # not cutting out the last two columns

  #EEMplot[EEMplot<0]=0 # have line so that if fluorescence is less than 0, becomes 0.
  explot = seq(240, 800, by = 2) 
  emplot = as.numeric(row.names(EEMplot))
  
  jpeg(file=plotpath)
  contour.plots(eems = as.matrix(EEMplot), Title = paste(sample.ID, "Delta", sep=""), ex = explot, em = emplot)  
  dev.off()
  
  ########## Calculating spectral indicies
  # Abs indicies
  # Calculate delta abs indicies
  # pre EEM files
  setwd(prechlor.files)
  preabs.samplename <- toString(all[i,5])
  preabs.temp <- as.data.frame(read.delim(preabs.samplename, header= TRUE, sep = ",", stringsAsFactors=FALSE))
  
  #post EEM files
  setwd(postchlor.files)
  postabs.samplename <- toString(all[i,3])
  postabs.temp <- as.data.frame(read.delim(postabs.samplename, header= TRUE, sep = ",", stringsAsFactors=FALSE))
  
  ############# trim
  # Ensure that absorbance files go from 800-240 nm
  #trim so that exitation and emission goes from the same
  # call function
  setwd("/Users/ashlee/SpecScripts") 
  source("Abstrim_DBP_function.R")
  
  pre.trim <- abstrim(abs = preabs.temp, minex = "X240")
  post.trim <- abstrim(abs = postabs.temp, minex = "X240")
  
  ############# calculate pre - post
  delta.abs <- pre.trim - post.trim
  
   ##########
  # Calculating fluorescence indicies
  # call function
  setwd("/Users/ashlee/SpecScripts") 
  source("Aqualog_Fluorindicies_v2.R")
  
  #below may need to be altered depending on the output of your scan
  # ex wavelengths
  ex370 <- as.numeric(grep(370, colnames(delta.temp)))
  ex254 <- as.numeric(grep(254, colnames(delta.temp)))
  ex310 <- as.numeric(grep(310, colnames(delta.temp)))
  ex274 <- as.numeric(grep(274, colnames(delta.temp)))
  ex276 <- as.numeric(grep(276, colnames(delta.temp)))
  ex320 <- as.numeric(grep(320, colnames(delta.temp)))
  ex340 <- as.numeric(grep(340, colnames(delta.temp)))
  
  # em wavelengths
  em470 <- as.numeric(grep(470.394, rownames(delta.temp)))
  em520 <- as.numeric(grep(520.522, rownames(delta.temp))) 
  em435 <- as.numeric(grep(435.609, rownames(delta.temp))) 
  em480 <- as.numeric(grep(479.697, rownames(delta.temp)))
  em300 <- as.numeric(grep(300.484, rownames(delta.temp)))
  em345 <- as.numeric(grep(344.826, rownames(delta.temp)))
  em380 <- as.numeric(grep(380.302, rownames(delta.temp)))
  em420 <- as.numeric(grep(420.587, rownames(delta.temp)))
  em436 <- as.numeric(grep(436.766, rownames(delta.temp)))
  em350 <- as.numeric(grep(350., rownames(delta.temp)))
  em410 <- as.numeric(grep(410., rownames(delta.temp)))
  em430 <- as.numeric(grep(430., rownames(delta.temp)))
  
  # call function that calculates fluorescent indicies
  Fluor.indtemp <- Fluor(eem = delta.temp)
  
  # bind sample ID
  Fluor.indtemp <- cbind(sample.ID, Fluor.indtemp)
  Spec.ind[i,] <- Fluor.indtemp
  colnames(Spec.ind) <- colnames(Fluor.indtemp)
}

#### End of loop!
#write file containing spectral indicies + sample IDs
#after loop is finished with all samples
corrpath <- file.path(delta.eems, paste(project, "FluorIndicies.csv", sep = ""))
write.table(Spec.ind, file = corrpath, row.names = FALSE, col.names = TRUE, sep = ",")

######## Prepping files for Cory McKnight modelling in Matlab
########

setwd(delta.eems) 
filelist_EEMScor <- list.files(pattern = "_DBPdelta.csv$")

# CM - prepping for CM PARAFAC model
# Function does three things: trims EEMS according to specified min eexitation wavlenegth
# take out row and column names in first column and row and save in CM folder with _i as per graph headins file
# Also creates ex and em files, as well as graph headings file as txt file and saves in CM file
# Lastly, creates and saves graph headings file as txt file, which the function returns to double check

ex.PARAFAC <- seq(240, 800, by = 2) #change if excitation wavlenegths are different

# call function
setwd("/Users/ashlee/SpecScripts") 
source("EEMCMtrimDBPdelta_function.R")
CMsave <- CMtrim(filedirectory = delta.eems, filelist = filelist_EEMScor, project = project, exmin = "X240",
                 directoryCM = directoryCM, ex = ex.PARAFAC)

# check graph headings file returned by function
CMsave

####################################################################################################################################
######### DOM Fluor
########
# Get files ready for DOMFLuor toolbox. 
# Need .csv file for ex, em and one csv file containing all of the fluorescence EEMS compiled
# Use Save Dr EEMS function
# Inputs include the filelist, the project, the vector containing sample names, and the excitation wavelength min you want to trim to
# File cuts EEMs from ex min that you want to

setwd("/Users/ashlee/SpecScripts") 
source("EEMSDrEEMsaveDBPDelta_function.R")

ex.DrEEMS = seq(240, 800, by = 2)
DrEEM.data = DrEEM(filelist = filelist_EEMScor, project = project, 
                   exmin = 'X240', filedirectory = delta.eems, ex = ex.DrEEMS)

# check DrEEM.data. This is the compiled EEMS for DrEEM PARAFAC modelling
head(DrEEM.data)
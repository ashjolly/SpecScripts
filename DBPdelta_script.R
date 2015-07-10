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

project <- "DBPdelta"
# create file that matches the pre and post files together

# prechlor files
setwd(prechlor.files)
filelist_pre <- list.files(pattern = "CorrectedNEW.csv$")
# Get sample ID from pre chlor files
sample.ID <- 0 #create sample ID variable
n = length(filelist_pre)
for (i in 1:n){
  sampleID.pre.temp <- strapplyc(filelist_pre[i], "DBP(.*)_DBPPre_CorrectedNEW", simplify = TRUE)
  sample.ID[i] <- sampleID.pre.temp 
}

filelist_pre <- as.data.frame(cbind(filelist_pre, sample.ID))

# Get sample ID from post chlor files
#postchlor
setwd(postchlor.files)
filelist_post <- list.files(pattern = "CorrectedNEW.csv$")

n = length(filelist_post)
for (i in 1:n){
  sampleID.post.temp <- strapplyc(filelist_post[i], "DBPChlor(.*)_DBPPost_CorrectedNEW", simplify = TRUE)
  sample.ID[i] <- sampleID.post.temp 
}

filelist_post <- as.data.frame(cbind(filelist_post, sample.ID))

# merge files by sample ID
all <- merge(filelist_post, filelist_pre, by = "sample.ID", all = FALSE)
Fluor.ind = data.frame(matrix(vector(), 0, 8))

# load each file and do the subtraction

n = dim(all)[1]
for (i in 1:n){
  # pre files
  setwd(prechlor.files)
  pre.samplename <- toString(all[i,3])
  pre.temp <- as.data.frame(read.delim(pre.samplename, header= TRUE, sep = ",", stringsAsFactors=FALSE))
  
  #post files
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
  sample.ID <- strapplyc(pre.samplename, "(.*)_DBPPre_CorrectedNEW", simplify = TRUE)
  
  # save delta eem file
  deltapath <- file.path(delta.eems, paste(sample.ID,"_", project,".csv", sep = ""))
  write.table(delta.temp, file = deltapath, row.names = TRUE,col.names = TRUE, sep = ",")

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
  Fluor.ind[i,] <- Fluor.indtemp
  colnames(Fluor.ind) <- colnames(Fluor.indtemp)
  
  ##########
  # last - plot delta eems as a contour plot
  #variables to change
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
}

#### End of loop!
#write file containing spectral indicies + sample IDs
#after loop is finished with all samples
corrpath <- file.path(delta.eems, paste(project, "FluorIndicies.csv", sep = ""))
write.table(Fluor.ind, file = corrpath, row.names = FALSE, col.names = TRUE, sep = ",")


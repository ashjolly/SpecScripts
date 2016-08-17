# file for correcting ABS/EEMS file associated with soil lysimeters

# File for correcting EEMS from Aqualog
# Note that this script is very similat to 'EEMS_Correctionandindicies_DBPpost" previously file (overwritten) 
# with the exception that you now specify whether you want to load corrected EEM files ("PEM.dat") files, or raw EEMS without any corrections 
# done in the software ("SYM.dat")

# The code will then correct SYM according to this order:
# (Note that the CDOM community typically corrects in the following order):
# Instrument Correct the Raman file - done in software
# Instrument Correct and Raman Normalize the Blank
# Instrument Correct the Sample - done in software
# Inner filter Correct the Sample
# Raman Normalize the Sample
# Blank Subtract

# PEM.dat files have been corrected by AJ for blank, IFE and Raleigh (just needs Raman)
# there was some question in summer 2015 about the order of this correction versus that accepted by the CDOM community
# (Aqualog software = blank, IFE, Raleigh and then Raman)
# According to Adam Gimore from Horiba, as well as from studies done by Rose Cory, the fact that Aqualog does the blank correctio
# prior to others shouldn't matter, and doesn't affect the calculation of inidicies or modelling
# Also, they showed that doing the blank correction first is more accurate if the quality of the blank is called into question
# Specifically, that doing the blank subtraction after Raman, etc can introduce negatives if there are elements in the blank
# and that this can amplify artifacts present within the blank...

# OK, so this code DOES THIS:
# 1. Take blank corrected (sample - blank) SYM files from Aqualog
# 2. Raman correct them
# 3. IFE correct them
# 4. account for dilution factor. 
# 5. Run eemscat.m file from matlab IN THE FUTURE. for now it doesn't
# 6. Calculate Abs and Fluor indicies from corrected files
# Trim and save files in correct locations for CM and Dr EEMS modelling

# 22july2015
# fudge I hope this is the last time I dos this... :0

# TO DOS:
# 1. INTERFACE WITH script for interpolating Raleigh data in Matlab. Call matlab from the script in correction loop and use to correct
# Note that Raleigh corrections are done in Matlab with "eemscat.m" file from Rasmus Bro et al. 
# This file interpolates in Raleigh and Raman regions, rather than
# Code allows user to specify which type of EEM is used. Ultimatley it would be great to run this from here, but not Raleigh correcting prior to 
# calculating indicies shouldn't matter

#2. Fe corrections?
##############################################################################################################

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

####
# directory with all of the fluorescence files
# CR Soil Lysimeters
directoryall <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Lysimeter/CRLys_Fluorescence/CRLys_EEMsCompiled"
#___________________________________________________________________________________
# directory for corrected EEMS - Raleigh corrected as well
#CR Soil lysimeters
directoryCorrectedEEMS <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Lysimeter/CRLys_Fluorescence/CRLys_CorrectedEEMS"
directoryCorrectedRaleigh <- directoryCorrectedEEMS

# Project name
project <- "CRLys"

#___________________________________________________________________________________
#######
# directory for saving EEMS for CM PARAFAC in 'Correct EEMS" file in the CM PARAFAC folder
# This is the same for all projects
directoryCM <-"/Users/user/Documents/MATLAB/CorrectedEEMS" 

#___________________________________________________________________________________
#######
# general directory
# CR Lysimeters 
directorygeneral <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Lysimeter/CRLys_Fluorescence"

## Absorbance - corrected and interpolated abs for soil lysimeter studies
directoryabs <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Lysimeter/CRLys_Fluorescence/CRLys_CorrInterpAbs"

######
#dilution file
top = c("sample.ID", "dilutionfactor")
# dilution file - lysimeters

dilution <-as.data.frame(read.delim("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Lysimeter/CRLys_Fluorescence/CRLys_dil.csv", 
                                  sep=",", header = TRUE, col.names = top))
dil <- dilution
######## create a graph headings file
# call function to create a graph headings file from abs, EEM and blank file
#setwd("/Users/user/SpecScripts") 
#source("EEMfilecompCR_function.R")
# select whether you want to work with sample-blank files (EEMfiletype ="SYM.dat"), 
# or files processed in Aqualog software (EEMfiletype ="PEM.dat")
#data.3 <- EEMfilecomp(workdir= directoryall, dil = dilution, EEMfiletype = "_SYM.dat")
#data <- EEMfilecompCR(workdir= directoryall, dil = dilution, EEMfiletype = "PEM.dat")
workdir= directoryall
setwd(workdir)
#above directory contains all blank, Abs and EEM files for correction and calculation from Aqualog

#########
#Blank files
filelist_Blank <- list.files(pattern = "_BEM.dat$")

#create column with sample ID - extracted from blank filename
y = length(filelist_Blank)

sample.ID <- 0 #create sample ID variable

for (i in 1:y){
  sample.ID.temp <- strapplyc(filelist_Blank[i], "(.*)_BEM", simplify = TRUE)
  sample.ID[i] <- sample.ID.temp
}
filelist_Blank <- cbind(filelist_Blank, sample.ID)

###########
#Abs 
filelist_Abs <- list.files(pattern = "_ABS.dat$")
#create column with sample ID - extracted from ABS filename
y = length(filelist_Abs)
sample.ID <- 0 #create sample ID variable
for (i in 1:y){
  sample.ID.temp <- strapplyc(filelist_Abs[i], "(.*)_ABS", simplify = TRUE)
  sample.ID[i] <- sample.ID.temp
}
filelist_Abs <- cbind(filelist_Abs, sample.ID)

#########
# raw EEMS files - note that these are IFM and RM
#filelist_EEMS <- list.files(pattern = "PEM.dat$")

#below is raw eem without any corrections
filelist_EEMS <- list.files(pattern = "_SYM.dat")

#create column with sample ID - extracted from EEMS filename
y = length(filelist_EEMS)
sample.ID <- 0
for (i in 1:y){
  sample.ID.temp <- strapplyc(filelist_EEMS[i], "(.*)_SYM.dat", simplify = TRUE)
  sample.ID[i] <- sample.ID.temp
}

filelist_EEMS <- cbind(filelist_EEMS, sample.ID)

# Get the eems that have been IFE/RM corrected
filelist_EEMSpem <- list.files(pattern = "_PEM.dat")

#create column with sample ID - extracted from EEMS filename
y = length(filelist_EEMSpem)
sample.ID <- 0
for (i in 1:y){
  sample.ID.temp <- strapplyc(filelist_EEMSpem[i], "(.*)_PEM", simplify = TRUE)
  sample.ID[i] <- sample.ID.temp
}
filelist_EEMSpem <- cbind(filelist_EEMSpem, sample.ID)

#######
# Merge blank, EEM, abs and dilution files according to sample ID
data.3 <- Reduce(function(x,y) merge(x, y, by = "sample.ID", all = TRUE), list(filelist_EEMS, filelist_Abs, filelist_Blank, dil, filelist_EEMSpem))

write.csv(data.3, file = paste0(directoryall, "/CRlys_data3.csv"))

data.3$filelist_Blank[is.na(data.3$filelist_Blank)] <- "CRLys_21Nov2013_WestDeep_BEM.dat"

#########################################################################################################
############### Correct files
# from 
#setwd("/Users/user/SpecScripts") 
#source("EEMcorrectionCR_function.R")
#EEMcorrect = EEMcorrection(data.3 = data.2, directoryall = directoryall, directoryCorrectedEEMS = directoryCorrectedEEMS) 

# set working directory where all uncorrected EEMS, blank and absorbance files are
setwd(directoryall)

# Run corrections on all files contained within the data.3 vector for the files without any corrections
# insert empty variables for populating with ex an em vectors
# Complete cases for PEM files only
data.4 <- data.3[!is.na(data.3["filelist_EEMS"]),]
data.5 <- unique(data.4)
n = nrow(data.4) #number of files you are going to correct in the file. Double check ths prior to proceeding

ex_all = data.frame(matrix(vector(), 5000, n))
em_all = data.frame(matrix(vector(), 5000, n))  

#ex and em from blank - create dataframe
ex_blank = data.frame(matrix(vector(), 5000, n))
em_blank = data.frame(matrix(vector(), 5000, n)) 

# create a abs blank dataframe
ave.blank = data.frame(matrix(vector(), 1, n))

# 5000 is a placeholder for the number of rows. Note that this speeds the loop in comparison to leaving it as zero.

### set up loop to correct files and calculate indicies from all files in folder

for (i in 1:n){
  #### load and trim files: functions to load, trim and correct EEMS data - blank, EEM and files
  #sample name
  samplename <- toString(data.4[i,1]) #column with the sample IDs
  
  #### load and trim EEMS
  #call function. Note that wd will change to sample WD in function
  setwd("/Users/user/SpecScripts") 
  source("EEMfileLoadTrim_function.R")
  EEM <- EEMtrim(graphheadings = data.4, samplewd = directoryall, loopnum = i)
  
  #ex and em wavelengths
  ex = as.numeric((sort(colnames(EEM), decreasing = TRUE)))
  em = as.numeric((sort(rownames(EEM), decreasing = FALSE)))
  
  # interpolate to make sure that the wavelengths go from the same wavelengths
  # interpolate for all excitation to regular emission from 248-831 nm, 1 nm incrementes
  em.seq <- seq(as.numeric(format(round(as.numeric(min(rownames(EEM))), 0), nsmall = 0)), as.numeric(format(round(as.numeric(max(rownames(EEM))), 0), nsmall = 0)), by = 1) 
  EEM.interp <- as.data.frame(apply(EEM, 2, function(g) approx(x = as.numeric(rownames(EEM)), y = g, xout = em.seq)$y))
  row.names(EEM.interp) <- em.seq  
  #colnames(EEM.interp) <- seq(239, 800, by = 3)
  
  # if the excitation doesn't go in 2 nm increments, interpolate so that it does 
  ex.seq <- seq(240,800, by = 2)
  EEM.interp <- na.omit(EEM.interp)
  EEM.interp2 <- as.data.frame(t(apply(EEM.interp, 1, function(g) approx(x = as.numeric(colnames(EEM.interp)), y = g, xout = ex.seq)$y)))
  colnames(EEM.interp2) <- ex.seq
  
  #### load and trim ABS
  # call function. Note that wd will change to sample WD in function
  setwd("/Users/user/SpecScripts") 
  source("EEMAbsLoadTrim_function.R")
  Abstrim <- ABStrim(graphheadings = data.4, samplewd = directoryall, loopnum = i, column = 3)
  
  # interpolate to make sure that the wavelengths go from the same wavelengths (2 nm increments)
  # interpolate for all excitation to regular emission from 248-831 nm, 1 nm incrementes
  # if the excitation doesn't go in 2 nm increments, interpolate so that it does 
  if (as.numeric(colnames(Abstrim)[1]) - as.numeric(colnames(Abstrim)[2]) != -2){
    ex.seq <- seq(800, 240, by = -2)
    Abs.inter <- as.data.frame(t(approx(x = as.numeric(colnames(Abstrim)), y = Abstrim[1,], xout = ex.seq)$y))
    colnames(Abs.inter) <- ex.seq
  } 
  
  #### load and trim blank
  setwd("/Users/user/SpecScripts")
  source("EEMBlankLoadTrim_function.R")
  Blktrim <- BLANKtrim(graphheadings = data.4, samplewd = directoryall, loopnum = i)
  
  # interpolate for all excitation to regular emission from 248-831 nm, 1 nm incrementes
  em.seq <- seq(as.numeric(format(round(as.numeric(min(rownames(Blktrim))), 0), nsmall = 0)), as.numeric(format(round(as.numeric(max(rownames(Blktrim))), 0), nsmall = 0)), by = 1) 
  Blk.interp <- as.data.frame(apply(Blktrim, 2, function(g) approx(x = as.numeric(rownames(Blktrim)), y = g, xout = em.seq)$y))
  row.names(Blk.interp) <- em.seq  

  # if the excitation doesn't go in 2 nm increments, interpolate so that it does 
  ex.seq <- seq(240,800, by = 2)
  Blk.interp <- na.omit(Blk.interp)
  Blk.interp2 <- as.data.frame(t(apply(Blk.interp, 1, function(g) approx(x = as.numeric(colnames(Blk.interp)), y = g, xout = ex.seq)$y)))
  colnames(Blk.interp2) <- ex.seq
  
  #### identify dilution factor in master file
  # Dilution = column 5 in data.3
  dil = data.4[i,5]
  
  ################################## Corrections
  ########### IFE: Correct raw EEM for IFE if sample has not been corrected for this
  # Determine if needs IFE done if EEM samplename is SYM.dat
  # Otherwise, if sample name is PEM, IFE has been done in the software
  EEMsampletype <- strapplyc(as.character(data.4[i,2]), paste(samplename, "(.*).dat", sep = ""), simplify = TRUE)
  
  if (EEMsampletype == "_SYM") {
    # call function
    setwd("/Users/user/SpecScripts") 
    source("EEMIFECorr_function.R")
    ex = as.numeric(colnames(EEM.interp2))
    em = as.numeric(row.names(EEM.interp2))
    EEM.IFC <- innerfilter(eem = EEM.interp2, abs = Abs.inter, em = em, ex = ex)
    # note that IFE should be between 0.4 and  0.98 according to McKnight 2001 (doi: 10.4319/lo.2001.46.1.0038)
  }
  
  # if sample type = PEM, do NOT do IFE correction - already done in file
  if (EEMsampletype == "PEM"){
    EEM.IFC <- EEM
  }
  
  ########### Normalize IFE EEM and blank file according to area under Raman peak
  # Need to do for all EEMs exported from Aqualog
  # call function
  setwd("/Users/user/SpecScripts") 
  source("Ramancorrect_v1.R")
  
  # wavelengths for Raman corrections
  # tell R where em = 375 nm, em = 430 nm; ex = 350 nm
  em375 <- as.numeric(grep("375", rownames(Blk.interp2)))
  em430 <- as.numeric(grep("430", rownames(Blk.interp2)))
  ex350 <- as.numeric(grep("350", colnames(Blk.interp2)))
  
  # get the Raman correction file from the Raman function stored
  Raman.area <- Ramancor(blank = Blk.interp2, em375, em430, ex350, em = as.numeric(row.names(Blk.interp2))) 
  
  # normalize the EEM file for the area underneath the raman curve calculated above 
  # Raman normalize the raw EEM
  EEM.ram = EEM.IFC/Raman.area 
  
  # Raman normalize the blank file (important for blank correction)
  blankram = Blk.interp2/Raman.area
  
  ###########################
  ##### Apply dilution factor to EEM and to Abs file
  EEM.dil = EEM.ram*dil 
  
  # blank correct Abs file
  # reference: 
  # Green, S. A., & Blough, N. V. (1994). Optical absorption and fluorescence properties of chromophoric dissolved organic matter in natural waters. Limnology and Oceanography, 39(8), 1903–1916. http://doi.org/10.4319/lo.1994.39.8.1903
  # Helms, J. R., Stubbins, A., & Ritchie, J. D. (2008). Absorption spectral slopes and slope ratios as indicators of molecular weight, source, and photobleaching of chromophoric dissolved organic matter. Limnology and ….
  blank <- rowMeans(Abs.inter[,grep("800", colnames(Abs.inter)):grep("700", colnames(Abs.inter))]) #find average value of abs between 800-700 nm
  ave.blank[i,1] <- blank
  abs.blank <- Abs.inter - blank # subtract the blank from the initial spectra
  
  # apply dilution factor to corrected absorbance
  Absdil = abs.blank*dil #second column = absorbance readings. Assumes Beer's Law applies (c~abs)
  
  ##################################
  ########### Correct for Raleigh Masking
  # Done via interpolation software avalaible for Matlab script - 'eemscat.m'
  # Access on 23July2015 from http://www.models.life.ku.dk/EEM_correction
  # Also see DOI http://dx.doi.org/10.1016/j.chemolab.2015.01.017 explaination of why interpolating
  # through the Raleigh scatter is the best option in terms of PARAFAC modelling
  
  # R script emulates the interpolation script found here (interpolation according to the matlab file above)
  # Don't know why I can;t get that gawd damned script to work. Oh well ;)
  
  # this portion of the script runs the eemscat.m function to get the interpolated spectra. Last correction before saving
  
  if (EEMsampletype == "_SYM"){
    # note that this will gap fill the second order Raleigh scatter with na.spline function in zoo
    # Use function to interpolate or add zeros for 1st order Raleigh ('yes' versus "no" for R1)
    setwd("/Users/user/SpecScripts") 
    source("EEMRaleighSoillys_function.R")
    EEM.rm <- raleigh(eem = na.omit(EEM.ram), slitwidth1 = 20, slitwidth2 = 20, ramanwidth = 5, R1 = "no")
  }

  # if Raleigh has already been done in Aqualog software (inserted 0's, not the best option)
   if(EEMsampletype == "_PEM"){
   EEM.rm <- EEM.dil
   }
  
  # change the corrected EEM so that excitation goes from 200-800, rather than 800-200
  EEMcorr <- EEM.rm[,sort(names(EEM.rm), decreasing = FALSE)]
  
  ###########
  ##### Save the corrected EEM
  corrpath <- file.path(directoryCorrectedEEMS, paste(samplename,"_", project,"_Corrected",".csv", sep = ""))
  write.table(EEMcorr, file = corrpath, row.names = TRUE,col.names = TRUE, sep = ",")
  
  ###########
  ##### Save the corrected absorbance file (trimmed and dilution factor accounted for). In same folder as EEMS
  abscorrpath <- file.path(directoryabs, paste(samplename,"_", project,"_AbsCorrected",".csv", sep = ""))
  write.table(Absdil, file = abscorrpath, row.names = TRUE,col.names = TRUE, sep = ",")
  
  ##########
  # last - plot corrected eems as a contour plot
  #variables to change
  zmax = max(EEMcorr,na.rm=TRUE) # put the max intensity of that you want to graph
  #EEMmax[i] <- zmax #to show the maximum fluorescence for all files
  xlimit <- range(300, 700, finite=TRUE)
  ylimit <- range(240, 450, finite = TRUE)
  
  numcont = 100 # number of contour levels you want: Change if you want
  
  ##### contour plotting function
  # call contour plot function
  setwd("/Users/user/SpecScripts") 
  source("EEM_contour_v1.R")
  
  #Plot contours and save in correction file
  plotpath <- file.path(directoryCorrectedEEMS, paste(samplename,"_", project, "_Contour",".jpeg", sep = ""))
  
  g <- length(EEMcorr)
  EEMplot <- EEMcorr # not cutting out the last two columns
  EEMplot[EEMplot<0]=0 # have line so that if fluorescence is less than 0, becomes 0.
  
  explot = as.numeric(colnames(as.matrix(EEMplot)))
  emplot = as.numeric(row.names(EEMplot))
  
  jpeg(file=plotpath)
  contour.plots(eems = as.matrix(EEMplot), Title = paste(samplename, project, sep = "_"), ex = explot, em = emplot, 
                zmax = zmax, zmin = 0, numcont = numcont)  
  dev.off()
  # note that the above is meant to be a crude graphing - better graphing done in matlab once
  # you figure out the max emission for your dataset (normalize all of the plots to this)
}

################
# Correct the PEM files - just for the Raman correct. Only in the instance where there isn't a SYM version
# Get subset with the PEM files that you will use
data.4 <- data.3[38:51,]
data.4 <- unique(data.4[!is.na(data.4["filelist_EEMSpem"]),]) # get complete cases of PEM only 
data.4[,2] <- data.4$filelist_EEMSpem # make sure the second column has PEM file names to make sure it works with functions

n = nrow(data.4) #number of files you are going to correct in the file. Double check ths prior to proceeding

ex_all = data.frame(matrix(vector(), 5000, n)) # 5000 is a placeholder for the number of rows. Note that this speeds the loop in comparison to leaving it as zero.
em_all = data.frame(matrix(vector(), 5000, n))  
#ex and em from blank - create dataframe
ex_blank = data.frame(matrix(vector(), 5000, n))
em_blank = data.frame(matrix(vector(), 5000, n)) 

# create a abs blank dataframe
ave.blank = data.frame(matrix(vector(), 1, n))

### set up loop to correct files and calculate indicies from all files in folder
for (i in 1:n){
  #### load and trim files: functions to load, trim and correct EEMS data - blank, EEM and files
  #sample name
  samplename <- toString(data.4[i,1]) #column with the sample IDs
  
  #### load and trim EEMS
  #call function. Note that wd will change to sample WD in function
  setwd("/Users/user/SpecScripts") 
  source("EEMfileLoadTrim_function.R")
  EEM <- EEMtrim(graphheadings = data.4, samplewd = directoryall, loopnum = i)
  
  #ex and em wavelengths
  ex = as.numeric((sort(colnames(EEM), decreasing = TRUE)))
  em = as.numeric((sort(rownames(EEM), decreasing = FALSE)))
  
  # interpolate to make sure that the wavelengths go from the same wavelengths
  # interpolate for all excitation to regular emission from 248-831 nm, 1 nm incrementes
  em.seq <- seq(as.numeric(format(round(as.numeric(min(rownames(EEM))), 0), nsmall = 0)), as.numeric(format(round(as.numeric(max(rownames(EEM))), 0), nsmall = 0)), by = 1) 
  EEM.interp <- as.data.frame(apply(EEM, 2, function(g) approx(x = as.numeric(rownames(EEM)), y = g, xout = em.seq)$y))
  row.names(EEM.interp) <- em.seq  
  #colnames(EEM.interp) <- seq(239, 800, by = 3)
  
  # if the excitation doesn't go in 2 nm increments, interpolate so that it does 
  ex.seq <- seq(240,800, by = 2)
  EEM.interp <- na.omit(EEM.interp)
  EEM.interp2 <- as.data.frame(t(apply(EEM.interp, 1, function(g) approx(x = as.numeric(colnames(EEM.interp)), y = g, xout = ex.seq)$y)))
  colnames(EEM.interp2) <- ex.seq
  
  #### load and trim ABS
  # call function. Note that wd will change to sample WD in function
  setwd("/Users/user/SpecScripts") 
  source("EEMAbsLoadTrim_function.R")
  Abstrim <- ABStrim(graphheadings = data.4, samplewd = directoryall, loopnum = i, column = 3)
  
  # interpolate to make sure that the wavelengths go from the same wavelengths (2 nm increments)
  # interpolate for all excitation to regular emission from 248-831 nm, 1 nm incrementes
  # if the excitation doesn't go in 2 nm increments, interpolate so that it does 
  if (as.numeric(colnames(Abstrim)[1]) - as.numeric(colnames(Abstrim)[2]) != -2){
    ex.seq <- seq(800, 240, by = -2)
    Abs.inter <- as.data.frame(t(approx(x = as.numeric(colnames(Abstrim)), y = Abstrim[1,], xout = ex.seq)$y))
    colnames(Abs.inter) <- ex.seq
  } 
  
  #### load and trim blank
  setwd("/Users/user/SpecScripts")
  source("EEMBlankLoadTrim_function.R")
  Blktrim <- BLANKtrim(graphheadings = data.4, samplewd = directoryall, loopnum = i)
  
  # interpolate for all excitation to regular emission from 248-831 nm, 1 nm incrementes
  em.seq <- seq(as.numeric(format(round(as.numeric(min(rownames(Blktrim))), 0), nsmall = 0)), as.numeric(format(round(as.numeric(max(rownames(Blktrim))), 0), nsmall = 0)), by = 1) 
  Blk.interp <- as.data.frame(apply(Blktrim, 2, function(g) approx(x = as.numeric(rownames(Blktrim)), y = g, xout = em.seq)$y))
  row.names(Blk.interp) <- em.seq  
  
  # if the excitation doesn't go in 2 nm increments, interpolate so that it does 
  ex.seq <- seq(240,800, by = 2)
  Blk.interp <- na.omit(Blk.interp)
  Blk.interp2 <- as.data.frame(t(apply(Blk.interp, 1, function(g) approx(x = as.numeric(colnames(Blk.interp)), y = g, xout = ex.seq)$y)))
  colnames(Blk.interp2) <- ex.seq
  
  #### identify dilution factor in master file
  # Dilution = column 5 in data.3
  dil = data.4[i,5]
  
  ################################## Corrections
  ########### IFE: Correct raw EEM for IFE if sample has not been corrected for this
  # Determine if needs IFE done if EEM samplename is SYM.dat
  # Otherwise, if sample name is PEM, IFE has been done in the software
  EEMsampletype <- strapplyc(as.character(data.4[i,2]), paste(samplename, "(.*).dat", sep = ""), simplify = TRUE)
  
  if (EEMsampletype == "_SYM") {
    # call function
    setwd("/Users/user/SpecScripts") 
    source("EEMIFECorr_function.R")
    ex = as.numeric(colnames(EEM.interp2))
    em = as.numeric(row.names(EEM.interp2))
    EEM.IFC <- innerfilter(eem = EEM.interp2, abs = Abs.inter, em = em, ex = ex)
    # note that IFE should be between 0.4 and  0.98 according to McKnight 2001 (doi: 10.4319/lo.2001.46.1.0038)
  }
  
  # if sample type = PEM, do NOT do IFE correction - already done in file
  if (EEMsampletype == "_PEM"){
    EEM.IFC <- EEM.interp2
  }
  
  ########### Normalize IFE EEM and blank file according to area under Raman peak
  # Need to do for all EEMs exported from Aqualog
  # call function
  setwd("/Users/user/SpecScripts") 
  source("Ramancorrect_v1.R")
  
  # wavelengths for Raman corrections
  # tell R where em = 375 nm, em = 430 nm; ex = 350 nm
  em375 <- as.numeric(grep("375", rownames(Blk.interp2)))
  em430 <- as.numeric(grep("430", rownames(Blk.interp2)))
  ex350 <- as.numeric(grep("350", colnames(Blk.interp2)))
  
  # get the Raman correction file from the Raman function stored
  Raman.area <- Ramancor(blank = Blk.interp2, em375, em430, ex350, em = as.numeric(row.names(Blk.interp2))) 
  
  # normalize the EEM file for the area underneath the raman curve calculated above 
  # Raman normalize the raw EEM
  EEM.ram = EEM.IFC/Raman.area 
  
  # Raman normalize the blank file (important for blank correction)
  blankram = Blk.interp2/Raman.area
  
  ###########################
  ##### Apply dilution factor to EEM and to Abs file
  EEM.dil = EEM.ram*dil 
  
  # blank correct Abs file
  # reference: 
  # Green, S. A., & Blough, N. V. (1994). Optical absorption and fluorescence properties of chromophoric dissolved organic matter in natural waters. Limnology and Oceanography, 39(8), 1903–1916. http://doi.org/10.4319/lo.1994.39.8.1903
  # Helms, J. R., Stubbins, A., & Ritchie, J. D. (2008). Absorption spectral slopes and slope ratios as indicators of molecular weight, source, and photobleaching of chromophoric dissolved organic matter. Limnology and ….
  blank <- rowMeans(Abs.inter[,grep("800", colnames(Abs.inter)):grep("700", colnames(Abs.inter))]) #find average value of abs between 800-700 nm
  ave.blank[i,1] <- blank
  abs.blank <- Abs.inter - blank # subtract the blank from the initial spectra
  
  # apply dilution factor to corrected absorbance
  Absdil = abs.blank*dil #second column = absorbance readings. Assumes Beer's Law applies (c~abs)
  
  ##################################
  ########### Correct for Raleigh Masking
  # Done via interpolation software avalaible for Matlab script - 'eemscat.m'
  # Access on 23July2015 from http://www.models.life.ku.dk/EEM_correction
  # Also see DOI http://dx.doi.org/10.1016/j.chemolab.2015.01.017 explaination of why interpolating
  # through the Raleigh scatter is the best option in terms of PARAFAC modelling
  
  # R script emulates the interpolation script found here (interpolation according to the matlab file above)
  # Don't know why I can;t get that gawd damned script to work. Oh well ;)
  
  # this portion of the script runs the eemscat.m function to get the interpolated spectra. Last correction before saving
  
  if (EEMsampletype == "_SYM"){
    # note that this will gap fill the second order Raleigh scatter with na.spline function in zoo
    # Use function to interpolate or add zeros for 1st order Raleigh ('yes' versus "no" for R1)
    setwd("/Users/user/SpecScripts") 
    source("EEMRaleighSoillys_function.R")
    EEM.rm <- raleigh(eem = na.omit(EEM.ram), slitwidth1 = 20, slitwidth2 = 20, ramanwidth = 5, R1 = "no")
  }
  
  # if Raleigh has already been done in Aqualog software (inserted 0's, not the best option)
  if(EEMsampletype == "_PEM"){
    EEM.rm <- EEM.dil
  }
  
  # change the corrected EEM so that excitation goes from 200-800, rather than 800-200
  EEMcorr <- EEM.rm[,sort(names(EEM.rm), decreasing = FALSE)]
  
  ###########
  ##### Save the corrected EEM
  corrpath <- file.path(directoryCorrectedEEMS, paste(samplename,"_", project,"_Corrected",".csv", sep = ""))
  write.table(EEMcorr, file = corrpath, row.names = TRUE,col.names = TRUE, sep = ",")
  
  ###########
  ##### Save the corrected absorbance file (trimmed and dilution factor accounted for). In same folder as EEMS
  abscorrpath <- file.path(directoryabs, paste(samplename,"_", project,"_AbsCorrected",".csv", sep = ""))
  write.table(Absdil, file = abscorrpath, row.names = TRUE,col.names = TRUE, sep = ",")
  
  ##########
  # last - plot corrected eems as a contour plot
  #variables to change
  zmax = max(EEMcorr,na.rm=TRUE) # put the max intensity of that you want to graph
  #EEMmax[i] <- zmax #to show the maximum fluorescence for all files
  xlimit <- range(300, 700, finite=TRUE)
  ylimit <- range(240, 450, finite = TRUE)
  
  numcont = 100 # number of contour levels you want: Change if you want
  
  ##### contour plotting function
  # call contour plot function
  setwd("/Users/user/SpecScripts") 
  source("EEM_contour_v1.R")
  
  #Plot contours and save in correction file
  plotpath <- file.path(directoryCorrectedEEMS, paste(samplename,"_", project, "_Contour",".jpeg", sep = ""))
  
  g <- length(EEMcorr)
  EEMplot <- EEMcorr # not cutting out the last two columns
  EEMplot[EEMplot<0]=0 # have line so that if fluorescence is less than 0, becomes 0.
  
  explot = as.numeric(colnames(as.matrix(EEMplot)))
  emplot = as.numeric(row.names(EEMplot))
  
  jpeg(file=plotpath)
  contour.plots(eems = as.matrix(EEMplot), Title = paste(samplename, project, sep = "_"), ex = explot, em = emplot, 
                zmax = zmax, zmin = 0, numcont = numcont)  
  dev.off()
}

########################################################################
# Calculate indicies from corrected EEMs and abs
########################################################################

# create master file with abs and EEMs corrected file names aligned according to sample ID
setwd(directoryCorrectedEEMS)
filelist_EEMS <- list.files(pattern = "_Corrected.csv$")
#create column with sample ID - extracted from corrected EEMS filename

y = length(filelist_EEMS) # length should be 50 (50 files)

sample.ID <- 0 #create sample ID variable

for (i in 1:y){
  sample.ID.temp <- strapplyc(filelist_EEMS[i], paste("(.*)", "_", project, "_Corrected.csv", sep = ""), simplify = TRUE)
  sample.ID[i] <- sample.ID.temp
}
filelist <- cbind(filelist_EEMS , sample.ID)

###########
#Abs - corrected
setwd(directoryabs)
filelist_Abs_corr <- unique(list.files(pattern = "AbsCorrected.csv"))
#create column with sample ID - extracted from ABS filename
y = length(filelist_Abs_corr)

sample.ID <- 0 #create sample ID variable
for (i in 1:y){
  sample.ID.temp <- strapplyc(filelist_Abs_corr[i], paste("(.*)", "_", project, "_AbsCorrected.csv", sep = ""), simplify = TRUE)
  sample.ID[i] <- sample.ID.temp
}
filelist_Abs_corr <- cbind(filelist_Abs_corr, sample.ID)

#######
# Merge EEM and Abs filenames by sample ID to create file with all of the filenames
filelist_EEMScor <- merge(filelist, filelist_Abs_corr,  by = "sample.ID", all = TRUE)

######## Calculate Abs and Fluor indicies
n = dim(filelist_EEMScor)[1]

Spectral.Indicies = data.frame(matrix(vector(), n,35)) #creating an empty vector

for (i in 1:n){
  
  ###########
  # Calculating absorbance indicies
  # load the Abs file
  
  ###########
  # Calculating absorbance indicies
  # load the Abs file
  setwd(directoryabs)
  abs.temp.corr <-as.data.frame(read.delim(as.character(filelist_EEMScor[i,3]), 
                                           header= TRUE, sep = ",", stringsAsFactors=FALSE))
  
  # Calculate absoprtion coefficients from blank corrected data
  # References:
  # Helms, J. R., Stubbins, A., & Ritchie, J. D. (2008). Absorption spectral slopes and slope ratios as indicators of molecular weight, source, and photobleaching of chromophoric dissolved organic matter. Limnology and ….
  
  l <- 10 / 1000 #10 mm (1cm) path length expressed in m
  Naperian = 2.303 * abs.temp.corr/l # naperian
  Decadic <- abs.temp.corr/l # decadic absorbance. Note that this is what spetrolyzer gives. USE THIS for the 
  
  # call function to calculate absorbance indicies
  setwd("/Users/user/SpecScripts") 
  source("Aqualog_Absindicies_v1.R")
  
  #call the function to calculate indicies
  Abs.ind.dec <- Abs(absorbance = Decadic)
  Abs.ind.Nap <- Abs(absorbance = Naperian)
  
  ###################################
  # Fluorescence calculation
  setwd(directoryCorrectedEEMS)
  EEMcorr <-as.data.frame(read.delim(as.character(filelist_EEMScor[i,2]), 
                                     header= TRUE, sep = ",", stringsAsFactors=FALSE))
  
  # wavelengths for Fluorescence indicies calculation
  # ex wavelengths
  ex370 <- as.numeric(grep("370", colnames(EEMcorr)))
  ex254 <- as.numeric(grep("254", colnames(EEMcorr)))
  ex310 <- as.numeric(grep("310", colnames(EEMcorr)))
  ex274 <- as.numeric(grep("274", colnames(EEMcorr)))
  ex276 <- as.numeric(grep("276", colnames(EEMcorr)))
  ex320 <- as.numeric(grep("320", colnames(EEMcorr)))
  ex340 <- as.numeric(grep("340", colnames(EEMcorr)))
  ex260 <- as.numeric(grep("260", colnames(EEMcorr)))
  ex290 <- as.numeric(grep("290", colnames(EEMcorr)))
  ex240 <- as.numeric(grep("240", colnames(EEMcorr)))
  ex270 <- as.numeric(grep("270", colnames(EEMcorr)))
  ex300 <- as.numeric(grep("300", colnames(EEMcorr)))
  
  # em wavelengths
  em470 <- as.numeric(grep("470", rownames(EEMcorr)))
  em520 <- as.numeric(grep("520", rownames(EEMcorr)))
  em435 <- as.numeric(grep("435", rownames(EEMcorr))) 
  em480 <- as.numeric(grep("480", rownames(EEMcorr)))
  em300 <- as.numeric(grep("300", rownames(EEMcorr)))
  em345 <- as.numeric(grep("345", rownames(EEMcorr)))
  em380 <- as.numeric(grep("380", rownames(EEMcorr)))
  em420 <- as.numeric(grep("420", rownames(EEMcorr)))
  em436 <- as.numeric(grep("436", rownames(EEMcorr)))
  em350 <- as.numeric(grep("350", rownames(EEMcorr)))
  em410 <- as.numeric(grep("410", rownames(EEMcorr)))
  em430 <- as.numeric(grep("430", rownames(EEMcorr)))
  em320 <- as.numeric(grep("320", rownames(EEMcorr)))
  em326 <- as.numeric(grep("326", rownames(EEMcorr)))
  em430 <- as.numeric(grep("430", rownames(EEMcorr)))
  em400 <- as.numeric(grep("400", rownames(EEMcorr)))
  em450 <- as.numeric(grep("450", rownames(EEMcorr)))
  
  # call function
  setwd("/Users/user/SpecScripts") 
  source("Aqualog_Fluorindicies_v2.R")
  
  # call function that calculates fluorescent indicies
  Fluor.ind <- Fluor(eem = EEMcorr)
  
  ##########
  # bind fluor indicies with abs indicies as well as the sample id
  samplename <- as.character(filelist_EEMScor[i,1]) # column name where sample ID is 
  
  Spectral.Ind <- cbind(samplename, Abs.ind.dec, Abs.ind.Nap, Fluor.ind) 
  top <- colnames(Spectral.Ind)
  Spectral.Indicies[i,]  <- Spectral.Ind
  colnames(Spectral.Indicies) <- top
}
# save the resultant file
# write file after function is complete
corrpath <- file.path(directoryCorrectedEEMS, paste(project, "AbsFluor_Indicies.csv", sep = ""))
write.table(Spectral.Indicies, file = corrpath, row.names = FALSE, col.names = TRUE, sep = ",")

####################################################################################################################################
############################ Cutting EEMS for Cory McKnight and DOM Fluor toolbox
########
# Ensure that EEMS are all the same size + works for both the CM code as well as the DOMFluor toolbox to code together
# For DBP, this means ex = 240-800 in 2 nm incrmenets, noting that the  

######## Prepping files for Cory McKnight modelling in Matlab
########
# CM - prepping for CM PARAFAC model
# Function does three things: trims EEMS according to specified min eexitation wavlenegth
# take out row and column names in first column and row and save in CM folder with _i as per graph headins file
# Also creates ex and em files, as well as graph headings file as txt file and saves in CM file
# Lastly, creates and saves graph headings file as txt file, which the function returns to double check

ex.PARAFAC <- seq(240, 800, by = 2) #change if excitation wavlenegths are different

# create file list for the soil extract files
setwd(directoryCorrectedEEMS)
filelist_EEMS <- list.files(pattern = "_Corrected.csv$")

setwd("/Users/user/SpecScripts") 
source("AbsEEMSfilecompsoil_function.R")

filelist_EEMScor <- abseemfilecomp(directoryRaleigh = directoryInterpEEM, 
                                   projectname = project,  
                                   directorynoncorabs = directoryabs,
                                   filelist_EEMScor = filelist_EEMS)

# call function
setwd("/Users/user/SpecScripts") 
source("EEMCMtrim_function.R")
CMsave <- CMtrim(filedirectory = directoryCorrectedEEMS, filelist = filelist_EEMS, project = project, exmin = "X240",
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

setwd("/Users/user/SpecScripts") 
source("EEMSDrEEMsave_function.R")

ex.DrEEMS = seq(240, 800, by = 2)

DrEEM.data = DrEEM(filelist = filelist_EEMS, project = project, 
                   exmin = 'X240', filedirectory = directoryCorrectedEEMS, ex = ex.DrEEMS, 
                   DrEEMfolder = "/Users/user/Documents/MATLAB/toolbox/CorrEEMS/CRLys")

# check DrEEM.data. This is the compiled EEMS for DrEEM PARAFAC modelling
head(DrEEM.data)


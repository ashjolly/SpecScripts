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
# pre DBP
#directoryall <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_all"
# post DBP
#directoryall <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_all"
# Waterlogged samples - file with all of the EEMS
directoryall <- "/Users/user/Dropbox/PhD Work/PhD Data/WL_data/WL_Fluorescence/WL_filesfromAqualog _codecorrected"
#___________________________________________________________________________________
# directory for corrected EEMS and corrected Abs files (multiplied by dilution file)
# pre DBP
#directoryCorrectedEEMS <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_correctedEEMs"
# post DBP
#directoryCorrectedEEMS <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_correctedEEMS"
# Waterlogged samples
directoryCorrectedEEMS <- "/Users/user/Dropbox/PhD Work/PhD Data/WL_data/WL_Fluorescence/WL_CorrectedEEMS"

#___________________________________________________________________________________
# directory for corrected EEMS - Raleigh corrected as well
# pre DBP
#directoryCorrectedRaleigh <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_correctedEEMSRaleigh"
# post DBP
#directoryCorrectedRaleigh <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_correctedEEMSRaleigh"
# Waterlogged
directoryCorrectedRaleigh <- "/Users/user/Dropbox/PhD Work/PhD Data/WL_data/WL_Fluorescence/WL_correctedEEMRaleigh"

#___________________________________________________________________________________
#######
# directory for saving EEMS for CM PARAFAC in 'Correct EEMS" file in the CM PARAFAC folder
# This is the same for all projects
directoryCM <-"/Users/user/Documents/MATLAB/CorrectedEEMS" 

#___________________________________________________________________________________
#######
# general directory
# pre DBP
#directorygeneral <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_prechlorination"
# post DBP
#directorygeneral <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_postchlorination"
# Waterlogged
directorygeneral <- "/Users/user/Dropbox/PhD Work/PhD Data/WL_data/WL_Fluorescence"

#___________________________________________________________________________________
# project
#project = "DBPPre"
#project = "DBPPost"
project = "WL" #waterlogged code

######
#dilution file
top = c("sample.ID", "dilutionfactor")
#DBP pre dilution file
#dilution <-as.data.frame(read.csv(paste(directorygeneral, "/", project, "_Aqualogdilution.csv", sep = ""), sep=",", header = TRUE, col.names = top))
# DBP post 
#dilution <-as.data.frame(read.csv(paste(directorygeneral, "/", "DBP_postchlor_Aqualogdilution.csv", sep = ""),sep=",", header = TRUE, col.names = top)) 
# Waterlogged
dilution <-as.data.frame(read.csv(paste(directorygeneral, "/", "WL_dilution_factorsAQUALOG_v2.csv", sep = ""),sep=",", header = FALSE, col.names = top)) 

#___________________________________________________________________________________
###########################
## ex and em positions within your eems. This is so the Fluor and Raman correction script can find the right columns in order to calculate Fluorescences indicies
# below may need to be altered depending on the output of your scan - double check exact emission wavelengths
# wavlengths for Raman corrections
# tell R where em = 375 nm, em = 430 nm; ex = 350 nm
#DBP
#em.375 = 375.7
#em.430 = 429.8
#WL
em.375 = 375.426
em.430 = 430.117
ex.350 = 350

# wavelengths for Fluorescence indicies calculation
# ex wavelengths
ex.370 <- 370
ex.254 <- 254
ex.310 <- 310
ex.274 <- 274
ex.276 <- 276
ex.320 <- 320
ex.340 <- 340
ex.260 <- 260
ex.290 <- 290
ex.240 <- 240
ex.270 <- 270
ex.300 <- 300

# em wavelengths - DBP
#em.470 <- 470.394
#em.520 <- 520.522
#em.435 <- 435.609
#em.480 <- 479.697
#em.300 <- 300.484
#em.345 <- 344.826
#em.380 <- 380.302
#em.420 <- 420.587
#em.436 <- 436.766
#em.350 <- 350.534
#em.410 <- 410.205
#em.430 <- 429.828
#em.320 <- 320.908
#em.326 <- 326.594
#em.400 <- 400.99
#em.450 <- 450.663

# em wavelengths - WL
em.470 <- 470.104
em.520 <- 520.23
em.435 <- 435.32
em.480 <- 479.988
em.300 <- 300.201
em.345 <- 345.111
em.380 <- 380.015
em.420 <- 420.298
em.436 <- 435.898
em.350 <- 350.249
em.410 <- 409.917
em.430 <- 430.117
em.320 <- 320.056
em.326 <- 326.31
em.400 <- 400.127
em.450 <- 450.373

ex.wavelengths <- data.frame(rbind(ex.350, ex.370, ex.254, ex.310,ex.274,ex.276,ex.320,ex.340,
                                   ex.260, ex.290, ex.240, ex.270, ex.300))

em.wavelengths <- data.frame(rbind(em.375, em.470, em.520,em.435,em.480,em.300,em.345,
                                   em.380,em.420,em.436,em.350,em.410,
                                   em.430, em.320,em.326,em.400, em.450))

#ex.wavelengths <- t(read.csv('/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_exemfile/DBP_exwave.csv', header = TRUE))
#rownames(ex.wavelengths) <- ex.wavelengths[,1]
#em.wavelengths <- read.csv('/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_exemfile/DBP_emwave.csv', header = TRUE)
#rownames(em.wavelengths) <- em.wavelengths[,1]
###########
# call function to create a graph headings file from abs, EEM and blank file
setwd("/Users/user/SpecScripts") 
source("EEMfilecomp_function.R")

# select whether you want to work with sample-blank files (EEMfiletype ="SYM.dat"), 
# or files processed in Aqualog software (EEMfiletype ="PEM.dat")

data.3 <- EEMfilecomp(workdir= directoryall, dil = dilution, EEMfiletype = "SYM.dat")

# read all EEMS, abs and blank files into three files?

# Should not have to change anything below this!

##########################################################################################################################################################
# Run corrections on all files contained within the data.3 vector
# This function will basically run different functions for running corrections
# Basically will take all uncorrected EEMS within a folder, and apply functions for corrections
# Avoids having the loop within the script - idea is to make corrections more adaptable between different projects and easier to apply

setwd("/Users/user/SpecScripts") 
source("EEMcorrection_function.R")

EEMcorrect = EEMcorrection(data.3 = data.3, directoryall = directoryall, directoryCorrectedEEMS = directoryCorrectedEEMS, 
                          slitwidth1 = 15, slitwidth2 = 15,
                          ex.wavelengths = ex.wavelengths, em.wavelengths = em.wavelengths)

# will return vector with ex and em wavlenegths for samplesa s well as for blank. Check?
EEMcorrect

########################################################################################
######## Loop over corrected EEMs within the corrected EEM folder to correct for Raleigh
# Prior to calculating indicies, etc.
# Note that corrected for Raleigh occurs within CM and DrEEM PARAFAC models, so no need to do this for files that will be used in PARAFAC

setwd(directoryCorrectedEEMS)
filelist_EEMS <- list.files(pattern = "_Corrected.csv$")

EEM.IFE.rm = lapply(filelist_EEMS, function(x) read.delim(x, header = TRUE, sep = ","))

for (i in 1:length(filelist_EEMS)){
  #sample name
  samplename <- strapplyc(filelist_EEMS[i], paste("(.*)", "_", project, "_Corrected.csv", sep = ""), simplify = TRUE)
  
  temp.EEMS <- as.data.frame(EEM.IFE.rm[i])
  
  # trim so that eem starts at "X240"
  setwd("/Users/user/SpecScripts") 
  source("EEMtrimstart_function.R")
  temp.EEMS <- eemtrimstart(eem = temp.EEMS, start = "X240")
  
  # note that this will gap fill the second order Raleigh scatter with na.spline function in zoo
  # Use function to interpolate or add zeros for 1st order Raleigh ('yes' versus "no" for R1)
  setwd("/Users/user/SpecScripts") 
  source("EEMRaleigh_function.R")
  EEM.rm.temp <- raleigh(eem = temp.EEMS, slitwidth1 = 20, slitwidth2 = 20,ramanwidth = 5, R1 = "no")
  
  ##### Save the Raleigh corrected EEM
  corrpath <- file.path(directoryCorrectedRaleigh, paste(samplename,"_",project,"_Raleighcorr",".csv", sep = ""))
  write.table(EEM.rm.temp, file = corrpath, row.names = TRUE,col.names = TRUE, sep = ",")
  
  ##########
  # last - plot corrected eems as a contour plot
  #variables to change
  zmax = max(EEM.rm.temp,na.rm=TRUE) # put the max intensity of that you want to graph
  #EEMmax[i] <- zmax #to show the maximum fluorescence for all files
  xlimit <- range(300, 700, finite=TRUE)
  ylimit <- range(240, 450, finite = TRUE)
  
  numcont = 100 # number of contour levels you want: Change if you want
  
  ##### contour plotting function
  # call contour plot function
  setwd("/Users/user/SpecScripts") 
  source("EEM_contour_v1.R")
  
  #Plot contours and save in correction file
  plotpath <- file.path(directoryCorrectedRaleigh, paste(samplename,"_", project, "_ContourRaleigh",".jpeg", sep = ""))
  
  EEMplot <- EEM.rm.temp # not cutting out the last two columns
  EEMplot[EEMplot<0]=0 # have line so that if fluorescence is less than 0, becomes 0.
  
  explot = colnames(as.matrix(EEMplot))
  explot = as.numeric(gsub("X","",explot))
  emplot = as.numeric(row.names(EEMplot))
  
  jpeg(file=plotpath)
  contour.plots(eems = as.matrix(EEMplot), Title = samplename, ex = explot, em = emplot, 
                zmax = zmax, zmin = 0, numcont = numcont)  
  dev.off()
}

########################################################################################
######### Loop over corrected files in the to calculate indicies
# create master file with abs and EEMs corrected file names aligned according to sample ID
# call function to create master file  with filenames
setwd(directoryCorrectedRaleigh)
filelist_EEMS <- list.files(pattern = "_Raleighcorr.csv$")

setwd("/Users/user/SpecScripts") 
source("AbsEEMSfilecomp_function.R")

filelist_EEMScor <- abseemfilecomp(directoryAbsEEMs = directoryCorrectedEEMS, directoryRaleigh = directoryCorrectedRaleigh, 
                                   projectname = project, 
                                   filelist_EEMScor = filelist_EEMS)

# set directory with EEMS that you corrected according to the loop above

# Call function that will loop over the files in this folder and apply abs and fluor correction indicies
setwd("/Users/user/SpecScripts") 
source("EEMSIndCalculation_function.R")

spec.indicies = calc.indicies(filelist_EEMScor = filelist_EEMScor, directoryCorrectedAbs = directoryCorrectedEEMS,
                              directoryEEMs = directoryCorrectedRaleigh, ex.wavelengths, em.wavelengths)

# Check spectral indicies
spec.indicies

######## 
#write file containing spectral indicies + sample IDs
#after loop is finished with all samples
corrpath <- file.path(directoryCorrectedRaleigh, paste(project, "SpectralIndicies.csv", sep = ""))
write.table(spec.indicies, file = corrpath, row.names = FALSE, col.names = TRUE, sep = ",")

####################################################################################################################################
############################ Cutting EEMS for Cory McKnight and DOM Fluor toolbox
########
# Ensure that EEMS are all the same size + works for both the CM code as well as the DOMFluor toolbox to code together
# For DBP, this means ex = 240-800 in 2 nm incrmenets, noting that the  

# filelist of corrected EEMS - for both types of PARAFAC modelling

setwd(directoryCorrectedEEMS) 
filelist_EEMScor <- list.files(pattern = "_Corrected.csv$")

setwd("/Users/user/SpecScripts") 
source("AbsEEMSfilecomp_function.R")

filelist_EEMScor_Raleigh <- abseemfilecomp(directoryAbsEEMs = directoryCorrectedEEMS, directoryRaleigh = directoryCorrectedRaleigh, 
                                   projectname = project, 
                                   filelist_EEMScor = filelist_EEMS)

######## Prepping files for Cory McKnight modelling in Matlab
########
# CM - prepping for CM PARAFAC model
# Function does three things: trims EEMS according to specified min eexitation wavlenegth
# take out row and column names in first column and row and save in CM folder with _i as per graph headins file
# Also creates ex and em files, as well as graph headings file as txt file and saves in CM file
# Lastly, creates and saves graph headings file as txt file, which the function returns to double check

ex.PARAFAC <- seq(240, 800, by = 2) #change if excitation wavlenegths are different

# call function
setwd("/Users/user/SpecScripts") 
source("EEMCMtrim_function.R")
CMsave <- CMtrim(filedirectory = directoryCorrectedEEMS, filelist = filelist_EEMScor, project = project, exmin = "X240",
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

# Take out the green roof and rainwater harvest sampples. These will be super-imposed on the analysis later
# samples = 56-63, 73-75, 91-93, 119+13
GR <- c("DBP0056_DBPPre_Corrected.csv", "DBP0057_DBPPre_Corrected.csv", "DBP0058_DBPPre_Corrected.csv", 
        "DBP0060_DBPPre_Corrected.csv", "DBP0061_DBPPre_Corrected.csv", "DBP0062_DBPPre_Corrected.csv", 
        "DBP0073_DBPPre_Corrected.csv", "DBP0074_DBPPre_Corrected.csv", "DBP0075_DBPPre_Corrected.csv", 
        "DBP0091_DBPPre_Corrected.csv", "DBP0092_DBPPre_Corrected.csv", "DBP0093_DBPPre_Corrected.csv", 
        "DBP0119_DBPPre_Corrected.csv", "DBP0013_DBPPre_Corrected.csv")

# chlorinated
GRc <- c("DBPChlor0056_DBPPost_Corrected.csv", "DBPChlor0057_DBPPost_Corrected.csv", "DBPChlor0058_DBPPost_Corrected.csv", 
        "DBPChlor0060_DBPPost_Corrected.csv", "DBPChlor0061_DBPPost_Corrected.csv", "DBPChlor0062_DBPPost_Corrected.csv", 
        "DBPChlor0073_DBPPost_Corrected.csv", "DBPChlor0074_DBPPost_Corrected.csv", "DBPChlor0075_DBPPost_Corrected.csv", 
        "DBPChlor0091_DBPPost_Corrected.csv", "DBPChlor0092_DBPPost_Corrected.csv", "DBPChlor0093_DBPPost_Corrected.csv", 
        "DBPChlor0119_DBPPost_Corrected.csv", "DBPChlor0013_DBPPost_Corrected.csv")


# remove the green rrof samples from the set prior to PARAFAC modelling
#filelist_EEMScor <- filelist_EEMScor[ !(filelist_EEMScor  %in% GR) ] 
filelist_EEMScor <- filelist_EEMScor[!(filelist_EEMScor  %in% GRc)]

setwd("/Users/user/SpecScripts") 
source("EEMSDrEEMsave_function.R")

ex.DrEEMS = seq(240, 800, by = 2)

DrEEM.data = DrEEM(filelist = filelist_EEMScor, project = project, 
                  exmin = 'X240', filedirectory = directoryCorrectedEEMS, ex = ex.DrEEMS, 
                  DrEEMfolder = "/Users/user/Documents/MATLAB/toolbox/CorrEEMS/DBPPost_noGR")

# check DrEEM.data. This is the compiled EEMS for DrEEM PARAFAC modelling
head(DrEEM.data)

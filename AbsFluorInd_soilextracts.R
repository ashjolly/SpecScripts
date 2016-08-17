# Calculting EEM and abs indicies from corrected files
# Taken from master file - too many changes necessary to change the master file
# 11aug2016 Ashlee J's PhD
###########

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

# For soil extracts
directoryCorrectedRaleigh <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_soilextracts/CR_soilextracts_fluorescence/CorrectedEEMSRaleigh"
directoryCM <-"/Users/user/Documents/MATLAB/CorrectedEEMS" 
directoryabs <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_soilextracts/CR_soilextracts_fluorescence/CorrInterpAbs"
project <- "SoilExtracts"
directoryCorrectedEEMS <- directoryCorrectedRaleigh
directoryInterpEEM <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_soilextracts/CR_soilextracts_fluorescence/CorrectedEEM_Interp"

#############
setwd(directoryCorrectedRaleigh)
filelist_EEMS <- list.files(pattern = "_Raleighcorr.csv$")

# for soil extracts (3 nm increments for absorbance) - interpolate abs to 200-800, 2 nm increments
# Only need to do once!
#setwd("/Users/user/SpecScripts") 
#source("CRsoilAbsinterp_function.R")
#abs.interpret(filelist_EEMScor, directoryCorrectedEEMS = directoryCorrectedRaleigh, directoryCorrectedAbs = directoryabs)

# create file list for the soil extract files
setwd("/Users/user/SpecScripts") 
source("AbsEEMSfilecompsoil_function.R")

filelist_EEMScor <- abseemfilecomp(directoryRaleigh = directoryCorrectedRaleigh, 
                                   projectname = project,  
                                   directorynoncorabs = directoryabs,
                                   filelist_EEMScor = filelist_EEMS)


# Call function that will loop over the files in this folder and apply abs and fluor correction indicies

n = dim(filelist_EEMScor)[1]

Spectral.Indicies = data.frame(matrix(vector(), n,35)) #creating an empty vector

for (i in 1:n){
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
  
  ##########
  # Calculating fluorescence indicies
  setwd(directoryCorrectedRaleigh)
  EEMcorr <-as.data.frame(read.delim(as.character(filelist_EEMScor[i,2]), 
                                     header= TRUE, sep = ",", stringsAsFactors=FALSE))
 
   # interpolate for all excitation to regular emission from 248-831 nm, 1 nm incrementes
  em.seq <- seq(as.numeric(format(round(as.numeric(min(rownames(EEMcorr))), 0), nsmall = 0)), as.numeric(format(round(as.numeric(max(rownames(EEMcorr))), 0), nsmall = 0)), by = 1) 
  EEM.interp <- as.data.frame(apply(EEMcorr, 2, function(g) approx(x = as.numeric(rownames(EEMcorr)), y = g, xout = em.seq)$y))
  row.names(EEM.interp) <- em.seq  
  colnames(EEM.interp) <- seq(239, 800, by = 3)
  
  # if the excitation doesn't go in 2 nm increments, interpolate so that it does 
  ex.seq <- seq(240,800, by = 2)
  EEM.interp <- na.omit(EEM.interp)
  EEM.interp2 <- as.data.frame(t(apply(EEM.interp, 1, function(g) approx(x = as.numeric(colnames(EEM.interp)), y = g, xout = ex.seq)$y)))
  colnames(EEM.interp2) <- ex.seq
  
  # save interpolated EEM
  samplename <- filelist_EEMScor[i,1]
  corrpath <- file.path(directoryInterpEEM, paste(samplename,"_", project,"_CorrInterp",".csv", sep = ""))
  write.table(EEM.interp2, file = corrpath, row.names = TRUE,col.names = TRUE, sep = ",")
  
  # wavelengths for Fluorescence indicies calculation
  # ex wavelengths
  ex370 <- as.numeric(grep("370", colnames(EEM.interp2)))
  ex254 <- as.numeric(grep("254", colnames(EEM.interp2)))
  ex310 <- as.numeric(grep("310", colnames(EEM.interp2)))
  ex274 <- as.numeric(grep("274", colnames(EEM.interp2)))
  ex276 <- as.numeric(grep("276", colnames(EEM.interp2)))
  ex320 <- as.numeric(grep("320", colnames(EEM.interp2)))
  ex340 <- as.numeric(grep("340", colnames(EEM.interp2)))
  ex260 <- as.numeric(grep("260", colnames(EEM.interp2)))
  ex290 <- as.numeric(grep("290", colnames(EEM.interp2)))
  ex240 <- as.numeric(grep("240", colnames(EEM.interp2)))
  ex270 <- as.numeric(grep("270", colnames(EEM.interp2)))
  ex300 <- as.numeric(grep("300", colnames(EEM.interp2)))
  
  # em wavelengths
  em470 <- as.numeric(grep("470", rownames(EEM.interp2)))
  em520 <- as.numeric(grep("520", rownames(EEM.interp2))) 
  em435 <- as.numeric(grep("435", rownames(EEM.interp2))) 
  em480 <- as.numeric(grep("480", rownames(EEM.interp2)))
  em300 <- as.numeric(grep("300", rownames(EEM.interp2)))
  em345 <- as.numeric(grep("345", rownames(EEM.interp2)))
  em380 <- as.numeric(grep("380", rownames(EEM.interp2)))
  em420 <- as.numeric(grep("420", rownames(EEM.interp2)))
  em436 <- as.numeric(grep("436", rownames(EEM.interp2)))
  em350 <- as.numeric(grep("350", rownames(EEM.interp2)))
  em410 <- as.numeric(grep("410", rownames(EEM.interp2)))
  em430 <- as.numeric(grep("430", rownames(EEM.interp2)))
  em320 <- as.numeric(grep("320", rownames(EEM.interp2)))
  em326 <- as.numeric(grep("326", rownames(EEM.interp2)))
  em430 <- as.numeric(grep("430", rownames(EEM.interp2)))
  em400 <- as.numeric(grep("400", rownames(EEM.interp2)))
  em450 <- as.numeric(grep("450", rownames(EEM.interp2)))
  
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

# write file after function is complete
corrpath <- file.path(directoryInterpEEM, paste(project, "Fluorescence_Indicies.csv", sep = ""))
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
setwd(directoryInterpEEM)
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
CMsave <- CMtrim(filedirectory = directoryInterpEEM, filelist = filelist_EEMS, project = project, exmin = "X240",
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
                   exmin = 'X240', filedirectory = directoryInterpEEM, ex = ex.DrEEMS, 
                   DrEEMfolder = "/Users/user/Documents/MATLAB/toolbox/CorrEEMS/CRSoilExtracts")

# check DrEEM.data. This is the compiled EEMS for DrEEM PARAFAC modelling
head(DrEEM.data)

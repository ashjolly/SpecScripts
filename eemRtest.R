# Script for double checking that indicies and corrections from your script and other EEM r packages is the same
# Same with indicies
# Created 24 NOv20145
# Ashlee's PhD
# Worse timing ever
#####################

# clean up list
rm(list = ls())
ls()

# Load packages
library("EEM")
devtools::install_github("PMassicotte/eemR")
library(eemR)

##### Test 1 - are corrections the same between your scripts and theirs?
# Plan - correct both ways, and minus spectra from each other.
# Difference should be minimal.

# file with EEMS
directoryall <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_all"

# general directory
directorygeneral <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_postchlorination"

# project
#project = "DBPPre"
project = "DBPPost"

#dilution file
top = c("sample.ID", "dilutionfactor")
#DBP post dilution file
#dilution <-as.data.frame(read.csv(paste(directorygeneral, "/", project, "_Aqualogdilution.csv", sep = ""),
#                                  sep=",", header = TRUE, col.names = top))

dilution <-as.data.frame(read.csv(paste(directorygeneral, "/", "DBP_postchlor_Aqualogdilution.csv", sep = ""),
                                  sep=",", header = TRUE, col.names = top)) #post chlorination

###########################
## ex an em positions within your eems. This is so the Fluor and Raman correction script can find the right columns in order to calculate Fluorescences indicies
# below may need to be altered depending on the output of your scan - double check exact emission wavelengths
# wavlengths for Raman corrections
# tell R where em = 375 nm, em = 430 nm; ex = 350 nm
em.375 = 375.7
em.430 = 429.8
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

# em wavelengths
em.470 <- 470.394
em.520 <- 520.522
em.435 <- 435.609
em.480 <- 479.697
em.300 <- 300.484
em.345 <- 344.826
em.380 <- 380.302
em.420 <- 420.587
em.436 <- 436.766
em.350 <- 350
em.410 <- 410
em.430 <- 430

ex.wavelengths <- data.frame(rbind(ex.350, ex.370, ex.254, ex.310,ex.274,ex.276,ex.320,ex.340))
em.wavelengths <- data.frame(rbind(em.375, em.430, em.470, em.520,em.435,em.480,em.300,em.345,
                                   em.380,em.420,em.436,em.350,em.410))

###########
# call function to create a graph headings file from abs, EEM and blank file
setwd("/Users/user/SpecScripts") 
source("EEMfilecomp_function.R")

# select whether you want to work with sample-blank files (EEMfiletype ="SYM.dat"), 
# or files processed in Aqualog software (EEMfiletype ="PEM.dat")

data.3 <- EEMfilecomp(workdir= directoryall, dil = dilution, EEMfiletype = "SYM.dat")

####### Load one sample file - EEM, abs, dilution, blank
setwd("/Users/user/SpecScripts") 
source("EEMfileLoadTrim_function.R")
EEM <- EEMtrim(graphheadings = data.3, samplewd = directoryall, loopnum = 2)
setwd("/Users/user/SpecScripts") 
source("EEMAbsLoadTrim_function.R")
Abs.AJ <- ABStrim(graphheadings = data.3, samplewd = directoryall, loopnum = 2)
setwd("/Users/user/SpecScripts") 
source("EEMBlankLoadTrim_function.R")
Blk.AJ <- BLANKtrim(graphheadings = data.3, samplewd = directoryall, loopnum = 2)

ex = as.numeric((sort(colnames(EEM), decreasing = TRUE)))
em = as.numeric((sort(rownames(EEM), decreasing = FALSE)))
#############################################################
# Load using eemR - blank, EEM and abs
file <- system.file(paste(directoryall, "/27feb2015001DBPChlor0002SYM.dat", sep = ""), package = "eemR")
eem <- eem_read(paste(directoryall, "/27feb2015001DBPChlor0002SYM.dat", sep = ""))
plot(eem)
abs <- eem_read(paste(directoryall, "/27feb2015001DBPChlor0002ABS.dat", sep = ""))
blank <- eem_read(paste(directoryall, "/27feb2015001DBPChlor0002BEM.dat", sep = ""))
##################### ##################### ##################### AJ corrections
# perform IFE - your correction
setwd("/Users/user/SpecScripts") 
source("EEMIFECorr_function.R")
EEM.IFC <- innerfilter(eem = EEM, abs = Abs.AJ, em = em, ex = ex)

# perform IFE - eemR
EEM.IFC.eemr <- eem_inner_filter_effect(eem,
                                 absorbance = abs,
                                 pathlength = 1) ## 1 cm fluo pathlenght)


########### Normalize IFE EEM and blank file according to area under Raman peak
# Need to do for all EEMs exported from Aqualog
# call function
setwd("/Users/user/SpecScripts") 
source("Ramancorrect_v1.R")

# wavelengths for Raman corrections
# tell R where em = 375 nm, em = 430 nm; ex = 350 nm
em375 <- as.numeric(grep(em.wavelengths[c("em.375"),], rownames(EEM)))
em430 <- as.numeric(grep(em.wavelengths[c("em.430"),], rownames(EEM)))
ex350 <- as.numeric(match(ex.wavelengths[c("ex.350"),], colnames(EEM)))

# get the Raman correction file from the Raman function stored
Raman.area <- Ramancor(blank = Blk.AJ) 

# normalize the EEM file for the area underneath the raman curve calculated above 
# Raman normalize the raw EEM
EEM.ram = EEM/Raman.area 

# Raman normalize the blank file (important for blank correction)
blankram = Blktrim/Raman.area

########### Try raman normalization via eemR package
res <- eem_raman_normalisation(eem, blank)

############### Perform Raleigh Masking
setwd("/Users/user/SpecScripts") 
source("EEMRaleigh_function.R")
EEM.rm <- raleigh(eem = EEM, slitwidth1 = 20, slitwidth2 = 20, ramanwidth = 10, R1 = 'no')  

# perform corrections on the EEMS using the eemR package

EEM.rm.eemr <- 
  
# take the difference between the two to see how the diferent appraches changed the 

##### Test 2 - check that indicies are the same 
# Calculated indicies from your own
# calculate indicies from eemR package

# Calculate the difference between the indicies
# correct EEM - your way
# IFE
setwd("/Users/user/SpecScripts") 
source("EEMIFECorr_function.R")
EEM.IFC <- innerfilter(eem = EEM, abs = Abs.AJ, em = em, ex = ex)
# Raman
setwd("/Users/user/SpecScripts") 
source("Ramancorrect_v1.R")
Raman.area <- Ramancor(blank = Blk.AJ) 

EEM.rm <- EEM.IFC/Raman.area

# normalize the EEM file for the area underneath the raman curve calculated above 
# Raman normalize the raw EEM
EEM.ram = EEM.rm/Raman.area 
#Raleigh 
setwd("/Users/user/SpecScripts") 
source("EEMRaleigh_function.R")
EEM.rm <- raleigh(eem = EEM.ram, slitwidth1 = 20, slitwidth2 = 20, ramanwidth = 10, R1 = 'no')  

### calculate indicies - AJ
# abs indicies
setwd("/Users/user/SpecScripts") 
source("Aqualog_Absindicies_v1.R")
#call the function to calculate indicies
Abs.ind <- Abs(absorbance = Abs.AJ)
# fluor indicies
# wavelengths for Fluorescence indicies calculation
# ex wavelengths
ex370 <- as.numeric(grep(ex.wavelengths[c("ex.370"),], colnames(EEM.rm)))
ex254 <- as.numeric(grep(ex.wavelengths[c("ex.254"),], colnames(EEM.rm)))
ex310 <- as.numeric(grep(ex.wavelengths[c("ex.310"),], colnames(EEM.rm)))
ex274 <- as.numeric(grep(ex.wavelengths[c("ex.274"),], colnames(EEM.rm)))
ex276 <- as.numeric(grep(ex.wavelengths[c("ex.276"),], colnames(EEM.rm)))
ex320 <- as.numeric(grep(ex.wavelengths[c("ex.320"),], colnames(EEM.rm)))
ex340 <- as.numeric(grep(ex.wavelengths[c("ex.340"),], colnames(EEM.rm)))

# em wavelengths
em470 <- as.numeric(grep(em.wavelengths[c("em.470"),], rownames(EEM.rm)))
em520 <- as.numeric(grep(em.wavelengths[c("em.520"),], rownames(EEM.rm))) 
em435 <- as.numeric(grep(em.wavelengths[c("em.435"),], rownames(EEM.rm))) 
em480 <- as.numeric(grep(em.wavelengths[c("em.480"),], rownames(EEM.rm)))
em300 <- as.numeric(grep(em.wavelengths[c("em.300"),], rownames(EEM.rm)))
em345 <- as.numeric(grep(em.wavelengths[c("em.345"),], rownames(EEM.rm)))
em380 <- as.numeric(grep(em.wavelengths[c("em.380"),], rownames(EEM.rm)))
em420 <- as.numeric(grep(em.wavelengths[c("em.420"),], rownames(EEM.rm)))
em436 <- as.numeric(grep(em.wavelengths[c("em.436"),], rownames(EEM.rm)))
em350 <- as.numeric(grep(em.wavelengths[c("em.350"),], rownames(EEM.rm)))
em410 <- as.numeric(grep(em.wavelengths[c("em.410"),], rownames(EEM.rm)))
em430 <- as.numeric(grep(em.wavelengths[c("em.430"),], rownames(EEM.rm)))

# call function
setwd("/Users/user/SpecScripts") 
source("Aqualog_Fluorindicies_v2.R")

# call function that calculates fluorescent indicies
Fluor.ind <- Fluor(eem = EEM.rm)

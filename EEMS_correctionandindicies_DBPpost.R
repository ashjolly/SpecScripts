# File for correcting EEMS from DBP post chlorination
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

######################

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

####
# DBP Post chlorination
# directory with all of the fluorescence files
directoryall <- "/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_all"

# directory for corrected EEMS and corrected Abs files (multiplied by dilution file)
directoryCorrectedEEMS <- "/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_correctedEEMs"

#######
# directory for saving EEMS for CM PARAFAC in 'Correct EEMS" file in the CM PARAFAC folder
# This is the same for all projects
directoryCM <-"/Users/ashlee/Documents/MATLAB/CorrectedEEMS" 

######
#dilution file
top = c("sample.ID", "dilutionfactor")
#DBP post dilution file
dilution <-as.data.frame(read.csv("/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_Aqualogdilution.csv", 
                                  sep=",", header = TRUE, col.names = top))

project = "DBPPost"

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

waves <- data.frame(rbind(ex.370, ex.254, ex.310,ex.274,ex.276,ex.320,ex.340,
                          em.470,em.520,em.435,em.480,em.300,em.345,em.380,em.420,em.436,em.350,em.410,em.430))

###########
# call function to create a graph headings file from abs, EEM and blank file
setwd("/Users/ashlee/SpecScripts") 
source("EEMfilecomp_function.R")

# select whether you want to work with sample-blank files (EEMfiletype ="SYM.dat"), 
# or files processed in Aqualog software (EEMfiletype ="PEM.dat")

data.3 <- EEMfilecomp(workdir= directoryall, dil = dilution, EEMfiletype = "SYM.dat")

# Should not have to change anything below this!
##########################################################################################################################################################
##########################################################################################################################################################
# Run corrections on all files contained within the data.3 vector
# This function will basically run different functions for running corrections
# Basically will take all uncorrected EEMS within a folder, and apply functions for corrections
# Avoids having the loop within the script - idea is to make corrections more adaptable between different projects and easier to apply

#setwd("/Users/ashlee/SpecScripts") 
#source("EEMcorrection_function.R")
#EEMcorrect = EEMcorrection(data.3 = data.3, directoryall = directoryall, directoryCorrectedEEMS = directoryCorrectedEEMS, 
 #                          slitwidth1 = 15, slitwidth2 = 15,
 #                          em.375 = em.375, em.430 = em.430)

# will return vector with ex and em wavlenegths for samplesa s well as for blank. Check?
#EEMcorrect
slitwidth1 = 15
slitwidth2 = 15

# set working directory where all uncorrected EEMS, blank and absorbance files are
setwd(directoryall)

# Run corrections on all files contained within the data.3 vector
# insert empty variables for populating with ex an em vectors
n = nrow(data.3) #number of files you are going to correct in the file. Double check ths prior to proceeding

ex_all = data.frame(matrix(vector(), 5000, n))
em_all = data.frame(matrix(vector(), 5000, n))  

#ex and em from blank - create dataframe
ex_blank = data.frame(matrix(vector(), 5000, n))
em_blank = data.frame(matrix(vector(), 5000, n)) 

# 5000 is a placeholder for the number of rows. Note that this speeds the loop in comparison to leaving it as zero.

### set up loop to correct files and calculate indicies from all files in folder

for (i in 1:n){
  
  #### load and trim files: functions to load, trim and correct EEMS data - blank, EEM and files
  #sample name
  samplename <- toString(data.3[i,1]) #column with the sample IDs
  
  #### load and trim EEMS
  #call function. Note that wd will change to sample WD in function
  setwd("/Users/ashlee/SpecScripts") 
  source("EEMfileLoadTrim_function.R")
  EEM <- EEMtrim(graphheadings = data.3, samplewd = directoryall, loopnum = i)
  
  #ex and em wavelengths
  ex = as.numeric((sort(colnames(EEM), decreasing = TRUE)))
  em = as.numeric((sort(rownames(EEM), decreasing = FALSE)))
  
  #em_all is a variable that holds all of the em wavelengths. It's a way of seeing if em wavelengths are different between samples
  em_all = rbind(em_all, rownames(EEM))
  #ex_all is a variable that holds all of the ex wavelengths. It's a way of seeing if ex wavelengths are different between samples
  ex_all= rbind(ex_all, colnames(EEM))
  
  #### load and trim ABS
  # call function. Note that wd will change to sample WD in function
  setwd("/Users/ashlee/SpecScripts") 
  source("EEMAbsLoadTrim_function.R")
  Abstrim <- ABStrim(graphheadings = data.3, samplewd = directoryall, loopnum = i)
  
  #### load and trim blank
  setwd("/Users/ashlee/SpecScripts")
  source("EEMBlankLoadTrim_function.R")
  Blktrim <- BLANKtrim(graphheadings = data.3, samplewd = directoryall, loopnum = i)
  
  #ex and em_all are variables that output the complete ex and em. Make sure that all of your samples are the same
  ex_blank = rbind(ex_blank, colnames(Blktrim))
  em_blank = rbind(em_blank, rownames(Blktrim))
  
  #### identify dilution factor in master file
  # Dilution = column 5 in data.3
  dil = data.3[i,5]
  
  ################################## Corrections
  ########### IFE: Correct raw EEM for IFE if sample has not been corrected for this
  # Determine if needs IFE done if EEM samplename is SYM.dat
  # Otherwise, if sample name is PEM, IFE has been done in the software
  EEMsampletype <- strapplyc(as.character(data.3[i,2]), paste(samplename, "(.*).dat", sep = ""), simplify = TRUE)
  
  if (EEMsampletype == "SYM") {
    # call function
    setwd("/Users/ashlee/SpecScripts") 
    source("EEMIFECorr_function.R")
    
    EEM.IFC <- innerfilter(eem = EEM, abs = Abstrim, em = em, ex = ex)
    # note that IFE should be between 0.4 and  0.98 according to McKnight 2001 (doi: 10.4319/lo.2001.46.1.0038)
  }
  
  # if sample type = PEM, do NOT do IFE correction - already done in file
  if (EEMsampletype == "PEM"){
    EEM.IFC <- EEM
  }
  
  ########### Normalize IFE EEM and blank file according to area under Raman peak
  # Need to do for all EEMs exported from Aqualog
  # call function
  setwd("/Users/ashlee/SpecScripts") 
  source("Ramancorrect_v1.R")
  
  # wavlengths for Raman corrections
  # tell R where em = 375 nm, em = 430 nm; ex = 350 nm
  em375 <-  as.numeric(grep(em.375, rownames(EEM)))
  em430 <-  as.numeric(grep(em.430, rownames(EEM)))
  ex350 <- as.numeric(match(350, colnames(EEM)))
  
  # get the Raman correction file from the Raman function stored
  Raman.area <- Ramancor(blank = Blktrim) 
  
  # normalize the EEM file for the area underneath the raman curve calculated above 
  # Raman normalize the raw EEM
  EEM.ram = EEM.IFC/Raman.area 
  
  # Raman normalize the blank file (important for blank correction)
  blankram = Blktrim/Raman.area
  
  ##################################
  ########### Correct for Blank
  # Only if have not been done in software
  # corrected EEM = Sample - blank. Both SYM and PEM have been instrument corrected
  
  #if (EEMsampletype == "") {
  #  EEM.blk <- EEM.ram - blankram
  #}
  
  # if sample type = SYM and PEM, do NOT do blank correction
  #elseif {
  #  EEM.blk <- EEM.ram
  #}
  ##### note that SYM + PEM files from Aqualog have already been blank corrected
  
  ###########################
  ##### Apply dilution factor to EEM and to Abs file
  EEM.dil = EEM.ram*dil 
  Absdil = Abstrim*dil #second column = absorbance readings. Assumes Beer's Law applies (c~abs)
  
  ##################################
  ########### Correct for Raleigh Masking
  # Done via interpolation software avalaible for Matlab script - 'eemscat.m'
  # Access on 23July2015 from http://www.models.life.ku.dk/EEM_correction
  # Also see DOI http://dx.doi.org/10.1016/j.chemolab.2015.01.017 explaination of why interpolating
  # through the Raleigh scatter is the best option in terms of PARAFAC modelling
  
  # R script emulates the interpolation script found here (interpolation according to the matlab file above)
  # Don't know why I can;t get that gawd damned script to work. Oh well ;)
  
  # this portion of the script runs the eemscat.m function to get the interpolated spectra. Last correction before saving
  
  if (EEMsampletype == "SYM"){
    
    # call function
    setwd("/Users/ashlee/SpecScripts") 
    source("EEMRaleigh_function.R")
    # note that this will gap fill the second order Raleigh scatter with na.spline function in zoo
    EEM.rm <- raleigh(eem = EEM.dil, slitwidth1, slitwidth2)
  }
  
  # if Raleigh has already been done in Aqualog software (inserted 0's, not the best option)
  if(EEMsampletype == "PEM"){
    EEM.rm <- EEM.dil
  }
  
  ##### Apply correction factor for Fe concentration
  ##### TO DO!!!!!
  
  # change the corrected EEM so that excitation goes from 200-800, rather than 800-200
  EEMcorr <- EEM.rm[,sort(names(EEM.rm), decreasing = FALSE)]
  
  ###########
  ##### Save the corrected EEM
  corrpath <- file.path(directoryCorrectedEEMS, paste(samplename,"_", project,"_Corrected",".csv", sep = ""))
  write.table(EEMcorr, file = corrpath, row.names = TRUE,col.names = TRUE, sep = ",")
  
  ###########
  ##### Save the corrected absorbance file (trimmed and dilution factor accounted for). In same folder as EEMS
  abscorrpath <- file.path(directoryCorrectedEEMS, paste(samplename,"_", project,"_AbsCorrected",".csv", sep = ""))
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
  setwd("/Users/ashlee/SpecScripts") 
  source("EEM_contour_v1.R")
  
  #Plot contours and save in correction file
  plotpath <- file.path(directoryCorrectedEEMS, paste(samplename,"_", project, "_Contour",".jpeg", sep = ""))
  
  g <- length(EEMcorr)
  EEMplot <- EEMcorr # not cutting out the last two columns
  EEMplot[EEMplot<0]=0 # have line so that if fluorescence is less than 0, becomes 0.
  
  explot = as.numeric(colnames(EEMplot)) 
  emplot = as.numeric(row.names(EEMplot))
  
  jpeg(file=plotpath)
  contour.plots(eems = as.matrix(EEMplot), Title = samplename, ex = explot, em = emplot)  
  dev.off()
}

#### End of corrections loop!

########################################################################################
######### Loop over corrected files in the to calculate indicies
Spectral.Indicies = data.frame(matrix(vector(), 5000, 17)) #creating an empty vector

# create master file with abs and EEMs corrected file names aligned according to sample ID
# call function to create master file  with filenames
setwd("/Users/ashlee/SpecScripts") 
source("AbsEEMSfilecomp_function.R")

filelist_EEMScor <- abseemfilecomp(directoryAbsEEMs = directoryCorrectedEEMS, projectname = project)

# set directory with EEMS that you corrected according to the loop above
# Call function that will loop over the files in this folder and apply abs and fluor correction indicies

setwd("/Users/ashlee/SpecScripts") 
source("EEMSIndCalculation_function.R")
Spectral.Indicies = calc.indicies(filelist_EEMScor = filelist_EEMScor, directoryCorrectedEEMS = directoryCorrectedEEMS, waves = waves)

######## end of loop!
#write file containing spectral indicies + sample IDs
#after loop is finished with all samples
corrpath <- file.path(directoryCorrectedEEMS, paste(project, "SpectralIndicies.csv", sep = ""))
write.table(Spectral.Indicies, file = corrpath, row.names = FALSE, col.names = TRUE, sep = ",")

# Check spectral indicies
Spectral.Indicies

####################################################################################################################################
############################ Cutting EEMS for Cory McKnight and DOM Fluor toolbox
########
# Ensure that EEMS are all the same size + works for both the CM code as well as the DOMFluor toolbox to code together
# For DBP, this means ex = 240-800 in 2 nm incrmenets, noting that the  

# filelist of corrected EEMS - for both types of PARAFAC modelling

setwd(directoryCorrectedEEMS) 
filelist_EEMScor <- list.files(pattern = "_Corrected.csv$")

######## Prepping files for Cory McKnight modelling in Matlab
########
# CM - prepping for CM PARAFAC model
# Function does three things: trims EEMS according to specified min eexitation wavlenegth
# take out row and column names in first column and row and save in CM folder with _i as per graph headins file
# Also creates ex and em files, as well as graph headings file as txt file and saves in CM file
# Lastly, creates and saves graph headings file as txt file, which the function returns to double check

ex.PARAFAC <- seq(240, 800, by = 2) #change if excitation wavlenegths are different

# call function
setwd("/Users/ashlee/SpecScripts") 
source("EEMCMtrim_function.R")
CMsave <- CMtrim(filedirectory = directoryCorrectedEEMS, filelist = filelist_EEMScor, project = project, exmin = "X240",
                 directoryCM = directoryCM, ex = ex.PARAFAC)

# check graph headings file returned by function
CMsave
####

####################################################################################################################################
######### DOM Fluor
########
# Get files ready for DOMFLuor toolbox. 
# Need .csv file for ex, em and one csv file containing all of the fluorescence EEMS compiled
# Use Save Dr EEMS function
# Inputs include the filelist, the project, the vector containing sample names, and the excitation wavelength min you want to trim to
# File cuts EEMs from ex min that you want to

setwd("/Users/ashlee/SpecScripts") 
source("EEMSDrEEMsave_function.R")

ex.DrEEMS = seq(240, 800, by = 2)
DrEEM.data = DrEEM(filelist = filelist_EEMScor, project = project, 
                   exmin = 'X240', filedirectory = directoryCorrectedEEMS, ex = ex.DrEEMS)

# check DrEEM.data. This is the compiled EEMS for DrEEM PARAFAC modelling
head(DrEEM.data)

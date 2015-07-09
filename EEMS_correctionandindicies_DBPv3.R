###
# things that Script needs to do

# 1. Correct EEMS for Raman
# 2. Correct EEMS and absorbance for dilution factor as well as iron interferences
# 2. Calculate absorbance indicies from absorbance file
# 4. Calculate fluorescence indicies (FI, Humification, freshness)
# 5. Save corrected file and file containing parameters
# 6. Plot contour and save as jpeg file

# Note that this script just calls multiple scripts in which the actual calculation is done. No
# actual calculations are done in this script, other than dilution factor. Everything else written 
# in another script

#   #Recommended progress of corrections from McKnight Lab:
# Instrument Correct the Raman file - done in software
# Instrument Correct and Raman Normalize the Blank
# Instrument Correct the Sample - done in software
# Inner filter Correct the Sample
# Raman Normalize the Sample
# Blank Subtract

## To do's - 
# Fe Corrections
# Investigate why it is cutting out the end of two columns
# also add in BIX to fluorescence indicies scrip
# Saves having to rename files as in matlab script! as well as copy them into matlab folder
###########

## set working directory
rm(list = ls())
ls()

#######
# get file list for the blank, absorbance and EEMS file

library(reshape)
library(plyr)
library(gsubfn)

######
# directories for different projects
# DBP Pre chlorination
#blank directory
#directoryblank <-"/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_blank" 
#abs directory
#directoryAbs <-"/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_Abs" 
# raw EEMS directory
#directoryEEMS <-"/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBPS_prechlor_EEM" 
# directory for corrected EEMS
#directoryCorrectedEEMS <- "/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_correctedEEMs"

####
# DBP Post chlorination
# directory with all of the fluorescence files
directoryall <- "/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_all"

# directory for corrected EEMS
directoryCorrectedEEMS <- "/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_correctedEEMs"

#######
# directory for saving EEMS for CM PARAFAC in 'Correct EEMS" file in the CM PARAFAC folder
# This is the same for all projects
directoryCM <-"/Users/ashlee/Documents/MATLAB/CorrectedEEMS" 

######
#dilution file
top = c("sample.ID", "dilutionfactor")

#DBP pre
#dilution <-as.data.frame(read.csv("/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_Aqualogdilution.csv", 
#                                  sep=",", header = TRUE, col.names = top))

#DBP post
dilution <-as.data.frame(read.csv("/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_Aqualogdilution.csv", 
                                  sep=",", header = TRUE, col.names = top))

#project -> "DBPPre"
project = "DBPPost"

### Should not have to change anything beyond this

###########
# call function to create a graph headings file from abs, EEM and blank file
setwd("/Users/ashlee/SpecScripts") 
source("EEMfilecomp_function.R")

data.3 <- EEMfilecomp(workdir= directoryall, dil = dilution)

#############################################################################
### set up loop for all files in folders
Spectral.Indicies = data.frame(matrix(vector(), 0, 17)) #creating an empty vector

n = nrow(data.3)

ex_all = data.frame(matrix(vector(), 0, n))
em_all = data.frame(matrix(vector(), 0, n))  

#blank - create dataframe
ex_blank = data.frame(matrix(vector(), 0, n))
em_blank = data.frame(matrix(vector(), 0, n))

#set working directory for scripts that have corrections etc
#note that column 1 in data must be = 

for (i in 1:n){
  
  # functions to load, trim and correct EEMS data - blank, EEM and files
  
  #### load and trim files
  #sample name
  samplename <- toString(data.3[i,1]) #column with the sample IDs
  
  #### load and trim EEMS
  #call function. Note that wd will change to sample WD in function
  setwd("/Users/ashlee/SpecScripts") 
  source("EEMfileLoadTrim_function.R")
  EEM <- EEMtrim(graphheadings = data.3, samplewd = directoryall, loopnum = i)

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
  
  #ex and em_all are variables that output the complete ex and em
  ex_blank = rbind(ex_blank, colnames(Blktrim))
  em_blank = rbind(em_blank, rownames(Blktrim))

  #### identify dilution factor in master file
  # Dilution = column 5 in data.3
  dil = data.3[i,5]
  
  #ex and em wavelengths
  ex = as.numeric((sort(colnames(EEM), decreasing = TRUE)))
  em = as.numeric((sort(rownames(EEM), decreasing = FALSE)))
  
  ################################## Corrections
  ########### Correct raw EEM for IFE
  # call function
  setwd("/Users/ashlee/SpecScripts") 
  source("EEMIFECorr_function.R")
  
  EEM.IFC <- innerfilter(eem = EEM, abs = Abstrim, em = em, ex = ex)
  # note that IFE should be between 0.4 and  0.98 according to McKnight 2001 (doi: 10.4319/lo.2001.46.1.0038)
  
  ########### Normalize IFE EEM and blank file according to area under Raman peak
  # call function
  #setwd("/Users/ashlee/Dropbox/R Scripts")
  setwd("/Users/ashlee/SpecScripts") 
  source("Ramancorrect_v1.R")

  # tell R where em = 375 nm, em = 430 nm; ex = 350 nm
  em375 <-  as.numeric(grep(375.7, rownames(EEM)))
  em430 <-  as.numeric(grep(429.8, rownames(EEM)))
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
  # Sample - blank = corrected EEM
  EEM.blk <- EEM.ram - blankram
  
  ###########################
  ##### Apply dilution factor to EEM and to Abs file
  EEM.dil = EEM.blk*dil 
  Absdil = Abstrim*dil #second column = absorbance readings. Assumes Beer's Law applies (c~abs)
  
  ##################################
  ########### Correct for Raleigh Masking
  # call function
  setwd("/Users/ashlee/SpecScripts") 
  source("EEMRaleigh_function.R")
  
  EEM.rm <- raleigh(eem = EEM.dil, slitwidth = 10)
  
  ##### Apply correction factor for Fe concentration
  ##### TO DO!!!!!
  
  # change the corrected EEM so that excitation goes from 200-800, rather than 800-200
  EEMcorr <- EEM.rm[,sort(names(EEM.rm), decreasing = FALSE)]
  
  ###########
  ##### Save the corrected EEM
  #x <- length(EEMcorr)
  #EEMcorr1 <- EEMcorr[,c(1:(x-4))] #cut out last three columns of Na data. Remove this line eventually ;)
  # Note that the file still contains two columns of Na's. Will save the original file, and remove this after corrections (saving for PARAFAC)
  
  corrpath <- file.path(directoryCorrectedEEMS, paste(samplename,"_", project,"_CorrectedNEW",".csv", sep = ""))
  write.table(EEMcorr, file = corrpath, row.names = TRUE,col.names = TRUE, sep = ",")
  
  ###########
  # Calculating absorbance indicies
  # call function
  setwd("/Users/ashlee/SpecScripts") 
  source("Aqualog_Absindicies_v1.R")
  
  #call the function to calculate indicies
  Abs.ind <- Abs(absorbance = Absdil)
  #Abs.all[i] <- cbind(samplename, Abs.ind) #Put sample number
  
  ##########
  # Calculating fluorescence indicies
  # call function
  setwd("/Users/ashlee/SpecScripts") 
  source("Aqualog_Fluorindicies_v2.R")
  
  #below may need to be altered depending on the output of your scan
  # ex wavelengths
  ex370 <- as.numeric(grep(370, colnames(EEMcorr)))
  ex254 <- as.numeric(grep(254, colnames(EEMcorr)))
  ex310 <- as.numeric(grep(310, colnames(EEMcorr)))
  ex274 <- as.numeric(grep(274, colnames(EEMcorr)))
  ex276 <- as.numeric(grep(276, colnames(EEMcorr)))
  ex320 <- as.numeric(grep(320, colnames(EEMcorr)))
  ex340 <- as.numeric(grep(340, colnames(EEMcorr)))
  
  # em wavelengths
  em470 <- as.numeric(grep(470.394, rownames(EEMcorr)))
  em520 <- as.numeric(grep(520.522, rownames(EEMcorr))) 
  em435 <- as.numeric(grep(435.609, rownames(EEMcorr))) 
  em480 <- as.numeric(grep(479.697, rownames(EEMcorr)))
  em300 <- as.numeric(grep(300.484, rownames(EEMcorr)))
  em345 <- as.numeric(grep(344.826, rownames(EEMcorr)))
  em380 <- as.numeric(grep(380.302, rownames(EEMcorr)))
  em420 <- as.numeric(grep(420.587, rownames(EEMcorr)))
  em436 <- as.numeric(grep(436.766, rownames(EEMcorr)))
  em350 <- as.numeric(grep(350., rownames(EEMcorr)))
  em410 <- as.numeric(grep(410., rownames(EEMcorr)))
  em430 <- as.numeric(grep(430., rownames(EEMcorr)))
  
  # call function that calculates fluorescent indicies
  Fluor.ind <- Fluor(eem = EEMcorr)
  
  ##########
  # bind fluor indicies with abs indicies as well as the sample id
  #Spectral.Indicies[i] <- cbind(samplename, Abs.ind, Fluor.ind) 
  
  Spectral.Ind <- cbind(samplename, Abs.ind, Fluor.ind) 
  top <- colnames(Spectral.Ind)
  Spectral.Indicies[i,]  <- Spectral.Ind
  colnames(Spectral.Indicies) <- top
  
  ##########
  # last - plot corrected eems as a contour plot
  #variables to change
  zmax = max(EEMcorr,na.rm=TRUE) # put the max intensity of that you want to graph
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
  plotpath <- file.path(directoryCorrectedEEMS, paste(samplename,"_", project, "_Contour",".jpeg", sep = ""))
  
  g <- length(EEMcorr)
  EEMplot <- EEMcorr # not cutting out the last two columns
  # EEMplot <- EEMcorr[,c(1:(g-4))] #cut out last three columns of Na data. Remove this line eventually ;)
  # EEMplot <- EEMcorr[,complete.cases(EEMcorr)] omits everything
  # create new dataset without missing data 
  
  EEMplot[EEMplot<0]=0 # have line so that if fluorescence is less than 0, becomes 0.
  explot = as.numeric(colnames(EEMplot)) #cut out the last two wavelengths.. these are filled with Nas and cut out before saving EEMS
  emplot = as.numeric(row.names(EEMplot))
  
  jpeg(file=plotpath)
  contour.plots(eems = as.matrix(EEMplot), Title = samplename, ex = explot, em = emplot)  
  dev.off()
  
  # note that the above is meant to be a crude graphing - better graphing done in matlab once
  # you figure out the max emission for your dataset (normalize all of the plots to this)
  
}

#### End of loop!
#write file containing spectral indicies + sample IDs
#after loop is finished with all samples
corrpath <- file.path(directoryCorrectedEEMS, paste("SpectralIndicies.csv", sep = ""))
write.table(Spectral.Indicies, file = corrpath, row.names = FALSE, col.names = TRUE, sep = ",")

############################ Cutting EEMS for Cory McKnight and DOM Fluor toolbox
########
# Ensure that EEMS are all the same size + works for both the CM code as well as the DOMFluor toolbox to code together
# For DBP, this means ex = 240-800 in 2 nm incrmenets, noting that the  
# filelist of corrected EEMS

# call function that trims EEMS according to the min ex wavelength that you specify.
# note that you have to look at teh files to see what the ex and em ranges of the samples are.
# These are held in em_all and ex_all

setwd("/Users/ashlee/SpecScripts") 
source("EEMCMtrim_function.R")

CMsave <- CMtrim(directory = directoryCorrectedEEMS, projectname = project, minex = "X240")

setwd(directoryCorrectedEEMS) 
filelist_EEMScor <- list.files(pattern = "_Corrected.csv$")

n = length(filelist_EEMScor)

######## Prepping files for Cory McKnight modelling in Matlab
########
# CM - take out row and column names in first column and row and save in CM folder

for (i in 1:n){
  temp.EEMS <- read.delim(filelist_EEMScor[i], header= TRUE, sep = ",")
  
  #trim so that exitation and emission goes from the same
  ex.temp <- colnames(temp.EEMS)
  
  if(ex.temp[1] != "X240") {
    # if first value in ex.temp is not 240, trim 
    ex.length <- length(ex.temp)
    # find column where the exitation wavelength is 240 to cut from
    x240 = as.numeric(match("X240",names(temp.EEMS)))
    temp.EEMS <- temp.EEMS[,c(x240:ex.length)]
  } 
  
  # cut out any columns containing Nas- this is 798 and 800 nm. Must cut last four rows of data from 20april2015
  #temp.EEMS.1 <- na.omit(temp.EEMS)
  g <- length(temp.EEMS)
  temp.EEMS.1 <- temp.EEMS[,c(1:(g-4))] #cut out the last four colomns manually
  
  #resave without the row and column names
  samplename <- strapplyc(filelist_EEMScor[i], "(.*)Prechlor", simplify = TRUE)
  corrpath <- file.path(directoryCM, paste(samplename,"PrechlorCorrCM_",i,".csv", sep = ""))
  write.table(temp.EEMS.1, file = corrpath, row.names = FALSE,col.names = FALSE, sep = ",")
  
}

# save ex and em in two separate files, to make it easier to read into CM PARAFAC files
corrpath <- file.path("/Users/ashlee/Documents/MATLAB/ExEmfiles", paste(project,"em",".csv", sep = ""))
write.table(em, file = corrpath, row.names = FALSE,col.names = FALSE, sep = ",")

corrpath <- file.path("/Users/ashlee/Documents/MATLAB/ExEmfiles", paste(project,"ex",".csv", sep = ""))
ex.PARAFAC <- seq(240, 796, by = 2)
write.table(ex.PARAFAC, file = corrpath, row.names = FALSE,col.names = FALSE, sep = ",")

######### DOM Fluor
########
# Get files ready for DOMFLuor toolbox. 
# Need .csv file for ex, em and one csv file containing all of the fluorescence compiled

#compiled file containing all fluorescence data
# first, cut files so that they go from 240 - 800 nm
#setwd(directoryCorrectedEEMS) 
#filelist_EEMScorr <- list.files(path = directoryCorrectedEEMS, pattern = "Corrected.csv$")
#x = length(filelist_EEMScorr)

n = length(filelist_EEMScor)
for (i in 1:n){
  temp.EEMS <- read.delim(filelist_EEMScor[i], header= TRUE, sep = ",")
  
  #trim so that exitation and emission goes from the same
  ex.temp <- colnames(temp.EEMS)
  
  if(ex.temp[1] != "X240") {
    # if first value in ex.temp is not 240, trim 
    ex.length <- length(ex.temp)
    # find column where the exitation wavelength is 240 to cut from
    x240 = as.numeric(match("X240",names(temp.EEMS)))
    temp.EEMS <- temp.EEMS[,c(x240:ex.length)]
  } 
  
  # cut out any columns containing Nas- this is 798 and 800 nm. Must cut last four rows of data from 20april2015
  #temp.EEMS.1 <- na.omit(temp.EEMS) #DOESN'T WORK!! omits everything
  g <- length(temp.EEMS)
  temp.EEMS.1 <- temp.EEMS[,c(1:(g-4))] #cut out the last four colomns manually
  
  # create a new dataset where the post-cut EEMS are compiled together by rows
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- temp.EEMS.1
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-temp.EEMS.1
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
}
#seems to have doubled first dataset, remove?
x <- length(em)
y <- dim(dataset)[1]
dataset.2 <- dataset[c((x+1):y),]
y <- dim(dataset.2)[1]
remove(x)
remove(y)

corrpath <- file.path("/Users/ashlee/Documents/MATLAB/ExEmfiles", paste(project,"ex",".csv", sep = ""))

write.table(dataset.2, file = file.path("/Users/ashlee/Documents/MATLAB/DOMFluor/", paste(project, "/fl.csv", sep = "")),
            row.names = FALSE, col.names = FALSE, sep = ",") #saved in matlab folder

#Ex file
write.table(ex.PARAFAC, file = file.path("/Users/ashlee/Documents/MATLAB/DOMFluor/", paste(project,"/Ex.csv", sep = "")),
            row.names = FALSE, col.names = FALSE, sep = ",") #saved in matlab folder

#Em
write.table(em, file = file.path("/Users/ashlee/Documents/MATLAB/DOMFluor/", paste(project,"/Em.csv", sep = "")), 
            row.names = FALSE, col.names = FALSE, sep = ",") #saved in matlab folder

#File containing sample names
write.table(sample.ID, file = file.path("/Users/ashlee/Documents/MATLAB/DOMFluor/", paste(project,"/01key.csv",sep = "")), 
            row.names = FALSE, col.names = FALSE, sep = ",") #saved in matlab folder

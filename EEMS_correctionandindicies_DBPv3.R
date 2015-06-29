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

#setwd(directoryall)

# Should not have to change anything beyond this

###########
# call function to create a graph headings file from abs, EEM and blank file
setwd("/Users/ashlee/SpecScripts") 
source("EEMfilecomp_function.R")

data.3 <- EEMfilecomp(workdir= directoryall, dil = dilution)

#############################################################################
### set up loop for all files in folders
Spectral.Indicies = data.frame(matrix(vector(), 0, 16)) #creating an empty vector

n = nrow(data.3)
ex_all = data.frame(matrix(vector(), 0, n))
em_all = data.frame(matrix(vector(), 0, n))

#set working directory for scripts that have corrections etc
#note that column 1 in data must be = 

for (i in 1:n){
  
  # functions to load and trim EEMS data - blank, EEM and files
  
  # load and trim EEMS
  #call function. Note that wd will change to sample WD in function
  setwd("/Users/ashlee/SpecScripts") 
  source("")
  
  #### load and trim ABS
  # call function. Note that wd will change to sample WD in function
  setwd("/Users/ashlee/SpecScripts") 
  source("")
  Abstrim <- ABStrim(graphheadings = data.3, samplewd = directoryall)
  
  # load and trim blank
  Blktrim <- BLANKtrim(graphheadings = data.3, samplewd = directoryall)
  
  # identify dilution factor in master file
  # Dilution = column 5 in data.3
  dil = data.3[i,5]
  

  # Read in the EEM file
  EEMSfilename <- toString(data.3[i,2]) # set EEMS file for the sample
  #setwd(directoryEEMS) 
  EEMSfile <- read.delim(EEMSfilename, header= FALSE, sep = "")
  #samplename <- strapplyc(EEMSfilename, "001(.*)PEM", simplify = TRUE)
  samplename <- toString(data.3[i,1]) 
  
  # Blankfilename <- test2[i,2] # set blank file for the sample
  Blankfilename <- toString(data.3[i,4]) 
  #setwd(directoryblank) 
  Blankfile <- read.delim(Blankfilename, header= FALSE, sep = "")
  
  #Absfile <- test2[i,3] # set Abs file for the sample
  Absfilename <- toString(data.3[i,3]) 
  #setwd(directoryAbs) 
  Absfile <- read.delim(Absfilename, header= FALSE, sep = "")
  
  
  ############## trim and arrange the abs, blank and eem files
  # EEM file
  #set em and ex based on the columns (emission) and rows (excitation) of the raw file
  x = ncol(EEMSfile)
  ex_initial = as.numeric(EEMSfile[1, c(4:x)])
  ex = as.numeric((sort(EEMSfile[1, c(4:x)], decreasing = FALSE)))
  ex_all[i,] = ex
  
  #emission = y axis
  y = nrow(EEMSfile)
  em = as.numeric(t(data.frame(EEMSfile[c(3:y), 1])))
  em_all[i,] = em
  
  # take out the first three rows and first column of data in EEMS. These just contain text
  EEMScut <- EEMSfile[c(3:y), c(2:x)]
  colnames(EEMScut) <- ex_initial
  rownames(EEMScut) <- em
  
  #### Blank
  # take out the first three rows and first column of data in Blank
  #emission = y axis
  y = nrow(Blankfile)
  #excitation - x axis
  x = ncol(Blankfile) 
  
  Blkcor <- Blankfile[c(2:y), c(2:x)]
  colnames(Blkcor) <- ex_initial
  rownames(Blkcor) <- em
  
  ######## trim abs file
  row.names(Absfile) <- Absfile[,1]
  #Abs <- data.frame(Absfile[,-1])
  Absfile1 = data.frame(t(Absfile))
  Absfile2 = Absfile1[-1,]
  
  ##################################
  ########### Correct for Raman
  # call function
  #setwd("/Users/ashlee/Dropbox/R Scripts")
  setwd("/Users/ashlee/SpecScripts") 
  
  source("Ramancorrect_v1.R")
  
  # tell R where em = 375 nm, em = 430 nm; ex = 350 nm
  em375 <-  as.numeric(grep(375.7, rownames(EEMScut)))
  em430 <-  as.numeric(grep(429.8, rownames(EEMScut)))
  ex350 <- as.numeric(match(350, colnames(EEMScut)))
  
  Raman.area <- Ramancor(blank = Blkcor) # get the Raman correction file from the Raman function stored
  
  # normalize all of the EEM file for the area underneath the raman file  
  EEMram = EEMScut/Raman.area 
  
  ###########################
  ##### Apply dilution factor to EEM and to Abs file
  EEMdil = EEMram*dil 
  Absdil = Absfile2*dil #second column = absorbance readings. Assumes Beer's Law applies (c~abs)
  
  ##### Apply correction factor for Fe concentration
  #####TO DO!!!!!
  
  # change the corrected EEM so that excitation goes from 200-800, rather than 800-200
  EEMcorr <- EEMdil[,sort(names(EEMdil), decreasing = FALSE)]
  # cutting out two columns - 200 and 202? ***************************
  
  ##### Save the corrected EEM
  #x <- length(EEMcorr)
  #EEMcorr1 <- EEMcorr[,c(1:(x-4))] #cut out last three columns of Na data. Remove this line eventually ;)
  # Note that the file still contains two columns of Na's. Will save the original file, and remove this after corrections (saving for PARAFAC)
  
  corrpath <- file.path(directoryCorrectedEEMS, paste(samplename,"_", project,"_Corrected",".csv", sep = ""))
  write.table(EEMcorr, file = corrpath, row.names = TRUE,col.names = TRUE, sep = ",")
  
  ###########
  # Calculating absorbance indicies
  # call function
  #setwd("/Users/ashlee/Dropbox/R Scripts") 
  source("Aqualog_Absindicies_v1.R")
  
  #call the function to calculate indicies
  Abs.ind <- Abs(absorbance = Absdil)
  #Abs.all[i] <- cbind(samplename, Abs.ind) #Put sample number
  
  ##########
  # Calculating fluorescence indicies
  # call function
  #setwd("/Users/ashlee/Dropbox/R Scripts") 
  source("Aqualog_Fluorindicies_v2.R")
  
  #below may need to be altered depending on the output of your scan
  ex370 <- as.numeric(grep(370, colnames(EEMcorr)))
  em470 <- as.numeric(grep(470.394, rownames(EEMcorr)))
  em520 <- as.numeric(grep(520.522, rownames(EEMcorr))) 
  em435 <- as.numeric(grep(435.609, rownames(EEMcorr))) 
  em480 <- as.numeric(grep(479.697, rownames(EEMcorr)))
  ex254 <- as.numeric(grep(254, colnames(EEMcorr)))
  em300 <- as.numeric(grep(300.484, rownames(EEMcorr)))
  em345 <- as.numeric(grep(344.826, rownames(EEMcorr)))
  ex310 <- as.numeric(grep(310, colnames(EEMcorr)))
  em380 <- as.numeric(grep(380.302, rownames(EEMcorr)))
  em420 <- as.numeric(grep(420.587, rownames(EEMcorr)))
  em436 <- as.numeric(grep(436.766, rownames(EEMcorr)))
  
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
  #setwd("/Users/ashlee/Dropbox/R Scripts") 
  source("EEM_contour_v1.R")
  
  #Plot contours and save in correction file
  plotpath <- file.path(directoryCorrectedEEMS, paste(samplename,"_", project, "_Contour",".jpeg", sep = ""))
  
  g <- length(EEMcorr)
  EEMplot <- EEMcorr[,c(1:(g-4))] #cut out last three columns of Na data. Remove this line eventually ;)
  #EEMplot <- EEMcorr[,complete.cases(EEMcorr)] omits everything
  # create new dataset without missing data 
  
  EEMplot[EEMplot<0]=0 # have line so that if fluorescence is less than 0, becomes 0.
  explot = ex[c(1:(length(ex)-2))] #cut out the last two wavelengths.. these are filled with Nas and cut out before saving EEMS
  
  jpeg(file=plotpath)
  contour.plots(eems = as.matrix(EEMplot), Title = samplename, ex = explot)  
  dev.off()
  
  # note that the above is meant to be a crude graphing - better graphing done in matlab once
  # you figure out the max emission for your dataset (normalize all of the plots to this)
  
}

####
#write file containing spectral indicies + sample IDs
#after loop is finished with all samples
corrpath <- file.path(directoryCorrectedEEMS, paste("SpectralIndicies.csv", sep = ""))
write.table(Spectral.Indicies, file = corrpath, row.names = FALSE, col.names = TRUE, sep = ",")

############################
########
# Ensure that EEMS are all the same size + works for both the CM code as well as the DOMFluor toolbox
# For DBP, this means ex = 240-800 in 2 nm incrmenets, noting that the  
#filelist of corrected EEMS

setwd(directoryCorrectedEEMS) 
filelist_EEMScor <- list.files(pattern = "_Corrected.csv$")

n = length(filelist_EEMScor)

######## 
# Prepping files for CM modelling in Matlab
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
corrpath <- file.path("/Users/ashlee/Documents/MATLAB/ExEmfiles", paste("DBPpre","em",".csv", sep = ""))
write.table(em, file = corrpath, row.names = FALSE,col.names = FALSE, sep = ",")

corrpath <- file.path("/Users/ashlee/Documents/MATLAB/ExEmfiles", paste("DBPpre","ex",".csv", sep = ""))
ex.PARAFAC <- seq(240, 796, by = 2)
write.table(ex.PARAFAC, file = corrpath, row.names = FALSE,col.names = FALSE, sep = ",")

######### DOM Fluor
# get files ready for DOMFLuor toolbox. Need .csv file for ex, em and one csv file containing all of the fluorescence compiled

#compiled file containing all fluorescence data
# first, cut files so that they go from 240 - 800 nm
#setwd(directoryCorrectedEEMS) 
#filelist_EEMScorr <- list.files(path = directoryCorrectedEEMS, pattern = "Corrected.csv$")
#x = length(filelist_EEMScorr)

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

write.table(dataset.2, file = "/Users/ashlee/Documents/MATLAB/DOMFluor/DBP_pre/fl.csv", 
            row.names = FALSE, col.names = FALSE, sep = ",") #saved in matlab folder

#Ex file
write.table(ex.PARAFAC, file = "/Users/ashlee/Documents/MATLAB/DOMFluor/DBP_pre/Ex.csv", 
            row.names = FALSE, col.names = FALSE, sep = ",") #saved in matlab folder

#Em
write.table(em, file = "/Users/ashlee/Documents/MATLAB/DOMFluor/DBP_pre/Em.csv", 
            row.names = FALSE, col.names = FALSE, sep = ",") #saved in matlab folder

#File containing sample names
write.table(sample.ID, file = "/Users/ashlee/Documents/MATLAB/DOMFluor/DBP_pre/01key.csv", 
            row.names = FALSE, col.names = FALSE, sep = ",") #saved in matlab folder

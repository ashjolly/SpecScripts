#
# 
# Script for compiling par and fp files and calculating DOC, TOC, NO3 averages from multiple measurements
# Jan 21 2015 AJ
# For DBP sample files/WL files
#######

# set working directory
rm(list = ls())
ls()

# set path - WL
spectro.direct <- "/Users/user/Dropbox/PhD Work/PhD Data/WL_data/WL_spectrodata"
fppar.directory <- "/Users/user/Dropbox/PhD Work/PhD Data/WL_data/WL_spectrodata/WL_fppar"
project <- "WL"

setwd(fppar.directory) 

# confirm path
getwd()

##########
## first compile all par and fp files within a folder

#install.packages("reshape")
library(reshape)
library(plyr)
library(gsubfn)

##########
#read in the fp and par files

#ensure the headings for all the files are in correct format: DateTime, Status_0, etc.. and samples have column
filelist.fp <- list.files(pattern = ".fp$")
filelist.par <- list.files(pattern = ".par") #note that one WL file is .par.txt

fp <- length(filelist.fp)
par <- length(filelist.par) #check to make sure that both are the same length

sample.ID <- 0 #create sample ID variable

for (i in 1:fp){
  sample.ID.temp <- strapplyc(filelist.fp[i], "_(.*)_.", simplify = TRUE)
  sample.ID[i] <- sample.ID.temp
}

sampleID.par <- 0
for (i in 1:par){
  sample.ID.temp <- strapplyc(filelist.par[i], "_(.*)_.", simplify = TRUE)
  sampleID.par[i] <- sample.ID.temp
}

filelist.fpsampleID <- as.data.frame(cbind(filelist.fp, sample.ID))
filelist.parsampleID <- as.data.frame(cbind(filelist.par, sampleID.par))
colnames(filelist.parsampleID) <- c("filelist.par", "sample.ID")

filelist.fppar <- as.data.frame(merge(filelist.fpsampleID, filelist.parsampleID, by = "sample.ID", all = TRUE))
#filelist.fppar <- as.data.frame(cbind(filelist.par, filelist.fpsampleID))
#remove(filelist.fpsampleID)

###########
#dilution file
dil.top = c("sample.ID", "dilutionfactor")
dilution <- as.data.frame(read.csv(paste(spectro.direct, "WL_dilutionfactors.csv", sep = "/"), sep=",", header = FALSE, col.names = dil.top))

filelist.all <- merge(filelist.fppar, dilution, by = "sample.ID", all = TRUE)
#make sure only unique entries are in
#filelist.all <- unique(filelist.pre)
remove(dilution)
remove(dil.top)

#############################################################################
# run loop to calculate averages, calculate absorbance parameters, account for dilution, correct DOC

n = nrow(filelist.all)
#abs.data <- data.frame(matrix(vector(), 0, 28))
for (i in 1:n){
  
  # par file
  parfilename <- toString(filelist.all[i,3])
  top.par = c('Date','Time', 'Status',  'Turb.FTU',  'Turb.status',  'NO3.N',  'NO3.status',	'TOC',	'TOC.status',	'DOC',	'DOC.status','SAC254',	'SAC254.status',	'Level [m]NaN-NaN_2',	'[Level_0.0_1.0_0.0_0.0]',	'Temp [°C]NaN-NaN_2',	'[Temp_0.0_1.0_0.0_0.0]',	'analogINNaN-NaN_2',	'[analogIN_0.0_1.0_0.0_0.0]')
  par.filepath <- file.path(fppar.directory, paste(parfilename, sep = ""))
  parfile <- read.csv(par.filepath, header= FALSE, sep = "", skip=2, na.strings = "NaN")
  t1 <- length(top.par)
  t2 <- dim(parfile)[2]
  
  if (t2 == t1){
    colnames(parfile) <- top.par
  }
  
  # fp file
  fpfilename <- toString(filelist.all[i,2])
  top.fp = c('Date','Time', 'Status_0', seq(200, 750, by = 2.5))
  fp.filepath <- file.path(fppar.directory, paste(fpfilename, sep = ""))
  fpfile <- read.csv(fp.filepath, header= FALSE, sep = "", skip=2, na.strings = "NaN", col.names = top.fp) 
                     
  #dilution factor
  dil = filelist.all[i,4]
  
  #merge fp and par files
  fppar <- merge(parfile, fpfile, by = 'Time', all = TRUE)
  
  ###########################
  # Calculating absorbance indicies
  # call function
  setwd("/Users/user/SpecScripts") 
  source("Spectro_SpecindiciesAbs_function.R")
  
  #call the function to calculate indicies
  Abs.indicies <- Abs.ind(spec = fppar) #path length = 33 mm
  
  # bind together calculated Abs inidcies - calculated from 
  Abs.all <- cbind(parfile[,c(4,6,8,10)], Abs.indicies)
  
  #take average of all columns
  Abs.ave <- t(data.frame(colMeans(Abs.all)))
  
  #take stdev of columns
  Abs.stdev <- data.frame(t(apply(Abs.all, 2, sd)))
  
  #account for dilution factor
  Abs.dil <- data.frame(Abs.ave*dil)
  
  #correct DOC according to Shimadzu calibration
  Abs.dil$DOCcorr <-  as.numeric(Abs.dil$DOC)*1.0271-0.1998
  
  # associate with sample id
  sample <- as.character(filelist.all[i,1])
  abs.data.temp <- cbind(sample, Abs.dil, Abs.stdev)
  
  #abs.data[i,] <- abs.data.temp
  #rm(abs.data.temp)
  
  # save file
  # if the  dataset doesn't exist, create it
  if (!exists("abs.data")){
    abs.data <- abs.data.temp
  }
  
  # if the merged dataset does exist, append to it
  if (exists("abs.data")){
    temp_dataset <- abs.data.temp
    abs.data<-rbind(abs.data, temp_dataset)
    rm(temp_dataset)
  }                  
}
###############
# need to remove duplicates?

# sort by sample ID
abs.data <- arrange(abs.data, sample) 

############### 
# write dataset
filepath <- file.path(spectro.direct, paste(project,"_calculatedspec.csv", sep = ""))
write.table(abs.data, file = filepath, sep = ",", qmethod = "double")

############### 
########### end
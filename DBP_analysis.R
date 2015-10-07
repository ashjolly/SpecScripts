#
#
# Script for analyzing DBP data
# 29sept2015
# For analyzxing data 
################################################################################
# clean up list
rm(list = ls())
ls()
################################################################################
# install libraries necessary for analysis
library('gsubfn')
library('abind')
library('zoo')

################################################################################
# Read in data used for analysis
## Water quality and location data

############ Fluorescence data
## Prechlorinated EEMS
pre.directory <- '/Users/ashlee/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_correctedEEMs'

## post Chlorinated EEMS

# prechlorination indicies
# postchlorination indicies

# Prechlorination CM Model
# Post clorination CM model

#Prechlorination PARAFAC model
# post chlorination PARAFAC model


################## DBP concentration
# THM data
# HAA data

################################################################################
# Pt 1 - how does chlorination change the spectral composition of EEMS?
# princip-al component analysis to look at how chlorination altered spectral characteristics (EEMS)
# References:
# http://www.r-bloggers.com/computing-and-visualizing-pca-in-r/
# http://planspace.org/2013/02/03/pca-3d-visualization-and-clustering-in-r/

## PCA on raw EEMS 
# Aim of this is to see regions that explain most of the variance within the pre chlorinated EEMS
# Examine all data lumped together to see any alterations 
# read in data - compiled data according to sample ID

# locate the prechlorinated corrected eems within the file
setwd(pre.directory) 
filelist_DBPpre <- list.files(pattern = "_Corrected.csv$")

# create graph heading variable
graphheadings = data.frame((0))

# compile all of the corrected pre chlorination EEMS and  correct from Raman and Raleigh scatter
n = length(filelist_DBPpre)
exmin = 'X240'
project = 'DBPPre'

# run loop over all files within the corrected file list
for (i in 1:n){
  # set working directory back to directory with sample ID + read in EEMs
  setwd(pre.directory) 
  temp.EEMS <- read.delim(filelist_DBPpre[i], header= TRUE, sep = ",")
  
  # ensure that the ex ranges are the same for all of the data - 240 to 200 nm in 2 nm increments
  ex.temp <- colnames(temp.EEMS)
  if(ex.temp[1] != exmin) {
    # if first value in ex.temp is not 240, trim 
    ex.length <- length(ex.temp)
    # find column where the exitation wavelength is 240 to cut from
    x240 = as.numeric(match(exmin,names(temp.EEMS)))
    temp.EEMS <- temp.EEMS[,c(x240:ex.length)]
  } 
  
  # Correct for Raleigh scatter using function - interpolates for first and second Raleigh
  setwd("/Users/ashlee/SpecScripts") 
  source("EEMRaleigh_function.R")
  temp.cut <- raleigh(eem = temp.EEMS, slitwidth1 = 15, slitwidth2 = 15)
  
  # get the emission variables from the EEM
  em = row.names(temp.cut)
  # get the excitation variables from the cut EEMS
  ex = colnames(temp.cut)
  
  #save as a array
  # if the merged dataset does exist, append to it by column
  if (exists("EEM.dataset")){
    #temp_dataset <- temp.cut
    EEM.dataset <- abind(EEM.dataset, temp.cut, along = 3)
    #rm(temp.cut)
  }
  
  # if the merged dataset doesn't exist, create it
  if (!exists("EEM.dataset")){
    EEM.dataset <- temp.cut
  }
  
  # Create graph headings variable to identify samples
  samplename <- strapplyc(filelist_DBPpre[i], paste("(.*)_", project, "_Corrected", sep = ""), simplify = TRUE)
  graphheadings[i,] <-samplename
  
}

# Compile and decompose EEM in array such that ex-em pairs are the columns and the sample ID is the row prior to pCA
# Use apply in order to speed this up? very slow as a loop?

n = dim(EEM.dataset)[3]

for (i in 1:n){
  temp.EEM <- as.data.frame(EEM.dataset[,,i])
  
  # Use function to unfold the EEM
  setwd("/Users/ashlee/SpecScripts") 
  source("EEMPCAcompile_function.R")
  PCA.data <- PCA.EEM(EEM = temp.EEM)
  
}

# add sample ID to first column of the PCA.data dataset
PCA.data <- cbind(graphheadings, PCA.data)

## Do PCA on the compiled pre-chlorinated data
pca.pre <- prcomp(PCA.data, center = TRUE, scale. = TRUE)

# Analyze PCA results
# print method
print(pca.pre)

# plot method - decide how many components to keep
plot(pca.pre, type = "l")

# summary method
summary(pca.pre)

# graph first principals of variation
library('devtools')
install_github("ggbiplot", "vqv")
library(ggbiplot)
g <- ggbiplot(pca.pre, obs.scale = 1, var.scale = 1, 
              groups = ir.species, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')

## PCA on the Pre and Post chlorinated EEMS
# Aim is to see where in the spectra you see the greatest changes between the pre and post chlorinated samples


###### PCA on PARAFAC results - CM model and DOMFluor model

###########
# Pt 2 - How does EEMS and Water Quality parameters predict DBP formation?
# thoughts:
# use machine learning algorythm to predict DBP parameters
# try multivraite linear model where
#[DBP] = theta0 + theta1x1+theta2x2
# where x = spectral and EEM characteristics, as well as water quality parameters (pH, anion concentrations, DOC)
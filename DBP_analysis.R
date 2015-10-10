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
pre.directory <- '/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_correctedEEMs'

## post Chlorinated EEMS
post.directory <- '/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_correctedEEMs'

# prechlorination indicies
# postchlorination indicies

# Prechlorination CM Model
# Post clorination CM model

# Prechlorination PARAFAC model
# post chlorination PARAFAC model

################## DBP concentration
# THM data
# HAA data

################################################################################
# Pt **** - how does chlorination change the spectral composition of EEMS?
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
  setwd("/Users/user/SpecScripts") 
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

#function for transposing and compiling EEMS according to column
PCA.EEM <- function(EEM){
  # Compile and decompose EEM such that ex-em pairs are the columns and the sample ID is the row
    temp.PCA = data.frame(t(EEM))
    colnames(temp.PCA) = (paste(colnames(EEM), row.names(EEM), sep = '_'))
    
    # if the merged dataset does exist, append to it by column
    if (exists("PCA.data")){
      #temp_dataset <- temp.cut
      PCA.data<-cbind(PCA.data, temp.PCA)
      #rm(temp.cut)
    }
    
    # if the merged dataset doesn't exist, create it
    if (!exists("PCA.data")){
      PCA.data <- temp.PCA
    }
}

# create empty vector
EEM.row = data.frame(matrix(vector(), 5000, 200000))

# assemble using vectorized solution
for (i in 1:n){
  # get the sample data frame
  temp.EEM <- EEM.dataset[,,i]
  
  # use function and apply to reorganize for PCA
  EEM.row[i,] <- as.data.frame(apply(temp.EEM, 2, PCA.EEM))
}

# Cut out the NAs inserted when you created a dataframe to put the PCA formatted data into
PCA.pre <- EEM.row[rowSums(is.na(EEM.row)) != ncol(EEM.row),]
PCA.pre <- EEM.row[colSums(is.na(EEM.row)) != nrow(EEM.row),]

# ensure that the column names are properly named
colnames(EEM.row) = (paste(colnames(EEM.row), row.names(EEM.dataset[,,1]), sep = '_'))
row.names(EEM.row) = graphheadings[,1]

#######################
## Do PCA on the compiled pre-chlorinated data
pca.pre <- prcomp(EEM.row, center = TRUE, scale. = TRUE)

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
g <- ggbiplot(pca.pre), obs.scale = 1, var.scale = 1) 
              #groups = ir.species, ellipse = TRUE, 
              #circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')

######################################## 
# Comparing pre to post chlorinated EEMS

#### 
# first add column to sample ID that specarates unchlorinated or chlorinated EEMS


######
# PCA on the Pre and Post chlorinated EEMS
# Aim is to see where in the spectra you see the greatest changes between the pre and post chlorinated samples




####
# Self organizing maps - pre and post chlorination EEMS
# Use to look at how chlorination affects the spectral characteristics of 



###### PCA on PARAFAC results - CM model and DOMFluor model

###########
# Pt 2 - How does EEMS and Water Quality parameters predict DBP formation?
# thoughts:
# use machine learning algorythm to predict DBP parameters
# try multivraite linear model where
#[DBP] = theta0 + theta1x1+theta2x2
# where x = spectral and EEM characteristics, as well as water quality parameters (pH, anion concentrations, DOC)
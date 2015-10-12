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
# directory for data used for analysis
save.directory <- '/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_analysisdata'

##############################################################
# functions used in script

################################################################################
# Analysis Script
# Pt **** - how does chlorination change the spectral composition of EEMS?
# princip-al component analysis to look at how chlorination altered spectral characteristics (EEMS)
# References:
# http://www.r-bloggers.com/computing-and-visualizing-pca-in-r/
# http://planspace.org/2013/02/03/pca-3d-visualization-and-clustering-in-r/

## PCA on raw EEMS 
# Aim of this is to see regions that explain most of the variance within the pre chlorinated EEMS
# Examine all data lumped together to see any alterations 

#####################################################################
## Do PCA  on the compiled pre-chlorinated data
# read in data compiled for PCA analysis - note that data is in array
PCA.pre <- readRDS(paste(save.directory, "/PCApre.rds", sep = ""))
# add in a variable that organizes according to which watersheds are drinking, which are protected, etc..
# To cluster.. see if there is a pattern within watersheds that are protected

# read in file containing pre chlor EEMs assembled for PCA analysis
# perofrm PCA analysis on all pre chlorinated EEMS
pca.pre <- prcomp(PCA.pre, center = TRUE, scale. = TRUE)

# Analyze PCA results
# print method
print(pca.pre)

# plot method - decide how many components to keep
plot(pca.pre, type = "l")
autoplot(prcomp(PCA.pre), data = PCA.pre)

# summary method
summary(pca.pre)
head(pca.pre$rotation)

# graph first principals of variation
library('devtools')
install_github("ggbiplot", "vqv")
library(ggbiplot)
g <- ggbiplot(pca.pre)
#, obs.scale = 1, var.scale = 1) 
              #groups = ir.species, ellipse = TRUE, 
              #circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
### Post chlor eems

PCA.post <- readRDS(paste(save.directory, "/PCApost.rds", sep = ""))


######################################## 
# Comparing pre to post chlorinated EEMS
#### 
# first add column to sample ID that specarates unchlorinated or chlorinated EEMS
# EEM.pre = prechlorinated EEM array
# EEM.post = post chlorinated EEM array
 #no need to do this? Go by sample ID?
# compile the pre and post chlorinated PCA data together



######
# PCA on the Pre and Post chlorinated EEMS
# Aim is to see where in the spectra you see the greatest changes between the pre and post chlorinated samples




####
# Self organizing maps - pre and post chlorination EEMS
# Use to look at how chlorination affects the spectral characteristics of 
# http://www.r-bloggers.com/self-organising-maps-for-customer-segmentation-using-r/
# Use Data compiled for PCA analysis - samples on the rows, variables on the columns



###### PCA on PARAFAC results - CM model and DOMFluor model

###########
# Pt 2 - How does EEMS and Water Quality parameters predict DBP formation?
# thoughts:
# use machine learning algorythm to predict DBP parameters
# try multivraite linear model where
#[DBP] = theta0 + theta1x1+theta2x2
# where x = spectral and EEM characteristics, as well as water quality parameters (pH, anion concentrations, DOC)
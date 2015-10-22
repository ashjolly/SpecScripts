# Script for analyzing DBP data
# 29sept2015
# For analysing data 
# Ashlee Jollymore's PhD
# DBP project
################################################################################

# clean up list
rm(list = ls())
ls()
################################################################################
# install libraries necessary for analysis
library('gsubfn')
library('abind')
library('zoo')
library('devtools')
install_github("ggbiplot", "vqv")
library(ggbiplot)
################################################################################
# Read in data used for analysis
## Water quality and location data
# directory for data used for analysis
save.directory <- '/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_analysisdata'

##############################################################
# functions used in script

################################################################################
# Pt *** DIfferences in water quailty parameters between sites
# Question - how different are sites in terms of water quality parameters?
# show in histogram of DOC, TOC, pH, temp, [DO], SUVA, Br, TN concentrations

# File with all of the water quality parameters

waterquality <- read.csv((paste(save.directory, '/DBP_master_v5.csv', sep = "")), header = TRUE)

# Histograms of various water quality parameters



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
pcapre.loadings <- pca.pre$rotation[,1:4]
plot(pcapre.loadings[,1:2], type = 'p')
plot(pcapre.loadings[,3:4], type = 'p')

# graph first principals of variation

g <- ggbiplot(pca.pre)
#, obs.scale = 1, var.scale = 1) 
              #groups = ir.species, ellipse = TRUE, 
              #circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')

###########################################################
### Post chlor eems
# read in compiled file
PCA.post <- readRDS(paste(save.directory, "/PCApost.rds", sep = ""))

# run PCA on the post chlorinated EEMS 
pca.post <- prcomp(PCA.post, center = TRUE, scale. = TRUE)

# plot method - decide how many components to keep
plot(pca.post, type = "l")
autoplot(prcomp(PCA.post), data = PCA.post)

# summary method
summary(pca.post)
head(pca.post$rotation)
pcapost.loadings <- pca.post$rotation[,1:4]
plot(pcapost.loadings[,1:2], type = 'p')
plot(pcapost.loadings[,3:4], type = 'p')

################################################################################  
# Comparing pre to post chlorinated EEMS
#### 
# first add column to sample ID that specarates unchlorinated or chlorinated EEMS
# EEM.pre = prechlorinated EEM array
# EEM.post = post chlorinated EEM array
 #no need to do this? Go by sample ID?
# compile the pre and post chlorinated PCA data together

PCA.all <- rbind(PCA.pre, PCA.post)

# create a column to deliniate pre versus post chlorinated EEMS
chlor = data.frame((0))
PCA.all <- cbind(chlor, PCA.all)

#insert chlor
PCA.all[1:117,1] <- "pre chlorination"
PCA.all[118:237,1] <- 'post chlorination'
#rewrite chlor
chlor <- PCA.all[,1]

# Do PCA on pre + post chlorinated EEMS - which wavelengths result in greatest difference
pca.all <- prcomp(PCA.all[,2:140501], center = TRUE, scale. = TRUE)

plot(pca.all, type = "l")

pcaall.loadings <- pca.all$rotation[,1:4]
plot(pcaall.loadings[,1:2], type = 'p')
plot(pcaall.loadings[,3:4], type = 'p')

g <- ggbiplot(pca.all, obs.scale = 1, var.scale = 1, 
  groups = PCA.all[,1], ellipse = TRUE, circle = TRUE)

####
# Self organizing maps - pre and post chlorination EEMS
# Use to look at how chlorination affects the spectral characteristics of 
# http://www.r-bloggers.com/self-organising-maps-for-customer-segmentation-using-r/
# Use Data compiled for PCA analysis - samples on the rows, variables on the columns

require('kohonen')

# fit pre + post to 5x4 SOM
DBPall.som <- som(as.matrix(PCA.all[,2:140501]), somgrid(5,4, 'hexagonal'))

#plot
plot(DBPall.som, type = 'mapping')

####################################################################################
# Analyzing the difference between CM models - pre and post chlorination EEMS
# load in pre and post CM results
CM.pre <- read.csv(paste(save.directory, "/DBPpre_componentsandloadings_CM.csv", sep = ""))
CM.post <- read.csv(paste(save.directory, "/DBPpost_componentsandloadings_CM.csv", sep = ""))

###### PCA on PARAFAC results - CM model and DOMFluor model 
CM.all <- rbind(CM.pre, CM.post)

# add in chlor column
CM.all <- cbind(chlor, CM.all)
row.names(CM.all) <- CM.all$sample.ID
#Remove Nas
CMPCA.all <- na.omit(CM.all)

# do PCA analysis on all dataset - look for components that explain largest differences between pre and post
PCA.CM <- prcomp(CMPCA.all[,3:15], center = TRUE, scale. = TRUE)

# choose top 2 components of variation, and plot against one another to see if you can cluster pre and post chlorination
PCA.CM$rotation

g <- ggbiplot(PCA.CM, obs.scale = 1, var.scale = 1, 
              groups = CMPCA.all$chlor, ellipse = FALSE, circle = FALSE)


###########
# Pt 2 - How does EEMS and Water Quality parameters predict DBP formation?
# Thoughts:
# use machine learning algorythm to predict DBP parameters
# try multivraite linear model where
#[DBP] = theta0 + theta1x1+theta2x2
# where x = spectral and EEM characteristics, as well as water quality parameters (pH, anion concentrations, DOC)
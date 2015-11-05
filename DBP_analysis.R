# Script for analyzing DBP data
# Created 29sept2015
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
library('plyr')
library('stringr')
library('gplots')

################################################################################
# Read in data used for analysis
## Water quality and location data
# directory for data used for analysis
save.directory <- '/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_analysisdata'

##############################################################
# functions used in script

################################################################################################################################################################
# Pt *** DIfferences in water quailty parameters between sites
# Question - how different are sites in terms of water quality parameters?
# show in histogram of DOC, TOC, pH, temp, [DO], SUVA, Br, TN concentrations

# File with all of the water quality parameters

waterquality <- read.csv((paste(save.directory, '/DBP_master_v6.csv', sep = "")), header = TRUE)

# read in absorbance indicies
spec.indicies <- read.csv("/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_correctedEEMs/DBPPreSpectralIndicies.csv",
                       header = TRUE)

#cut out Nas from spectral indicies file
ind <- apply(spec.indicies, 1, function(x) all(is.na(x)))   # function for removing Nas.. from all data
spec.indicies <- spec.indicies[ !ind, ]

############### Creating table of water quality parameters
# desire pH, DO, Water temp, Br, F, TN, DOC, SUVA, EC, nitrate, TSS: max, min, average, stdev

### Data cleaning and manipulation
# Create Sample ID column so that the absorbance and water quality data can be compared
# create vector with just sample number from water quality vector
sample <- sapply(strsplit(as.character(waterquality$DBPCode), split='_', fixed=TRUE), function(x) (x[2]))
samplename <- paste("DBP", sample, sep = "")
samplename <- str_pad(sample, 4, pad = "0")

waterquality$samplename <- paste("DBP", samplename, sep = "")
remove(samplename, sample)

# Variables that don't need data cleaning: pH, DO, Water temp, TN, DOC, EC

# Variables that do need cleaning
# - Br - need to change entries that are below detection levels 
waterquality$Br[waterquality$Br < 0.02] <- 0 # replace all values below 0.02 mg/L with 0 = below detection limits

# - F  - need to change entries below detection limits
waterquality$F[waterquality$F < 0.02] <- 0 # replace all values below 0.02 mg/L with 0 = below detection limits

# - NO3 - need to change entries before detection limits
waterquality$NO3[waterquality$NO3 < 0.02] <- 0 # replace all values below 0.02 mg/L with 0 = below detection limits

# - SUVA - need to calculate SUVA
# find the specific absorbance at 254 nm from absorbance spectra
sac254 <- cbind(as.character(spec.indicies$samplename), spec.indicies$abs254) 
NPOC <- cbind(waterquality$samplename, waterquality$NPOC_DOC_corrected)

column <- c("samplename", "abs254")
colnames(sac254) <- column
column <- c("samplename", "NPOC")
colnames(NPOC) <- column

sac254 <- data.frame(merge(sac254, NPOC, by = "samplename", all = TRUE))
remove(column, NPOC) #remove extra variables

# calculate SUVA as sac254/[DOC]
sac254$SUVA <- as.numeric(as.character(sac254$abs254))/as.numeric(as.character(sac254$NPOC))

# bind back inro the waterquality dataframe, according to sample ID
waterquality <- merge(waterquality, sac254, by = "samplename", all = TRUE)

# - TSS - Ensure that data that is negative for TSS is 0 (where 0 = below detection limit)
waterquality$TSS_mgL[waterquality$TSS_mgL< 0] <- 0  

##############
# Find stats (average, max, min, stdev) for water quality parameters
# Make a subgrouping just with the parameters that you want to simplify
quality.stats <- data.frame(waterquality$samplename, waterquality$pH, waterquality$EC_mScm,
                       waterquality$DO_mgL, waterquality$Water.Temp, waterquality$TSS_mgL,
                       waterquality$NPOC_DOC_corrected, waterquality$TN, waterquality$SUVA, waterquality$F, 
                       waterquality$Br, waterquality$NO3
                       )

# Run the functions length, mean, and sd on the value of "change" for each group, 
# broken down by sex + condition
wq.mean <- apply(quality.stats[,2:12], 2, mean, na.rm = TRUE)
wq.max <- apply(quality.stats[,2:12], 2, max, na.rm = TRUE)
wq.min <- apply(quality.stats[,2:12], 2, min, na.rm = TRUE)
wq.sd <- apply(quality.stats[,2:12], 2, sd, na.rm = TRUE)
#wq.length <- apply(quality.stats[,2:12], 2, length, na.rm = TRUE)
#wq.se <- wq.sd/wq.length

wq.stat <- rbind(wq.mean, wq.max, wq.min, wq.sd)  # bind stats back together
write.table(wq.stat, file = paste(save.directory, "waterqualitystats.csv", sep = ""),  sep = ",")

########################
# Heat map - do heat map to look at clustering of data according to water quality and spectral parameters
# Create a dataframe with the above water quality parameters and spectral parameters
colnames(quality.stats) <- c("samplename", "pH", "EC", "DO", "Temperature", "TSS", "DOC", "TN", "SUVA", "F", "Br", "NO3")
wq.heat <- merge(spec.indicies, quality.stats, by = 'samplename', all = TRUE)

# calculate the percent above/below the mean value that the value is for that parameter
# Calculated as (value - ave)/ave *100

remove(wq.permean)

n = dim(wq.heat)[2]
for (i in 2:n){
  m <- mean(wq.heat[,i], na.rm = TRUE)
  percent <- data.frame(apply(data.frame(wq.heat[,i]), 1, function(x) ((x-m)/m)*100))
  
  # if the merged dataset  exists, append to it by row
  if (exists("wq.permean")){
    wq.permean <- cbind(wq.permean, percent)
  }
  
  # if the merged dataset doesn't exist, create it
  if (!exists("wq.permean")){
    wq.permean <- percent
  }
}

colnames(wq.permean) <- colnames(wq.heat)[2:28] # assign column names 
row.names(wq.permean) <- wq.heat$samplename # assign row names as sample ID

# Arrange data for heat map
wq.permean[wq.permean > 200] <- NA # replace values that are greater than 200% with Na's - doesn't make sense
wq.permean[wq.permean < -200] <- NA # replace values that are less than -200% with Na's - doesn't make sense

# get rid of Zlonoskia area column... very littel data!
wq.permean$HIX_Zsonlay_area = NULL
wq.permean$HIX_Zsonlay_sum = NULL
wq.permean$F = NULL
wq.permean$Br = NULL

rnames <- row.names(wq.permean)
mat.data <- data.matrix(t(wq.permean[1:117,]))          # convert data to matrix

# do heat map
# create colour palette
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

# creates a 5 x 5 inch image
png(paste(save.directory, "DBP_WQ_heatmap.png", sep = ""),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          # Change the data within the heat map boxes
          #cellnote = mat.data,  # same data set for cell labels
          #notecex = 0.8,          # Change the font size of the data labels
          #notecol="black",      # change font color of cell labels to black
          
          # labels
          main = "Patterns in Water Quality Parameters - DBP", # heat map title
          
          # dendorgram and groupings
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=c("row"),     # only draw a row dendrogram
          density.info="histogram",  # turns on density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          
          # appearance
          margins =c(8,15),     # widens margins around plot
          col= my_palette,       # use on color palette defined earlier 
          cexCol=1.5, 
          cexRow = 1.5,          # decrease row font size to fit
          srtCol=45,           # rotate the x labels at 45 deg so they fit
          #axisnames = FALSE,
          
          na.color = 'white',   # colour of NA blocks
          keysize = 1,          # size of the colour key
          Rowv = TRUE,
          Colv = TRUE)            # turn on column clustering
dev.off()

################################################################################################################################################################
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
# See eaxmple at http://www.r-bloggers.com/computing-and-visualizing-pca-in-r/
library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)

g <- ggbiplot(pca.pre, obs.scale = 1, var.scale = 1)
#, obs.scale = 1, var.scale = 1) # From initial code
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
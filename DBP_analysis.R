# Script for analyzing DBP data
# Created 29sept2015
# For analysing data 
# Ashlee Jollymore's PhD
# DBP project

# References:
# http://www.sthda.com/english/wiki/principal-component-analysis-how-to-reveal-the-most-important-variables-in-your-data-r-software-and-data-mining#at_pco=smlwn-1.0&at_si=563bc26c64fc73c8&at_ab=per-2&at_pos=0&at_tot=1
# https://cran.r-project.org/web/packages/EEM/vignettes/vignette.html
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
library('FactoMineR')
library('nlme')
library(devtools)
install_github("ggbiplot", "vqv")
library(ggbiplot)
library("factoextra")
library("EEM")
devtools::install_github("PMassicotte/eemR")

################################################################################
# Read in data used for analysis
## Water quality and location data
# directory for data used for analysis
save.directory <- '/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_analysisdata'
pre.directory <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_correctedEEMSRaleigh"
post.directory <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_correctedEEMSRaleigh"
##############################################################
# functions used in script

make.eem <- function (ex, PCAcomponents){
  n <- length(unique(ex))
  PCA.var <-  data.frame(matrix(vector(), as.numeric(length(unique(PCAcomponents$em))), length(unique(PCAcomponents$ex))))
  
  for (i in 1:n){
    ex.temp <- unique(ex)[i]
    test <- data.frame(subset(PCAcomponents, PCAcomponents$ex == ex.temp))
    row.names(test) <- unique(PCAcomponents$em)
    colnames(test)[1] <- ex.temp
    
    PCA.var[,i] <- test[,1]
  }
  
  row.names(PCA.var) <- unique(PCAcomponents$em)
  colnames(PCA.var) <- unique(PCAcomponents$ex)
  return(PCA.var)
}

################################################################################################################################################################
# Pt 1 Differences in water quailty parameters between sites
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

# Exclude specific variables (due to little data, etc)
wq.permean$HIX_Zsonlay_area = NULL
wq.permean$HIX_Zsonlay_sum = NULL
wq.permean$F = NULL
wq.permean$Br = NULL

rnames <- row.names(wq.permean)
mat.data <- data.matrix(t(wq.permean[1:117,]))          # convert data to matrix

# do heat map
# create colour palette
my_palette <- colorRampPalette(c("red", "yellow", "blue"))(n = 299)

# creates a 5 x 5 inch image
png(paste(save.directory, "/DBP_WQ_heatmap.png", sep = ""),    # create PNG for the heat map        
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

########################## Heat map of quartile data
# Aim is to do the same heat map as above, except on where the value falls from the mean
# show where the individual value falls on quartile
# reference - http://www.r-bloggers.com/quartiles-deciles-and-percentiles/
# Cumulative distribution - ecdf function in R

wq.cdf <- as.data.frame(apply(wq.heat[,2:28], 2, function(x) ecdf(x)(x))) # calculate CDF for variables

# Do a heat map of the CDF from water quality parameters
# Exclude specific variables (due to little data, etc)
wq.cdf$HIX_Zsonlay_area = NULL
wq.cdf$HIX_Zsonlay_sum = NULL
wq.cdf$F = NULL
wq.cdf$Br = NULL

mat.data <- data.matrix(t(wq.cdf[1:117,]))          # convert data to matrix
colnames(mat.data) <- wq.heat[1:117,1]              # add column names - sample ID
  
# do heat map
# create colour palette
my_palette <- colorRampPalette(c("light blue", "dark blue"))(n = 299)

# creates a 5 x 5 inch image
png(paste(save.directory, "/DBP_WQ_heatmap_CDF.png", sep = ""),    # create PNG for the heat map        
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
          main = "Patterns in Water Quality Parameters - CDF DBP", # heat map title
          
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

######################## Examining EEM Data
## PCA on raw EEMS 
# Aim of this is to see regions that explain most of the variance within the pre chlorinated EEMS
# Examine all data lumped together to see any alterations 
#####################################################################
## Do PCA  on the compiled pre-chlorinated data
###############################
# Prechlorinated EEMS
# locate the prechlorinated corrected eems within the file
setwd(pre.directory) 
filelist_DBPpre <- list.files(pattern = "_Raleighcorr.csv$")

# call function to compile as a PCA object
setwd("/Users/user/SpecScripts") 
source("PCAfilecomp_function.R")
PCA.EEMpre <- PCA.eem(filelist = filelist_DBPpre, directory = pre.directory) 
PCA.pre <- PCA.EEMpre
# add in a variable that organizes according to which watersheds are drinking, which are protected, etc..
# To cluster.. see if there is a pattern within watersheds that are protected

# Preprocessing - 
# mean centering - done by center = TRUE in pca
# mean scaling - done by scale = TRUE in pca

# Normalization - according to 'eem' package in r
normalize


# read in file containing pre chlor EEMs assembled for PCA analysis
# perofrm PCA analysis on all pre chlorinated EEMS
pca.pre <- prcomp(PCA.pre, center = TRUE, scale. = TRUE)

# plot the loadings
plotLoading(pca.pre, ncomp = 2)

# plot the scores from the PCA - EEM package
plotScore

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

g <- ggbiplot(pca.pre, obs.scale = 1, var.scale = 1)
#, obs.scale = 1, var.scale = 1) # From initial code
#groups = ir.species, ellipse = TRUE, 
#circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')

############ find the areas that explain the most variance
pca.pre.Facto <-  PCA(PCA.pre, graph = FALSE)

# plot the contribution of each wavelength to the first 5 components to the PCA
PCA.pre.c <- data.frame(pca.pre.Facto$var$contrib)

# extract em and ex wavelengths from dataset (rownames)
PCA.pre.c$ex = sapply(strsplit(as.character(row.names(PCA.pre.c)), split='_', fixed=TRUE), function(x) (x[1])) # excitation wavelengths
PCA.pre.c$em <- as.numeric(sapply(strsplit(as.character(row.names(PCA.pre.c)), split='_', fixed=TRUE), function(x) (x[2]))) #emission wavelenghts

# use function to make a 3-D eem out of the dimentions
dim1.PCApre <- make.eem(ex = unique(PCA.pre.c$ex), PCAcomponents = data.frame(PCA.pre.c[,c(1,6,7)]))
dim2.PCApre <- make.eem(ex = unique(PCA.pre.c$ex), PCAcomponents = data.frame(PCA.pre.c[,c(2,6,7)]))
dim3.PCApre <- make.eem(ex = unique(PCA.pre.c$ex), PCAcomponents = data.frame(PCA.pre.c[,c(3,6,7)]))
dim4.PCApre <- make.eem(ex = unique(PCA.pre.c$ex), PCAcomponents = data.frame(PCA.pre.c[,c(4,6,7)]))
dim5.PCApre <- make.eem(ex = unique(PCA.pre.c$ex), PCAcomponents = data.frame(PCA.pre.c[,c(5,6,7)]))

# plot as contour plots using function
setwd("/Users/user/SpecScripts") 
source("EEM_contour_v1.R")

#variables to change
xlimit <- range(300, 700, finite=TRUE)
ylimit <- range(240, 800, finite = TRUE)

numcont = 100 # number of contour levels you want: Change if you want
#Plot contours and save in correction file

explot = seq(240, 800, by = 2)
emplot = as.numeric(row.names(dim1.PCApre))

########### specific to first 
plotpath <- file.path(save.directory, "PCApre_dim1_Contour.jpeg")
zmax = max(dim1.PCApre,na.rm=TRUE) # put the max intensity of that you want to graph
zmin = min(dim1.PCApre,na.rm=TRUE)

jpeg(file=plotpath)
contour.plots(eems = as.matrix(dim1.PCApre), Title = "PCA- Pre First Component", ex = explot, em = emplot,
              zmax,zmin,numcont)  
dev.off()

########## specific to second 
plotpath <- file.path(save.directory, "PCApre_dim2_Contour.jpeg")
zmax = max(dim2.PCApre,na.rm=TRUE) # put the max intensity of that you want to graph
zmin = min(dim2.PCApre,na.rm=TRUE)

jpeg(file=plotpath)
contour.plots(eems = as.matrix(dim2.PCApre), Title = "PCA- Pre Component2", ex = explot, em = emplot,
              zmax,zmin,numcont)  
dev.off()
########## specific to third 
plotpath <- file.path(save.directory, "PCApre_dim3_Contour.jpeg")
zmax = max(dim3.PCApre,na.rm=TRUE) # put the max intensity of that you want to graph
zmin = min(dim3.PCApre,na.rm=TRUE)

jpeg(file=plotpath)
contour.plots(eems = as.matrix(dim3.PCApre), Title = "PCA- Pre Component 3", ex = explot, em = emplot,
              zmax,zmin,numcont)  
dev.off()
########## specific to second 
plotpath <- file.path(save.directory, "PCApre_dim4_Contour.jpeg")
zmax = max(dim4.PCApre,na.rm=TRUE) # put the max intensity of that you want to graph
zmin = min(dim4.PCApre,na.rm=TRUE)

jpeg(file=plotpath)
contour.plots(eems = as.matrix(dim4.PCApre), Title = "PCA- Pre Component 4", ex = explot, em = emplot,
              zmax,zmin,numcont)  
dev.off()
########## specific to 5 
plotpath <- file.path(save.directory, "PCApre_dim5_Contour.jpeg")
zmax = max(dim5.PCApre,na.rm=TRUE) # put the max intensity of that you want to graph
zmin = min(dim5.PCApre,na.rm=TRUE)

jpeg(file=plotpath)
contour.plots(eems = as.matrix(dim5.PCApre), Title = "PCA- Pre Component 5", ex = explot, em = emplot,
              zmax,zmin,numcont)  
dev.off()

######### Examining 6 - Component PARAFAC fit to Prechlorinated EEMS
# Compare PCA model for the custom 6 component model
# Import Fmax values for 6 component model 
Pre.Fmax <- read.csv("/Users/user/Documents/MATLAB/toolbox/PARAFACresults/DOMFluor_DREEMSHydbrid/DBPPre/percentloadings.csv", header = FALSE, sep = ",")

colnames(Pre.Fmax) <- c("samplename", "Pre_C1", "Pre_C2", "Pre_C3", "Pre_C4", "Pre_C5", "Pre_C6")

# change the sample name column to reflect the pre chlronation format - important for later merging
sample <- sapply(strsplit(as.character(Pre.Fmax$samplename), split='_', fixed=TRUE), function(x) (x[1]))
samplename <- paste("DBP", sample, sep = "")
samplename <- str_pad(sample, 4, pad = "0")

Pre.Fmax$samplename <- samplename

# principal component analysis on PARAFAC results
PCA.pre.6 <- prcomp(Pre.Fmax[,2:7], center = TRUE, scale. = TRUE)

# plot PCA
png(paste(save.directory, "DBP_PCA6comppre.png", sep = ""),    # create graphic for the         
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

ggbiplot(PCA.pre.6, obs.scale = 1, var.scale = 1, 
         ellipse = FALSE, circle = FALSE)+ggtitle("PCA Results- 6 Component PARAFAC Model Pre-Chlorination")
dev.off()

summary(PCA.pre.6)

##### Variation within Components - prechlorination boxplots
# assemble data with component in one column and FMax in another

remove(PCApre6)

for (i in 1:6){
  temp.C <- data.frame(Pre.Fmax[,i+1])
  temp.component <- colnames(Pre.Fmax)[i+1]
  temp.C$component <- temp.component
  
  # if the merged dataset  exists, append to it by row
  if (exists("PCApre6")){
    PCApre6 <- rbind(PCApre6, temp.C)
  }
  
  # if the merged dataset doesn't exist, create it
  if (!exists("PCApre6")){
    PCApre6 <- temp.C
  }
}

colnames(PCApre6)[1] <- "values"

ggplot(PCApre6, aes(x=component, y=values, fill=component)) + geom_boxplot() 

##### Looking at overall variation within all water quality parameters
# Which parameters explain the greatest degree of variation within the dataset?
# Compile PARAFAC model fits with water quality parameters
wq.all <- merge(Pre.Fmax, wq.heat, by = "samplename", all = TRUE) 

# do PCA on all variables - which variables explain the greatest degree of variation?
wq.all.pca <- prcomp(na.omit(wq.all[,2:34]), center = TRUE, scale. = TRUE)
summary(wq.all.pca)

# plot PCA results
# plot PCA
png(paste(save.directory, "/DBP_PCAWQdata.png", sep = ""),    # create graphic for the         
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

ggbiplot(wq.all.pca, obs.scale = 1, var.scale = 1, 
         ellipse = FALSE, circle = FALSE)+ggtitle("PCA Results- WQ Parameters")
dev.off()

# show relative contribution of each variable to the first 5 components of PCA
wq.pca <- PCA(na.omit(wq.all[,2:34]), graph = FALSE)
head(wq.pca$var$contrib)

# plot contribution to first 2 PCA components
fviz_pca_contrib(wq.pca, choice = "var", axes = 1)
fviz_pca_contrib(wq.pca, choice = "var", axes = 2)

################################################################################################################################################################
# Pt 2 - how does chlorination change the spectral composition of EEMS?
# princip-al component analysis to look at how chlorination altered spectral characteristics (EEMS)
# References:
# http://www.r-bloggers.com/computing-and-visualizing-pca-in-r/
# http://planspace.org/2013/02/03/pca-3d-visualization-and-clustering-in-r/
###########################################################
# Do PCA on all of the PARAFAC models results (custom model) - pre, post, pre+post, delta
# post chlorination EEMs - 2 component model
DBPpost <- read.csv("/Users/user/Documents/MATLAB/toolbox/PARAFACresults/DOMFluor_DREEMSHydbrid/DBPPost/percentloadings.csv", header = FALSE, sep = ",")
colnames(DBPpost) <- c("samplename","","Post_C1", "Post_C2")
DBPpost[,2] <- NULL # get rid of empty second column

# get rid of chlor in sample name
DBPpost$samplename <- gsub("Chlor","",DBPpost$samplename)

# 3 component PARAFAC delta EEMS model
DBPdelta <- read.csv("/Users/user/Documents/MATLAB/toolbox/PARAFACresults/DOMFluor_DREEMSHydbrid/DBPdelta/percentloadings.csv", header = FALSE, sep = ",")
colnames(DBPdelta) <- c("samplename","","Delta_C1", "Delta_C2", "Delta_C3")

DBPdelta[,2] <- NULL # get rid of empty second column

#preandpost - 6 component model
DBPprepost <- read.csv("/Users/user/Documents/MATLAB/toolbox/PARAFACresults/DOMFluor_DREEMSHydbrid/DBPprepost/percentloadings.csv", header = FALSE, sep = ",")

colnames(DBPprepost) <- c("samplename","Post6_C1", "Post6_C2", 
                          "Post6_C3", "Post6_C4",
                          "Post6_C5", "Post6_C6")

# take only the post chlorination model out
DBPprepost$prepost <- sapply(strsplit(as.character(DBPprepost$samplename), split='_', fixed=TRUE), function(x) (x[2]))

# take only post chlorinated eems
DBPprepost.post <- subset(DBPprepost, DBPprepost$prepost == "DBPPost")
DBPprepost.post$prepost <- NULL # get rid of column

DBPprepost.post$samplename <- sapply(strsplit(as.character(DBPprepost.post$samplename), split='_', fixed=TRUE), function(x) (x[1]))
DBPprepost.post$samplename <-  gsub("Chlor","",DBPprepost.post$samplename)

# compile all together according to sample ID
PCA.PARAFAC.all <- data.frame(Reduce(function(x,y) merge(x,y, by = "samplename", all = FALSE), 
                          list(Pre.Fmax,DBPpost,DBPdelta, DBPprepost.post)))

# do PCA on all variables - which variables explain the greatest degree of variation?
PCA.PARAFAC <- prcomp(na.omit(PCA.PARAFAC.all[,2:18]), center = TRUE, scale. = TRUE)
summary(PCA.PARAFAC)

# plot PCA results on all PARAFAC data
png(paste(save.directory, "/DBP_PCAPARAFACdata.png", sep = ""),    # create graphic for the         
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

ggbiplot(PCA.PARAFAC, obs.scale = 1, var.scale = 1, 
         ellipse = FALSE, circle = FALSE)+ggtitle("PCA Results- PARAFAC Parameters")
dev.off()

# show relative contribution of each variable to the first 5 components of PCA
PCA.PARAFAC <- PCA(na.omit(PCA.PARAFAC.all[,2:18]), graph = TRUE)
head(PCA.PARAFAC$var$contrib)

# plot contribution to first 2 PCA components
fviz_pca_contrib(PCA.PARAFAC, choice = "var", axes = 1)
fviz_pca_contrib(PCA.PARAFAC, choice = "var", axes = 2)
#####################################
# PCA on pre and post - 6 Component model
# add in column for pre versus post chlorination
DBPprepost$prepost <- gsub("DBP", "", DBPprepost$prepost)
  
# Do PCA
PCA.Prepost6 <- prcomp(na.omit(DBPprepost[,2:7]), center = TRUE, scale. = TRUE)
summary(PCA.Prepost6)

# plot PCA results on all PARAFAC data
png(paste(save.directory, "/PCA_prepost6comp.png", sep = ""),    # create graphic for the         
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

ggbiplot(PCA.Prepost6, obs.scale = 1, var.scale = 1, 
         groups = DBPprepost$prepost, ellipse = TRUE, circle = FALSE)+ggtitle("PCA - 6 Comp PARAFAC Pre vs Post")
dev.off()

#####################################
# PCA on corrected EEMS
### Post chlor eems - compile Raleigh, IFE and Raman corrected files for PCA
setwd(post.directory) 
filelist_DBPpost <- list.files(pattern = "_Raleighcorr.csv$")

# call function to compile as a PCA object
setwd("/Users/user/SpecScripts") 
source("PCAfilecomp_function.R")
PCA.EEMpost <- PCA.eem(filelist = filelist_DBPpost, directory = post.directory) 


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

################################################################################################################################################################  
# Comparing pre to post chlorinated EEMS
#### 
# first add column to sample ID that separates unchlorinated or chlorinated EEMS
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
# Will try doing PCA using FactoMineR package - lets you choose the number of variables..

pca.all <- prcomp(PCA.all[,2:140501], center = TRUE, scale. = TRUE)
pca.all.Facto <-  PCA(PCA.all[,2:140501], graph = FALSE)
summary(pca.all)

plot(pca.all, type = "l")

pcaall.loadings <- pca.all$rotation[,1:4]
plot(pcaall.loadings[,1:2], type = 'p')
plot(pcaall.loadings[,3:4], type = 'p')

# PLot PCA. This plot is super complex as it is taking all of the wavelengths.
g <- ggbiplot(pca.all, obs.scale = 1, var.scale = 1, 
              groups = PCA.all[,1], ellipse = FALSE, circle = FALSE)

# To make plotting nicer, choose the variables that have the greatest contribution to the PCA (otherwise your plot is a mess)
head(pca.all.Facto$var$contrib)   # The larger the value of the contribution, the more the variable contributes to the component.

# write to a csv file for examination
write.table(pca.all.Facto$var$contrib, file = paste(save.directory, "/PrePost_PCAvariables.csv", sep = ""), sep = ",")

############### Plot the contribution of each wavelength for most important components from the PCA as a contour plot
# plot the contribution of each wavelength to the first *** components to the PCA
PCA.c <- data.frame(pca.all.Facto$var$contrib)

#extract em and ex wavelengths from dataset (rownames)
PCA.c$ex = unique(sapply(strsplit(as.character(row.names(PCA.c)), split='_', fixed=TRUE), function(x) (x[1]))) # excitation wavelengths
PCA.c$em <- unique(as.numeric(sapply(strsplit(as.character(row.names(PCA.c)), split='_', fixed=TRUE), function(x) (x[2])))) #emission wavelenghts

####
dim1.PCAprepost <- make.eem(ex = unique(PCA.c$ex), PCAcomponents = data.frame(PCA.c[,c(1,6,7)]))
dim2.PCAprepost <- make.eem(ex = unique(PCA.c$ex), PCAcomponents = data.frame(PCA.c[,c(2,6,7)]))
dim3.PCAprepost <- make.eem(ex = unique(PCA.c$ex), PCAcomponents = data.frame(PCA.c[,c(3,6,7)]))
dim4.PCAprepost <- make.eem(ex = unique(PCA.c$ex), PCAcomponents = data.frame(PCA.c[,c(4,6,7)]))
dim5.PCAprepost <- make.eem(ex = unique(PCA.c$ex), PCAcomponents = data.frame(PCA.c[,c(5,6,7)]))

# plot as contour plots using function
setwd("/Users/user/SpecScripts") 
source("EEM_contour_v1.R")

#variables to change
xlimit <- range(300, 700, finite=TRUE)
ylimit <- range(240, 800, finite = TRUE)

numcont = 20 # number of contour levels you want: Change if you want
#Plot contours and save in correction file

explot = seq(240, 800, by = 2)
emplot = as.numeric(row.names(dim1.PCAprepost))

########### specific to first 
plotpath <- file.path(save.directory, "PCAprepost_dim1_Contour.jpeg")
zmax = max(dim1.PCAprepost,na.rm=TRUE) # put the max intensity of that you want to graph
zmin = min(dim1.PCAprepost,na.rm=TRUE)

jpeg(file=plotpath)
contour.plots(eems = as.matrix(dim1.PCAprepost), Title = "PCA- Pre and Post First Component", ex = explot, em = emplot,
              zmax,zmin,numcont)  
dev.off()

########## specific to second 
plotpath <- file.path(save.directory, "PCAprepost_dim2_Contour.jpeg")
zmax = max(dim2.PCAprepost,na.rm=TRUE) # put the max intensity of that you want to graph
zmin = min(dim2.PCAprepost,na.rm=TRUE)

jpeg(file=plotpath)
contour.plots(eems = as.matrix(dim2.PCAprepost), Title = "PCA- Pre and Post Component2", ex = explot, em = emplot,
              zmax,zmin,numcont)  
dev.off()
########## specific to third 
plotpath <- file.path(save.directory, "PCAprepost_dim3_Contour.jpeg")
zmax = max(dim3.PCAprepost,na.rm=TRUE) # put the max intensity of that you want to graph
zmin = min(dim3.PCAprepost,na.rm=TRUE)

jpeg(file=plotpath)
contour.plots(eems = as.matrix(dim3.PCAprepost), Title = "PCA- Pre and Post Component 3", ex = explot, em = emplot,
              zmax,zmin,numcont)  
dev.off()
########## specific to second 
plotpath <- file.path(save.directory, "PCAprepost_dim4_Contour.jpeg")
zmax = max(dim4.PCAprepost,na.rm=TRUE) # put the max intensity of that you want to graph
zmin = min(dim4.PCAprepost,na.rm=TRUE)

jpeg(file=plotpath)
contour.plots(eems = as.matrix(dim4.PCAprepost), Title = "PCA- Pre and Post Component 4", ex = explot, em = emplot,
              zmax,zmin,numcont)  
dev.off()
########## specific to 5 
plotpath <- file.path(save.directory, "PCAprepost_dim5_Contour.jpeg")
zmax = max(dim5.PCAprepost,na.rm=TRUE) # put the max intensity of that you want to graph
zmin = min(dim5.PCAprepost,na.rm=TRUE)

jpeg(file=plotpath)
contour.plots(eems = as.matrix(dim5.PCAprepost), Title = "PCA- Pre and Post Component 5", ex = explot, em = emplot,
              zmax,zmin,numcont)  
dev.off()
################################################################################################
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

# change the sample name column to reflect the pre chlorination format - important for later merging
sample <- sapply(strsplit(as.character(CM.pre$sample.ID), split='DBPPre', fixed=TRUE), function(x) (x[1]))
samplename <- str_pad(sample, 4, pad = "0")

CM.pre$samplename <- samplename

####### post chlorination
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

summary(PCA.CM) # get stats for PCA on CM - pre and post chlorination

# choose top 2 components of variation, and plot against one another to see if you can cluster pre and post chlorination
PCA.CM$rotation[,1]

plot(PCA.CM$rotation[,1], xlab = "CM Component", ylab = "PCA Rotation", main = 'CM Pre and Post PCA - Component 1')

png(paste(save.directory, "DBP_CM_PCA.png", sep = ""),    # create graphic for the         
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

g <- ggbiplot(PCA.CM, obs.scale = 1, var.scale = 1, 
              groups = CMPCA.all$chlor, ellipse = FALSE, circle = FALSE)+ggtitle("PCA Results- 13 Component PARAFAC Model")

g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'vertical', 
               legend.position = 'right')
dev.off()

##########################
# Analyzing pre and post chlorinated eems with 6 component model
# load post chlorinated EEM fit to 6 component model
prepost.Fmax <- read.csv("/Users/user/Documents/MATLAB/toolbox/PARAFACresults/DOMFluor_DREEMSHydbrid/DBPprepost/percentloadings.csv", header = FALSE, sep = ",")

colnames(prepost.Fmax) <- c("samplename", "C1", "C2", "C3", "C4", "C5", "C6")

# bind column with chlorination status
prepost.Fmax <- cbind(prepost.Fmax, chlor)

# principal component analysis on PARAFAC results
PCA.prepost.6 <- prcomp(prepost.Fmax[,2:7], center = TRUE, scale. = TRUE)

# plot PCA
g <- ggbiplot(PCA.prepost.6, obs.scale = 1, var.scale = 1, 
              groups = prepost.Fmax$chlor, ellipse = FALSE, circle = FALSE)+ggtitle("PCA Results- 6 Component PARAFAC Model Pre+Post Chlorination")

summary(PCA.prepost.6)

# find the components that explain the most variation
pca.prepost6com.Facto <-  PCA(prepost.Fmax[,2:7], graph = FALSE)
PCA.comp.6c <- data.frame(pca.prepost6com.Facto$var$contrib)
x = c("C1", "C2", "C3", "C4", "C5", "C6")
x = seq(1,6, by = 1)
# plot variable contributions of each component
g <- ggplot(PCA.comp.6c, aes(x, y = value, col = 'PCA Dim')) + 
  geom_line(aes(y = PCA.comp.6c$Dim.1, col = "Dim 1")) + 
  geom_line(aes(y = PCA.comp.6c$Dim.2, col = "Dim 2")) + 
  geom_line(aes(y = PCA.comp.6c$Dim.3, col = "Dim 3")) +
  geom_line(aes(y = PCA.comp.6c$Dim.4, col = "Dim 4")) +
  geom_line(aes(y = PCA.comp.6c$Dim.5, col = "Dim 5")) +
  xlab("PARAFAC Components")+
  ylab("Percent Contribution")+
  ggtitle("Contribution to PCA")

######################
# Exploring variation in the pre and post chlorinated EEMS - heat map (clustering)
# prepare prepost.Fmax dataframe for heat map - clustering chlorination status on x axis, and the components on the y

# add column = 1 for prechlor, 2 for post chlor

prepost.Fmax$numchlor <- ifelse(prepost.Fmax$chlor == "pre chlorination",1, 2)

rnames <- prepost.Fmax$numchlor
mat.data <- data.matrix(t(prepost.Fmax[,2:7]))          # convert data to matrix
colnames(mat.data) <- rnames

# do heat map
# create colour palette
my_palette <- colorRampPalette(c("lightblue", "darkblue"))(n = 299)

# creates a 5 x 5 inch image
png(paste(save.directory, "/DBP_PARAFAC6com_heatmap.png", sep = ""),    # create PNG for the heat map        
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
          main = "6 Component PARAFAC - Pre versus Post Chlorination", # heat map title
          
          # dendorgram and groupings
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=c("row"),     # only draw a row dendrogram
          density.info="none",  # turns on density plot inside color legend
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
          Colv = FALSE)            # turn of column clustering
dev.off()

######################
# take 2 on heat map - take average of each of the components within the pre component, and calculate the difference in the post chlorination EEMS

# calculate the average of each of the 6 components within the pre chloration dataset
Pre.ave <- t(data.frame(apply(Pre.Fmax[,2:7], 2, mean)))

# Calculate difference from average component for the post chlorination EEM
# (post - preave)/preave*100
post.6comp <- subset(prepost.Fmax, prepost.Fmax$numchlor==2)

n = dim(post.6comp)[2]-2
for (i in 2:n){
  m <- Pre.ave[,(i-1)]
  percent <- data.frame(apply(data.frame(post.6comp[,i]), 1, function(x) (x-m)))
  
  # add to post.6comp dataframe
  post.6comp[,(i+n+1)] <- percent
  }

## do heat map

mat.data <- data.matrix(t(post.6comp[,10:15]))          # convert data to matrix
row.names(mat.data) <- x
colnames(mat.data) <- post.6comp$samplename

# do heat map
# create colour palette
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

# creates a 5 x 5 inch image
png(paste(save.directory, "/DBP_prepostdiff6comp_heatmap.png", sep = ""),    # create PNG for the heat map        
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
          main = "6 Component PARAFAC - Ave Pre - Post Chlorination", # heat map title
          
          # dendorgram and groupings
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=c("row"),     # only draw a row dendrogram
          density.info="none",  # turns on density plot inside color legend
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
          Colv = FALSE)            # turn of column clustering
dev.off()

######################
# using clustering model to predict which PARAFAC components best predict chlorination status
# Include PCA wavelengths from fit to pre and post models?

# http://stats.stackexchange.com/questions/30406/how-to-assess-predictive-power-of-set-of-categorical-predictors-of-a-binary-outc

# logistic regression

mod <- gls(prepost.Fmax$numchlor ~ prepost.Fmax$C1 +prepost.Fmax$C2+prepost.Fmax$C3+prepost.Fmax$C4+prepost.Fmax$C5+prepost.Fmax$C6, 
            correlation = corARMA(p = 1, q = 2))

#########################################################################################################################
# Pt 3 - How does EEMS and Water Quality parameters predict DBP formation?
# Thoughts:
# use machine learning algorythm to predict DBP parameters
# try multivraite linear model where
#[DBP] = theta0 + theta1x1+theta2x2
# where x = spectral and EEM characteristics, as well as water quality parameters (pH, anion concentrations, DOC)

# First, mutiple linear regression to rank the contribution of water quality variables to 

########### compile the DBP concentrations + clean the data
HAA <- as.data.frame(read.csv("/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_HAAData/HAAdata_analysis.csv", header = TRUE))
THM <- read.csv("/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_THMData/THMdata_analysis.csv", header = TRUE)

# calculate THM and HAA total (sum of all of the different species for that sample)
HAA$total <- apply(HAA[,2:10], 1, sum)
THM$total <- apply(THM[,2:5], 1, sum)

######## apply linear model, where y = HAA and THMs, and x = spectral variables...
# try #1: 
# y = total HAA and THMs (two models)
# x = PARAFAC components from the prechlorinated EEMS, water quality parameters, spectral parameters

# TO DO - add PCA components for pre chlorinated EEMs; add in the delta eems for the 6 component fit?

####### Bind spectral and water quality parameters with DBP data
waterquality.mod <- waterquality[,c(1,9,10,11,12,15,16:20,22,25,35:37)]
HAA.waterq <- Reduce(function(x, y) merge(x, y, by = 'samplename', all=FALSE), list(HAA, waterquality.mod[,c(1,15)],CM.pre, spec.indicies, Pre.Fmax))
THM.waterq <- Reduce(function(x, y) merge(x, y, by = 'samplename', all=FALSE), list(THM, waterquality.mod[,c(1,15)],CM.pre, spec.indicies, Pre.Fmax))

# get rid of columns you don't need
HAA.waterq$sample.ID <- NULL
THM.waterq$sample.ID <- NULL

# Do linear models for total HAA and total THMs
# Note that expect some of the variables to co-relate, thus use gls linear fit model
HAA.waterq <- lapply(HAA.waterq, as.numeric)
THM.waterq <- lapply(THM.waterq, as.numeric)

# use tree to look at the correlation between variables - first, CM model
library("tree")
HAA.total <- CM.model.HAA$total
HAA.CM <-lapply(CM.model.HAA[,2:17], as.numeric) # convert to numeric prior to running model

model.HAAtotal <- tree(HAA.total ~ ., data = HAA.CM) # run tree model on CM data
plot(model.HAAtotal) # plot tree model
text(model.HAAtotal) # put in text labels into the tree label

model.HAAtotal <- lm(HAA.waterq$total ~ HAA.waterq$NPOC + HAA.waterq$C1.x + HAA.waterq$C2.x +
                       HAA.waterq$C3.x + HAA.waterq$C4.x + HAA.waterq$C5.x +
                       HAA.waterq$C6.x + HAA.waterq$C7 + HAA.waterq$C8 + HAA.waterq$C9 + 
                       HAA.waterq$C10 + HAA.waterq$C11 + HAA.waterq$C12 + HAA.waterq$C13)
summary(model.HAAtotal)

# plot all of the 
plot(HAA.waterq$total ~ HAA.waterq$NPOC)
plot(HAA.waterq$total ~ HAA.waterq$C1.x)
plot(HAA.waterq$total ~ HAA.waterq$C2.x)
plot(HAA.waterq$total ~ HAA.waterq$C3.x)
plot(HAA.waterq$total ~ HAA.waterq$C4.x)
plot(HAA.waterq$total ~ HAA.waterq$C5.x)
plot(HAA.waterq$total ~ HAA.waterq$C6.x)
plot(HAA.waterq$total ~ HAA.waterq$C7)
plot(HAA.waterq$total ~ HAA.waterq$C8)
plot(HAA.waterq$total ~ HAA.waterq$C9)
plot(HAA.waterq$total ~ HAA.waterq$C10)
plot(HAA.waterq$total ~ HAA.waterq$C11)
plot(HAA.waterq$total ~ HAA.waterq$C12)
plot(HAA.waterq$total ~ HAA.waterq$C13)

## try custom PARAFAC components - linear model

model.HAAtotal.6comp <- lm(HAA.waterq$total ~ HAA.waterq$NPOC + HAA.waterq$C1.y + HAA.waterq$C2.y +
                       HAA.waterq$C3.y + HAA.waterq$C4.y + HAA.waterq$C5.y +
                       HAA.waterq$C6.y+ HAA.waterq$C1.x + HAA.waterq$C2.x +
                         HAA.waterq$C3.x + HAA.waterq$C4.x + HAA.waterq$C5.x +
                         HAA.waterq$C6.x + HAA.waterq$C7 + HAA.waterq$C8 + HAA.waterq$C9 + 
                         HAA.waterq$C10 + HAA.waterq$C11 + HAA.waterq$C12 + HAA.waterq$C13)

summary(model.HAAtotal.6comp)

NPOC.model <- lm(HAA.waterq$total ~ HAA.waterq$NPOC)
summary(NPOC.model)

####################### Total THMs
#NPOC and THMs
NPOC.model.THM <- lm(THM.waterq$total ~ THM.waterq$NPOC)
summary(NPOC.model.THM)
plot(THM.waterq$total ~ THM.waterq$NPOC) #note high outliers! are these green roofs? may have to remove


model.THMtotal<- lm(THM.waterq$total ~ THM.waterq$NPOC + THM.waterq$C1.y + THM.waterq$C2.y +
                             THM.waterq$C3.y + THM.waterq$C4.y + THM.waterq$C5.y +
                             THM.waterq$C6.y+ THM.waterq$C1.x + THM.waterq$C2.x +
                             THM.waterq$C3.x + THM.waterq$C4.x + THM.waterq$C5.x +
                             THM.waterq$C6.x + THM.waterq$C7 + THM.waterq$C8 + THM.waterq$C9 + 
                             THM.waterq$C10 + THM.waterq$C11 + THM.waterq$C12 + THM.waterq$C13)
summary(model.THMtotal)

# high NPOC is an issue? exponential model better?

#########################################
# try #2: 
# y = specific HAA and THMs (multiple models)
# x = PARAFAC components from the prechlorinated EEMS, water quality parameters, spectral parameters
# TO DO - add PCA components for pre chlorinated EEMs; add in the delta eems for the 6 component fit?


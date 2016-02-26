# Script for analyzing DBP data
# Created 29sept2015
# For analysing data 
# Ashlee Jollymore's PhD
# DBP project
#
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
library(ggplot2)
library(plyr)
library(eeptools)
library(MASS)
library(Hmisc)
library('corrplot') #package corrplot
################################################################################
# Read in data used for analysis
## Water quality and location data
# directory for data used for analysis
save.directory <- '/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_analysisdata'
pre.directory <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_correctedEEMSRaleigh"
post.directory <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_correctedEEMSRaleigh"
Delta.directory <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_delta"

# Cory McKnight fits
CM.pre.directory <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_pre_CM_PARAFAC"

########### colour blind colour palettes
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##############################################################
# functions used in script
make.eem <- function (ex, PCAcomponents){
  n <- length(unique(ex))
  PCA.var <-  data.frame(matrix(vector(), length(unique(PCAcomponents$em)), length(unique(ex))))
  row.names(PCA.var) <- unique(PCAcomponents$em)
  colnames(PCA.var) <- unique(PCAcomponents$ex)
  
  for (i in 1:n){
    ex.temp <- unique(ex)[i]
    test <- data.frame(subset(PCAcomponents, PCAcomponents$ex == ex.temp))
    row.names(test) <- unique(test$em)
    colnames(test)[1] <- ex.temp
    PCA.var[,i] <- test[,1]
  }
  row.names(PCA.var) <- unique(PCAcomponents$em)
  colnames(PCA.var) <- unique(PCAcomponents$ex)
  return(PCA.var)
}


####################################################################
# Reading in data - ALL DATA
# File with all of the water quality parameters
waterquality <- read.csv((paste(save.directory, '/DBP_master_v6.csv', sep = "")), header = TRUE)
# csv file that shows all of the watershed groupings
groupings <- read.csv((paste(save.directory, '/Watersheds_codes.csv', sep = "")), header = TRUE)

# merge the watershed codes to the water quality tables
waterquality <- merge(waterquality, groupings, by = "DBPCode")

# Create Sample ID column so that the absorbance and water quality data can be compared
# create vector with just sample number from water quality vector
sample <- sapply(strsplit(as.character(waterquality$DBPCode), split='_', fixed=TRUE), function(x) (x[2]))
samplename <- paste("DBP", sample, sep = "")
samplename <- str_pad(sample, 4, pad = "0")
waterquality$samplename <- paste("DBP", samplename, sep = "")
remove(samplename, sample)

# read in absorbance indicies - pre chlorination
spec.indicies <- read.csv(paste(pre.directory,'/DBPPreSpectralIndicies_corrABS.csv', sep = ""),header = TRUE)

# read in absorbance indicies - post chlorination
spec.indicies.post <- read.csv(paste(post.directory,'/DBPPostSpectralIndicies.csv', sep = ""),header = TRUE)

##########
# pre chlorination EEMs
# pre chlorination CM fits - CM directory
CM.pre <- read.csv(paste(CM.pre.directory, "/DBPpre_componentsandloadings_CM.csv", sep = "")) #percent from CM components
CM.pre.Fmax <- read.csv(paste(CM.pre.directory, "/DBPpre_componentsandloadings_CM_Fmax.csv", sep = "")) #percent from CM components

# change the sample name column to reflect the pre chlorination format - important for later merging
sample <- sapply(strsplit(as.character(CM.pre$sample.ID), split='DBPPre', fixed=TRUE), function(x) (x[1]))
samplename <- str_pad(sample, 4, pad = "0")
CM.pre$samplename <- samplename
CM.pre.Fmax$samplename <- samplename
remove(sample)

# pre chlorination PARAFAC - Fmax values for 6 component model 
Pre.Fmax <- read.csv("/Users/user/Documents/MATLAB/toolbox/PARAFACresults/DOMFluor_DREEMSHydbrid_withoutGR/DBPPre_withoutGR/Fmax.csv", header = FALSE, sep = ",")
Pre.Fmax.key <- t(read.csv("/Users/user/Documents/MATLAB/toolbox/CorrEEMS/DBPPre_noGR/01key.csv", header = FALSE, sep = ","))
Pre.Fmax$samplename <- Pre.Fmax.key
colnames(Pre.Fmax) <- c("DR_C1", "DR_C2", "DR_C3", "DR_C4", "DR_C5", "DR_C6", "samplename")

# calculate Fmax Percent - pre chlorination
Pre.Fmax.per <- cbind(Pre.Fmax[,1:6]/rowSums(Pre.Fmax[,1:6])*100, Pre.Fmax[,7])

#############
# pre and post chlorination CM data
prepost.CM.Fmax <- read.csv("/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_preandpost/DBP_prepost_CMresults/DBPprepost_componentsandloadings_CM_Fmax.csv", 
                         header = TRUE, sep = ",")
prepost.CM.per <- prepost.CM.Fmax <- read.csv("/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_preandpost/DBP_prepost_CMresults/DBPprepost_componentsandloadings_CM.csv", 
                                              header = TRUE, sep = ",")
##############
# Pre and post chlorination PARAFAC fits
prepost.6comp.Fmax <- read.csv("/Users/user/Documents/MATLAB/toolbox/PARAFACresults/DOMFluor_DREEMSHydbrid_withoutGR/DBPprepost_withoutGR/Fmax.csv", header = FALSE, sep = ",")
# get the sample ID from the key
prepost.key <- t(read.csv("/Users/user/Documents/MATLAB/toolbox/CorrEEMS/DBPprepost_noGR/01key.csv", header = FALSE, sep = ","))
prepost.6comp.Fmax$samplename <- prepost.key # add in the sample ID to the prepost Fmax file

# add in column names
colnames(prepost.6comp.Fmax) <- c("DR_C1", "DR_C2", "DR_C3", "DR_C4", "DR_C5", "DR_C6", "samplename")

# calculate Fmax Percent
prepost.Fmax.per <- cbind(prepost.6comp.Fmax[,1:6]/rowSums(prepost.6comp.Fmax[,1:6])*100, prepost.6comp.Fmax[,7])

#########
# DBP concentrations
HAA <- as.data.frame(read.csv("/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_HAAData/HAAdata_analysis.csv", header = TRUE))
THM <- read.csv("/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_THMData/THMdata_analysis.csv", header = TRUE)

#########
# Take out the green roof and rainwater harvest samples. These will be super-imposed on the analysis later
# samples = 56-63, 73-75, 91-93, 119+13
GR <- c("DBP0056", "DBP0057", "DBP0058", "DBP0060", "DBP0061", "DBP0062", 
        "DBP0073", "DBP0074", "DBP0075", "DBP0091", "DBP0092", "DBP0093", "DBP0119", "DBP0013")
# chlorinated
GRc <- c("DBPChlor0056", "DBPChlor0057", "DBPChlor0058", "DBPChlor0060", "DBPChlor0061", "DBPChlor0062", 
        "DBPChlor0073", "DBPChlor0074", "DBPChlor0075", "DBPChlor0091", "DBPChlor0092", "DBPChlor0093", "DBPChlor0119", "DBPChlor0013")

#remove from pre sample sets
waterquality <- waterquality[ !(waterquality$samplename  %in% GR), ] # from wq
spec.indicies <- spec.indicies[ !(spec.indicies$samplename  %in% GR), ] # 
Pre.Fmax <- Pre.Fmax[ !(Pre.Fmax$samplename  %in% GR), ] #
CM.pre <- CM.pre[ !(CM.pre$samplename  %in% GR), ] # 
PCA.EEMpre <- PCA.EEMpre[ !(rownames(PCA.EEMpre)  %in% GR), ] # 
# post
spec.indicies.post <- spec.indicies.post[ !(spec.indicies.post$samplename  %in% GRc), ] # 
CM.post  <- CM.post[ !(CM.post$samplename  %in% GRc), ] # 
DBPpost<- DBPpost[ !(DBPpost$samplename  %in% GRc), ] # 
PCA.EEMpost <- PCA.EEMpost[ !(rownames(PCA.EEMpost)  %in% GRc), ] # 
# Delta
PCA.EEMdelta <- PCA.EEMdelta[ !(PCA.EEMdelta$samplename  %in% GR), ] # 
DBPdelta <- DBPdelta[ !(DBPdelta$samplename  %in% GR), ] # 
spec.indicies.delta <- spec.indicies.delta[ !(spec.indicies.delta$samplename  %in% GR), ]
#pre and post chlor
DBPprepost <- DBPprepost[ !(DBPprepost$samplename  %in% rbind(GR, GRc)), ] # 
prepost.Fmax  <- prepost.Fmax[ !(prepost.Fmax$samplename  %in% GR), ] # remove prechlor
prepost.Fmax  <- prepost.Fmax[ !(prepost.Fmax$samplename  %in% GRc), ] # remove postchlor
HAA <- HAA[ !(HAA$samplename  %in% GR), ] # 
THM <- THM[ !(THM$samplename  %in% GR), ] # 
  
################################################################################################################################################################
# Pt 1 Differences in water quality parameters between sites
# Question - how different are sites in terms of water quality parameters?

# cut out Nas from spectral indicies file
ind <- apply(spec.indicies, 1, function(x) all(is.na(x)))   # function for removing Nas.. from all data
spec.indicies <- spec.indicies[ !ind, ]

ind <- apply(spec.indicies.post, 1, function(x) all(is.na(x)))   # function for removing Nas.. from all data
spec.indicies.post <- spec.indicies.post[ !ind, ]

############### 
# Table 1- Creating table of water quality parameters
# desire pH, DO, Water temp, Br, F, TN, DOC, SUVA, EC, nitrate, TSS: max, min, average, stdev

### Data cleaning and manipulation
# Variables that don't need data cleaning: pH, DO, Water temp, TN, DOC, EC
# Variables that do need cleaning
# - Br - need to change entries that are below detection levels 
waterquality$Br[waterquality$Br < 0.02] <- 0 # replace all values below 0.02 mg/L with 0 = below detection limits

# - F  - need to change entries below detection limits
waterquality$F[waterquality$F < 0.02] <- 0 # replace all values below 0.02 mg/L with 0 = below detection limits

# - NO3 - need to change entries before detection limits
waterquality$NO3[waterquality$NO3 < 0.02] <- 0 # replace all values below 0.02 mg/L with 0 = below detection limits

# - TSS - Ensure that data that is negative for TSS is 0 (where 0 = below detection limit)
waterquality$TSS_mgL[waterquality$TSS_mgL< 0] <- 0  

# - SUVA - calculate SUVA
# find the specific absorbance at 254 nm from absorbance spectra
waterquality <- merge(waterquality, spec.indicies, by = "samplename") # merge in the spectral indicies with the water quality parameters
waterquality$SUVA <- waterquality$abs254.dec/waterquality$NPOC_DOC_corrected

# Normalize peaks a,c,b,t by DOC
waterquality$peakA.norm <- waterquality$peakA/waterquality$NPOC_DOC_corrected
waterquality$peakC.norm <- waterquality$peakC/waterquality$NPOC_DOC_corrected
waterquality$peakB.norm <- waterquality$peakB/waterquality$NPOC_DOC_corrected
waterquality$peakT.norm <- waterquality$peakT/waterquality$NPOC_DOC_corrected

# merge CM fits into the waterquality dataframe to get the redox index
waterquality <- merge(waterquality, CM.pre, by = "samplename")

##############
# Find stats (average, max, min, stdev) for water quality parameters - grouping by region, and then watershed

# Make a subgrouping just with the parameters that you want to simplify
quality.stats <- data.frame(waterquality$samplename, waterquality$pH, waterquality$EC_mScm,
                       waterquality$DO_mgL, waterquality$Water.Temp, waterquality$TSS_mgL,
                       waterquality$NPOC_DOC_corrected, waterquality$TN, waterquality$SUVA, waterquality$F, 
                       waterquality$Br, waterquality$NO3,
                       waterquality$abs254.dec, waterquality$SUVA, waterquality$e2e3.dec, waterquality$e4e6.dec, waterquality$SR.dec,
                       waterquality$HIX_ohno_area, waterquality$FI, waterquality$FrI, 
                       waterquality$peakA.norm, waterquality$peakC.norm, waterquality$peakB.norm,
                       waterquality$peakT.norm, waterquality$OFI,
                       waterquality$redox, waterquality$perprotein,
                       waterquality$Region, waterquality$Watershed
                       )
# write.csv(quality.stats, file= paste(save.directory, 'teststats.csv', sep = "/"))
# plot to see outliers and odd data
# plot(waterquality$samplename, waterquality$perprotein)

# group first by region
region.mean <- aggregate(quality.stats[,2:27],list(region = quality.stats$waterquality.Region),mean, na.rm = TRUE)
region.max <- aggregate(quality.stats[,2:27],list(region = quality.stats$waterquality.Region),max, na.rm = TRUE)
region.min <- aggregate(quality.stats[,2:27],list(region = quality.stats$waterquality.Region),min, na.rm = TRUE)
region.sd <- aggregate(quality.stats[,2:27],list(region = quality.stats$waterquality.Region),sd, na.rm = TRUE)
region <- rbind(region.mean, region.max, region.min, region.sd)

# group secondly by watershed
Watershed.mean <- aggregate(quality.stats[,2:27],list(Watershed = quality.stats$waterquality.Watershed),mean, na.rm = TRUE)
Watershed.max <- aggregate(quality.stats[,2:27],list(Watershed = quality.stats$waterquality.Watershed),max, na.rm = TRUE)
Watershed.min <- aggregate(quality.stats[,2:27],list(Watershed = quality.stats$waterquality.Watershed),min, na.rm = TRUE)
Watershed.sd <- aggregate(quality.stats[,2:27],list(Watershed = quality.stats$waterquality.Watershed),sd, na.rm = TRUE)
watershed <- rbind(Watershed.mean, Watershed.max, Watershed.min, Watershed.sd)

# Do stats for all of the 
wq.mean <- apply(quality.stats[,2:27], 2, mean, na.rm = TRUE)
wq.max <- apply(quality.stats[,2:27], 2, max, na.rm = TRUE)
wq.min <- apply(quality.stats[,2:27], 2, min, na.rm = TRUE)
wq.sd <- apply(quality.stats[,2:27], 2, sd, na.rm = TRUE)
#wq.length <- apply(quality.stats[,2:12], 2, length, na.rm = TRUE)
#wq.se <- wq.sd/wq.length
wq.stat <- rbind(wq.mean, wq.max, wq.min, wq.sd)  # bind stats back together
write.table(wq.stat, file = paste(save.directory, "waterqualitystats.csv", sep = ""),  sep = ",")

# Do t-tests to see if certain groups are significantly different.


######################## 
# Part 1 Figure 2 - Boxplot of Fmax values from PARAFAC fits to the 6-component model
# Partitioned according to region.
###############   Variation within Components - prechlorination boxplots
# assemble data with component in one column and FMax in another
# merge the watershed codes to the water quality tables
Pre.Fmax <- merge(Pre.Fmax, waterquality[,c(1,36:38)], by = "samplename")

# unfold the dataframe into a form that can do histograms easily
remove(PCApre6)
for (i in 1:6){
  temp.C <- data.frame(Pre.Fmax[,i+1])
  temp.component <- colnames(Pre.Fmax)[i+1]
  temp.C$component <- temp.component
  temp.C$Region <- Pre.Fmax$Region
  temp.C$samplename <- Pre.Fmax$samplename
  
  # if the merged dataset  exists, append to it by row
  if (exists("PCApre6")){
    PCApre6 <- rbind(PCApre6, temp.C)
  }
  
  # if the merged dataset doesn't exist, create it
  if (!exists("PCApre6")){
    PCApre6 <- temp.C
  }
}

colnames(PCApre6)[1] <- "values" #rename first column that contains Fmax values

# save the boxplot as a figure in file
png(paste(save.directory, "/PrePARAFACboxplot.png", sep = ""),    # create graphic for the         
    width = 5*300,        # 5 x 300 pixels
    height = 3*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

  ggplot(PCApre6, aes(x=component, y=values, fill=Region)) + 
    geom_boxplot(outlier.shape = NA)  +#remove extreme values
    labs(title="PARAFAC Fmax by Region",x="PARAFAC Components", y = "Fmax Values") 
dev.off()

# Express as the percentgae, rather than the Fmax value:
colnames(Pre.Fmax.per)[7] <- "samplename"
Pre.Fmax.per <- merge(Pre.Fmax.per, waterquality[,c(1,36:38)], by = "samplename")

# unfold the dataframe into a form that can do boxplots easily
remove(PCApre6.per)
for (i in 1:6){
  temp.C <- data.frame(Pre.Fmax.per[,i+1])
  temp.component <- colnames(Pre.Fmax.per)[i+1]
  temp.C$component <- temp.component
  temp.C$Region <- Pre.Fmax.per$Region
  temp.C$samplename <- Pre.Fmax.per$samplename
  
  # if the merged dataset  exists, append to it by row
  if (exists("PCApre6.per")){
    PCApre6.per <- rbind(PCApre6.per, temp.C)
  }
  
  # if the merged dataset doesn't exist, create it
  if (!exists("PCApre6.per")){
    PCApre6.per <- temp.C
  }
}

colnames(PCApre6.per)[1] <- "values" #rename first column that contains Fmax values

# save the boxplot as a figure in file
png(paste(save.directory, "/PrePARAFAC_per_boxplot.png", sep = ""),    # create graphic for the         
    width = 5*300,        # 5 x 300 pixels
    height = 3*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

ggplot(PCApre6.per, aes(x=component, y=values, fill=Region)) + 
  geom_boxplot(outlier.shape = NA)  +#remove extreme values
  labs(title="PARAFAC Fmax by Region",x="PARAFAC Components", y = "Fmax (% of total fluorescence)") 
dev.off()

########################
# Part 1 Figure 3 - PCA on the water quality, CM, spectral proxies, and PARAFAC Fmax values.
# Can we see groupings according to the region?
# What variables account for the most variation between sites?
# Compile: water quality proxies, spectral proxies, PARAFAC Fmax, CMFmax
colnames(quality.stats)[1] <- "samplename"
wq.all <- data.frame(Reduce(function(x,y) merge(x,y, by = "samplename", all = TRUE), 
                            list(Pre.Fmax, quality.stats)))

# First, PCA with all of the variables (including the water quality parameters)
# Take out columns of data that you donn't want to include in PCA
wq.all.select <- data.frame(wq.all$waterquality.Region, wq.all$DR_C1, wq.all$DR_C2,wq.all$DR_C3,wq.all$DR_C4,wq.all$DR_C5,wq.all$DR_C6,
                       wq.all$waterquality.NPOC_DOC_corrected,
                       wq.all$waterquality.SUVA, wq.all$waterquality.NO3,wq.all$waterquality.abs254.dec,
                       wq.all$waterquality.e2e3.dec, wq.all$waterquality.e4e6.dec, wq.all$waterquality.SR.dec,
                       wq.all$waterquality.HIX_ohno_area, wq.all$waterquality.FI, wq.all$waterquality.FrI,
                       wq.all$waterquality.peakA.norm,wq.all$waterquality.peakC.norm,wq.all$waterquality.peakB.norm,wq.all$waterquality.peakT.norm,
                       wq.all$waterquality.OFI,wq.all$waterquality.redox
                       )
colnames(wq.all.select) <- gsub("\\wq.all.", "", colnames(wq.all.select))
colnames(wq.all.select) <- gsub("\\waterquality.", "", colnames(wq.all.select))
colnames(wq.all.select) <- gsub("\\DR_", "", colnames(wq.all.select))
colnames(wq.all.select) <- gsub("\\.norm", "", colnames(wq.all.select))
colnames(wq.all.select) <- gsub("\\_ohno_area", "", colnames(wq.all.select))
colnames(wq.all.select) <- gsub("\\.dec", "", colnames(wq.all.select))
colnames(wq.all.select)[8] <- "DOC"

# replace infinities in data with Nans
is.na(wq.all.select) <- sapply(wq.all.select, is.infinite)
wq.all.select <- na.omit(wq.all.select)

# do PCA on all variables - which variables explain the greatest degree of variation?
wq.all.pca <- prcomp(wq.all.select[,2:23], center = TRUE, scale. = TRUE, na.action=na.omit)
summary(wq.all.pca)

# plot PCA results
png(paste(save.directory, "/DBP_PCAWQdata.png", sep = ""),    # create graphic for the         
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

ggbiplot(wq.all.pca, obs.scale = 1, var.scale = 1, groups = na.omit(wq.all.select[,1]),
         ellipse = TRUE, circle = FALSE) +ggtitle("PCA Results- WQ Parameters") #+geom_point(colour = cbPalette[2:7]) 
dev.off()

# show relative contribution of each variable to the first 5 components of PCA
wq.pca <- PCA(na.omit(wq.all.select[,2:23]), graph = TRUE)
head(wq.pca$var$contrib)
PCA.contrib <- data.frame(wq.pca$var$contrib)
# sort variables by incresing contribution across the first 5 component
wq.sortvar <- PCA.contrib[order(-PCA.contrib$Dim.1,-PCA.contrib$Dim.2,-PCA.contrib$Dim.3,-PCA.contrib$Dim.4,-PCA.contrib$Dim.5), ]
write.csv(wq.sortvar, file = paste(save.directory, "PCAcontributions.csv", sep ="/")) #write contributions to a csv file that are sorted
# scree plot
fviz_screeplot(wq.pca, ncp=6) # first 6 components
# plot contribution to first 2 PCA components
fviz_contrib(wq.pca, choice = "var", axes = 1)
fviz_contrib(wq.pca, choice = "var", axes = 2)

#################### 
# Figure ** Correlation between DOC concentration and spectral parameters
# Question - which spectral parameters correlate best to DOC concentration?

pairs(wq.all.select[,2:23], panel = panel.smooth) #figure marigins too large

model1 <- lm(wq.all.select$DOC ~ wq.all.select$C1+wq.all.select$C2+wq.all.select$C3 + 
    wq.all.select$C4 + wq.all.select$C5 +wq.all.select$C6 +
    wq.all.select$SUVA +
    wq.all.select$abs254+
    wq.all.select$e2e3+ wq.all.select$e4e6 + wq.all.select$SR +
    wq.all.select$HIX + wq.all.select$FI+ wq.all.select$FrI+
    wq.all.select$peakA +wq.all.select$peakC+wq.all.select$peakB+wq.all.select$peakT +
    wq.all.select$OFI+wq.all.select$redox)
# Investigate linear model fits
summary(model1)
anova(model1)
coefficients(model1)
fitted(model1)
residuals(model1)
vcov(model1)

# Use step function to remove variables with less correlation to DOC concentration
step = stepAIC(model1, direction = 'both')

# logistic regression function to look at DOC - do this by region?
# http://www.stat.columbia.edu/~martin/W2024/R11.pdf
DOC.logr <- glm(wq.all.select$DOC ~ wq.all.select$C1+wq.all.select$C2+wq.all.select$C3 + 
                wq.all.select$C4 + wq.all.select$C5 +wq.all.select$C6 +
                 wq.all.select$SUVA +
                 wq.all.select$abs254+
                 wq.all.select$e2e3+ wq.all.select$e4e6 + wq.all.select$SR +
                 wq.all.select$HIX + wq.all.select$FI+ wq.all.select$FrI+
                 wq.all.select$peakA +wq.all.select$peakC+wq.all.select$peakB+wq.all.select$peakT +
                 wq.all.select$OFI+wq.all.select$redox,
                data = wq.all.select, family = "gaussian")
summary(DOC.logr)
beta =coef(DOC.logr)
stepAIC(DOC.logr, direction = 'both')

DOC.logr.select <- glm(DOC ~ C1 + C3 
                      + abs254 + peakC + peakT +OFI,
                      data = wq.all.select, family = "gaussian")
summary(DOC.logr.select)
#get error terms for model
results.reduced =glm(DOC ~ 1, data = wq.all.select, family = "gaussian")
anova(results.reduced,DOC.logr.select , test="Chisq")

################################################################################################################################################################
# Pt 2 - how does chlorination change the spectral composition of EEMS?
# princip-al component analysis to look at how chlorination altered spectral characteristics (EEMS)
# References:
# http://www.r-bloggers.com/computing-and-visualizing-pca-in-r/
# http://planspace.org/2013/02/03/pca-3d-visualization-and-clustering-in-r/
###########################################################
# calculate the percent difference in water quality and spectral parameters
# Percent difference = (pre-post)/pre*100

# read in the csv file with the CM fits for the pre and post chlorinated EEMs
CM.prepost.fmax <- read.csv("/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_preandpost/DBP_prepost_CMresults/DBPprepost_componentsandloadings_CM_Fmax.csv", header = TRUE, sep = ",")
CM.prepost.percent <- read.csv("/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_preandpost/DBP_prepost_CMresults/DBPprepost_componentsandloadings_CM.csv", header = TRUE, sep = ",")

# split the CM file into the pre and post file
CM.pre.percent <- CM.prepost.percent[grep("DBP0", CM.prepost.percent$sample.ID), ] # prechlorinated EEMS
CM.post.percent <- CM.prepost.percent[grep("DBPChlor", CM.prepost.percent$sample.ID), ] # postchlorinated EEMS

# merge with the spectral indicies
CM.pre.percent$samplename <- sub(' *\\DBPprepost.*$', '', CM.pre.percent$sample.ID) # creTe sample name with correct format
spec.indicies <- merge(spec.indicies, CM.pre.percent, by = "samplename", all = FALSE)

CM.post.percent$samplename <- sub(' *\\DBPprepost.*$', '', CM.post.percent$sample.ID) # create sample name with correct format
spec.indicies.post <- merge(spec.indicies.post, CM.post.percent, by = "samplename", all = FALSE)

## merge the 6 component PARAFAC model into both the pre and post chlorinated spec indicies
spec.indicies.PARAFAC <- merge(spec.indicies, Pre.Fmax[,1:7], by = "samplename", all = TRUE)
# post chlorination
PARAFAC.post.Fmax <- prepost.6comp.Fmax[grep("DBPChlor", prepost.6comp.Fmax$samplename), ] # postchlorinated EEMS PARAFAC Fmax
spec.indicies.post.PARAFAC <- merge(spec.indicies.post, PARAFAC.post.Fmax, by = "samplename")

## Merge the pre and post chlorinated EEMS and calculate the difference between them
# Get the sample ID on the post chlorinated EEM
spec.indicies.post.PARAFAC$samplename.post <- spec.indicies.post.PARAFAC$samplename
sample <- sapply(strsplit(as.character(spec.indicies.post.PARAFAC$samplename.post), split='Chlor', fixed=TRUE), function(x) (paste(x[1],x[2], sep = "")))
samplename <- str_pad(sample, 4, pad = "0")
spec.indicies.post.PARAFAC$samplename <- samplename
remove(sample, samplename)
spec.indicies.post.PARAFAC$samplename.post <- NULL
colnames(spec.indicies.post.PARAFAC) <- colnames(spec.indicies.PARAFAC)

# merge the two together
spec.all <- merge(spec.indicies.PARAFAC,spec.indicies.post.PARAFAC,by="samplename", all = FALSE)

#calculate the percent change (pre-post)/pre*100
delta.spec <- (spec.all[,grepl("*\\.x$",names(spec.all))] - spec.all[,grepl("*\\.y$",names(spec.all))])/spec.all[,grepl("*\\.x$",names(spec.all))]*100
# calculate delat abs 272
delta.abs272 <- spec.all$abs272.dec.x - spec.all$abs272.dec.y

delta.spec <- cbind(spec.all[,1,drop=FALSE],delta.abs272,delta.spec) # add in sample names

## boxplots to show percent changes
# choose the proxies that are important
delta.spec.select <- data.frame(delta.spec$samplename, 
                                delta.spec$abs254.dec.x, delta.spec$abs272.dec.x,
                                delta.spec$e2e3.dec.x, delta.spec$e4e6.dec.x,
                                #delta.spec$CDOM.total.int.dec.x, 
                                delta.spec$SR.dec.x,
                                delta.spec$FI.x,delta.spec$HIX_ohno_area.x,
                                delta.spec$FrI.x,
                                delta.spec$peakA.x,delta.spec$peakC.x, delta.spec$peakB.x,delta.spec$peakT.x,
                                delta.spec$OFI.x,
                                delta.spec$perprotein.x, delta.spec$redox.x,
                                delta.spec$DR_C1.x, delta.spec$DR_C2.x,delta.spec$DR_C3.x,delta.spec$DR_C4.x,delta.spec$DR_C5.x,delta.spec$DR_C6.x
                                )
colnames(delta.spec.select) <- gsub("\\.x", "", colnames(delta.spec.select))
colnames(delta.spec.select) <- gsub("\\.dec", "", colnames(delta.spec.select))
colnames(delta.spec.select) <- gsub("\\DR_", "", colnames(delta.spec.select))
colnames(delta.spec.select) <- gsub("delta.spec.", "", colnames(delta.spec.select))
colnames(delta.spec.select) <- gsub("_ohno_area", "", colnames(delta.spec.select))

# unfold the dataframe into a form that can do histograms easily
remove(perc.spec)
for (i in 1:(dim(delta.spec.select)[2]-1)){
  temp.C <- data.frame(delta.spec.select[,i+1])
  temp.component <- colnames(delta.spec.select)[i+1]
  temp.C$proxy <- temp.component
  # if the merged dataset  exists, append to it by row
  if (exists("perc.spec")){
    perc.spec <- rbind(perc.spec, temp.C)
  }
  # if the merged dataset doesn't exist, create it
  if (!exists("perc.spec")){
    perc.spec <- temp.C
  }
  remove(temp.C)
}

# do box plot
colnames(perc.spec)[1] <- "Percent" #rename first column that contains Fmax values
p <- ggplot(perc.spec, aes(x=proxy, y=Percent)) +
  geom_boxplot(outlier.shape = NA) + #remove extreme values
  coord_cartesian(ylim = c(-500, 500)) + # change y limits on boxplot
  labs(title="Percent Change Upon Chlorination",x="Spectral Proxy", y = "Percent Change (Upon Chlorination)")
p + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# do t tests to look at which variables are significantly changed by chlorination (pre versus post)
t.abs254 <- t.test(spec.indicies.PARAFAC$abs254.dec, spec.indicies.post.PARAFAC$abs254.dec)
t.abs272 <- t.test(spec.indicies.PARAFAC$abs272.dec, spec.indicies.post.PARAFAC$abs272.dec)
t.e2e3 <- t.test(spec.indicies.PARAFAC$e2e3.dec, spec.indicies.post.PARAFAC$e2e3.dec)
t.e4e6 <- t.test(spec.indicies.PARAFAC$e4e6.dec, spec.indicies.post.PARAFAC$e4e6.dec)
t.CDOM <- t.test(spec.indicies.PARAFAC$CDOM.total.int.dec, spec.indicies.post.PARAFAC$CDOM.total.int.dec)
t.SR <- t.test(spec.indicies.PARAFAC$SR.dec,spec.indicies.post.PARAFAC$SR.dec)
# For FI, get rid of Inf values
spec.indicies.post$FI[!is.finite(spec.indicies.post.PARAFAC$FI)] <- NaN
FI.post <- na.omit(spec.indicies.post$FI)
t.FI <- t.test(spec.indicies.PARAFAC$FI,FI.post, na.rm = TRUE)

t.HIX_ohno_area <- t.test(spec.indicies.PARAFAC$HIX_ohno_area, spec.indicies.post.PARAFAC$HIX_ohno_area)
t.FrI <- t.test(spec.indicies.PARAFAC$FrI, spec.indicies.post.PARAFAC$FrI)
t.peakA <- t.test(spec.indicies.PARAFAC$peakA, spec.indicies.post.PARAFAC$peakA)
t.peakC <- t.test(spec.indicies.PARAFAC$peakC, spec.indicies.post.PARAFAC$peakC)
t.peakB <- t.test(spec.indicies.PARAFAC$peakB, spec.indicies.post.PARAFAC$peakB)
t.peakT <- t.test(spec.indicies.PARAFAC$peakT, spec.indicies.post.PARAFAC$peakT)
t.OFI <- t.test(spec.indicies.PARAFAC$OFI, spec.indicies.post.PARAFAC$OFI)
t.perprotein <- t.test(spec.indicies.PARAFAC$perprotein, spec.indicies.post.PARAFAC$perprotein)
t.redox <- t.test(spec.indicies.PARAFAC$redox, spec.indicies.post.PARAFAC$redox)
t.DR_C1 <- t.test(spec.indicies.PARAFAC$DR_C1, spec.indicies.post.PARAFAC$DR_C1)
t.DR_C2 <- t.test(spec.indicies.PARAFAC$DR_C2, spec.indicies.post.PARAFAC$DR_C2)
t.DR_C3 <- t.test(spec.indicies.PARAFAC$DR_C3, spec.indicies.post.PARAFAC$DR_C3)
t.DR_C4 <- t.test(spec.indicies.PARAFAC$DR_C4, spec.indicies.post.PARAFAC$DR_C4)
t.DR_C5 <- t.test(spec.indicies.PARAFAC$DR_C5, spec.indicies.post.PARAFAC$DR_C5)
t.DR_C6 <- t.test(spec.indicies.PARAFAC$DR_C6, spec.indicies.post.PARAFAC$DR_C6)

# do correlation matrix between the percent change upon chlorination
delta.spec.select <- do.call(data.frame,lapply(delta.spec.select[,2:22], function(x) replace(x, is.infinite(x),NA)))
corr.matrix.delta <- cor(na.omit(delta.spec.select)) # correlation matric between variables

png(paste(save.directory, "DBPdelta_correlations.png", sep = "/"),    # create graphic for the correlations plot       
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

  corrplot(corr.matrix.delta, method = "circle") #plot matrix
dev.off()

# write correlation matrix to csv file
write.csv(corr.matrix.delta, file = paste(save.directory, "DBPDelta_Corrmatrix.csv", sep ="/")) #write correlation matrix to a csv file

#########################################################################################################################
# Pt 3 - How does EEMS and Water Quality parameters predict DBP formation?
# Thoughts:
# use machine learning algorythm to predict DBP parameters
# try multivraite linear model where
#[DBP] = theta0 + theta1x1+theta2x2
# where x = spectral and EEM characteristics, as well as water quality parameters (pH, anion concentrations, DOC)

# First, mutiple linear regression to rank the contribution of water quality variables to 

########### compile the DBP concentrations + clean the data

# calculate THM and HAA total (sum of all of the different species for that sample)
HAA$total <- apply(HAA[,2:10], 1, sum)
THM$total <- apply(THM[,2:5], 1, sum)

######## apply linear model, where y = HAA and THMs, and x = spectral variables...
# try #1: 
# y = total HAA and THMs (two models)
# x = PARAFAC components from the prechlorinated EEMS, water quality parameters, spectral parameters

# TO DO - add PCA components for pre chlorinated EEMs; add in the delta eems for the 6 component fit?

####### Bind spectral and water quality parameters with DBP data
waterquality.mod <- data.frame(waterquality$samplename, waterquality$NPOC_DOC_corrected, waterquality$SUVA, 
                               waterquality$abs254.dec, waterquality$abs272.dec, delta.spec$delta.abs272 , 
                               waterquality$e2e3.dec, 
                               waterquality$e4e6.dec, waterquality$SR.dec, 
                               waterquality$FI, waterquality$HIX_ohno_area, waterquality$FrI, 
                               waterquality$peakA.norm, waterquality$peakC.norm, waterquality$peakB.norm, waterquality$peakT.norm, 
                               waterquality$OFI, waterquality$perprotein, waterquality$redox)

colnames(waterquality.mod) <- gsub("\\waterquality.", "", colnames(waterquality.mod))
colnames(waterquality.mod) <- gsub("\\.norm", "", colnames(waterquality.mod))
colnames(waterquality.mod)[6] <- "delta.abs272"

# add in the PARAFAC fits
waterquality.mod.1 <- merge(waterquality.mod, Pre.Fmax, by = 'samplename', all = TRUE)

HAA.waterq <- Reduce(function(x, y) merge(x, y, by = 'samplename', all=FALSE), list(HAA, waterquality.mod, Pre.Fmax))
THM.waterq <- Reduce(function(x, y) merge(x, y, by = 'samplename', all=FALSE), list(THM, waterquality.mod, Pre.Fmax))

# get rid of columns you don't need
HAA.waterq$sample.ID <- NULL
THM.waterq$sample.ID <- NULL

# calculate THM and HAA yield
HAA.waterq$total.yield <- HAA.waterq$total/HAA.waterq$NPOC_DOC_corrected
THM.waterq$total.yield <- THM.waterq$total/THM.waterq$NPOC_DOC_corrected

# Do linear models for total HAA and total THMs
# Note that expect some of the variables to co-relate, thus use gls linear fit model
HAA.waterq <- lapply(HAA.waterq, as.numeric)
THM.waterq <- lapply(THM.waterq, as.numeric)

############ total HAAs
# use tree to look at the correlation between variables - first, CM model
#library("tree")
#HAA.total <- CM.model.HAA$total
#HAA.CM <-lapply(CM.model.HAA[,2:17], as.numeric) # convert to numeric prior to running model
# Note this reference about trying to fit linear models to highly correalated data
# https://www.researchgate.net/post/R_Squared_Value_is_high_about_070_however_the_p_value_for_all_my_independent_value_is_over_005_How_could_this_happen
###############

model.HAAtotal <- tree(HAA.total ~ ., data = HAA.CM) # run tree model on CM data
plot(model.HAAtotal) # plot tree model
text(model.HAAtotal) # put in text labels into the tree label

# initial model - linear model
lmmodel.HAAtotal <- lm(HAA.waterq$total ~ HAA.waterq$NPOC_DOC_corrected*HAA.waterq$SUVA 
                         + HAA.waterq$DR_C1 + HAA.waterq$DR_C2 +
                           HAA.waterq$DR_C3 + HAA.waterq$DR_C4 + HAA.waterq$DR_C5 +
                           HAA.waterq$DR_C6 + HAA.waterq$perprotein + HAA.waterq$redox 
                         + HAA.waterq$abs254.dec + HAA.waterq$abs272.dec + HAA.waterq$delta.abs272
                         + HAA.waterq$e2e3.dec + HAA.waterq$e4e6.dec + HAA.waterq$SR.dec)
summary(lmmodel.HAAtotal)

# Use step function to remove variables with less correlation to DOC concentration
step.lm.haa = summary(step(lmmodel.HAAtotal, direction = 'both'))

# generalized least squres model - specific variables
glmmodel.HAAtotal.select <- glm(HAA.waterq$total ~ HAA.waterq$NPOC_DOC_corrected 
                      #* HAA.waterq$SUVA #+ HAA.waterq$abs254.dec 
                      + HAA.waterq$abs272.dec 
                      + HAA.waterq$DR_C1 + HAA.waterq$DR_C2 + HAA.waterq$DR_C3 + HAA.waterq$DR_C4 + HAA.waterq$DR_C5 + HAA.waterq$DR_C6,
                      family = "gaussian")

summary(glmmodel.HAAtotal.select)

# Use step function to remove variables with less correlation to DOC concentration
step.glm.haa.select = summary(step(glmmodel.HAAtotal.select, direction = 'both'))
test <- anova(glmmodel.HAAtotal.select)

# generalized least squres model
glmmodel.HAAtotal <- glm(HAA.waterq$total ~ HAA.waterq$NPOC_DOC_corrected 
                         + HAA.waterq$SUVA + HAA.waterq$abs254.dec + HAA.waterq$abs272.dec + HAA.waterq$delta.abs272
                         + HAA.waterq$DR_C1 + HAA.waterq$DR_C2 + HAA.waterq$DR_C3 + HAA.waterq$DR_C4 + HAA.waterq$DR_C5 + HAA.waterq$DR_C6  
                         + HAA.waterq$perprotein + HAA.waterq$redox 
                         + HAA.waterq$e2e3.dec + HAA.waterq$e4e6.dec + HAA.waterq$SR.dec,
                         family = "gaussian")

summary(glmmodel.HAAtotal)

# Use step function to remove variables with less correlation to DOC concentration
step.glm.haa = summary(step(glmmodel.HAAtotal, direction = 'both'))

# Do individual lm fits
summary(lm(HAA.waterq$total~ HAA.waterq$NPOC_DOC_corrected))$r.squared
summary(lm(HAA.waterq$total ~ HAA.waterq$SUVA))$r.squared
summary(lm(HAA.waterq$total ~ HAA.waterq$abs254.dec))$r.squared
summary(lm(HAA.waterq$total ~ HAA.waterq$abs272.dec))$r.squared
summary(lm(HAA.waterq$total ~ HAA.waterq$delta.abs272))$r.squared
summary(lm(HAA.waterq$total ~ HAA.waterq$DR_C1))$r.squared
summary(lm(HAA.waterq$total ~ HAA.waterq$DR_C2))$r.squared
summary(lm(HAA.waterq$total ~ HAA.waterq$DR_C3))$r.squared
summary(lm(HAA.waterq$total ~ HAA.waterq$DR_C4))$r.squared
summary(lm(HAA.waterq$total~ HAA.waterq$DR_C5))$r.squared
summary(lm(HAA.waterq$total~ HAA.waterq$DR_C6))$r.squared
summary(lm(HAA.waterq$total~ HAA.waterq$perprotein))$r.squared
summary(lm(HAA.waterq$total ~ HAA.waterq$redox))$r.squared
summary(lm(HAA.waterq$total~ HAA.waterq$e2e3))$r.squared
summary(lm(HAA.waterq$total ~ HAA.waterq$e4e6))$r.squared
summary(lm(HAA.waterq$total ~ HAA.waterq$CDOM.total))$r.squared
summary(lm(HAA.waterq$total ~ HAA.waterq$slope_ratio))$r.squared
summary(lm(HAA.waterq$total ~ HAA.waterq$SR))$r.squared
summary(lm(HAA.waterq$total ~ HAA.waterq$FI))$r.squared
summary(lm(HAA.waterq$total ~ HAA.waterq$FrI))$r.squared
summary(lm(HAA.waterq$total ~ HAA.waterq$peakt.peakC))$r.squared
summary(lm(HAA.waterq$total ~ HAA.waterq$HIX_ohno_area))$r.squared

# export the linear model results
#write.table(HAA.total.lm, file = paste(save.directory, "HAATotal_lmresults.csv", sep = ","))
# plot the total HAA versus all variables above
plot(HAA.waterq$total ~ HAA.waterq$NPOC_DOC_corrected)
plot(HAA.waterq$total ~ HAA.waterq$SUVA)
plot(HAA.waterq$total ~ HAA.waterq$abs254.dec)
plot(HAA.waterq$total ~ HAA.waterq$abs272.dec)
plot(HAA.waterq$total ~ HAA.waterq$delta.abs272)
plot(HAA.waterq$total ~ HAA.waterq$DR_C1)
plot(HAA.waterq$total ~ HAA.waterq$DR_C2)
plot(HAA.waterq$total ~ HAA.waterq$DR_C3)
plot(HAA.waterq$total ~ HAA.waterq$DR_C4)
plot(HAA.waterq$total ~ HAA.waterq$DR_C5)
plot(HAA.waterq$total ~ HAA.waterq$DR_C6)
plot(HAA.waterq$total ~ HAA.waterq$perprotein)
plot(HAA.waterq$total ~ HAA.waterq$redox)
plot(HAA.waterq$total ~ HAA.waterq$e2e3)
plot(HAA.waterq$total ~ HAA.waterq$e4e6)
plot(HAA.waterq$total ~ HAA.waterq$CDOM.total)
plot(THM.waterq$total ~ HAA.waterq$slope_ratio)
plot(HAA.waterq$total ~ HAA.waterq$SR)
plot(HAA.waterq$total ~ HAA.waterq$FI)
plot(HAA.waterq$total ~ HAA.waterq$FrI)
plot(HAA.waterq$total ~ HAA.waterq$peakt.peakC)
plot(HAA.waterq$total ~ HAA.waterq$HIX_ohno_area)

####################### plot the significant variables nicely using ggplot
library(ggpmisc)

HAA.waterq.data <- as.data.frame(HAA.waterq) # convert to data frame
HAA.waterq.data<- HAA.waterq.data[-22,]
my.formula = y~x
# HAA versus DOC
ggplot(HAA.waterq.data, aes(x=NPOC_DOC_corrected, y=total)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  ggtitle("[THAAs] versus [DOC]") +
  labs(x="DOC Concentration (mg/L)",y="THAAs Concentration (ug/L)")

# HAA versus abs254
ggplot(HAA.waterq.data, aes(x=abs254.dec, y=total)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  ggtitle("[THAAs] versus a254") +
  labs(x="a254 (1/m)",y="THAAs Concentration (ug/L)")

# HAA versus abs272
ggplot(HAA.waterq.data, aes(x=abs272.dec, y=total)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  ggtitle("[THAAs] versus a272") +
  labs(x="a272 (1/m)",y="THAAs Concentration (ug/L)")

# HAA versus C1
ggplot(HAA.waterq.data, aes(x=DR_C1, y=total)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  ggtitle("[THAAs] versus C1") +
  labs(x="C1 Fmax (R.U)",y="THAAs Concentration (ug/L)")

# HAA versus C2
ggplot(HAA.waterq.data, aes(x=DR_C2, y=total)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  ggtitle("[THAAs] versus C2") +
  labs(x="C2 Fmax (R.U)",y="THAAs Concentration (ug/L)")

# HAA versus C3
ggplot(HAA.waterq.data, aes(x=DR_C3, y=total)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  ggtitle("[THAAs] versus C3") +
  labs(x="C3 Fmax (R.U)",y="THAAs Concentration (ug/L)")

####################### Total THMs
#NPOC and THMs
NPOC.model.THM <- lm(THM.waterq$total ~ THM.waterq$NPOC)
summary(NPOC.model.THM)
plot(THM.waterq$total ~ THM.waterq$NPOC) #note high outliers! are these green roofs? may have to remove
test <- lm(THM.waterq$total ~ THM.waterq$NPOC)

# initial linear model - lm model
model.THMtotal.lm <- lm(THM.waterq$total ~ THM.waterq$NPOC_DOC_corrected 
                      + THM.waterq$DR_C1 + THM.waterq$DR_C2 +
                        THM.waterq$DR_C3 + THM.waterq$DR_C4 + THM.waterq$DR_C5 +
                        THM.waterq$DR_C6 + THM.waterq$perprotein + THM.waterq$redox 
                      + THM.waterq$SUVA + THM.waterq$abs254.dec + THM.waterq$abs272.dec + THMwaterq$delta
                      + THM.waterq$e2e3.dec + THM.waterq$e4e6.dec + THM.waterq$SR.dec
                    )
summary(model.THMtotal.lm)

# try step function
lm.step.THM <- summary(step(model.THMtotal.lm, direction="both")) # summary of linear model

# initial linear model - glm model
glmmodel.THMtotal <- glm(THM.waterq$total ~ THM.waterq$NPOC_DOC_corrected 
                         + THM.waterq$SUVA + THM.waterq$abs254.dec + THM.waterq$abs272.dec + THM.waterq$delta.abs272
                         + THM.waterq$DR_C1 + THM.waterq$DR_C2 + THM.waterq$DR_C3 + THM.waterq$DR_C4 + THM.waterq$DR_C5 + THM.waterq$DR_C6  
                         + THM.waterq$perprotein + THM.waterq$redox 
                         + THM.waterq$e2e3.dec + THM.waterq$e4e6.dec + THM.waterq$SR.dec,
                         family = "gaussian")
summary(glmmodel.THMtotal)

# try step function
glm.step.THM <- summary(step(glmmodel.THMtotal, direction="both")) # summary of linear model

# Do individual lm fits
summary(lm(THM.waterq$total~ THM.waterq$NPOC_DOC_corrected))$r.squared
summary(lm(THM.waterq$total ~ THM.waterq$SUVA))$r.squared
summary(lm(THM.waterq$total ~ THM.waterq$abs254.dec))$r.squared
summary(lm(THM.waterq$total ~ THM.waterq$abs272.dec))$r.squared
summary(lm(THM.waterq$total ~ THM.waterq$delta.abs272))$r.squared
summary(lm(THM.waterq$total ~ THM.waterq$DR_C1))$r.squared
summary(lm(THM.waterq$total ~ THM.waterq$DR_C2))$r.squared
summary(lm(THM.waterq$total ~ THM.waterq$DR_C3))$r.squared
summary(lm(THM.waterq$total ~ THM.waterq$DR_C4))$r.squared
summary(lm(THM.waterq$total~ THM.waterq$DR_C5))$r.squared
summary(lm(THM.waterq$total~ THM.waterq$DR_C6))$r.squared
summary(lm(THM.waterq$total~ THM.waterq$perprotein))$r.squared
summary(lm(THM.waterq$total ~ THM.waterq$redox))$r.squared
summary(lm(THM.waterq$total~ THM.waterq$e2e3))$r.squared
summary(lm(THM.waterq$total ~ THM.waterq$e4e6))$r.squared
summary(lm(THM.waterq$total ~ THM.waterq$CDOM.total))$r.squared
summary(lm(THM.waterq$total ~ THM.waterq$slope_ratio))$r.squared
summary(lm(THM.waterq$total ~ THM.waterq$SR))$r.squared
summary(lm(THM.waterq$total ~ THM.waterq$FI))$r.squared
summary(lm(THM.waterq$total ~ THM.waterq$FrI))$r.squared
summary(lm(THM.waterq$total ~ THM.waterq$peakt.peakC))$r.squared
summary(lm(THM.waterq$total ~ THM.waterq$HIX_ohno_area))$r.squared


#write.table(THM.total.lm, file = paste(save.directory, "THMTotal_lmresults.csv", sep = ","))
# plot the total THM versus all variables above
plot(THM.waterq$total ~ THM.waterq$NPOC_DOC_corrected)
plot(THM.waterq$total ~ THM.waterq$SUVA)
plot(THM.waterq$total ~ THM.waterq$abs254.dec)
plot(THM.waterq$total ~ THM.waterq$abs272.dec)
plot(THM.waterq$total ~ THM.waterq$delta.abs272)
plot(THM.waterq$total ~ THM.waterq$DR_C1)
plot(THM.waterq$total ~ THM.waterq$DR_C2)
plot(THM.waterq$total ~ THM.waterq$DR_C3)
plot(THM.waterq$total ~ THM.waterq$DR_C4)
plot(THM.waterq$total ~ THM.waterq$DR_C5)
plot(THM.waterq$total ~ THM.waterq$DR_C6)
plot(THM.waterq$total ~ THM.waterq$perprotein)
plot(THM.waterq$total ~ THM.waterq$redox)
plot(THM.waterq$total ~ THM.waterq$e2e3)
plot(THM.waterq$total ~ THM.waterq$e4e6)
plot(THM.waterq$total ~ THM.waterq$CDOM.total)
plot(THM.waterq$total ~ THM.waterq$slope_ratio)
plot(THM.waterq$total ~ THM.waterq$SR)
plot(THM.waterq$total ~ THM.waterq$FI)
plot(THM.waterq$total ~ THM.waterq$FrI)
plot(THM.waterq$total ~ THM.waterq$peakt.peakC)
plot(THM.waterq$total ~ THM.waterq$HIX_ohno_area)

####### plot the siginifcant variables 
THM.waterq.data <- as.data.frame(THM.waterq) # convert to data frame
my.formula = y~x

# HAA versus C3
ggplot(THM.waterq.data, aes(x=DR_C4, y=total)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  ggtitle("[TTHMs] versus C4") +
  labs(x="Fmax (R.U)",y="TTHMs Concentration (ug/L)")
#C5
ggplot(THM.waterq.data, aes(x=DR_C5, y=total)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  ggtitle("[TTHMs] versus C5") +
  labs(x="Fmax (R.U)",y="TTHMs Concentration (ug/L)")
#C6
ggplot(THM.waterq.data, aes(x=DR_C6, y=total)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  ggtitle("[TTHMs] versus C6") +
  labs(x="Fmax (R.U)",y="TTHMs Concentration (ug/L)")



#################################################################################
######### correlation matrix - spectral parameters and total THM/HAA
#Figure - look at the correlation between different variables. Eventually add in THM and HAA formation potential
# merge THM/HAA with wq parameters
wq.all.select.cor <- data.frame(wq.all$samplename, wq.all$waterquality.Region, wq.all$DR_C1, wq.all$DR_C2,wq.all$DR_C3,wq.all$DR_C4,wq.all$DR_C5,wq.all$DR_C6,
                            wq.all$waterquality.NPOC_DOC_corrected,
                            wq.all$waterquality.SUVA, wq.all$waterquality.abs254.dec, waterquality$abs272.dec, delta.spec$delta.abs272,
                            wq.all$waterquality.e2e3.dec, wq.all$waterquality.e4e6.dec, wq.all$waterquality.SR.dec,
                            wq.all$waterquality.HIX_ohno_area, wq.all$waterquality.FI, wq.all$waterquality.FrI,
                            wq.all$waterquality.peakA.norm,wq.all$waterquality.peakC.norm,wq.all$waterquality.peakB.norm,wq.all$waterquality.peakT.norm,
                            wq.all$waterquality.OFI,wq.all$waterquality.redox
)

colnames(wq.all.select.cor) <- gsub("\\wq.all.", "", colnames(wq.all.select.cor))
corr.data1 <- merge(wq.all.select.cor, HAA, by = "samplename", all = TRUE)
corr.data <- merge(corr.data1, THM, by = "samplename", all = TRUE)

corr.data2 <- cbind(corr.data, corr.data$total.x, corr.data$total.y)
corr.data2 <- na.omit(corr.data2)

corr.matrix <- cor(corr.data[,3:25]) # correlation matric between variables
corr.matrix.dbps <- cor(corr.data2[,c(3:25, 41:42)])

# bind two rows of DBPS correlations to corr matrix of other variables
test <- corr.matrix.dbps[24:25,-24:-25]
test2 <- corr.matrix.dbps[,24:25]

corr.all.2 <- rbind(corr.matrix, test)
corr.all <- as.matrix(cbind(corr.all.2, test2))

colnames(corr.all) <- c('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'DOC', 'SUVA', 'abs254', 'abs272', 'Delta abs272', 'e2e3', 'e4e6', 'SR', 'HIX', 'FI', 'FrI', 'Peak A', 'Peak C', 'Peak B', 'Peak T', 'OFI', 'Redox Index', "Total HAAs", "Total THMs")
row.names(corr.all) <- c('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'DOC', 'SUVA', 'abs254', 'abs272', 'Delta abs272', 'e2e3', 'e4e6', 'SR', 'HIX', 'FI', 'FrI', 'Peak A', 'Peak C', 'Peak B', 'Peak T', 'OFI', 'Redox Index',"Total HAAs", "Total THMs")

write.csv(corr.all, file = paste(save.directory, "WQpre_Corrmatrix.csv", sep ="/")) #write correlation matrix to a csv file
#corr.all <- na.omit(corr.all)

# plot the correlation matrix of the spectral parameters
png(paste(save.directory, "DBPtotalDBPs_correlations.png", sep = "/"),    # create graphic for the         
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

corrplot(na.omit(corr.all), method = "circle") #plot matrix
dev.off()

#########################################
# try #2: 
# y = specific HAA and THMs (multiple models)
# x = PARAFAC components from the prechlorinated EEMS, water quality parameters, spectral parameters
# TO DO - add PCA components for pre chlorinated EEMs; add in the delta eems for the 6 component fit?

# Try cca from vegan package
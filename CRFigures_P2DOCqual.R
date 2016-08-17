# Script for figures associated with CR Paper 2 - Quality dynamics
# Ashlee J's PhD
#
# 27 June2016
############################################################################################

# Clear environment
rm(list = ls())
ls()

###### Necessary toolboxes
library(pacman)
p_load(lubridate,plyr,ggplot2,reshape2,car,grid,gridBase,gridExtra, taRifx,magrittr, ggpmisc,
       stringi,stringr,plyr,dplyr,tidyr,EcoHydRology,openair,reshape,lmomco, zoo, hydroTSM, rowr, 
       reshape2, pgirmess, nlme, chron, repmis, ggbiplot, FactoMineR, factoextra)

library('devtools')
install_github("ggbiplot", "vqv")
library(ggbiplot)
library('corrplot') #package corrplot
library(reshape2)
library(gridExtra)
library(grid)
library(gtable)
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
library(reshape2)
library(gridExtra)
library(grid)
library(gtable)

############ functions
mean.function <- function(filename, interval){
  mean <- ddply(filename,.(interval),summarise,mean.turb=mean(Turb.FTU,na.rm = TRUE),
                mean.DOC.=mean(DOCcorr,na.rm = TRUE),mean.SUVA.=mean(SUVA,na.rm = TRUE),
                mean.e2e3.=mean(e2e3,na.rm = TRUE),
                mean.e4e6.=mean(e4e6,na.rm = TRUE),mean.slope295275.=mean(slope295_275,na.rm = TRUE),
                mean.slope400350.=mean(slope400_350,na.rm = TRUE),
                mean.sloperatio.=mean(slope_ratio,na.rm = TRUE),
                mean.CDOMtot1.=mean(CDOM.total1,na.rm = TRUE),
                mean.CDOMtot2.=mean(CDOM.total2,na.rm = TRUE),mean.S1.=mean(S1,na.rm = TRUE),
                mean.S2.=mean(S2,na.rm = TRUE),mean.S3.=mean(S3,na.rm = TRUE),
                mean.SR.=mean(SR,na.rm = TRUE),mean.cgm2.=mean(cflux.gm2,na.rm = TRUE),
                mean.cmgm2=mean(cflux.mgm2,na.rm = TRUE),mean.Q=mean(Q.L.s,na.rm = TRUE), mean.bt.=mean(bt,na.rm = TRUE),mean.qf.=mean(qft.x,na.rm = TRUE),
                mean.bfper=mean(bfper,na.rm = TRUE),mean.qfper=mean(qfper,na.rm = TRUE))  
  return(mean)
}

# function for adding column of pre/post time to data
logstatus.f <- function(filename){
  filename$logstatus = ifelse(as.Date(filename$date) <= "2010-10-10 00:00:00", "pre", "post")
  return(filename)
}

# function for partitioning between wet and dry seasons
wetdry.f <- function(filename){
  filename$hydro <- ifelse(month(filename$date) %in% 4:9, "dry", "wet") 
  return(filename)
}


########################################################################################
###### Set paths for data
fig.dir <- "/Users/user/Dropbox/PhD Work/Thesis chapters/CR Chapter/CRDOCconcP1/Figures" # directory for saving figures
spectro.path <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Spectrodata/parfp_all"
discharge.path <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Stream/CR_discharge"
climate.path <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Climate"
groundwater.path <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Stream/CR_trutracks/SN_0904203"

# set path
#setwd("/Users/user/Dropbox/par and fp compilation")

##### Read in data + ensure format OK
#### DOC data from spectro - compiled file. No calculated spectral parameters
spectro.all<- read.csv(file = "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Spectrodata/SpectralOutput_2009t2014_2_noparameters.csv", head=TRUE,sep=",")
attach(spectro.all)
names(spectro.all)
spectro.all$date <- as.POSIXct(strptime(spectro.all$date, format = "%y-%m-%d %H:%M", tz = "America/Los_Angeles"))

##CHECK AS THE DATA FROM THE COMPILATION IS IN GMT!!!
spectro.all <- data.frame(spectro.all)

#### Discharge from stream
# Function for compiling discharge from all
# Note - took out the open air function to convert to 30 min. Not working?
source("/Users/user/SpecScripts/CRDischargeCompilation_function.R")
discharge <- discharge.comp(discharge.dir = discharge.path)
discharge$date <- as.POSIXct(strptime(discharge$date2, format = "%Y-%m-%d %H:%M"), tz = "America/Los_Angeles")
discharge.ave <- timeAverage(discharge, avg.time = "30 min", start.date = "2007-11-15", fill = TRUE)
attr(discharge.ave$date, "tzone") <- "America/Los_Angeles"
attr(discharge$date, "tzone") <- "America/Los_Angeles"

# read in discharge file from Mark
#mark.discharge <- source_data("https://github.com/UBCecohydro/ecohydro.datasets/blob/master/Campbell.River/Q.WQ.July2016.RDS?raw=true")
mark.discharge <- readRDS("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Stream/CR_discharge/Q.WQ.July2016.RDS")

# add column for pre/post harvest
mark.discharge <- logstatus.f(mark.discharge)
# add column for wet/dry  period
mark.discharge <- wetdry.f(mark.discharge)

#### Climate variables from Biomet
# Function for compiling climate from biomet
source("/Users/user/SpecScripts/CRClimateCompilation_function.R")
climate <- as.data.frame(climate.comp(climate.dir = climate.path))


#### Groundwater depth from groundwater well - trutrack closest to the weir
source("/Users/user/SpecScripts/CRGroundwaterComp_function.R")
groundwater <- groundwater.comp(ground.dir = groundwater.path)
# convert groundwater trace to 30 minute using open air, and convert time zone from GMT to Pacific
groundwater$date <- as.POSIXct(strptime(groundwater$date, format = "%Y-%m-%d %H:%M"), tz = "GMT")
groundwater.ave <- timeAverage(groundwater, avg.time = "30 min", start.date = "2010-09-22", fill = TRUE)
attr(groundwater.ave$date, "tzone") <- "America/Los_Angeles"

## Seven Anion Data from IC 
anions <- read.csv("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Stream/CR_IC/Compiled_ICSevenAnions_CR_Discharge_Sampledates_v4.csv", header = TRUE, sep = ",")
anions$sample <- anions$Sample.ID
# pad the sample ID for merging
anions$sample <- paste0("CR", str_pad(sapply(strsplit(as.character(anions$sample), "R"), "[", 2), 4, pad = "0"))
########################################################################################################################
# details for the plotting
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

### Default theme
theme = theme_set(theme_bw() + 
                    theme(strip.text = element_text(size=14),
                          axis.text=element_text(size=14, color = "black"), 
                          axis.title = element_text(size=14), legend.text = element_text(size=14),
                          plot.title = element_text(size=14)))
##############################
# Constants etc
area.hect <- 91 #watershed area in hectares

#######################
# Data preparation
########################

# calculate absorbance indicies from the spectro (timeseries) data
# first, correct abs spectra
source("/Users/user/SpecScripts/SpecAbsCorr_function.R")
abspec.corr <- abscorrspec(spectra = spectro.all)
# put back in date/DOCcorr/NO3
abspec.corr <- cbind(abspec.corr, spectro.all$date, spectro.all$DOCcorr, spectro.all$NO3.N)
colnames(abspec.corr)[222] <- "date"
colnames(abspec.corr)[223] <- "DOCcorr"
colnames(abspec.corr)[224] <- "NO3"

# calculate the abs spectral indicies from the corrected spectra
#source("/Users/user/SpecScripts/CRSpectralIndAbs_function.R")
#abs.ind.1 <- as.data.frame(Abs.ind(spec = abspec.corr[1:10000,], pathlength = 35))
#abs.ind.2 <- as.data.frame(Abs.ind(spec = abspec.corr[10001:20000,], pathlength = 35))
#abs.ind.3 <- as.data.frame(Abs.ind(spec = abspec.corr[20001:30000,], pathlength = 35))
#abs.ind.4 <- as.data.frame(Abs.ind(spec = abspec.corr[30001:39046,], pathlength = 35))
#spec.39047 <- Abs.ind(spec = abspec.corr[39047,], pathlength = 35)
#spec.39048 <- Abs.ind(spec = abspec.corr[39048,], pathlength = 35)
#spec.39049 <- Abs.ind(spec = abspec.corr[39049,], pathlength = 35)
#abs.ind.5 <- as.data.frame(Abs.ind(spec = abspec.corr[39050:dim(abspec.corr)[1],], pathlength = 35))

#spec.39048 <- c("NA", "NA","NA", "NA", "NA", "NA","NA", "NA", "NA" )
#spec.39049 <- c("NA", "NA","NA", "NA", "NA", "NA","NA", "NA", "NA" )

# bind together
#abs.all <- rbind(as.data.frame(abs.ind.1), abs.ind.2, abs.ind.3, abs.ind.4, 
#                 spec.39047, spec.39048, spec.39049, abs.ind.5)

# add in date
#abs.all$date <- abspec.corr$date
#abs.all$DOCcorr <- abspec.corr$DOCcorr
#abs.all$NO3 <- abspec.corr$NO3
#colnames(abs.all)[6:9] <- c("S1", "S2", "S3", "SlopeRatio")

# this takes forever! Save data to the absorbance folder...
#write.csv(abs.all, 
#          file = paste0("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Spectrodata", "/CRabsind_all.csv"))
#abs.all$date <- as.POSIXct(strptime(abs.all$date, format = "%Y-%m-%d %H:%M"))
#attr(abs.all$date, "tzone") <- "America/Los_Angeles"

## Load the pre-calculated spectral indicies
abs.all <- read.csv(file = paste0("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Spectrodata", "/CRabsind_all.csv"))
abs.all$date <- as.POSIXct(strptime(abs.all$date, format = "%Y-%m-%d %H:%M"))
attr(abs.all$date, "tzone") <- "America/Los_Angeles"

# merge the abs with the discharge/precip
abs.Q <- Reduce(function(x, y) merge(x, y, by = "date", all=TRUE), list(abs.all, mark.discharge, climate))
# add in column with pre/post, wet/dry
abs.Q$hydrolog <- paste(abs.Q$logstatus, abs.Q$hydro, sep = "/")
# add in month column
abs.Q$month <- format(abs.Q$date, "%m")

# Timeseries plots of spectral indicies
time.SUVA <- ggplot(abs.Q, aes(date, SUVA)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("SUVA") + theme()

time.e2e3 <- ggplot(abs.Q, aes(date, e2e3)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("e2e3") + theme()

time.e4e6 <- ggplot(abs.Q, aes(date, e4e6)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("e4e6") + theme()
  
time.SR <-   ggplot(abs.Q, aes(date, SlopeRatio)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("e4e6") + theme() + scale_y_continuous(limits = c(-5, 5))

## PCA on data - which variables change
## Correlation matrix - abs spectral/DOC conc/Q/bf/qf/pH/EC/Precip/soilT
CR.spec.abs <- na.omit(abs.Q[,c(3:6,11:16,20:26,29:30,31:32,35:36)])
corr.matrix <- cor(as.data.frame(CR.spec.abs[,c(1:6, 8:17, 20:23)]), use="pairwise.complete.obs") # correlation matric between variables
# plot the correlation matrix of the spectral parameters
pdf(file=paste0(fig.dir,"/CR_absindicies_correlations.pdf"), width = 11, height = 8.5)
  corrplot(na.omit(corr.matrix), method = "circle") #plot matrix
dev.off()
write.csv(corr.matrix, paste0(fig.dir, "/CRabsind_corr.csv"))

# Linear relationships between DOC; e2e3, e4e6, SUVA, Slope Ratio - linear model

## ANOVA/Box plots on pre/post for significant changes

######################################################################
####### Fluorescence data
# get indicies
CR.Flindicies <- read.csv("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR EEMS/CR_fluorescence/CREEMS_corrIFE_RM_Ram/CRFluorescence_Indicies.csv", header = TRUE)
CR.Flindicies$sample <- CR.Flindicies$samplename

# get the CM fits
CR.CM.Fmax <- read.csv("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR EEMS/CR_CMPARAFACResults/CR_componentsandloadings_CM_Fmax.csv")
CR.CM.per <- read.csv("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR EEMS/CR_CMPARAFACResults/CR_componentsandloadings_CM.csv")
# get sample ID column for both
CR.CM.Fmax$sample <- sapply(strsplit(as.character(CR.CM.Fmax$sample.ID), "CRCorr"), "[", 1)
CR.CM.per$sample <- sapply(strsplit(as.character(CR.CM.per$sample.ID), "CRCorr"), "[", 1)

# Get PARAFAC fits 
CR.4comp <- read.csv("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR EEMS/CR_PARAFAC_4comp/FMax_4.csv", header = FALSE)
# add in the title C1-4
colnames(CR.4comp) <- c("PARAFAC_C1", "PARAFAC_C2", "PARAFAC_C3", "PARAFAC_C4")
# add in sample ID.
CRkey <- t(as.data.frame(read.csv("/Users/user/Documents/MATLAB/toolbox/CorrEEMS/CR/CR_01key.csv", header = FALSE)))
CR.4comp$sample <- CRkey
# calculate sum, and percentage of each as the total fluorescence
CR.4comp$sum <- rowSums(CR.4comp[,1:4])
CR.4comp.per <- CR.4comp[,1:4]/CR.4comp$sum*100
colnames(CR.4comp.per) <- c("C1_per", "C2_per", "C3_per", "C4_per")
CR.4comp <- cbind(CR.4comp, CR.4comp.per)

# merge all fluorescence data together by sample ID
CR.fl <- Reduce(function(x, y) merge(x, y, by = "sample", all=TRUE), list(CR.Flindicies, CR.CM.per, CR.CM.Fmax, CR.4comp))
# pad the sample ID to four digits
CR.fl$sample <- paste0("CR", str_pad(sapply(strsplit(as.character(CR.fl$sample), "R"), "[", 2), 4, pad = "0"))

# merge flruoescenc data with the anion data from the IC
CR.grab <- merge(CR.fl, anions, by = "sample", all = TRUE)

# for grab sample data, compile by sample ID, and ensure sample date/time correctly associated.
# Fluorescence (CM), Fluorescence (PARAFAC), Fluorescence (Indicies)
# get file with the sample ID and date
sample.date <- read.csv("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Stream/CR_sampledate_grabsample/CR_discharge_sampledate.csv", header = TRUE)
sample.date$sample <- paste0("CR", str_pad(sapply(strsplit(as.character(sample.date$Sample.ID), "R"), "[", 2), 4, pad = "0"))
sample.date$date <- as.POSIXct(strptime(sample.date$Sample.Date, format = "%Y-%m-%d %H:%M", tz = "GMT"))
attr(sample.date$date, "tzone") <- "America/Los_Angeles" # change timezeon from GMT to PST

# change the sample.date to 30 min
#sample.date.ave <-  timeAverage(sample.date, avg.time = "30 min", start.date = "2007-11-14 12:30", fill = TRUE)
#attr(sample.date.ave$date, "tzone") <- "America/Los_Angeles"
sample.date$datenew = sample.date$date - 60*10
sample.date$dateoriginal <- sample.date$date
sample.date$date <- sample.date$datenew

# merge with discharge
#sample.date.Q <- merge(sample.date, discharge, by = "date", all = TRUE)
#write.csv(sample.date.Q, file = "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/sampledateQTRUE.csv")

#sample.date.Qave <- merge(sample.date, mark.discharge, by = "date", all = TRUE)
#write.csv(sample.date.Qave, file = "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/sampledateQmark.csv")
sample.date.QaveFALSE <- merge(sample.date, mark.discharge, by = "date", all = FALSE)
#write.csv(sample.date.QaveFALSE, file = "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/sampledateQmarkFALSE.csv")
# looks good! use this
discharge.true <- sample.date.QaveFALSE
# Merge into the anion and fluorescence data
CR.grab <- merge(CR.grab, discharge.true, by = "sample", all = TRUE)
CR.grab$date <- as.POSIXct(strptime(CR.grab$Sample.Date.y, format = "%Y-%m-%d %H:%M", tz = "America/Los_Angeles"))

# Add in the pre/post logging
CR.grab <- logstatus.f(CR.grab)
# Add in wet/dry
CR.grab <- wetdry.f(CR.grab)
# add in the time of day (diurnal cycles)
# add in the month
CR.grab$month <- format(CR.grab$date,"%m")
# Add in pre/post + wet/dry
CR.grab$hydrolog <- paste(CR.grab$logstatus, CR.grab$hydro, sep = "/")

# clean up the dataframe 
CR.grab$Sample.ID <- NULL
CR.grab$Sample.ID.x <- NULL
CR.grab$Sample.ID.y <- NULL
CR.grab$datenew.x <- NULL
CR.grab$dateoriginal.x <- NULL
CR.grab$date.x <- NULL
CR.grab$datenew.y <- NULL
CR.grab$date.y <- NULL
CR.grab$dateoriginal.y <- NULL
CR.grab$samplename <- NULL
CR.grab$Batch <- NULL
CR.grab$OK.Not <- NULL
CR.grab$Notes <- NULL

# add in the discharge/quickflow/stormflow to calculate flow weighted mean
CR.grab$hour <- format(CR.grab$date,"%H")

# calculated the discharge weighted mean for grab samples - c/Q
CR.grab.DWM.Ls <- as.data.frame(cbind(CR.grab[,c(2:13,15:29,44:47,52:58,69,70,71,72)]/CR.grab[,63], CR.grab$logstatus))
CR.grab.DWM.mmday <- CR.grab[,c(2:13,15:29,44:47,52:58,69,70,71,72)]/CR.grab[,67]
CR.grab.DWM.Ls$Br_mgL <- NULL
CR.grab.DWM.Ls$OFI <- NULL
CR.grab.DWM.Ls$F_mgL <- NULL

################## PCA Analysis
# Figure 1 PCA analysis
# Question: Which varaibles explain the most variance within the dataset (partitioned according to pre/post, wet/dry)?
# Do PCA on Flow-Weighted Fluorescence Parameters
CR.grab.fl <- na.omit(CR.grab[,c(2:12,15:29, 44:47, 53,54,56,58,63, 78)])
CR.grab.fl.FWM <- CR.grab.fl[,c(1:34)]/CR.grab.fl[,35]
CR.grab.fl.FWM$logstatus <- CR.grab.fl$logstatus

CR.grab.DWM.PCA <- prcomp(na.omit(CR.grab.fl.FWM[,c(1:31, 33:34)]), center = TRUE, scale. = TRUE, na.action=na.omit)
summary(CR.grab.DWM.PCA)

# plot PCA results
pdf(file=paste0(fig.dir,"/CR_PCADWMLs.pdf"), width = 11, height = 8.5)
ggbiplot(CR.grab.DWM.PCA, obs.scale = 1, var.scale = 1, groups = na.omit(CR.grab.fl.FWM$logstatus),
         ellipse = TRUE, circle = FALSE, varname.abbrev = FALSE) +
  #scale_colour_manual(values=cbPalette[1:6], name="Region") +
  theme(legend.direction = 'vertical', legend.position = 'right') 
dev.off()

# show relative contribution of each variable to the first 5 components of PCA
fl.FWM.pca <- PCA(CR.grab.fl.FWM[,c(1:31, 33:34)], graph = TRUE)
head(fl.FWM.pca$var$contrib)
PCA.contrib <- data.frame(fl.FWM.pca$var$contrib)
# sort variables by incresing contribution across the first 5 component
wq.sortvar <- PCA.contrib[order(-PCA.contrib$Dim.1,-PCA.contrib$Dim.2,-PCA.contrib$Dim.3,-PCA.contrib$Dim.4,-PCA.contrib$Dim.5), ]
write.csv(wq.sortvar, file = paste(fig.dir, "/PCAcontributions_FWM.csv", sep ="/")) #write contributions to a csv file that are sorted
# scree plot
fviz_screeplot(fl.FWM.pca, ncp=6) # first 6 components
# plot contribution to first 2 PCA components
fviz_contrib(fl.FWM.pca, choice = "var", axes = 1)
fviz_contrib(fl.FWM.pca, choice = "var", axes = 2)

############### Do PCA on non-flowweighted means - only fluorescence parameters
CR.grab.fl <- na.omit(CR.grab[,c(2,3,7:11,28:29,48:52,58,63,68:70,83:84,86)])
CR.grab.PCA <- prcomp(na.omit(CR.grab.fl[,1:19]), center = TRUE, scale. = TRUE, na.action=na.omit)
summary(CR.grab.PCA)

# plot PCA results
pdf(file=paste0(fig.dir,"/CR_PCA_grabfl.pdf"), width = 11, height = 8.5)
ggbiplot(CR.grab.PCA, obs.scale = 1, var.scale = 1, groups = na.omit(CR.grab.fl$hydrolog),
         ellipse = TRUE, circle = TRUE, varname.abbrev = FALSE) +
  scale_colour_manual(values=cbPalette[1:6], name="logstatus") +
  theme(legend.direction = 'vertical', legend.position = 'right') 
dev.off()

# show relative contribution of each variable to the first 5 components of PCA
grabfl.pca <- PCA(CR.grab.fl[,1:19], graph = TRUE)
head(grabfl.pca$var$contrib)
PCA.contrib <- data.frame(grabfl.pca$var$contrib)
# sort variables by incresing contribution across the first 5 component
wq.sortvar <- PCA.contrib[order(-PCA.contrib$Dim.1,-PCA.contrib$Dim.2,-PCA.contrib$Dim.3,-PCA.contrib$Dim.4,-PCA.contrib$Dim.5), ]
write.csv(wq.sortvar, file = paste(fig.dir, "/PCAcontributions_flgrab.csv", sep ="/")) #write contributions to a csv file that are sorted
# scree plot
fviz_screeplot(grabfl.pca, ncp=6) # first 6 components
# plot contribution to first 2 PCA components
fviz_contrib(grabfl.pca, choice = "var", axes = 1)
fviz_contrib(grabfl.pca, choice = "var", axes = 2)

## PCA Only on fluorescence and absorbance parameters
CR.grab.fl <- na.omit(CR.grab[,c(2,3,7:11,28:29,49:52,83:84,86)])
CR.grab.PCA <- prcomp(na.omit(CR.grab.fl[,1:13]), center = TRUE, scale. = TRUE, na.action=na.omit)
summary(CR.grab.PCA)

# plot PCA results
pdf(file=paste0(fig.dir,"/CR_PCA_grabfl.pdf"), width = 11, height = 8.5)
ggbiplot(CR.grab.PCA, obs.scale = 1, var.scale = 1, groups = na.omit(CR.grab.fl$hydrolog),
         ellipse = TRUE, circle = TRUE, varname.abbrev = FALSE) +
  scale_colour_manual(values=cbPalette[1:6], name="logstatus") +
  theme(legend.direction = 'vertical', legend.position = 'right') 
dev.off()

# show relative contribution of each variable to the first 5 components of PCA
grabfl.pca <- PCA(CR.grab.fl[,1:13], graph = TRUE)
head(grabfl.pca$var$contrib)
PCA.contrib <- data.frame(grabfl.pca$var$contrib)
# sort variables by incresing contribution across the first 5 component
wq.sortvar <- PCA.contrib[order(-PCA.contrib$Dim.1,-PCA.contrib$Dim.2,-PCA.contrib$Dim.3,-PCA.contrib$Dim.4,-PCA.contrib$Dim.5), ]
write.csv(wq.sortvar, file = paste(fig.dir, "/PCAcontributions_flgrab.csv", sep ="/")) #write contributions to a csv file that are sorted
# scree plot
fviz_screeplot(grabfl.pca, ncp=6) # first 6 components
# plot contribution to first 2 PCA components
fviz_contrib(grabfl.pca, choice = "var", axes = 1)
fviz_contrib(grabfl.pca, choice = "var", axes = 2)

##############################################################################################################
# Figure Correlation Matrix - Water Quality variables and others (discharge, climate)
CR.grab.fl <- na.omit(CR.grab[,c(2:12,15:29, 44:47, 63:65,78)])
corr.matrix <- cor(CR.grab.fl[,1:33]) # correlation matric between variables
# plot the correlation matrix of the spectral parameters
pdf(file=paste0(fig.dir,"/CR_flindicies_correlations.pdf"), width = 11, height = 8.5)
corrplot(na.omit(corr.matrix), method = "circle") #plot matrix
dev.off()

# Correlation matrix with the high frequency data
CR.abs.cor <- na.omit(abs.Q[,c(3,11,12,13,14:16,20:26,31:32,35:36)])
corr.matrix <- cor(CR.abs.cor) # correlation matric between variables
# plot the correlation matrix of the spectral parameters
pdf(file=paste0(fig.dir,"/CR_abs_correlations.pdf"), width = 11, height = 8.5)
corrplot(na.omit(corr.matrix), method = "circle") #plot matrix
dev.off()
write.csv(corr.matrix, paste0(fig.dir, ("/CRcorrmatrix_abs.csv"))) #write correlation matrix
#########################################################################################################
# Figure - ANOVA analysis to look at significance of changes in pre and post 
# Boxplot of flow weighted means for the pre/post period per month for different parameters
# do anova test on each of the pairs?

# SUVA
SUVA.box <- ggplot(abs.Q, aes(x=month, y=SUVA, fill=logstatus)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('SUVA')) +
  ylab(expression('SUVA')) +
  xlab(expression('Month')) +
  coord_cartesian(xlim = c(1,12), ylim=c(2,5)) 
SR.box <- ggplot(abs.Q, aes(x=month, y=SlopeRatio, fill=logstatus)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('Slope Ratio - by month')) +
  ylab(expression('Slope Ratio')) +
  xlab(expression('Month')) +
  coord_cartesian(xlim = c(1,12), ylim=c(0,2.5))
FI.box <- ggplot(CR.grab, aes(x=month, y=FI, fill=logstatus)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('FI')) +
  ylab(expression('FI')) +
  xlab(expression('Month')) 
HIX.box <- ggplot(subset(CR.grab, HIX_ohno_area >0), aes(x=month, y=HIX_ohno_area, fill=logstatus)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('HIX')) +
  ylab(expression('HIX')) +
  xlab(expression('Month')) +
  coord_cartesian(xlim = c(1,12), ylim=c(0,2))
FrI.box <- ggplot(subset(CR.grab, FrI >0), aes(x=month, y=FrI, fill=logstatus)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('Freshness Index')) +
  ylab(expression('FrI (A.U)')) +
  xlab(expression('Month')) +
  coord_cartesian(xlim = c(1,12), ylim=c(0,1.5))
Cl.box <- ggplot(subset(CR.grab, Cl_mgL >0), aes(x=month, y=Cl_mgL, fill=logstatus)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('Cl')) +
  ylab(expression('Cl (mg/L)')) +
  xlab(expression('Month')) +
  coord_cartesian(xlim = c(1,12), ylim=c(0,4))
SO4.box <- ggplot(subset(CR.grab, SO4_S_mgL >0), aes(x=month, y=SO4_S_mgL, fill=logstatus)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('SO4')) +
  ylab(expression('SO4 (mg/L)')) +
  xlab(expression('Month')) +
  coord_cartesian(xlim = c(1,12), ylim=c(0,1.5))
Redox.box <- ggplot(subset(CR.grab, redox >0), aes(x=month, y=redox, fill=logstatus)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('redox Index')) +
  ylab(expression('redox Index (A.U)')) +
  xlab(expression('Month')) +
  coord_cartesian(xlim = c(1,12), ylim=c(0.2,0.6))
PerProtein.box <- ggplot(subset(CR.grab, perprotein >0), aes(x=month, y=perprotein, fill=logstatus)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('Percent Protein')) +
  ylab(expression('Percent Protein (%)')) +
  xlab(expression('Month')) +
  coord_cartesian(xlim = c(1,12), ylim=c(0,25))
C1.box <- ggplot(subset(CR.grab, C1_per >0), aes(x=month, y=C1_per, fill=logstatus)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('C1 %')) +
  ylab(expression('C1 Fmax (%)')) +
  xlab(expression('Month')) +
  coord_cartesian(xlim = c(1,12), ylim=c(10,60))
C2.box <- ggplot(subset(CR.grab, C2_per >0), aes(x=month, y=C2_per, fill=logstatus)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('C2 %')) +
  ylab(expression('C2 Fmax (%)')) +
  xlab(expression('Month')) +
  coord_cartesian(xlim = c(1,12), ylim=c(15,22))
C3.box <- ggplot(subset(CR.grab, C3_per >0), aes(x=month, y=C3_per, fill=logstatus)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('C3 %')) +
  ylab(expression('C3 Fmax (%)')) +
  xlab(expression('Month')) +
  coord_cartesian(xlim = c(1,12), ylim=c(0,25))
C4.box <- ggplot(subset(CR.grab, C4_per >0), aes(x=month, y=C4_per, fill=logstatus)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('C4 %')) +
  ylab(expression('C4 Fmax (%)')) +
  xlab(expression('Month')) +
  coord_cartesian(xlim = c(1,12), ylim=c(25,35))
###### Save as PDF files
pdf(file=paste0(fig.dir,"/CRFigures_Boxind.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(SUVA.box, SR.box, FI.box, FrI.box, HIX.box, Cl.box, SO4.box, ncol = 2)
dev.off()
# save fluorescence parameters in another
pdf(file=paste0(fig.dir,"/CRFigures_Boxfluorescence.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(Redox.box, PerProtein.box, C1.box, C2.box, C3.box, C4.box, ncol = 2)
dev.off()

## ANOVA on pre/post by month

#################################################################
## Conc/Q relationships: pre/post; wet/dry. Look for significant relationship
# DOC
DOC.cQ <- ggplot(abs.Q, aes(x=Q.mm.d, y=DOCcorr, color = hydrolog)) +
  scale_color_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(method="lm",aes(fill = hydrolog)) + 
  labs(x="Q (mm/day)", y=expression(paste("DOC (mg L"[-1]))) +
  labs(title="cQ plot of Spectral DOC vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) +
  coord_cartesian(xlim = c(0,13), ylim=c(0,10)) 
# SUVA
SUVA.cQ <- ggplot(subset(abs.Q,Q.mm.d > 1), aes(x=Q.mm.d, y=SUVA, color = hydrolog)) +
  scale_color_manual(values=cbPalette[1:4], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:4], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(method="lm",aes(fill = hydrolog)) + 
  labs(x="Q (mm/day)", y=expression(paste("SUVA"[254]))) +
  labs(title="cQ plot of Spectral SUVA vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) +
  coord_cartesian(xlim = c(0,13), ylim=c(3,4)) 
# SR 
SR.cQ <- ggplot(subset(abs.Q,Q.mm.d > 1), aes(x=Q.mm.d, y=SlopeRatio, color = hydrolog)) +
  scale_color_manual(values=cbPalette[1:4], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:4], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(method="lm",aes(fill = hydrolog)) + 
  labs(x="Q (mm/day)", y=expression(paste("S"[R]))) +
  labs(title="cQ plot of Spectral Slope Ratio vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) + 
  coord_cartesian(xlim = c(0,13), ylim=c(0.9,1.5)) 
## Fluorescence - FI
FI.cQ <- ggplot(CR.grab, aes(x=Q.mm.d, y=FI, color = hydrolog)) +
  scale_color_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(method="lm",aes(fill = hydrolog)) + 
  labs(x="Q (mm/day)", y=expression(paste("FI"))) +
  labs(title="cQ plot of Spectral FI vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) + 
  coord_cartesian(xlim = c(0,13), ylim=c(1.25,2)) 

HIX.cQ <- ggplot(CR.grab, aes(x=Q.mm.d, y=HIX_ohno_area, color = hydrolog)) +
  scale_color_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(method="lm",aes(fill = hydrolog)) + 
  labs(x="Q (mm/day)", y=expression(paste("FI"))) +
  labs(title="cQ plot of Spectral HIX vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) + 
  coord_cartesian(xlim = c(0,13), ylim=c(-0.5,2)) 

FrI.cQ <- ggplot(CR.grab, aes(x=Q.mm.d, y=FrI, color = hydrolog)) +
  scale_color_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(method="lm",aes(fill = hydrolog)) + 
  labs(x="Q (mm/day)", y=expression(paste("Freshness Index"))) +
  labs(title="cQ plot of Spectral FrI vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) + 
  coord_cartesian(xlim = c(0,13), ylim=c(0.1,1)) 
# anions
Cl.cQ <- ggplot(subset(CR.grab, Cl_mgL >0.1), aes(x=Q.mm.d, y=Cl_mgL, color = hydrolog)) +
  scale_color_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(method="lm",aes(fill = hydrolog)) + 
  labs(x="Q (mm/day)", y=expression(paste("Cl (mg/L)"))) +
  labs(title="cQ plot of Spectral Cl vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) + 
  coord_cartesian(xlim = c(0,13), ylim=c(0,20)) 

SO4.cQ <- ggplot(subset(CR.grab, SO4_S_mgL >0.1), aes(x=Q.mm.d, y=SO4_S_mgL, color = hydrolog)) +
  scale_color_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(method="lm",aes(fill = hydrolog)) + 
  labs(x="Q (mm/day)", y=expression(paste("SO4 (mg/L)"))) +
  labs(title="cQ plot of Spectral SO4 vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) + 
  coord_cartesian(xlim = c(0,13), ylim=c(0,1)) 

## Flourescence - Redox
Redox.cQ <- ggplot(CR.grab, aes(x=Q.mm.d, y=redox, color = hydrolog)) +
  scale_color_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(method="lm",aes(fill = hydrolog)) + 
  labs(x="Q (mm/day)", y=expression(paste("Redox Index (A.U)"))) +
  labs(title="cQ plot of Spectral Redox Index (A.U) vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) + 
  coord_cartesian(xlim = c(0,13), ylim=c(0.25,0.4)) 

## FLruorescence - Per Protein
PerProtein.cQ <- ggplot(CR.grab, aes(x=Q.mm.d, y=perprotein, color = hydrolog)) +
  scale_color_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(method="lm",aes(fill = hydrolog)) + 
  labs(x="Q (mm/day)", y=expression(paste("Percent Protein (%)"))) +
  labs(title="cQ plot of Spectral Percent Protein (%) vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) + 
  coord_cartesian(xlim = c(0,13), ylim=c(0,10)) 
## Fluorescence - P1
C1.cQ <- ggplot(subset(CR.grab, PARAFAC_C1 >0), aes(x=Q.mm.d, y=C1_per, color = hydrolog)) +
  scale_color_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(method="lm",aes(fill = hydrolog)) + 
  labs(x="Q (mm/day)", y=expression(paste("% C1"))) +
  labs(title="cQ plot of Spectral % C1 vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) + 
  coord_cartesian(xlim = c(0,13), ylim=c(20,45)) 

## Fluorescence - P2
C2.cQ <- ggplot(subset(CR.grab, PARAFAC_C2 >0), aes(x=Q.mm.d, y=C2_per, color = hydrolog)) +
  scale_color_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(method="lm",aes(fill = hydrolog)) + 
  labs(x="Q (mm/day)", y=expression(paste("C2 %"))) +
  labs(title="cQ plot of Spectral C2% vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) + 
  coord_cartesian(xlim = c(0,13), ylim=c(18,21)) 
## Fluorescence - P3
C3.cQ <- ggplot(subset(CR.grab, PARAFAC_C3 >0), aes(x=Q.mm.d, y=C3_per, color = hydrolog)) +
  scale_color_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(method="lm",aes(fill = hydrolog)) + 
  labs(x="Q (mm/day)", y=expression(paste("C3%"))) +
  labs(title="cQ plot of Spectral C3% vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) + 
  coord_cartesian(xlim = c(0,13), ylim=c(15,24)) 

## Fluorescence - P4
C4.cQ <- ggplot(subset(CR.grab, PARAFAC_C4 >0), aes(x=Q.mm.d, y=C4_per, color = hydrolog)) +
  scale_color_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(method="lm",aes(fill = hydrolog)) + 
  labs(x="Q (mm/day)", y=expression(paste("C4%"))) +
  labs(title="cQ plot of Spectral C4% vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) + 
  coord_cartesian(xlim = c(0,13), ylim=c(24,31)) 

## save all as a pdf together
pdf(file=paste0(fig.dir,"/CRFigures_CQrel.pdf"), width = 8.5, height = 11) #save figure
DOC.cQ 
SR.cQ 
SUVA.cQ 
FI.cQ 
HIX.cQ 
FrI.cQ 
Cl.cQ
SO4.cQ 
Redox.cQ 
PerProtein.cQ
C1.cQ
C2.cQ
C3.cQ
C4.cQ
dev.off()

# Linear models? 

############################################################################################
# Soil Extracts
# linear model between soil depth and quality/DOC variables. Which variables best explain changes in soil depth?


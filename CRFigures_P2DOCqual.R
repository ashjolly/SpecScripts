# Script for figures associated with CR Paper 2 - Quality dynamics
# Ashlee J's PhD
# 27 June2016
############################################################################################

# Clear environment
rm(list = ls())
ls()

###### Necessary toolboxes
library(pacman)
p_load(lubridate,plyr,ggplot2,reshape2,car,grid,gridBase,gridExtra, taRifx,magrittr, ggpmisc,
       stringi,stringr,plyr,dplyr,tidyr,EcoHydRology,openair,reshape,lmomco, zoo, hydroTSM, rowr, 
       reshape2, pgirmess, nlme, chron, repmis, ggbiplot, FactoMineR, factoextra, stringr, rms, splines, gtable,
       'gsubfn', 'abind', 'zoo', 'corrplot', 'gplots', "EEM", eeptools, MASS, Hmisc)

library('devtools')
install_github("ggbiplot", "vqv")
library(ggbiplot)
devtools::install_github("PMassicotte/eemR")

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
fig.dir <- "/Users/user/Dropbox/PhD Work/Thesis chapters/CR Chapter/Figures" # directory for saving figures
spectro.path <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Spectrodata/parfp_all"
discharge.path <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Stream/CR_discharge"
climate.path <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Climate"
groundwater.path <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Stream/CR_trutracks/SN_0904203"
soilextract.path <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_soilextracts"
soillys.path <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Lysimeter"
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
discharge1 <- discharge.comp(discharge.dir = discharge.path)
discharge1$date <- as.POSIXct(strptime(discharge$date2, format = "%Y-%m-%d %H:%M"), tz = "America/Los_Angeles")
discharge.ave <- timeAverage(discharge, avg.time = "30 min", start.date = "2007-11-15", fill = TRUE)
attr(discharge.ave$date, "tzone") <- "America/Los_Angeles"
attr(discharge1$date, "tzone") <- "America/Los_Angeles"

# read in discharge file from Mark
#mark.discharge <- source_data("https://github.com/UBCecohydro/ecohydro.datasets/blob/master/Campbell.River/Q.WQ.July2016.RDS?raw=true")
mark.discharge <- readRDS("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Stream/CR_discharge/Q.WQ.July2016.RDS")
mark.discharge$date <- as.POSIXct(strptime(mark.discharge$date, format = "%Y-%m-%d %H:%M"), tz = "America/Los_Angeles")
#discharge.ave <- timeAverage(discharge, avg.time = "30 min", start.date = "2007-11-15", fill = TRUE)
attr(mark.discharge$date, "tzone") <- "America/Los_Angeles"
discharge <- mark.discharge #rename the discharge from Mark

# add column for pre/post harvest
discharge <- logstatus.f(discharge)
# add column for wet/dry  period
discharge <- wetdry.f(discharge)

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
abspec.corr <- cbind(abspec.corr, spectro.all$date, spectro.all$DOCcorr, spectro.all$NO3.N, spectro.all$SAC254)
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
abs.all$abs254 <- spectro.all$SAC254

# merge the abs with the discharge/precip
abs.Q <- Reduce(function(x, y) merge(x, y, by = "date", all=TRUE), list(abs.all, mark.discharge, climate))
# add in column with pre/post, wet/dry
abs.Q <- logstatus.f(abs.Q)
abs.Q <- wetdry.f(abs.Q)
abs.Q$hydrolog <- paste(abs.Q$logstatus, abs.Q$hydro, sep = "/")
# add in month column
abs.Q$month <- format(abs.Q$date, "%m")

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

############################################################
# Figure 1 - timeseries of in situ parameters - SUVA, slope ratio, doc concentration etc.
# Plot with Q, [DOC], SUVA, e2e3, e4e5, SR, and PARAFAC components
############################################################

# Timeseries plots of spectral indicies + DOC + Q
time.Q <- ggplot(subset(abs.Q, DOCcorr >=1), aes(date, Q.L.s)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("Q (L/s)") + theme()

time.DOC <- ggplot(subset(abs.Q, DOCcorr >=1), aes(date, DOCcorr)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("[DOC] (mg/L)") + theme()

time.abs254 <- ggplot(subset(abs.Q, DOCcorr >=1), aes(date, abs254)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("abs254") + theme()

time.SUVA <- ggplot(subset(abs.Q, DOCcorr >=1), aes(date, SUVA)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("SUVA") + theme()

time.e2e3 <- ggplot(abs.Q, aes(date, e2e3)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("e2e3") + theme()

time.e4e6 <- ggplot(abs.Q, aes(date, e4e6)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("e4e6") + theme()+ scale_y_continuous(limits = c(-5, 1000))
  
time.SR <-   ggplot(subset(abs.Q, DOCcorr >=1), aes(date, SlopeRatio)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("Slope Ratio") + theme() + scale_y_continuous(limits = c(0, 5))

# timeseries look really odd?!! Jump in data.
# show as boxplots?

PARA.perprotein <- melt(CR.grab[,c(64, 28)], id = "date")
pd <- position_dodge(.65)
time.perprotein <- ggplot(PARA.perprotein, aes(date, value, colour = variable, shape = variable)) + 
  geom_point(position = pd, size = 2) +
  scale_color_manual(values=cbPalette[1]) +
  xlab("Date") + ylab("% Protein (CM PARAFAC Model)") + theme(legend.position="top")  +
  ylim(0, 50) 

# redox index
PARA.redox <- melt(CR.grab[,c(64, 29)], id = "date")
pd <- position_dodge(.65)
time.redox <- ggplot(PARA.redox, aes(date, value, colour = variable, shape = variable)) + 
  geom_point(position = pd, size = 2) +
  scale_color_manual(values=cbPalette[2]) +
  xlab("Date") + ylab("Redox Index (CM PARAFAC Model)") + theme(legend.position="top")  +
  ylim(0, 0.8) 

# PARAFAC - 4 component %
PARA.per <- melt(CR.grab[,c(64, 49:52)], id = "date")
pd <- position_dodge(.65)
time.PARAFACPer <- ggplot(PARA.per, aes(date, value, colour = variable, shape = variable)) + 
  geom_point(position = pd, size = 2) +
  scale_color_manual(values=cbPalette[1:4], 
                     name="PARAFAC Component",
                     breaks=c("C1_per", "C2_per", "C3_per", "C4_per"),
                     labels=c("C1", "C2", "C3", "C4")) +
  xlab("Date") + ylab("PARAFAC Component %") + theme(legend.position="top")  +
  ylim(0, 70) +
  scale_shape_manual(values=c(17,18,3,4), 
                     name="PARAFAC Component",
                     breaks=c("C1_per", "C2_per", "C3_per", "C4_per"),
                     labels=c("C1", "C2", "C3", "C4")) 
# PARAFAC fmax
PARA.fmax <- melt(CR.grab[,c(64, 44:47)], id = "date")
time.PARAFACfmax <- ggplot(PARA.fmax, aes(date, value, color = variable, shape = variable)) + 
  geom_point(position = pd, size = 2) +
  xlab("Date") + ylab("PARAFAC Component Fmax") + 
  theme(legend.position="top")  +
  ylim(0, 1) +
  scale_color_manual(values=cbPalette[1:4], 
                     name="PARAFAC Component",
                     breaks=c("PARAFAC_C1", "PARAFAC_C2", "PARAFAC_C3", "PARAFAC_C4"),
                     labels=c("C1", "C2", "C3", "C4")) +
  scale_shape_manual(values=c(17,18,3,4), 
                     name="PARAFAC Component",
                     breaks=c("PARAFAC_C1", "PARAFAC_C2", "PARAFAC_C3", "PARAFAC_C4"),
                     labels=c("C1", "C2", "C3", "C4")) 

# Arrange the timeseries on top of each other
pdf(file=paste0(fig.dir,"/CRFigures_QualityF1_TimeseriesAbs.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(time.Q, time.DOC, time.abs254, time.SUVA, time.SR, ncol = 1)
dev.off()

############################################# 
# Figure 2 - Timeseries of fluorescence data - PARAFAC components over time.

# Arrange timeseries for PARAFAC data
pdf(file=paste0(fig.dir,"/CRFigures_QualityF1_Timeseriesfluor.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(time.Q, time.DOC, time.perprotein, time.redox, time.PARAFACPer, time.PARAFACfmax, ncol = 1)
dev.off()

############################################# 
# plot example EEMto show OM fluorescence characteristics
# possible nice EEMs to use: 1017, 1018, 1021, 1027, 1034
EEMdir <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR EEMS/CR_fluorescence/CREEMS_corrIFE_RM_Ram"
EEM <- read.csv(file = paste0(EEMdir, "/CR1031_CR_Corrected.csv"), header = TRUE)
samplename  <- "CR1031"

#variables to change
zmax = max(EEM,na.rm=TRUE) # put the max intensity of that you want to graph
#EEMmax[i] <- zmax #to show the maximum fluorescence for all files
xlimit <- range(300, 500, finite=TRUE)
ylimit <- range(240, 450, finite = TRUE)

numcont = 20 # number of contour levels you want: Change if you want

##### contour plotting function
# call contour plot function
setwd("/Users/user/SpecScripts") 
source("EEM_contour_v1.R")

#Plot contours and save in correction file
plotpath <- file.path(fig.dir, paste(samplename,"_ContourRaleigh.jpeg", sep = ""))

EEMplot <- EEM # not cutting out the last two columns
EEMplot[EEMplot<0]=0 # have line so that if fluorescence is less than 0, becomes 0.

explot = colnames(as.matrix(EEMplot))
explot = as.numeric(gsub("X","",explot))
emplot = as.numeric(row.names(EEMplot))

jpeg(file=plotpath)
contour.plots(eems = as.matrix(EEMplot), Title = samplename, ex = explot, em = emplot, 
              zmax = zmax, zmin = 0, numcont = numcont)  
dev.off()

############################################# 
# Figure 3 - Diurnal and seasonal changes
# Figure 3A - Possibility of diurnal changes
# Timeseries of abs within dry summer months

# plot July and August for the time
julyaug <- subset(abs.Q, format.Date(abs.Q$date, "%m") == "07" | format.Date(abs.Q$date, "%m")=="08")
# add in air temp from biomet
julyaug$day <-  julyaug$date
year(julyaug$day) <- 0
julyaug.abs254 <- ggplot(subset(julyaug, julyaug$DOCcorr >=1), aes(day, abs254)) + 
  geom_point(size = 0.4) +
  xlab("Date") + ylab("abs254") + 
  #scale_color_manual(values=cbPalette[1:2]) + # colours
  theme(legend.position="top") + theme() + ggtitle("Diurnal Cycles: July - August abs254 ")

# SUVA
julyaug.SUVA <- ggplot(subset(julyaug, julyaug$DOCcorr >=1), aes(day, SUVA)) + 
  geom_point(size = 0.4) +
  xlab("Date") + ylab("SUVA") + 
  #scale_color_manual(values=cbPalette[1:2]) + # colours
  theme(legend.position="top") + theme() + ggtitle("Diurnal Cycles: July - August SUVA ")

# SR
julyaug.SR <- ggplot(subset(julyaug, julyaug$DOCcorr >=1), aes(day, SR)) + 
  geom_point(size = 0.4) +
  xlab("Date") + ylab("Slope Ratio") + 
  #scale_color_manual(values=cbPalette[1:2]) + # colours
  theme(legend.position="top") + theme() + ggtitle("Diurnal Cycles: July - August Slope Ratio ")
# abs 350 - need to calculate
#julyaug.SR <- ggplot(subset(julyaug, julyaug$DOCcorr >=1), aes(day, SR)) + 
#  geom_point(size = 0.4) +
 # xlab("Date") + ylab("Slope Ratio") + 
#  #scale_color_manual(values=cbPalette[1:2]) + # colours
#  theme(legend.position="top") + theme() + ggtitle("Diurnal Cycles: July - August Slope Ratio ")

# Do by hour - last week of july
# Show the DOC/Q/airT over one day period - take last week of july
julDOC <- subset(julyaug, format.Date(julyaug$date, "%m") == "07" & format.Date(julyaug$date, "%d") >= "25" & format.Date(julyaug$date, "%d")<="31")
# take the mean hourly DOC concentration by year (3 different years)
meanJulydoc <- ddply(julDOC,.(format(julDOC$date, format='%y')),
                     summarise, meanhourDOC= mean(DOCcorr, na.rm = TRUE))
julDOC$hour <-  julDOC$date
year(julDOC$hour) <- 0
month(julDOC$hour) <- 0
day(julDOC$hour) <- 0
#julDOC$hour <- as.POSIXct(strptime(julDOC$hour, format="%H:%M"))
jul.abs254<- ggplot(subset(julDOC, julDOC$DOCcorr >=1), aes(hour, abs254)) + 
  geom_point(size = 0.4) +
  xlab("Date") + ylab("abs254") + 
  #scale_color_manual(values=cbPalette[1:2]) + # colours
  theme(legend.position="top") + theme() + ggtitle("Diurnal Cycles: July 25 - 31 abs254 ")

jul.SUVA<- ggplot(subset(julDOC, julDOC$DOCcorr >=1), aes(hour, SUVA)) + 
  geom_point(size = 0.4) +
  xlab("Date") + ylab("SUVA") + 
  #scale_color_manual(values=cbPalette[1:2]) + # colours
  theme(legend.position="top") + theme() + ggtitle("Diurnal Cycles: July 25 - 31 SUVA ") +
  scale_y_continuous(limits = c(3.5, 4))

jul.SR<- ggplot(subset(julDOC, julDOC$DOCcorr >=1), aes(hour, SR)) + 
  geom_point(size = 0.4) +
  xlab("Date") + ylab("Slope Ratio") + 
  #scale_color_manual(values=cbPalette[1:2]) + # colours
  theme(legend.position="top") + theme() + ggtitle("Diurnal Cycles: July 25 - 31 Slope Ratio ") +
  scale_y_continuous(limits = c(1.4, 2.2))

## Save plot - hourly data as a subplot of the whole week
#A viewport taking up a fraction of the plot area
vp <- viewport(width = 0.4, height = 0.4, x = 0.8, y = 0.8)
#Just draw the plot twice
png(file=paste0(fig.dir,"/CRFigures_diurnal_abs254.png"))
print(julyaug.abs254)
print(jul.abs254, vp = vp)
dev.off()
# SUVA - plot
vp <- viewport(width = 0.4, height = 0.4, x = 0.8, y = 0.2)
png(file=paste0(fig.dir,"/CRFigures_diurnal_SUVA.png"))
print(julyaug.SUVA)
print(jul.SUVA, vp = vp)
dev.off()

# SR - plot
vp <- viewport(width = 0.4, height = 0.4, x = 0.8, y = 0.8)
png(file=paste0(fig.dir,"/CRFigures_diurnal_SR.png"))
print(julyaug.SR)
print(jul.SR, vp = vp)
dev.off()

# Try plotting all - DOESN'T WORK _______________
#pdf(file=paste0(fig.dir,"/CRFigures_Quality_DiurnalAbs.pdf"), width = 8.5, height = 11) #save figure
#grid.arrange(print(jul.SR, vp = vp),print(jul.SUVA, vp = vp), print(jul.abs254, vp = vp), ncol = 1)
#dev.off()

########### Boxplots of PARAFAC (4 comp and redox and percent protein) by day/night - diurnal cycles
# do a column for d/n
daynight <- function(filename){
  filename$daynight <- ifelse(hour(filename$date) %in% 6:19, "day", "night") # 6 in the morning to 7 at night
  return(filename)
}

abs.Q <- daynight(abs.Q)
CR.grab <- daynight(CR.grab)

# Add in wet/dry + daynight
CR.grab$hydrodn <- paste(CR.grab$hydro, CR.grab$daynight, sep = "/")

# melt dataset to do boxplot - absorbance
box.crgrab <- melt(CR.grab[,c(2,3,7,28,29,49:52,44:47,89)], id = "hydrodn")
# boxplot - FI
diurnal.FI.box <- ggplot(subset(box.crgrab, variable  == "FI"), aes(x=hydrodn, y=value, fill=hydrodn)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('Diurnal - FI')) +
  ylab(expression('FI')) +
  xlab(expression('Day/Night')) 

# boxplot - HIX
diurnal.HIX.box <- ggplot(subset(box.crgrab, variable  == "HIX_ohno_area"), aes(x=hydrodn, y=value, fill=hydrodn)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('Diurnal - HIX')) +
  ylab(expression('HIX')) +
  xlab(expression('Day/Night')) + 
  coord_cartesian(ylim = c(0, 1.5))

diurnal.FrI.box <- ggplot(subset(box.crgrab, variable  == "FrI"), aes(x=hydrodn, y=value, fill=hydrodn)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('Diurnal - FrI')) +
  ylab(expression('FrI')) +
  xlab(expression('Day/Night')) + 
  coord_cartesian(ylim = c(0, 1.5))

diurnal.perprotein.box <- ggplot(subset(box.crgrab, variable  == "perprotein"), aes(x=hydrodn, y=value, fill=hydrodn)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('Diurnal - % Protein')) +
  ylab(expression('% Protein (%)')) +
  xlab(expression('Day/Night')) + 
  coord_cartesian(ylim = c(0, 15))

diurnal.redox.box <- ggplot(subset(box.crgrab, variable  == "redox"), aes(x=hydrodn, y=value, fill=hydrodn)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('Diurnal - Redox Index')) +
  ylab(expression('Redox Index')) +
  xlab(expression('Day/Night')) + 
  coord_cartesian(ylim = c(0.2, 0.5))
# PARAFAC Fmax
diurnal.C1.box <- ggplot(subset(box.crgrab, variable  == "PARAFAC_C1"), aes(x=hydrodn, y=value, fill=hydrodn)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('Diurnal - C1')) +
  ylab(expression('C1 Fmax')) +
  xlab(expression('Day/Night')) + 
  coord_cartesian(ylim = c(0, 1))
diurnal.C2.box <- ggplot(subset(box.crgrab, variable  == "PARAFAC_C2"), aes(x=hydrodn, y=value, fill=hydrodn)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('Diurnal - C2')) +
  ylab(expression('C2 Fmax')) +
  xlab(expression('Day/Night')) + 
  coord_cartesian(ylim = c(0, 1))
diurnal.C3.box <- ggplot(subset(box.crgrab, variable  == "PARAFAC_C3"), aes(x=hydrodn, y=value, fill=hydrodn)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('Diurnal - C3')) +
  ylab(expression('C3 Fmax')) +
  xlab(expression('Day/Night')) + 
  coord_cartesian(ylim = c(0, 1))
diurnal.C4.box <- ggplot(subset(box.crgrab, variable  == "PARAFAC_C4"), aes(x=hydrodn, y=value, fill=hydrodn)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)) +       # change colour to colour blind
  ggtitle(expression('Diurnal - C4')) +
  ylab(expression('C4 Fmax')) +
  xlab(expression('Day/Night')) + 
  coord_cartesian(ylim = c(0, 1))

# save boxplots - 6 in all
pdf(file=paste0(fig.dir,"/CRFigures_Quality_DiurnalFluor.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(diurnal.FI.box,diurnal.HIX.box, diurnal.FrI.box, diurnal.perprotein.box, diurnal.redox.box,
             diurnal.C1.box, diurnal.C2.box, diurnal.C3.box, diurnal.C4.box, ncol = 2)
dev.off()

#### Figure 3B - Seasonal Changes - boxplot of variables according to month - both in situ and PARAFAC components
# Need to show? Already show the seasonal boxplots for pre/post harvest? 

#############################################################################
# Figure 4 Correlation Matrix
## Correlation matrix - abs spectral/DOC conc/Q/bf/qf/pH/EC/Precip/soilT 
CR.spec.abs <- na.omit(abs.Q[,c(3:6,11:16,20:26,29:30,31:32,35:36)])
corr.matrix <- cor(as.data.frame(CR.spec.abs[,c(1:6, 8:17, 20:23)]), use="pairwise.complete.obs") # correlation matric between variables
# plot the correlation matrix of the spectral parameters
pdf(file=paste0(fig.dir,"/CR_absindicies_correlations.pdf"), width = 11, height = 8.5)
  corrplot(na.omit(corr.matrix), method = "circle") #plot matrix
dev.off()
write.csv(corr.matrix, paste0(fig.dir, "/CRabsind_corr.csv"))

####
# Figure Correlation Matrix - Water Quality variables and others (discharge, climate)
# clustering as per https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
# order="hclust", addrect=2) for clustering 
CR.grab.fl <- as.data.frame(na.omit(CR.grab[,c(2:12,15:29, 44:47, 63:65,68:70,74:78, 81, 58)]))
corr.matrix <- cor(CR.grab.fl[,c(1:2,6:10,12:31,34:43)]) # correlation matric between variables
# plot the correlation matrix of the spectral parameters
pdf(file=paste0(fig.dir,"/CR_flindicies_correlations.pdf"), width = 11, height = 8.5)
corrplot(na.omit(corr.matrix), method = "circle", order="hclust", addrect=2) #plot matrix
dev.off()
write.csv(corr.matrix, paste0(fig.dir, ("/CRcorrmatrix_abs.csv"))) #write correlation matrix

# show very unclear trends!!

# Correlation matrix with the high frequency data
CR.abs.cor <- na.omit(abs.Q[,c(3,11,12,13,14:17,21:27,30:31,34:35)])
corr.matrix <- cor(CR.abs.cor) # correlation matric between variables
# plot the correlation matrix of the spectral parameters
pdf(file=paste0(fig.dir,"/CR_abs_correlations.pdf"), width = 11, height = 8.5)
corrplot(na.omit(corr.matrix), method = "circle", order="hclust", addrect=2) #plot matrix
dev.off()
write.csv(corr.matrix, paste0(fig.dir, ("/CRcorrmatrix_abs.csv"))) #write correlation matrix

# do correlation matrix for the absorbance and fluorescence parameters together along with DOC and Q (baseflow, quickflow)
# make sure that the absorbance parameters 

# Linear relationships between DOC; e2e3, e4e6, SUVA, Slope Ratio - linear model - future work?

#########################################################################################################
# Part 2 - What is the effect of forest harvest on DOC characteristics?
#########################################################################################################
# Figure - ANOVA analysis to look at significance of changes in pre and post 
# Boxplot of flow weighted means for the pre/post period per month for different parameters
# do anova test on each of the pairs?
# ANOVA/Box plots on pre/post for significant changes
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
##################################
## Test for significant differences
# ANOVA on pre/post by month
# test for normality

# First, test whether data is normal or not.
## Have a look at the densities - SUVA
plot(density(abs.Q$SUVA, na.rm = TRUE))
## Perform the test - too many observations
#shapiro.test(abs.Q$SUVA)

## Plot using a qqplot - SUVA
qqnorm(abs.Q$SUVA);qqline(abs.Q$SUVA, col = 2)
# abs254 - not in the abs.Q data!! FUUCCCKKKK
qqnorm(abs.Q$SR);qqline(abs.Q$SR, col = 2)
# 

#

kruskal.test(jan.DOC.sub[,c(8,10,12)]) 
library(dunn.test)
dunn.test(jan.DOC.sub[,c(6,8,10,12)], g=DOCcorr, kw=TRUE)
kruskal.test(feb.DOC.sub[,c(6,8,10,12)]) 
kruskal.test(jan.DOC.sub[,c(4,6,8,10,12)]) 
kruskalmc(jan.DOC.sub$date, jan.DOC.sub[,c(4,6,8,10,12)]) 

#################################################################
# C12 and C1 (%) from cory mkcknight - by pre/post, wet dry.
# compare the two components as per 
# Strohmeier, S, K H Knorr, M Reichert, and S Frei. 2013. “Concentrations and Fluxes of Dissolved Organic Carbon in Runoff From a Forested Catchment: Insights From High Frequency Measurements.” …. doi:10.5194/bg-10-905-2013.
# http://stackoverflow.com/questions/2397097/how-can-a-data-ellipse-be-superimposed-on-a-ggplot2-scatterplot
library(ggplot2)
library(devtools)
library(digest)
source_url("https://raw.github.com/low-decarie/FAAV/master/r/stat-ellipse.R")    
# omit data where log status = NA
CR.grab.elipse <- CR.grab[!is.na(CR.grab$logstatus),]
C1C12.13comp <- qplot(data=CR.grab.elipse, x=C1.x, y=C12.x, colour=hydrolog, shape = hydrolog)+stat_ellipse() +
  xlim(0, 20) + ylim(0,15) + 
  scale_color_manual(values=cbPalette[1:4]) +
  labs(x="CM C1 (%)", y="CM C12 (%)", title = "DOC Origin - Pre and Post Logging") +
  coord_cartesian(xlim = c(5,15), ylim=c(0,15)) 

# save figure
pdf(file=paste0(fig.dir,"/CRFigures_13compC1vC13.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(C1C12.13comp, ncol = 1)
dev.off()
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

## Plot the raw data to see
ggplot(abs.Q, aes(x=Q.mm.d, y=DOCcorr, color = hydrolog)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Q.mm.day", y=expression(paste("[DOC] (mg/L)"))) +
  labs(title="Soil Extracts - [DOC] versus Depth")
# SUVA
SUVA.cQ <- ggplot(subset(abs.Q,Q.mm.d > 1), aes(x=Q.mm.d, y=SUVA, color = hydrolog)) +
  geom_point() +
  scale_color_manual(values=cbPalette[1:4], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:4], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(method="lm",aes(fill = hydrolog)) + 
  labs(x="Q (mm/day)", y=expression(paste("SUVA"[254]))) +
  labs(title="cQ plot of Spectral SUVA vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) +
  coord_cartesian(xlim = c(0,13), ylim=c(3,4)) 
# SR 
SR.cQ <- ggplot(subset(abs.Q,Q.mm.d > 1), aes(x=Q.mm.d, y=SlopeRatio, color = hydrolog)) +
  geom_point() +
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

## Fluorescence - Per Protein
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

### Conc-Q relationships for select variables: DOC, C1, C4, CM1, CM12 (%)
CR.grab.omit <- CR.grab[!is.na(CR.grab$logstatus),]
C1.cQ <- ggplot(subset(CR.grab.omit, Q.mm.d >=0.3), aes(x=Q.mm.d, y=C1_per, color = hydrolog)) +
  geom_point() +
  scale_color_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(data = subset(CR.grab.omit, Q.mm.d >=0.3), aes(x=Q.mm.d, y=C1_per, fill = hydrolog), method = 'lm', 
              formula = y ~ ns(x, 2), size = 1, se= FALSE)+
  labs(x="Q (mm/day)", y=expression(paste("C1%"))) +
  labs(title="cQ plot of Spectral C1% vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) 

C4.cQ <- ggplot(subset(CR.grab.omit, Q.mm.d >=0.3), aes(x=Q.mm.d, y=C4_per, color = hydrolog)) +
  geom_point() +
  scale_color_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(data = subset(CR.grab.omit, Q.mm.d >=0.3), aes(x=Q.mm.d, y=C4_per, fill = hydrolog), method = 'lm', 
              formula = y ~ ns(x, 2), size = 1, se= FALSE)+
  labs(x="Q (mm/day)", y=expression(paste("C4%"))) +
  labs(title="cQ plot of Spectral C4% vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) 

DOC.cQ <- ggplot(subset(abs.Q, Q.mm.d >=0.3), aes(x=Q.mm.d, y=DOCcorr, color = hydrolog)) +
  geom_point() +
  scale_color_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  # stat_smooth(method="lm",aes(fill = hydrolog)) + 
  labs(x="Q (mm/day)", y=expression(paste("[DOC] (mg/L)"))) +
  labs(title="cQ plot of Spectral [DOC] vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) + 
  coord_cartesian(xlim = c(0,40), ylim=c(0,20))

CMC1.cQ <- ggplot(subset(CR.grab.omit, Q.mm.d >=0.3), aes(x=Q.mm.d, y=C1.x, color = hydrolog)) +
  geom_point() +
  scale_color_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(data = subset(CR.grab.omit, Q.mm.d >=0.3), aes(x=Q.mm.d, y=C1.x, fill = hydrolog), method = 'lm', 
              formula = y ~ ns(x, 2), size = 1, se= FALSE)+
  labs(x="Q (mm/day)", y=expression(paste("CM C1%"))) +
  labs(title="cQ plot of Spectral CM C1% vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) + 
  coord_cartesian(xlim = c(0,30), ylim=c(0,20)) 

CMC12.cQ <- ggplot(subset(CR.grab.omit, Q.mm.d >=0.3), aes(x=Q.mm.d, y=C12.x, color = hydrolog)) +
  geom_point() +
  scale_color_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=cbPalette[1:5], guide=guide_legend(reverse=TRUE)) +
  stat_smooth(data = subset(CR.grab.omit, Q.mm.d >=0.3), aes(x=Q.mm.d, y=C12.x, fill = hydrolog), method = 'lm', 
              formula = y ~ ns(x, 2), size = 1, se= FALSE)+
  labs(x="Q (mm/day)", y=expression(paste("CM C12%"))) +
  labs(title="cQ plot of Spectral CM C12% vs. Q") +
  theme(legend.position="bottom", legend.title=element_blank()) #

# save 
pdf(file=paste0(fig.dir,"/CRFigures_CQrel_select.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(DOC.cQ, C1.cQ, C4.cQ, CMC1.cQ, CMC12.cQ, ncol = 2)
dev.off()
############
#DO linear models and save the coefficients (m,b, r2) in a table
# DOC- Q
abs.Q.pre <- subset(abs.Q, abs.Q$logstatus == "pre")
abs.Q.post <- subset(abs.Q, abs.Q$logstatus == "post")

spectrolm.DOC.pre <- c(summary(lm(as.numeric(abs.Q.pre$DOCcorr) ~ as.numeric(abs.Q.pre$Q.mm.d)))$coefficients[1:2], 
                   summary(lm(as.numeric(abs.Q.pre$DOCcorr) ~ as.numeric(abs.Q.pre$Q.mm.d)))$r.squared)
spectrolm.DOC.post <- c(summary(lm(as.numeric(abs.Q.post$DOCcorr) ~ as.numeric(abs.Q.post$Q.mm.d)))$coefficients[1:2], 
                       summary(lm(as.numeric(abs.Q.post$DOCcorr) ~ as.numeric(abs.Q.post$Q.mm.d)))$r.squared)
spectrolm.SUVA.pre <- c(summary(lm(as.numeric(abs.Q.pre$SUVA) ~ as.numeric(abs.Q.pre$Q.mm.d)))$coefficients[1:2], 
                       summary(lm(as.numeric(abs.Q.pre$SUVA) ~ as.numeric(abs.Q.pre$Q.mm.d)))$r.squared)
spectrolm.SUVA.post <- c(summary(lm(as.numeric(abs.Q.post$SUVA) ~ log(as.numeric(abs.Q.post$Q.mm.d))))$coefficients[1:2], 
                        summary(lm(as.numeric(abs.Q.post$SUVA) ~ log(as.numeric(abs.Q.post$Q.mm.d))))$r.squared)

## linear model
# Do linear models (glm) to look at relationships between DOC and various parameters, as well as hour, pre/post, wet/dry, month..

################## ################## ################## ################## ################## 
# Supplemental Figures
################## ################## ################## ################## ################## 

# PCA Analysis
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


############################################################################################
##### CR Soil Characteristics
# get 5 component PARAFAC fit - extracts and lysimeters
soil.5comp <- read.delim("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_SoilCorrEEMS/CR_Soil_PARAFAC/FMax_5.csv",sep = ",", header = FALSE)
# combine the 5 component with the sample key
soil.5comp$sample.ID <- t(read.delim("/Users/user/Documents/MATLAB/toolbox/CorrEEMS/CRSoil/CRSoil_01key.csv",sep = ",", header = FALSE))
colnames(soil.5comp)[1:5] <- c("C1", "C2", "C3", "C4", "C5")
#express as percent
soil.5comp.per <- cbind(soil.5comp[,1:5]/rowSums(soil.5comp[,1:5])*100, soil.5comp[,6])
colnames(soil.5comp.per)[6] <- "sample.ID"

############################################################################################
# Soil Extracts
## Compile all of the parameters you need - DOC (shimadzu) and aqualog (abs and fluor)
# linear model between soil depth and quality/DOC variables. Which variables best explain changes in soil depth?
extractwd <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_soilextracts"
### Get data
extract.DOC <- read.delim("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_soilextracts/CRSoil_DOC.csv", sep = ",")
extract.absfluor <- read.delim("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_soilextracts/CR_soilextracts_fluorescence/CorrectedEEM_Interp/SoilExtractsAbsFluor_Indicies.csv", sep = ",", header = TRUE)
extract.absfluor$sample.ID <- extract.absfluor$samplename
extract.CMFmax <- read.delim("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_soilextracts/CR_soilextracts_fluorescence/CR_PARAFACCMresults/SoilExtracts_componentsandloadings_CM_Fmax.csv", sep = ",", header = TRUE)
extract.CMper <- read.delim("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_soilextracts/CR_soilextracts_fluorescence/CR_PARAFACCMresults/SoilExtracts_componentsandloadings_CM.csv",sep = ",", header = TRUE)
# Get the 5-cmponent PARAFAC fit for the soil extracts only
extract.5comp.per <- soil.5comp.per[32:67,]
extract.5comp.fmax <- soil.5comp[32:67,]
extract.5comp.per$sample.ID <- sapply(strsplit(as.character(extract.5comp.per$sample.ID), split='_', fixed=TRUE), function(x) (x[1]))
extract.5comp.fmax$sample.ID <- sapply(strsplit(as.character(extract.5comp.fmax$sample.ID), split='_', fixed=TRUE), function(x) (x[1]))

# merge all together by the sample ID
extract.fluor <- merge(extract.CMFmax, extract.CMper, by = "sample.ID") #merged Cory McKnight Results
# create a new sample ID file as extract01
extract.fluor$sample.ID <- sub("(.*?)Soil.*", "\\1", as.character(extract.fluor$sample.ID))
#extract.fluor <- merge(extract.fluor, extract.5comp.per, by = "sample.ID") # combine the 4 component, and the CM results
#extract.fluor <- merge(extract.fluor, extract.absfluor, by = "sample.ID")
# merge with the DOC concentrations
#extract.all <- merge(extract.DOC, extract.fluor, by = "sample.ID")
# merge all together
extract.all <- Reduce(function(x,y) {merge(x,y, by = "sample.ID")}, list(extract.fluor, extract.5comp.per, extract.5comp.fmax, extract.absfluor, extract.DOC))

# associate the sample ID with the different levels, and repetitions within soil pits
extract.code <- read.delim("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_soilextracts/CRSoilExtracts_DepthCode.csv", sep = ",", header = TRUE)
# pad the numbers to two digits
extract.code <- apply(extract.code, 2,  function(x) str_pad(x, 2, pad = "0"))
# sample code is last two digits of the 
extract.all$samplecode  <- stri_sub(extract.all$sample.ID,-2,-1)
extract.all$extract <- c("dH2O", "dH2O", "dH2O", "dH2O", "dH2O", "dH2O", "dH2O", "dH2O", "dH2O", "dH2O", "dH2O", "dH2O",
                              "K2SO4", "K2SO4", "K2SO4", "K2SO4", "K2SO4", "K2SO4", "K2SO4", "K2SO4", "K2SO4", "K2SO4", "K2SO4", "K2SO4", 
                              "KCl", "KCl", "KCl", "KCl", "KCl", "KCl", "KCl", "KCl", "KCl", "KCl", "KCl", "KCl") 

# Associating sample depth to aggregate by depth 
extract.all$depth[extract.all$samplecode == c("01")] <- "01"
extract.all$depth[extract.all$samplecode == c("05")] <- "01"
extract.all$depth[extract.all$samplecode == c("09")] <- "01"
extract.all$depth[extract.all$samplecode == c("02")] <- "02"
extract.all$depth[extract.all$samplecode == c("06")] <- "02"
extract.all$depth[extract.all$samplecode == c("10")] <- "02"
extract.all$depth[extract.all$samplecode == c("03")] <- "03"
extract.all$depth[extract.all$samplecode == c("07")] <- "03"
extract.all$depth[extract.all$samplecode == c("11")] <- "03"
extract.all$depth[extract.all$samplecode == c("04")] <- "04"
extract.all$depth[extract.all$samplecode == c("08")] <- "04"
extract.all$depth[extract.all$samplecode == c("12")] <- "04"

## associate with the real depth in cm
extract.all$depth_cm[extract.all$depth == c("01")] <- 05
extract.all$depth_cm[extract.all$depth == c("02")] <- 26.5
extract.all$depth_cm[extract.all$depth == c("03")] <- 42.5
extract.all$depth_cm[extract.all$depth == c("04")] <- 69.5

# calculate SUVA
extract.all$SUVA <- extract.all$abs254/extract.all$NPOC_mgL

# do linear fits for all of the parameters, and get table with R2 and with all of the slopes. For all three types of soil extractions
extract.h2o <- subset(extract.all, extract.all$extract == "dH2O")
extract.kcl <- subset(extract.all, extract.all$extract == "KCl")
extract.k2so4 <- subset(extract.all, extract.all$extract == "K2SO4")

# Do linear fits - must copy and paste the right extract so I don't have to copy and paste all of the code below..
extractlm.DOC <- lm(as.numeric(extract.h2o$NPOC_mgL) ~ as.numeric(extract.h2o$depth_cm))
lm.DOC <- c(summary(extractlm.DOC)$coefficients[1:2], summary(extractlm.DOC)$r.squared)
extractlm.TN <- lm(as.numeric(extract.h2o$TN_mgL) ~ as.numeric(extract.h2o$depth_cm))
lm.TN <- c(summary(extractlm.TN)$coefficients[1:2], summary(extractlm.TN)$r.squared)
# per protein
extractlm.pp <- lm(as.numeric(extract.h2o$perprotein) ~ as.numeric(extract.h2o$depth_cm))
lm.pp <- c(summary(extractlm.pp)$coefficients[1:2], summary(extractlm.pp)$r.squared)
# redox
extractlm.redox <- lm(as.numeric(extract.h2o$redox) ~ as.numeric(extract.h2o$depth_cm))
lm.redox <- c(summary(extractlm.redox)$coefficients[1:2], summary(extractlm.redox)$r.squared)
# C1
extractlm.C1 <- lm(as.numeric(extract.h2o$C1) ~ as.numeric(extract.h2o$depth_cm))
lm.C1 <- c(summary(extractlm.C1)$coefficients[1:2], summary(extractlm.C1)$r.squared)
# C2
extractlm.C2 <- lm(as.numeric(extract.h2o$C2) ~ as.numeric(extract.h2o$depth_cm))
lm.C2 <- c(summary(extractlm.C2)$coefficients[1:2], summary(extractlm.C2)$r.squared)
# C3
extractlm.C3 <- lm(as.numeric(extract.h2o$C3) ~ as.numeric(extract.h2o$depth_cm))
lm.C3 <- c(summary(extractlm.C3)$coefficients[1:2], summary(extractlm.C3)$r.squared)
# C4
extractlm.C4 <- lm(as.numeric(extract.h2o$C4) ~ as.numeric(extract.h2o$depth_cm))
lm.C4 <- c(summary(extractlm.C4)$coefficients[1:2], summary(extractlm.C4)$r.squared)
# C5
extractlm.C5 <- lm(as.numeric(extract.h2o$C5) ~ as.numeric(extract.h2o$depth_cm))
lm.C5 <- c(summary(extractlm.C5)$coefficients[1:2], summary(extractlm.C5)$r.squared)
# abs254
extractlm.abs254 <- lm(as.numeric(extract.h2o$abs254) ~ as.numeric(extract.h2o$depth_cm))
lm.abs254 <- c(summary(extractlm.abs254)$coefficients[1:2], summary(extractlm.abs254)$r.squared)
# SUVA
extractlm.SUVA <- lm(as.numeric(extract.h2o$SUVA) ~ as.numeric(extract.h2o$depth_cm))
lm.SUVA <- c(summary(extractlm.SUVA)$coefficients[1:2], summary(extractlm.abs254)$r.squared)
# abs272
extractlm.abs272 <- lm(as.numeric(extract.h2o$abs272) ~ as.numeric(extract.h2o$depth_cm))
lm.abs272 <- c(summary(extractlm.abs272)$coefficients[1:2], summary(extractlm.abs272)$r.squared)
# e2e3
extractlm.e2e3 <- lm(as.numeric(extract.h2o$e2e3) ~ as.numeric(extract.h2o$depth_cm))
lm.e2e3 <- c(summary(extractlm.e2e3)$coefficients[1:2], summary(extractlm.e2e3)$r.squared)
# e4e6
extractlm.e4e6 <- lm(as.numeric(extract.h2o$e4e6) ~ as.numeric(extract.h2o$depth_cm))
lm.e4e6 <- c(summary(extractlm.e4e6)$coefficients[1:2], summary(extractlm.e4e6)$r.squared)
# SR
extractlm.SR <- lm(as.numeric(extract.h2o$SR) ~ as.numeric(extract.h2o$depth_cm))
lm.SR <- c(summary(extractlm.SR)$coefficients[1:2], summary(extractlm.SR)$r.squared)
# FI
extractlm.FI <- lm(as.numeric(extract.h2o$FI) ~ as.numeric(extract.h2o$depth_cm))
lm.FI <- c(summary(extractlm.FI)$coefficients[1:2], summary(extractlm.FI)$r.squared)
# HIX_ohno_area
extractlm.HIX_ohno_area <- lm(as.numeric(extract.h2o$HIX_ohno_area) ~ as.numeric(extract.h2o$depth_cm))
lm.HIX_ohno_area <- c(summary(extractlm.HIX_ohno_area)$coefficients[1:2], summary(extractlm.HIX_ohno_area)$r.squared)
# FrI
extractlm.FrI <- lm(as.numeric(extract.h2o$FrI) ~ as.numeric(extract.h2o$depth_cm))
lm.FrI <- c(summary(extractlm.FrI)$coefficients[1:2], summary(extractlm.FrI)$r.squared)
# peakA
extractlm.peakA <- lm(as.numeric(extract.h2o$peakA) ~ as.numeric(extract.h2o$depth_cm))
lm.peakA <- c(summary(extractlm.peakA)$coefficients[1:2], summary(extractlm.peakA)$r.squared)
# peakC
extractlm.peakC <- lm(as.numeric(extract.h2o$peakC) ~ as.numeric(extract.h2o$depth_cm))
lm.peakC <- c(summary(extractlm.peakC)$coefficients[1:2], summary(extractlm.peakC)$r.squared)
# peakB
extractlm.peakB <- lm(as.numeric(extract.h2o$peakB) ~ as.numeric(extract.h2o$depth_cm))
lm.peakB <- c(summary(extractlm.peakB)$coefficients[1:2], summary(extractlm.peakB)$r.squared)
# peakT
extractlm.peakT <- lm(as.numeric(extract.h2o$peakT) ~ as.numeric(extract.h2o$depth_cm))
lm.peakT <- c(summary(extractlm.peakT)$coefficients[1:2], summary(extractlm.peakT)$r.squared)
# OFI
extractlm.OFI <- lm(as.numeric(extract.h2o$OFI) ~ as.numeric(extract.h2o$depth_cm))
lm.OFI <- c(summary(extractlm.OFI)$coefficients[1:2], summary(extractlm.OFI)$r.squared)
# CoryMcKnightC1
extractlm.CM.C1 <- lm(as.numeric(extract.h2o$C1.y) ~ as.numeric(extract.h2o$depth_cm))
lm.CM.C1 <- c(summary(extractlm.CM.C1)$coefficients[1:2], summary(extractlm.CM.C1)$r.squared)
# CoryMcKnightC2
extractlm.CM.C2 <- lm(as.numeric(extract.h2o$C2.y) ~ as.numeric(extract.h2o$depth_cm))
lm.CM.C2 <- c(summary(extractlm.CM.C2)$coefficients[1:2], summary(extractlm.CM.C2)$r.squared)
# CoryMcKnightC3
extractlm.CM.C3 <- lm(as.numeric(extract.h2o$C3.y) ~ as.numeric(extract.h2o$depth_cm))
lm.CM.C3 <- c(summary(extractlm.CM.C3)$coefficients[1:2], summary(extractlm.CM.C3)$r.squared)
# CoryMcKnightC5
extractlm.CM.C4 <- lm(as.numeric(extract.h2o$C4.y) ~ as.numeric(extract.h2o$depth_cm))
lm.CM.C4 <- c(summary(extractlm.CM.C4)$coefficients[1:2], summary(extractlm.CM.C4)$r.squared)
# CoryMcKnightC5
extractlm.CM.C5 <- lm(as.numeric(extract.h2o$C5.y) ~ as.numeric(extract.h2o$depth_cm))
lm.CM.C5 <- c(summary(extractlm.CM.C5)$coefficients[1:2], summary(extractlm.CM.C5)$r.squared)
# CoryMcKnightC6
extractlm.CM.C6 <- lm(as.numeric(extract.h2o$C6.y) ~ as.numeric(extract.h2o$depth_cm))
lm.CM.C6 <- c(summary(extractlm.CM.C6)$coefficients[1:2], summary(extractlm.CM.C6)$r.squared)
# CoryMcKnightC7
extractlm.CM.C7 <- lm(as.numeric(extract.h2o$C7.y) ~ as.numeric(extract.h2o$depth_cm))
lm.CM.C7 <- c(summary(extractlm.CM.C7)$coefficients[1:2], summary(extractlm.CM.C7)$r.squared)
# CoryMcKnightC8
extractlm.CM.C8 <- lm(as.numeric(extract.h2o$C8.y) ~ as.numeric(extract.h2o$depth_cm))
lm.CM.C8 <- c(summary(extractlm.CM.C8)$coefficients[1:2], summary(extractlm.CM.C8)$r.squared)
# CoryMcKnightC9
extractlm.CM.C9 <- lm(as.numeric(extract.h2o$C9.y) ~ as.numeric(extract.h2o$depth_cm))
lm.CM.C9 <- c(summary(extractlm.CM.C9)$coefficients[1:2], summary(extractlm.CM.C9)$r.squared)
# CoryMcKnightC10
extractlm.CM.C10 <- lm(as.numeric(extract.h2o$C10.y) ~ as.numeric(extract.h2o$depth_cm))
lm.CM.C10 <- c(summary(extractlm.CM.C10)$coefficients[1:2], summary(extractlm.CM.C10)$r.squared)
# CoryMcKnightC11
extractlm.CM.C11 <- lm(as.numeric(extract.h2o$C11.y) ~ as.numeric(extract.h2o$depth_cm))
lm.CM.C11 <- c(summary(extractlm.CM.C11)$coefficients[1:2], summary(extractlm.CM.C11)$r.squared)
# CoryMcKnightC12
extractlm.CM.C12 <- lm(as.numeric(extract.h2o$C12.y) ~ as.numeric(extract.h2o$depth_cm))
lm.CM.C12 <- c(summary(extractlm.CM.C12)$coefficients[1:2], summary(extractlm.CM.C12)$r.squared)
# CoryMcKnightC13
extractlm.CM.C13 <- lm(as.numeric(extract.h2o$C13.y) ~ as.numeric(extract.h2o$depth_cm))
lm.CM.C13 <- c(summary(extractlm.CM.C13)$coefficients[1:2], summary(extractlm.CM.C13)$r.squared)
# Percent protein
extractlm.perprotein <- lm(as.numeric(extract.h2o$perprotein) ~ as.numeric(extract.h2o$depth_cm))
lm.perprotein <- c(summary(extractlm.perprotein)$coefficients[1:2], summary(extractlm.perprotein)$r.squared)
# Redox
extractlm.redox <- lm(as.numeric(extract.h2o$redox) ~ as.numeric(extract.h2o$depth_cm))
lm.redox <- c(summary(extractlm.redox)$coefficients[1:2], summary(extractlm.redox)$r.squared)

# bind together to write to csv
extract.linear.h2o <- rbind(lm.DOC, lm.TN, lm.abs254, lm.SUVA, lm.abs272, lm.C1, lm.C2, lm.C3, lm.C4, lm.C5,
                      lm.e2e3, lm.e4e6, lm.SR, lm.FI, lm.HIX_ohno_area, lm.FrI, lm.peakA, 
                      lm.perprotein, lm.redox,
                      lm.peakC, lm.peakB, lm.peakT, lm.OFI, 
                      lm.CM.C1, lm.CM.C2, lm.CM.C3,
                      lm.CM.C4, lm.CM.C5, lm.CM.C6, lm.CM.C7, lm.CM.C8, 
                      lm.CM.C9, lm.CM.C10, lm.CM.C11, lm.CM.C12, lm.CM.C13)
colnames(extract.linear.h2o) <- c("intercept.h2o", "slope.h2o", "R2.h2o")
# write the table with slopes/r2 values for all of the extracts
write.csv(extract.linear.h2o, file = paste0(soilextract.path, "/lmfitsh2o.csv"))

# Do scatter plot with the parameters versus soil depth 0 aggregating by soil depth - for depth
# for DOC
mean.dh2o <- ddply(extract.h2o,~depth_cm,summarise,
      mean.DOC=mean(NPOC_mgL),sd.DOC=sd(NPOC_mgL),
      mean.TN=mean(TN_mgL),sd.TN=sd(TN_mgL),
      mean.perprotein=mean(perprotein),sd.perprotein=sd(perprotein),
      mean.redox=mean(redox),sd.redox=sd(redox),
      mean.C1 = mean(C1),sd.C1 = sd(C1),
      mean.C2 = mean(C2),sd.C2 = sd(C2),
      mean.C3 = mean(C3),sd.C3 = sd(C3),
      mean.C4 = mean(C4),sd.C4 = sd(C4),
      mean.C5 = mean(C5),sd.C5 = sd(C5),
      mean.abs254=mean(abs254),sd.abs254=sd(abs254),
      mean.SUVA=mean(SUVA),sd.SUVA=sd(SUVA),
      mean.abs272=mean(abs272),sd.abs272=sd(abs272),
      mean.e2e3=mean(e2e3),sd.e2e3=sd(e2e3), mean.e4e6=mean(e4e6),sd.e4e6=sd(e4e6),
      mean.SR=mean(SR),sd.SR=sd(SR), mean.FI=mean(FI),sd.FI=sd(FI),
      mean.HIX=mean(HIX_ohno_area),sd.HIX=sd(HIX_ohno_area), mean.FrI=mean(FrI),sd.FrI=sd(FrI),
      mean.peakA=mean(peakA),sd.peakA=sd(peakA), mean.peakC=mean(peakC),sd.peakC=sd(peakC),
      mean.peakB=mean(peakB),sd.peakB=sd(peakB), mean.peakT=mean(peakT),sd.peakT=sd(peakT),
      mean.peakOFI=mean(OFI),sd.peakOFI=sd(OFI),
      mean.CM.C1 = mean(C1.y),sd.C1 = sd(C1.y),
      mean.CM.C2 = mean(C2.y),sd.C2 = sd(C2.y),
      mean.CM.C3 = mean(C3.y),sd.C3 = sd(C3.y),
      mean.CM.C4 = mean(C4.y),sd.C4 = sd(C4.y),
      mean.CM.C5 = mean(C5.y),sd.C5 = sd(C5.y),
      mean.CM.C6 = mean(C6.y),sd.C6 = sd(C6.y),
      mean.CM.C7 = mean(C7.y),sd.C7 = sd(C7.y),
      mean.CM.C8 = mean(C8.y),sd.C8 = sd(C8.y),
      mean.CM.C9 = mean(C9.y),sd.C9 = sd(C9.y),
      mean.CM.C10 = mean(C10.y),sd.C10 = sd(C10.y),
      mean.CM.C11 = mean(C11.y),sd.C11 = sd(C11.y),
      mean.CM.C12 = mean(C12.y),sd.C12 = sd(C12.y),
      mean.CM.C13 = mean(C13.y),sd.C13 = sd(C13.y))
mean.dh2o$extract <- "dh2o"
mean.dh2o <- melt(mean.dh2o, id = c("extract", "depth_cm"))

mean.kcl <- ddply(extract.kcl,~depth_cm,summarise,
                   mean.DOC=mean(NPOC_mgL),sd.DOC=sd(NPOC_mgL),
                   mean.TN=mean(TN_mgL),sd.TN=sd(TN_mgL),
                   mean.perprotein=mean(perprotein),sd.perprotein=sd(perprotein),
                   mean.redox=mean(redox),sd.redox=sd(redox),
                   mean.C1 = mean(C1),sd.C1 = sd(C1),
                   mean.C2 = mean(C2),sd.C2 = sd(C2),
                   mean.C3 = mean(C3),sd.C3 = sd(C3),
                   mean.C4 = mean(C4),sd.C4 = sd(C4),
                   mean.C5 = mean(C5),sd.C5 = sd(C5),
                   mean.abs254=mean(abs254),sd.abs254=sd(abs254),
                   mean.SUVA=mean(SUVA),sd.SUVA=sd(SUVA),
                   mean.abs272=mean(abs272),sd.abs272=sd(abs272),
                   mean.e2e3=mean(e2e3),sd.e2e3=sd(e2e3), mean.e4e6=mean(e4e6),sd.e4e6=sd(e4e6),
                   mean.SR=mean(SR),sd.SR=sd(SR), mean.FI=mean(FI),sd.FI=sd(FI),
                   mean.HIX=mean(HIX_ohno_area),sd.HIX=sd(HIX_ohno_area), mean.FrI=mean(FrI),sd.FrI=sd(FrI),
                   mean.peakA=mean(peakA),sd.peakA=sd(peakA), mean.peakC=mean(peakC),sd.peakC=sd(peakC),
                   mean.peakB=mean(peakB),sd.peakB=sd(peakB), mean.peakT=mean(peakT),sd.peakT=sd(peakT),
                   mean.peakOFI=mean(OFI),sd.peakOFI=sd(OFI),
                  mean.CM.C1 = mean(C1.y),sd.C1 = sd(C1.y),
                  mean.CM.C2 = mean(C2.y),sd.C2 = sd(C2.y),
                  mean.CM.C3 = mean(C3.y),sd.C3 = sd(C3.y),
                  mean.CM.C4 = mean(C4.y),sd.C4 = sd(C4.y),
                  mean.CM.C5 = mean(C5.y),sd.C5 = sd(C5.y),
                  mean.CM.C6 = mean(C6.y),sd.C6 = sd(C6.y),
                  mean.CM.C7 = mean(C7.y),sd.C7 = sd(C7.y),
                  mean.CM.C8 = mean(C8.y),sd.C8 = sd(C8.y),
                  mean.CM.C9 = mean(C9.y),sd.C9 = sd(C9.y),
                  mean.CM.C10 = mean(C10.y),sd.C10 = sd(C10.y),
                  mean.CM.C11 = mean(C11.y),sd.C11 = sd(C11.y),
                  mean.CM.C12 = mean(C12.y),sd.C12 = sd(C12.y),
                  mean.CM.C13 = mean(C13.y),sd.C13 = sd(C13.y))
mean.kcl$extract <- "kcl"
mean.kcl <- melt(mean.kcl, id = c("extract", "depth_cm"))

mean.k2s04 <- ddply(extract.k2so4,~depth_cm,summarise,
                   mean.DOC=mean(NPOC_mgL),sd.DOC=sd(NPOC_mgL),
                   mean.TN=mean(TN_mgL),sd.TN=sd(TN_mgL),
                   mean.perprotein=mean(perprotein),sd.perprotein=sd(perprotein),
                   mean.redox=mean(redox),sd.redox=sd(redox),
                   mean.C1 = mean(C1),sd.C1 = sd(C1),
                   mean.C2 = mean(C2),sd.C2 = sd(C2),
                   mean.C3 = mean(C3),sd.C3 = sd(C3),
                   mean.C4 = mean(C4),sd.C4 = sd(C4),
                   mean.C5 = mean(C5),sd.C5 = sd(C5),
                   mean.abs254=mean(abs254),sd.abs254=sd(abs254),
                   mean.SUVA=mean(SUVA),sd.SUVA=sd(SUVA),
                   mean.abs272=mean(abs272),sd.abs272=sd(abs272),
                   mean.e2e3=mean(e2e3),sd.e2e3=sd(e2e3), mean.e4e6=mean(e4e6),sd.e4e6=sd(e4e6),
                   mean.SR=mean(SR),sd.SR=sd(SR), mean.FI=mean(FI),sd.FI=sd(FI),
                   mean.HIX=mean(HIX_ohno_area),sd.HIX=sd(HIX_ohno_area), mean.FrI=mean(FrI),sd.FrI=sd(FrI),
                   mean.peakA=mean(peakA),sd.peakA=sd(peakA), mean.peakC=mean(peakC),sd.peakC=sd(peakC),
                   mean.peakB=mean(peakB),sd.peakB=sd(peakB), mean.peakT=mean(peakT),sd.peakT=sd(peakT),
                   mean.peakOFI=mean(OFI),sd.peakOFI=sd(OFI),
                   mean.CM.C1 = mean(C1.y),sd.C1 = sd(C1.y),
                   mean.CM.C2 = mean(C2.y),sd.C2 = sd(C2.y),
                   mean.CM.C3 = mean(C3.y),sd.C3 = sd(C3.y),
                   mean.CM.C4 = mean(C4.y),sd.C4 = sd(C4.y),
                   mean.CM.C5 = mean(C5.y),sd.C5 = sd(C5.y),
                   mean.CM.C6 = mean(C6.y),sd.C6 = sd(C6.y),
                   mean.CM.C7 = mean(C7.y),sd.C7 = sd(C7.y),
                   mean.CM.C8 = mean(C8.y),sd.C8 = sd(C8.y),
                   mean.CM.C9 = mean(C9.y),sd.C9 = sd(C9.y),
                   mean.CM.C10 = mean(C10.y),sd.C10 = sd(C10.y),
                   mean.CM.C11 = mean(C11.y),sd.C11 = sd(C11.y),
                   mean.CM.C12 = mean(C12.y),sd.C12 = sd(C12.y),
                   mean.CM.C13 = mean(C13.y),sd.C13 = sd(C13.y))
#colnames(mean.k2s04) <- paste(colnames(mean.k2s04), "k2s04", sep = "")
#colnames(mean.k2s04)[1] <- 'depth'
mean.k2s04$extract <- "k2so4"
mean.k2s04 <- melt(mean.k2s04, id = c("extract", "depth_cm"))
mean.all <-rbind(mean.dh2o, mean.kcl, mean.k2s04)

# DOC by depth for all extract types
extract.DOCplot <- ggplot(subset(mean.all, mean.all$variable == "mean.DOC"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("[DOC] (mg/L)"))) +
  labs(title="Soil Extracts - [DOC] versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.DOC")$value - subset(mean.all, mean.all$variable == "sd.DOC")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.DOC")$value + subset(mean.all, mean.all$variable == "sd.DOC")$value,
                    position = "dodge", width=0.2))
extract.TNplot <- ggplot(subset(mean.all, mean.all$variable == "mean.TN"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("[TN] (mg/L)"))) +
  labs(title="Soil Extracts - [TN] versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.TN")$value - subset(mean.all, mean.all$variable == "sd.TN")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.TN")$value + subset(mean.all, mean.all$variable == "sd.TN")$value,
                    position = "dodge", width=0.2))
extract.perproteinplot <- ggplot(subset(mean.all, mean.all$variable == "mean.perprotein"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("perprotein (%)"))) +
  labs(title="Soil Extracts - perprotein versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.perprotein")$value - subset(mean.all, mean.all$variable == "sd.perprotein")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.perprotein")$value + subset(mean.all, mean.all$variable == "sd.perprotein")$value,
                    position = "dodge", width=0.2))
extract.redox <- ggplot(subset(mean.all, mean.all$variable == "mean.redox"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("Redox (A.U)"))) +
  labs(title="Soil Extracts - Redox versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.redox")$value - subset(mean.all, mean.all$variable == "sd.redox")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.redox")$value + subset(mean.all, mean.all$variable == "sd.redox")$value,
                    position = "dodge", width=0.2))
extract.C1 <- ggplot(subset(mean.all, mean.all$variable == "mean.C1"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("C1 (A.U)"))) +
  labs(title="Soil Extracts - C1 versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.C1")$value - subset(mean.all, mean.all$variable == "sd.C1")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.C1")$value + subset(mean.all, mean.all$variable == "sd.C1")$value,
                    position = "dodge", width=0.2))
extract.C2 <- ggplot(subset(mean.all, mean.all$variable == "mean.C2"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("C2 (A.U)"))) +
  labs(title="Soil Extracts - C2 versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.C2")$value - subset(mean.all, mean.all$variable == "sd.C2")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.C2")$value + subset(mean.all, mean.all$variable == "sd.C2")$value,
                    position = "dodge", width=0.2))
extract.C3 <- ggplot(subset(mean.all, mean.all$variable == "mean.C3"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("C3 (A.U)"))) +
  labs(title="Soil Extracts - C3 versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.C3")$value - subset(mean.all, mean.all$variable == "sd.C3")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.C3")$value + subset(mean.all, mean.all$variable == "sd.C3")$value,
                    position = "dodge", width=0.2))
extract.C4 <- ggplot(subset(mean.all, mean.all$variable == "mean.C4"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("C4 (A.U)"))) +
  labs(title="Soil Extracts - C4 versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.C4")$value - subset(mean.all, mean.all$variable == "sd.C4")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.C4")$value + subset(mean.all, mean.all$variable == "sd.C4")$value,
                    position = "dodge", width=0.2))
extract.C5 <- ggplot(subset(mean.all, mean.all$variable == "mean.C5"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("C5 (A.U)"))) +
  labs(title="Soil Extracts - C5 versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.C5")$value - subset(mean.all, mean.all$variable == "sd.C5")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.C5")$value + subset(mean.all, mean.all$variable == "sd.C5")$value,
                    position = "dodge", width=0.2))
extract.abs254 <- ggplot(subset(mean.all, mean.all$variable == "mean.abs254"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("abs254 (A.U)"))) +
  labs(title="Soil Extracts - abs254 versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.abs254")$value - subset(mean.all, mean.all$variable == "sd.abs254")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.abs254")$value + subset(mean.all, mean.all$variable == "sd.abs254")$value,
                    position = "dodge", width=0.2))
extract.SUVA<- ggplot(subset(mean.all, mean.all$variable == "mean.SUVA"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("SUVA (A.U)"))) +
  labs(title="Soil Extracts - SUVA versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.SUVA")$value - subset(mean.all, mean.all$variable == "sd.SUVA")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.SUVA")$value + subset(mean.all, mean.all$variable == "sd.SUVA")$value,
                    position = "dodge", width=0.2))

extract.abs272 <- ggplot(subset(mean.all, mean.all$variable == "mean.abs272"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("abs272 (A.U)"))) +
  labs(title="Soil Extracts - abs272 versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.abs272")$value - subset(mean.all, mean.all$variable == "sd.abs272")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.abs272")$value + subset(mean.all, mean.all$variable == "sd.abs272")$value,
                    position = "dodge", width=0.2))
extract.e2e3 <- ggplot(subset(mean.all, mean.all$variable == "mean.e2e3"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("e2e3 (A.U)"))) +
  labs(title="Soil Extracts - e2e3 versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.e2e3")$value - subset(mean.all, mean.all$variable == "sd.e2e3")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.e2e3")$value + subset(mean.all, mean.all$variable == "sd.e2e3")$value,
                    position = "dodge", width=0.2))
extract.e4e6 <- ggplot(subset(mean.all, mean.all$variable == "mean.e4e6"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("e4e6 (A.U)"))) +
  labs(title="Soil Extracts - e4e6 versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.e4e6")$value - subset(mean.all, mean.all$variable == "sd.e4e6")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.e4e6")$value + subset(mean.all, mean.all$variable == "sd.e4e6")$value,
                    position = "dodge", width=0.2))
extract.SR <- ggplot(subset(mean.all, mean.all$variable == "mean.SR"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("SR (A.U)"))) +
  labs(title="Soil Extracts - SR versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.SR")$value - subset(mean.all, mean.all$variable == "sd.SR")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.SR")$value + subset(mean.all, mean.all$variable == "sd.SR")$value,
                    position = "dodge", width=0.2))
extract.FI <- ggplot(subset(mean.all, mean.all$variable == "mean.FI"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("FI (A.U)"))) +
  labs(title="Soil Extracts - FI versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.FI")$value - subset(mean.all, mean.all$variable == "sd.FI")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.FI")$value + subset(mean.all, mean.all$variable == "sd.FI")$value,
                    position = "dodge", width=0.2))
extract.HIX_ohno_area <- ggplot(subset(mean.all, mean.all$variable == "mean.HIX"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("HIX (A.U)"))) +
  labs(title="Soil Extracts - HIX versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.HIX")$value - subset(mean.all, mean.all$variable == "sd.HIX")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.HIX")$value + subset(mean.all, mean.all$variable == "sd.HIX")$value,
                    position = "dodge", width=0.2))

extract.FrI <- ggplot(subset(mean.all, mean.all$variable == "mean.FrI"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("FrI (A.U)"))) +
  labs(title="Soil Extracts - FrI versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.FrI")$value - subset(mean.all, mean.all$variable == "sd.FrI")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.FrI")$value + subset(mean.all, mean.all$variable == "sd.FrI")$value,
                    position = "dodge", width=0.2))
extract.peakA <- ggplot(subset(mean.all, mean.all$variable == "mean.peakA"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("peakA (A.U)"))) +
  labs(title="Soil Extracts - peakA versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.peakA")$value - subset(mean.all, mean.all$variable == "sd.peakA")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.peakA")$value + subset(mean.all, mean.all$variable == "sd.peakA")$value,
                    position = "dodge", width=0.2))
extract.peakC <- ggplot(subset(mean.all, mean.all$variable == "mean.peakC"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("peakC (A.U)"))) +
  labs(title="Soil Extracts - peakC versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.peakC")$value - subset(mean.all, mean.all$variable == "sd.peakC")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.peakC")$value + subset(mean.all, mean.all$variable == "sd.peakC")$value,
                    position = "dodge", width=0.2))
extract.peakB <- ggplot(subset(mean.all, mean.all$variable == "mean.peakB"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("peakB (A.U)"))) +
  labs(title="Soil Extracts - peakB versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.peakB")$value - subset(mean.all, mean.all$variable == "sd.peakB")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.peakB")$value + subset(mean.all, mean.all$variable == "sd.peakB")$value,
                    position = "dodge", width=0.2))
extract.peakT <- ggplot(subset(mean.all, mean.all$variable == "mean.peakT"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("peakT (A.U)"))) +
  labs(title="Soil Extracts - peakT versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.peakT")$value - subset(mean.all, mean.all$variable == "sd.peakT")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.peakT")$value + subset(mean.all, mean.all$variable == "sd.peakT")$value,
                    position = "dodge", width=0.2))
extract.OFI <- ggplot(subset(mean.all, mean.all$variable == "mean.peakOFI"), aes(x = as.numeric(depth_cm), y = as.numeric(value), colour = extract, shape=extract)) + 
  geom_point() +    
  geom_smooth(method="lm", se = FALSE) +
  labs(x="Depth (cm)", y=expression(paste("OFI (A.U)"))) +
  labs(title="Soil Extracts - OFI versus Depth") + 
  geom_errorbar(aes(ymin = subset(mean.all, mean.all$variable == "mean.peakOFI")$value - subset(mean.all, mean.all$variable == "sd.peakOFI")$value, 
                    ymax = subset(mean.all, mean.all$variable == "mean.peakOFI")$value + subset(mean.all, mean.all$variable == "sd.peakOFI")$value,
                    position = "dodge", width=0.2))

## Correlation plot - correlation plot for slope trends? How to do?
# save figures - DOC, abs254, C1-4
pdf(file=paste0(fig.dir,"/CRFigures_SoilExtracts.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(extract.DOCplot, extract.abs254, extract.C1, extract.C2, extract.C3, extract.C4, extract.C5, ncol = 2)
dev.off()

############################################################################################
# Soil Lysimeters
# Same type of analysis as for the siol extracts - change of parameters with different depths
lys.icshim <- read.delim("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Lysimeter/CR_Lys_v3.csv", header = TRUE, sep = ",")
lys.icshim$Date <- gsub("-", "", paste(lys.icshim$Date))
#lys.icshim$sample.ID <- apply(lys.icshim[ , c("Date", "Location", "Depth") ] , 1 , paste , collapse = "" )

# get spectral parameters - abs and fluor
lys.abdflu <- read.delim("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Lysimeter/CRLys_Fluorescence/CRLys_CorrectedEEMS/CRLysAbsFluor_Indicies.csv", header = TRUE, sep = ",")
lys.abdflu$sample.ID <- gsub("_", "", paste(lys.abdflu$samplename))
lys.abdflu$sample.ID <- gsub("CRLys", "", paste(lys.abdflu$sample.ID))
# get the CM fits = percent and Fmax
lys.CM.Fmax <- read.delim("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Lysimeter/CRLys_Fluorescence/CRLys_CMPARAFACResults/CRLys_componentsandloadings_CM_Fmax.csv", header = TRUE, sep = ",")
lys.CM.per <- read.delim("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Lysimeter/CRLys_Fluorescence/CRLys_CMPARAFACResults/CRLys_componentsandloadings_CM.csv", header = TRUE, sep = ",")
lys.CM <- merge(lys.CM.Fmax, lys.CM.per, by = "sample.ID")
lys.CM$sample.ID <- gsub("CRLys_", "", paste(lys.CM$sample.ID))
lys.CM$sample.ID <- gsub("CRLysCorrCM", "", paste(lys.CM$sample.ID))
lys.CM$sample.ID <- gsub("\\_[0-9]*$", "",  paste(lys.CM$sample.ID))
lys.CM$sample.ID <- gsub("_", "", paste(lys.CM$sample.ID))

# get the 5 component fits
lys.5comp.per <- soil.5comp.per[1:31,]
lys.5comp.per$sample.ID <- gsub("_CRLys", "", lys.5comp.per$sample.ID)
lys.5comp.per$sample.ID <- gsub("CRLys_", "", lys.5comp.per$sample.ID)
lys.5comp.per$sample.ID <- gsub("_", "", lys.5comp.per$sample.ID)
colnames(lys.5comp.per)[1:5] <- c("C1", "C2", "C3", "C4", "C5")

lys.5comp.fmax <- soil.5comp[1:31,]
lys.5comp.fmax$sample.ID <- gsub("_CRLys", "", lys.5comp.fmax$sample.ID)
lys.5comp.fmax$sample.ID <- gsub("CRLys_", "", lys.5comp.fmax$sample.ID)
lys.5comp.fmax$sample.ID <- gsub("_", "", lys.5comp.fmax$sample.ID)
colnames(lys.5comp.fmax)[1:5] <- c("C1", "C2", "C3", "C4", "C5")
#####
# get file with the actual depths of the lysimeters
lys.depth <- read.delim("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Lysimeter/CRLys_depth.csv",sep = ",", header = TRUE)

# merge all together
lys.all <- Reduce(function(x, y) merge(x, y, by = "sample.ID", all=TRUE), list(lys.icshim, lys.abdflu, lys.CM, lys.5comp.per, lys.5comp.fmax))

# associate with the right depth
lys.all$locdepth <- paste(lys.all$Location, lys.all$Depth, sep = "")
lys.all$depth_cm[lys.all$locdepth == c("EastShallow")]  <- lys.depth$Depth[lys.depth$Lysimeter == c("EastShallow")] 
lys.all$depth_cm[lys.all$locdepth == c("EastDeep")]  <- lys.depth$Depth[lys.depth$Lysimeter == c("EastDeep")] 
lys.all$depth_cm[lys.all$locdepth == c("WestShallow")]  <- lys.depth$Depth[lys.depth$Lysimeter == c("WestShallow")] 
lys.all$depth_cm[lys.all$locdepth == c("WestDeep")]  <- lys.depth$Depth[lys.depth$Lysimeter == c("WestDeep")] 
lys.all$depth_cm <- as.numeric(lys.all$depth_cm)

# calculate SUVA
lys.all$SUVA <- lys.all$abs254/lys.all$NPOC

# plot by location/depth - boxplots
CRLys.box.DOC <- ggplot(lys.all, aes(x=depth_cm, y=NPOC, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +        #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - DOC Trends')) +
  ylab(expression('[DOC] (mg/L)')) +
  xlab(expression('Lysimeter Location'))

CRLys.box.TN <- ggplot(lys.all, aes(x=depth_cm, y=TN, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - TN Trends')) +
  ylab(expression('[TN] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.Cl <- ggplot(lys.all, aes(x=depth_cm, y=cl...mg.L., group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - Cl Trends')) +
  ylab(expression('[Cl] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.no3 <- ggplot(lys.all, aes(x=depth_cm, y=no3.mgL, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - NO3 Trends')) +
  ylab(expression('[NO3] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.so4 <- ggplot(lys.all, aes(x=depth_cm, y=so4.mg.L, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - SO4 Trends')) +
  ylab(expression('[SO4] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.abs254 <- ggplot(lys.all, aes(x=depth_cm, y=abs254, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - abs254 Trends')) +
  ylab(expression('[abs254] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.SUVA <- ggplot(lys.all, aes(x=depth_cm, y=SUVA, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - SUVA Trends')) +
  ylab(expression('[SUVA]254 ')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.abs272 <- ggplot(lys.all, aes(x=depth_cm, y=abs272, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - abs272 Trends')) +
  ylab(expression('[abs272] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.e2e3 <- ggplot(lys.all, aes(x=depth_cm, y=e2e3, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - e2e3 Trends')) +
  ylab(expression('[e2e3] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.e4e6 <- ggplot(lys.all, aes(x=depth_cm, y=e4e6, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - e4e6 Trends')) +
  ylab(expression('[e4e6] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.SR <- ggplot(lys.all, aes(x=depth_cm, y=SR, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - SR Trends')) +
  ylab(expression('[SR] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.FI <- ggplot(lys.all, aes(x=depth_cm, y=FI, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - FI Trends')) +
  ylab(expression('[FI] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.HIX_ohno_area <- ggplot(lys.all, aes(x=depth_cm, y=HIX_ohno_area, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - HIX_ohno_area Trends')) +
  ylab(expression('[HIX_ohno_area] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.FrI <- ggplot(lys.all, aes(x=depth_cm, y=FrI, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - FrI Trends')) +
  ylab(expression('[FrI] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.peakA <- ggplot(lys.all, aes(x=depth_cm, y=peakA, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - peakA Trends')) +
  ylab(expression('[peakA] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.peakC <- ggplot(lys.all, aes(x=depth_cm, y=peakC, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - peakC Trends')) +
  ylab(expression('[peakC] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.peakB <- ggplot(lys.all, aes(x=depth_cm, y=peakB, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - peakB Trends')) +
  ylab(expression('[peakB] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.peakT <- ggplot(lys.all, aes(x=depth_cm, y=peakT, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - peakT Trends')) +
  ylab(expression('[peakT] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.OFI <- ggplot(lys.all, aes(x=depth_cm, y=OFI, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - OFI Trends')) +
  ylab(expression('[OFI] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.perprotein <- ggplot(lys.all, aes(x=depth_cm, y=perprotein, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - perprotein Trends')) +
  ylab(expression('[perprotein] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.redox <- ggplot(lys.all, aes(x=depth_cm, y=redox, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - redox Trends')) +
  ylab(expression('[redox] (mg/L)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.C1 <- ggplot(lys.all, aes(x=depth_cm, y=C1, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - C1 Trends')) +
  ylab(expression('[C1] (%)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.C2 <- ggplot(lys.all, aes(x=depth_cm, y=C2, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - C2 Trends')) +
  ylab(expression('[C2] (%)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.C3 <- ggplot(lys.all, aes(x=depth_cm, y=C3, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - C3 Trends')) +
  ylab(expression('[C3] (%)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.C4 <- ggplot(lys.all, aes(x=depth_cm, y=C4, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - C4 Trends')) +
  ylab(expression('[C4] (%)')) +
  xlab(expression('Lysimeter Location'))
CRLys.box.C5 <- ggplot(lys.all, aes(x=depth_cm, y=C5, group = depth_cm, fill = Location)) + 
  geom_boxplot(outlier.shape = NA)  +          #remove extreme values
  scale_fill_manual(values=c(cbPalette)[4:5]) +       # change colour to colour blind
  ggtitle(expression('CR Lysimeter - C5 Trends')) +
  ylab(expression('[C5] (%)')) +
  xlab(expression('Lysimeter Location'))
# save the figures - DOC, abs254 and the C1-C5
pdf(file=paste0(fig.dir,"/CRFigures_SoilLys.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(CRLys.box.DOC, CRLys.box.abs254, CRLys.box.C1, CRLys.box.C2, CRLys.box.C3, CRLys.box.C4, CRLys.box.C5, ncol = 2)
dev.off()

###################################
# do linear fits - need to do??? Is it fair/ok to compare the two, considering the different positions? etter to simply compare at the two different locations?
lm.DOC <- c(summary(lm(lys.all$NPOC ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$NPOC ~ lys.all$depth_cm))$r.squared)
lm.TN <- c(summary(lm(lys.all$TN ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$TN ~ lys.all$depth_cm))$r.squared)
lm.Cl <- c(summary(lm(lys.all$cl...mg.L. ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$cl...mg.L. ~ lys.all$depth_cm))$r.squared)
lm.no3 <- c(summary(lm(lys.all$no3.mgL ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$no3.mgL ~ lys.all$depth_cm))$r.squared)
lm.so4 <- c(summary(lm(lys.all$so4.mg.L ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$so4.mg.L ~ lys.all$depth_cm))$r.squared)
lm.abs254 <- c(summary(lm(lys.all$abs254 ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$abs254 ~ lys.all$depth_cm))$r.squared)
lm.abs272 <- c(summary(lm(lys.all$abs272 ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$abs272 ~ lys.all$depth_cm))$r.squared)
lm.e2e3 <- c(summary(lm(lys.all$e2e3 ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$e2e3 ~ lys.all$depth_cm))$r.squared)
lm.e4e6 <- c(summary(lm(lys.all$e4e6 ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$e4e6 ~ lys.all$depth_cm))$r.squared)
lm.SR <- c(summary(lm(lys.all$SR ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$SR ~ lys.all$depth_cm))$r.squared)
lm.FI <- c(summary(lm(lys.all$FI ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$FI ~ lys.all$depth_cm))$r.squared)
lm.HIX_ohno_area <- c(summary(lm(lys.all$HIX_ohno_area ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$HIX_ohno_area ~ lys.all$depth_cm))$r.squared)
lm.FrI <- c(summary(lm(lys.all$FrI ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$FrI ~ lys.all$depth_cm))$r.squared)
lm.peakA <- c(summary(lm(lys.all$peakA ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$peakA ~ lys.all$depth_cm))$r.squared)
lm.peakC <- c(summary(lm(lys.all$peakC ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$peakC ~ lys.all$depth_cm))$r.squared)
lm.peakB <- c(summary(lm(lys.all$peakB ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$peakB ~ lys.all$depth_cm))$r.squared)
lm.peakT <- c(summary(lm(lys.all$peakT ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$peakT ~ lys.all$depth_cm))$r.squared)
lm.OFI <- c(summary(lm(lys.all$OFI ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$OFI ~ lys.all$depth_cm))$r.squared)
lm.perprotein <- c(summary(lm(lys.all$perprotein ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$perprotein ~ lys.all$depth_cm))$r.squared)
lm.redox<- c(summary(lm(lys.all$redox ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$redox ~ lys.all$depth_cm))$r.squared)
lm.C1<- c(summary(lm(lys.all$C1 ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$C1 ~ lys.all$depth_cm))$r.squared)
lm.C2<- c(summary(lm(lys.all$C2 ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$C2 ~ lys.all$depth_cm))$r.squared)
lm.C3<- c(summary(lm(lys.all$C3 ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$C3 ~ lys.all$depth_cm))$r.squared)
lm.C4<- c(summary(lm(lys.all$C4 ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$C4 ~ lys.all$depth_cm))$r.squared)
lm.C5<- c(summary(lm(lys.all$C5 ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$C5 ~ lys.all$depth_cm))$r.squared)
lm.C1.y<- c(summary(lm(lys.all$C1.y ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$C1.y ~ lys.all$depth_cm))$r.squared)
lm.C2.y<- c(summary(lm(lys.all$C2.y ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$C2.y ~ lys.all$depth_cm))$r.squared)
lm.C3.y<- c(summary(lm(lys.all$C3.y ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$C3.y ~ lys.all$depth_cm))$r.squared)
lm.C4.y<- c(summary(lm(lys.all$C4.y ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$C4.y ~ lys.all$depth_cm))$r.squared)
lm.C5.y<- c(summary(lm(lys.all$C5.y ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$C5.y ~ lys.all$depth_cm))$r.squared)
lm.C6.y<- c(summary(lm(lys.all$C6.y ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$C6.y ~ lys.all$depth_cm))$r.squared)
lm.C7.y<- c(summary(lm(lys.all$C7.y ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$C7.y ~ lys.all$depth_cm))$r.squared)
lm.C8.y<- c(summary(lm(lys.all$C8.y ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$C8.y ~ lys.all$depth_cm))$r.squared)
lm.C9.y<- c(summary(lm(lys.all$C9.y ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$C9.y ~ lys.all$depth_cm))$r.squared)
lm.C10.y<- c(summary(lm(lys.all$C10.y ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$C10.y ~ lys.all$depth_cm))$r.squared)
lm.C11.y<- c(summary(lm(lys.all$C11.y ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$C11.y ~ lys.all$depth_cm))$r.squared)
lm.C12.y<- c(summary(lm(lys.all$C12.y ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$C12.y ~ lys.all$depth_cm))$r.squared)
lm.C13.y<- c(summary(lm(lys.all$C13.y ~ lys.all$depth_cm))$coefficients[1:2], summary(lm(lys.all$C13.y ~ lys.all$depth_cm))$r.squared)

# bind together to write to csv
lys.linear <- rbind(lm.DOC, lm.TN, lm.abs254, lm.SUVA, lm.abs272,
                    lm.Cl, lm.no3, lm.so4,
                    lm.e2e3, lm.e4e6, lm.SR, lm.FI, lm.HIX_ohno_area, lm.FrI, lm.peakA, 
                    lm.peakC, lm.peakB, lm.peakT, lm.OFI,
                    lm.perprotein, lm.redox, 
                    lm.C1, lm.C2, lm.C3, lm.C4, lm.C5,
                    lm.C1.y, lm.C2.y, lm.C3.y, lm.C4.y, lm.C5.y, lm.C6.y, lm.C7.y, lm.C8.y,
                    lm.C9.y, lm.C10.y, lm.C11.y, lm.C12.y, lm.C13.y)

colnames(lys.linear) <- c("intercept", "slope", "R2")
# write the table with slopes/r2 values for all of the extracts
write.csv(lys.linear, file = paste0(soillys.path, "/lmfits_lys.csv"))

###############################################################################
# Show soil characteristics alongside timeseries data - show discharge, alongside DOC and DOC from each of the lysimeters
# get discharge over the year time over which lysimeters and soils were sampled.
# Show discharge versus DOC/abs254/SR; also show discharge with the 4 component PARAFAC fits (grab samples)
# Soil extract date January 22, 2014.
# Soil lysimeters: from July 2013 - July 2014

# timeseries of discharge - july 2013 - july 2014
stream.dischsub <- subset(abs.Q, date > "2013-07-01" & date < "2014-07-01")
CR.grab.sub <- subset(CR.grab, date > "2013-07-01" & date < "2014-07-01")

# plot discharge for the year
time.disc <- ggplot(subset(stream.dischsub, stream.dischsub$DOCcorr >=2), aes(date, Q.L.s)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("Q (L/s)") + theme()

# plot DOC for the year
time.DOC <- ggplot(subset(stream.dischsub, stream.dischsub$DOCcorr >=2), aes(date, DOCcorr)) + 
  geom_point(size = 0.4, color = cbPalette[6]) +
  xlab("Date") + ylab("DOC") + theme()

# plot abs254 for the year
time.abs254 <- ggplot(subset(stream.dischsub, stream.dischsub$DOCcorr >=2), aes(date, abs254)) + 
  geom_point(size = 0.4, color = cbPalette[7]) +
  xlab("Date") + ylab("abs254") + theme()

time.SUVA <- ggplot(subset(stream.dischsub, stream.dischsub$DOCcorr >=2), aes(date, SUVA)) + 
  geom_point(size = 0.4, color = cbPalette[7]) +
  xlab("Date") + ylab("SUVA254") + theme()

time.e2e3 <- ggplot(subset(stream.dischsub, stream.dischsub$DOCcorr >=2), aes(date, e2e3)) + 
  geom_point(size = 0.4, color = cbPalette[4]) +
  xlab("Date") + ylab("E2E3") + 
  scale_y_continuous(limits = c(0, 8)) + theme()

time.e4e6 <- ggplot(subset(stream.dischsub, stream.dischsub$DOCcorr >=2), aes(date, e4e6)) + 
  geom_point(size = 0.4, color = cbPalette[3]) +
  xlab("Date") + ylab("E4E6") + 
  scale_y_continuous(limits = c(0, 1000)) + theme()

# plot slope ratio for the year
time.SR<- ggplot(subset(stream.dischsub, stream.dischsub$DOCcorr >=2), aes(date, SR)) + 
  geom_point(size = 0.4, color = cbPalette[8]) +
  xlab("Date") + ylab("Slope Ratio") + theme()

# plot discharge and PARAFAC points (4 component fits) for the year
# massage the component percents to have all points on one graph
#PARAFAC FMAX
Sub.PARA.per <- melt(CR.grab.sub[,c(64, 49:52)], id = "date")
pd <- position_dodge(.65)
time.PARAFACPer <- ggplot(Sub.PARA.per, aes(date, value, colour = variable, shape = variable)) + 
  geom_point(position = pd, size = 2) +
  scale_color_manual(values=cbPalette[1:4], 
                    name="PARAFAC Component",
                    breaks=c("C1_per", "C2_per", "C3_per", "C4_per"),
                    labels=c("C1", "C2", "C3", "C4")) +
  xlab("Date") + ylab("PARAFAC Component %") + theme(legend.position="top")  +
  ylim(0, 70) +
  scale_shape_manual(values=c(17,18,3,4), 
                     name="PARAFAC Component",
                     breaks=c("C1_per", "C2_per", "C3_per", "C4_per"),
                     labels=c("C1", "C2", "C3", "C4")) 
# PARAFAC fmax
Sub.PARA.fmax <- melt(CR.grab.sub[,c(64, 44:47)], id = "date")
time.PARAFACfmax <- ggplot(Sub.PARA.fmax, aes(date, value, color = variable, shape = variable)) + 
  geom_point(position = pd, size = 2) +
  xlab("Date") + ylab("PARAFAC Component Fmax") + 
  theme(legend.position="top")  +
  ylim(0, 1) +
  scale_color_manual(values=cbPalette[1:4], 
                      name="PARAFAC Component",
                      breaks=c("PARAFAC_C1", "PARAFAC_C2", "PARAFAC_C3", "PARAFAC_C4"),
                      labels=c("C1", "C2", "C3", "C4")) +
  scale_shape_manual(values=c(17,18,3,4), 
                     name="PARAFAC Component",
                     breaks=c("PARAFAC_C1", "PARAFAC_C2", "PARAFAC_C3", "PARAFAC_C4"),
                     labels=c("C1", "C2", "C3", "C4")) 

# Arrange the timeseries on top of each other
pdf(file=paste0(fig.dir,"/CRFigures_SoilTimeseries.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(time.disc, time.DOC, time.abs254, time.SUVA, time.PARAFACPer,time.PARAFACfmax, ncol = 1)
dev.off()

# Supplementary figure
pdf(file=paste0(fig.dir,"/CRFigures_SoilTimeseries_supplemental.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(time.disc, time.DOC, time.SR, time.e2e3, time.e4e6, time.PARAFACfmax, ncol = 1)
dev.off()

#######################
# Figure - Scatter plot of characteristics by different method
# DOC, parafac components
# get the lysimeter + extracts - DOC, fmax and percent PARAFAC
lys.sub <- lys.all[,c(12,16,83:88)]
lys.sub$depth[lys.sub$locdept == c("WestDeep")] <- "04"
lys.sub$depth[lys.sub$locdept == c("EastDeep")] <- "03"
lys.sub$depth[lys.sub$locdept == c("EastShallow")] <- "02"
lys.sub$depth[lys.sub$locdept == c("WestShallow")] <- "01"
lys.sub$type <- "lys"

lys.sub.melt <- melt(lys.sub, id = c("depth", "type"))

extract.sub <- extract.all[,c(41,75,80,81,82, 35:39)]
colnames(extract.sub)[3] <- "type"
colnames(extract.sub)[2] <- "NPOC"
extract.sub$depth_cm <- NULL
extract.sub.melt <- melt(extract.sub, id = c("depth", "type"))

melt.all <- rbind(lys.sub.melt, extract.sub.melt) # bind all together

# do scatter plots
DOC.scat <- ggplot(subset(melt.all, melt.all$variable == "NPOC"), aes(type, as.numeric(value), colour = type, shape = depth)) + 
  geom_point(position = pd, size = 2) +
  scale_color_manual(values=cbPalette[1:4], 
                     name="Method",
                     breaks=c("lys", "dH2O", "K2SO4", "KCl"),
                     labels=c("Lysimeter", "dH2O", "K2SO4", "KCl")) +
  xlab("Method") + ylab("[DOC] (mg/L)") + 
  ylim(0, 25) +
  scale_shape_manual(values=c(17,18,3,4), 
                     name="Depth Code",
                     breaks=c("01", "02", "03", "04"),
                     labels=c("1", "2", "3", "4")) 

C1.scat <- ggplot(subset(melt.all, melt.all$variable == "C1.y"), aes(type, as.numeric(value), colour = type, shape = depth)) + 
  geom_point(position = pd, size = 2) +
  scale_color_manual(values=cbPalette[1:4], 
                     name="Method",
                     breaks=c("lys", "dH2O", "K2SO4", "KCl"),
                     labels=c("Lysimeter", "dH2O", "K2SO4", "KCl")) +
  xlab("Method") + ylab("C1 Fmax") + 
  #ylim(0, 25) +
  scale_shape_manual(values=c(17,18,3,4), 
                     name="Depth Code",
                     breaks=c("01", "02", "03", "04"),
                     labels=c("1", "2", "3", "4")) 

C2.scat <- ggplot(subset(melt.all, melt.all$variable == "C2.y"), aes(type, as.numeric(value), colour = type, shape = depth)) + 
  geom_point(position = pd, size = 2) +
  scale_color_manual(values=cbPalette[1:4], 
                     name="Method",
                     breaks=c("lys", "dH2O", "K2SO4", "KCl"),
                     labels=c("Lysimeter", "dH2O", "K2SO4", "KCl")) +
  xlab("Method") + ylab("C2 Fmax") + 
  #ylim(0, 25) +
  scale_shape_manual(values=c(17,18,3,4), 
                     name="Depth Code",
                     breaks=c("01", "02", "03", "04"),
                     labels=c("1", "2", "3", "4")) 

C3.scat <- ggplot(subset(melt.all, melt.all$variable == "C3.y"), aes(type, as.numeric(value), colour = type, shape = depth)) + 
  geom_point(position = pd, size = 2) +
  scale_color_manual(values=cbPalette[1:4], 
                     name="Method",
                     breaks=c("lys", "dH2O", "K2SO4", "KCl"),
                     labels=c("Lysimeter", "dH2O", "K2SO4", "KCl")) +
  xlab("Method") + ylab("C3 Fmax") + 
  #ylim(0, 25) +
  scale_shape_manual(values=c(17,18,3,4), 
                     name="Depth Code",
                     breaks=c("01", "02", "03", "04"),
                     labels=c("1", "2", "3", "4")) 

C4.scat <- ggplot(subset(melt.all, melt.all$variable == "C4.y"), aes(type, as.numeric(value), colour = type, shape = depth)) + 
  geom_point(position = pd, size = 2) +
  scale_color_manual(values=cbPalette[1:4], 
                     name="Method",
                     breaks=c("lys", "dH2O", "K2SO4", "KCl"),
                     labels=c("Lysimeter", "dH2O", "K2SO4", "KCl")) +
  xlab("Method") + ylab("C4 Fmax") + 
  #ylim(0, 25) +
  scale_shape_manual(values=c(17,18,3,4), 
                     name="Depth Code",
                     breaks=c("01", "02", "03", "04"),
                     labels=c("1", "2", "3", "4")) 

C5.scat <- ggplot(subset(melt.all, melt.all$variable == "C5.y"), aes(type, as.numeric(value), colour = type, shape = depth)) + 
  geom_point(position = pd, size = 2) +
  scale_color_manual(values=cbPalette[1:4], 
                     name="Method",
                     breaks=c("lys", "dH2O", "K2SO4", "KCl"),
                     labels=c("Lysimeter", "dH2O", "K2SO4", "KCl")) +
  xlab("Method") + ylab("C5 Fmax") + 
  #ylim(0, 25) +
  scale_shape_manual(values=c(17,18,3,4), 
                     name="Depth Code",
                     breaks=c("01", "02", "03", "04"),
                     labels=c("1", "2", "3", "4")) 
# Arrange the timeseries on top of each other
pdf(file=paste0(fig.dir,"/CRFigures_SoilComparisons.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(DOC.scat, C1.scat, C2.scat, C3.scat,C4.scat, C5.scat, ncol = 2)
dev.off()

## get mean concentrations for eeach of the indicies - by lys, extract method
lys.mean <- as.data.frame(colMeans(lys.all[,c(6:14, 16:87, 89:90)], na.rm = TRUE)) # get mean
# get stdev
lys.sd <- as.data.frame(apply(as.data.frame(lys.all[,c(6:14, 16:87, 89:90)]), 2, function(x) sd(x, na.rm = TRUE)))

lys.stats <- cbind(lys.mean, lys.sd)
write.csv(lys.stats,  paste0(fig.dir, "/CRlysstats.csv"))

require(plyr)
 
extract.mean <- ddply(extract.all[,c(2:39, 41:78, 80, 82, 83)],~extract,colwise(mean))
extract.sd <- ddply(extract.all[,c(2:39, 41:78, 80, 82, 83)],~extract,colwise(sd))
extract.sd$type <-'sd'
extract.mean$type <- 'mean'
extract.stats <- rbind (extract.mean, extract.sd)
write.csv(extract.stats,  paste0(fig.dir, "/CRextractstats.csv"))

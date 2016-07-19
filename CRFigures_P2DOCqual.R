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
       reshape2, pgirmess, nlme)

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


#####
# add pH/EC/DO


## Seven Anion Data from IC 
anions <- read.csv("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Stream/CR_IC/Compiled_ICSevenAnions_CR_Discharge_Sampledates_v4.csv", header = TRUE, sep = ",")

## Fluorescence Data


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
source("/Users/user/SpecScripts/CRSpectralIndAbs_function.R")
abs.ind.1 <- as.data.frame(Abs.ind(spec = abspec.corr[1:10000,], pathlength = 35))
abs.ind.2 <- as.data.frame(Abs.ind(spec = abspec.corr[10001:20000,], pathlength = 35))
abs.ind.3 <- as.data.frame(Abs.ind(spec = abspec.corr[20001:30000,], pathlength = 35))
abs.ind.4 <- as.data.frame(Abs.ind(spec = abspec.corr[30001:39046,], pathlength = 35))
spec.39047 <- Abs.ind(spec = abspec.corr[39047,], pathlength = 35)
spec.39048 <- Abs.ind(spec = abspec.corr[39048,], pathlength = 35)
spec.39049 <- Abs.ind(spec = abspec.corr[39049,], pathlength = 35)
abs.ind.5 <- as.data.frame(Abs.ind(spec = abspec.corr[39050:dim(abspec.corr)[1],], pathlength = 35))

spec.39048 <- c("NA", "NA","NA", "NA", "NA", "NA","NA", "NA" )
spec.39049 <- c("NA", "NA","NA", "NA", "NA", "NA","NA", "NA" )

# bind together
abs.all <- rbind(as.data.frame(abs.ind.1), abs.ind.2, abs.ind.3, abs.ind.4, 
                 spec.39047, spec.39048, spec.39049, abs.ind.5)

# add in date
abs.all$date <- abspec.corr$date
abs.all$DOCcorr <- abspec.corr$DOCcorr
abs.all$NO3 <- abspec.corr$NO3

# this takes forever! Save data to the absorbance folder...
write.csv(abs.all, 
          file = paste0("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Spectrodata", "/CRabsind_all.csv"))

# Make sure date is formatted

# Do timeseries plots of spectral indicies

time.SUVA <- ggplot(abs.all, aes(date, SR)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("") + theme()

time.e2e3 <- ggplot(abs.ind, aes(date, Q.L.s)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("") + theme()

time.e4e6 <- 
  
time.SR <-   

# get the CM fits
# Get PARAFAC fits 



# for grab sample data, compile by sample ID, and ensure sample date/time correctly associated.


# Figure 1 PCA analysis
# Question: Which varaibles explain the most variance within the dataset (partitioned according to pre/post, wet/dry)?

# Correlation Matrix - Water Quality variables and others (discharge, climate)

############################################################################################
# Soil Extracts
# linear model between soil depth and quality/DOC variables. Which variables best explain changes in soil depth?


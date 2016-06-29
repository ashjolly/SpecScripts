# Script for figures associated with CR Paper 1 - Concentration dynamics
# Ashlee J's PhD
#
# 27 June2016
############################################################################################

# Clear environment
rm(list = ls())
ls()

###### Necessary toolboxes
library(reshape)
library(plyr)
library(openair)
library(EcoHydRology)
library(ggplot2)
library(tidyr)
library(dplyr)
library(plyr)
library(stringr)
library(stringi)
library(ggplot2)
library(pacman)
library(ggplot2)
p_load(lubridate,plyr,ggplot2,reshape2,car,grid,gridBase,gridExtra, taRifx,magrittr, ggpmisc)

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



###### Set paths for data
fig.dir <- "/Users/user/Dropbox/PhD Work/Thesis chapters/CR Chapter/CRDOCconcP1/Figures" # directory for saving figures
spectro.path <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Spectrodata/parfp_all"
discharge.path <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Stream/CR_discharge"
climate.path <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Climate"
groundwater.path <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Stream/CR_trutracks/SN_0904203"

# set path
#setwd("/Users/user/Dropbox/par and fp compilation")

##### Read in data
#### DOC data from spectro - compiled file. No calculated spectral parameters
spectro.all<- read.csv(file = "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Spectrodata/SpectralOutput_2009t2014_2_noparameters.csv", head=TRUE,sep=",")
attach(spectro.all)
names(spectro.all)
spectro.all$date <- as.POSIXct(strptime(spectro.all$date, format = "%y-%m-%d %H:%M", tz = "America/Los_Angeles"))
spectro.all <- data.frame(spectro.all)

#### Discharge from stream
# Function for compiling discharge from all
# Note - took out the open air function to convert to 30 min. Not working?
source("/Users/user/SpecScripts/CRDischargeCompilation_function.R")
discharge <- discharge.comp(discharge.dir = discharge.path)

#### Climate variables from Biomet
# Function for compiling climate from biomet
source("/Users/user/SpecScripts/CRClimateCompilation_function.R")
climate <- climate.comp(climate.dir = climate.path)

#### Groundwater depth from groundwater well - trutrack closest to the weir
source("/Users/user/SpecScripts/CRGroundwaterComp_function.R")
groundwater <- groundwater.comp(ground.dir = groundwater.path)

#####
# add pH/EC/DO


## Seven Anion Data from IC 
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

########################################################################################################################
#### Figure 2
# timeseries of discharge, precip, groundwater, DOC arranged on top of one another
# Table of max, min sum of climatic and discharge variables
# As per figure 1 in Strohmeier, S, K H Knorr, M Reichert, and S Frei. 2013. 
# “Concentrations and Fluxes of Dissolved Organic Carbon in Runoff From a Forested Catchment: Insights From High Frequency Measurements.” …. doi:10.5194/bg-10-905-2013.

time.disc <- ggplot(discharge, aes(date, Q.L.s)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("Q (L/s)") + theme()

time.precip <- ggplot(climate, aes(date, Precip)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("Precipitation (mm)") + theme()

time.soilT <- ggplot(climate, aes(date, SoilT)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("Soil Temperature (degC)") + theme()

time.airT <- ggplot(climate, aes(date, Tair)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("Air Temperature (degC)") + theme()

time.groundwater <- ggplot(groundwater, aes(date, water.height.average..mm.0904203)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("Groundwater Height (mm)") + theme()

time.DOC <- ggplot(spectro.all, aes(date, DOCcorr)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("[DOC] (mg/L)") + theme()

# Save as stacked figure
pdf(file=paste0(fig.dir,"/CRFigures_timeseries.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(time.DOC, ncol=1)
#grid.text("A", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
#grid.text("B", x=unit(0.5, "npc"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
dev.off()

##############
## Effect of harvest on discharge
# 8.64 is the conversion factor from L s 1 ha 1 to mm d 1.From: Sørensen, Rasmus, Eva Ring, Markus Meili, Lars Högbom, Jan Seibert, Thomas Grabs, Hjalmar Laudon, and Kevin Bishop. 2009. “Forest Harvest Increases Runoff Most During Low Flows in Two Boreal Streams.” AMBIO: a Journal of the Human Environment 38 (7): 357–63. doi:10.1579/0044-7447-38.7.357.

## Cumulative distribution function of discharge - by year
CDF.runoff <- ggplot(discharge, aes(x = Q.L.s)) + 
  stat_ecdf(aes(group = format(discharge$date, format='%Y'), colour = format(discharge$date, format='%Y')))
pdf(file=paste0(fig.dir,"/CRFigures_CDFrunoff.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(CDF.runoff, ncol=1)
dev.off()            
            
## Boxplot of mean daily Q (L/s) - presented per month
# calculate the mean daily Q
meandaily.Q <- ddply(discharge,.(format(discharge$date, format='%Y-%m-%d')),
                     summarise, mean.Q=mean(Q.L.s, na.rm = TRUE))
colnames(meandaily.Q)[1] <- "date"
meandaily.Q$date <- as.POSIXct(strptime(meandaily.Q$date, format = "%Y-%m-%d"))

# boxplot the mean daily Q by month
boxplot.discharge <- ggplot(meandaily.Q, aes(x=format(meandaily.Q$date, format='%Y-%m'), y=mean.Q)) + 
  geom_boxplot() + xlab("Date") + 
  scale_x_date(breaks = "1 year", minor_breaks = "1 month") +
  ylab("Mean Daily Q (L/s)") + theme()

# save
pdf(file=paste0(fig.dir,"/CRFigures_boxplotQ.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(boxplot.discharge, ncol=1)
dev.off()

# do monthly mean flow  
meanmonthly.Q <- ddply(discharge,.(format(discharge$date, format='%Y-%m')),
                       summarise, mean.Q=mean(Q.L.s, na.rm = TRUE))
colnames(meanmonthly.Q)[1] <- "date"
meanmonthly.Q$date1 <- as.character(meanmonthly.Q$date)
meanmonthly.Q$date1 <- as.POSIXct(strptime(meanmonthly.Q$date1, format = '%Y-%m')) # not working!

meanmonth.Q <- ggplot(data=meanmonthly.Q, aes(x=date, y=mean.Q, group=1)) +
  geom_line() +
  geom_point() +
  scale_x_date(breaks = "1 year", minor_breaks = "1 month") +
  expand_limits(y=0) +
  xlab("Date") + ylab("Monthly Mean Q (L/s)") +
  ggtitle("Monthly Mean Q")

# save
pdf(file=paste0(fig.dir,"/CRFigures_monthlymeanQ.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(meanmonth.Q, ncol=1)
dev.off()

## Percent change in discharge table 
# Sørensen, Rasmus, Eva Ring, Markus Meili, Lars Högbom, Jan Seibert, Thomas Grabs, Hjalmar Laudon, and Kevin Bishop. 2009. “Forest Harvest Increases Runoff Most During Low Flows in Two Boreal Streams.” AMBIO: a Journal of the Human Environment 38 (7): 357–63. doi:10.1579/0044-7447-38.7.357.
# Laudon, Hjalmar, Johannes Hedtjärn, Jakob Schelker, Kevin Bishop, Rasmus Sørensen, and Anneli Ågren. 2009. “Response of Dissolved Organic Carbon Following Forest Harvesting in a Boreal Forest.” AMBIO: a Journal of the Human Environment 38 (7): 381–86. doi:10.1579/0044-7447-38.7.381.

## Calculate the runoff from discharge - how much water exported by month

# Wetting ratio Precip/Q - changes over time

## Calculation of DOC Flux
# Method 1: daily average C = daily average DOC x average Q x 24 hours
#calculation of DOC flux

#calculate carbon flux - every 30 minutes
spectro.all$cflux.mgs = (spectro.all$DOCcorr * spectro.all$Q.L.s )
q = as.numeric(match("DOCcorr",names(spectro.all)))
c.flux= spectro.all[complete.cases(spectro.all[,q]),]

n <- dim(c.flux)[1]
for (i in 2:n) {
  timeinterval <- c.flux$date[i]-c.flux$date[(i-1)]
  c.flux$timeinterval <- timeinterval
}

for (i in 1:n){
  cflux.mg = (c.flux$DOCcorr[i] * c.flux$Q.L.s[i] * c.flux$timeinterval[i]*60)
  c.flux$cflux.mg[i] <- cflux.mg
}
c.flux$cflux.g = c.flux$cflux.mg/1000

#select only columns with cflux in mg and g and merge back into spectro.all
date = as.numeric(match("date",names(c.flux)))
cflux.mg.column = as.numeric(match("cflux.mg",names(c.flux)))
cflux.g.column = as.numeric(match("cflux.g",names(c.flux)))
column <- c(date, cflux.mg.column,cflux.g.column)
c.flux <- c.flux[,column]
spectro.all <- merge(spectro.all, c.flux, by = "date", all = TRUE)

#### c flux per hectares
# in hectares
# 1 hectare = 10000 m2.  91 hectares total
area.h = 91
area.m2 = area.h* 10000 
spectro.all$cflux.gm2 = spectro.all$cflux.g/area.m2
spectro.all$cflux.mgm2 = spectro.all$cflux.mg/area.m2



#### Figure 3
# How does Forest Harvest Affect DOC Concentration?
# Figure: Boxplot of mean DOC flux; table of ranked ANOVA scored to look for significant differences

#### Figure 4
# Figure 4: Hysterisis loops for event pre/post harvest. 
# Percent of DOC from top precipitation events pre and post harvest. 
# Show the amount of discharge generated and the DOC flux for an event pre and post harvest. Is this relationship changed?

#### Figure 5
# Timing of DOC and precipitation/groundwater levels: is DOC being mobilized into the stream faster after harvest?
# Figure 5 How to show this? Average time between peak precip to peak discharge, DOC?


#### Figure 6
# What are the drivers of DOC concentration? Model of DOC and climatic/discharge/stream chemistry variables.
# Figure 6: Relationship between [DOC] and variables. Outline model of DOC concentration within the pre and post harvest period. Correlation table? Also, include soil temp?

# Script for figures associated with CR Paper 1 - Concentration dynamics
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
       reshape2, pgirmess, nlme, scales, dunn.test, 'corrplot', lubridate) 

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
fig.dir <- "/Users/user/Dropbox/PhD Work/Thesis chapters/CR Chapter/Figures" # directory for saving figures
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
spectro.all <- logstatus.f(spectro.all)

##CHECK AS THE DATA FROM THE COMPILATION IS IN GMT!!!
spectro.all <- data.frame(spectro.all)

#### Discharge from stream
# Function for compiling discharge from all
# Note - took out the open air function to convert to 30 min. Not working?
#source("/Users/user/SpecScripts/CRDischargeCompilation_function.R")
#discharge <- discharge.comp(discharge.dir = discharge.path)

# add column for pre/post harvest
#discharge <- logstatus.f(discharge)
# add column for wet/dry  period
#discharge <- wetdry.f(discharge)

#mark.discharge <- source_data("https://github.com/UBCecohydro/ecohydro.datasets/blob/master/Campbell.River/Q.WQ.July2016.RDS?raw=true")
mark.discharge <- readRDS("/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Stream/CR_discharge/Q.WQ.July2016.RDS")
mark.discharge$date <- as.POSIXct(strptime(mark.discharge$date, format = "%Y-%m-%d %H:%M"), tz = "America/Los_Angeles")
#discharge.ave <- timeAverage(discharge, avg.time = "30 min", start.date = "2007-11-15", fill = TRUE)
attr(mark.discharge$date, "tzone") <- "America/Los_Angeles"
discharge <- mark.discharge

# add column for pre/post harvest
discharge <- logstatus.f(discharge)
# add column for wet/dry  period
discharge <- wetdry.f(discharge)
# cut off the discharge so that it ends at Dec 2014
discharge <- subset(discharge, discharge$date <= "2014-12-31")

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
# in mark discharge

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
##############################
# Constants etc
area.hect <- 91 #watershed area in hectares

####################################################################################
# Precipitation sum per month - when does precipitation come?
# figure 1 - weather over the year
# calculate the sum of precip over each month
#climate$date1<- as.POSIXct(strptime(climate$date, format = "%Y-%m-%d HH:MM:SS", tz = ""))
climate$month <- format(climate$date, format = "%m")
meanmonth.precip <- ddply(climate,.(format(climate$date, format = "%Y-%m")),
                          summarise, sum.P = sum(as.numeric(Precip), na.rm = TRUE))

colnames(meanmonth.precip)[1] <- "date"
meanmonth.precip$date <- as.Date(paste(meanmonth.precip$date, "-01", sep = ""), format = "%Y-%m-%d")
meanmonth.precip$month <- format(meanmonth.precip$date, format = "%m")

meanmonth.precip <- wetdry.f(meanmonth.precip)
# average precip by month
mean.P <- ddply(meanmonth.precip,.(format(meanmonth.precip$date, format='%m')),
                summarise, ave.P = mean(sum.P, na.rm = TRUE), sd.p = sd(sum.P, na.rm = TRUE))
mean.P <- as.data.frame(mean.P)
colnames(mean.P) <- c("month", "meanP", "sdP")

# make boxplot of average monthly P and daily temp - not working! Boxplot instead
meanP.month <- ggplot(meanmonth.precip,
                      aes(y = sum.P, x = month, fill = hydro)) +
  geom_boxplot() +  
  ggtitle("Mean Monthly  Precipitation") + 
  xlab("Month") +
  ylab(expression("Mean Monthly Precipitation (mm/month)")) + 
  scale_fill_manual(name = "", 
                    values = cbPalette[5:6]) +
  theme(legend.position="top")

## get average temperature by month
mean.Tdaily <- ddply(climate,.(format(climate$date, format='%Y-%m-%d')),
                     summarise, ave.T = mean(Tair, na.rm = TRUE), sd.T = sd(Tair, na.rm = TRUE))
colnames(mean.Tdaily)[1] <- "date"
mean.Tdaily$date <- as.Date(mean.Tdaily$date, format = "%Y-%m-%d")
mean.Tdaily$month <- format(mean.Tdaily$date, format = "%m")
mean.Tdaily <- wetdry.f(mean.Tdaily)

# daily mean temp by month
meanT.plot <- ggplot(mean.Tdaily,
                     aes(y = ave.T, x = month, fill = hydro)) +
  geom_boxplot() + 
  ggtitle("Monthly Average Temperature") + 
  xlab("Month") +
  ylab(expression("Mean Temperature ("~degree~"C)")) +
  scale_fill_manual(name = "", 
                    values = cbPalette[5:6]) +
  theme(legend.position="none")

# save the precip and averahe daily temp boxplots as one plot
pdf(file=paste0(fig.dir,"/CRFigures_climate.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(meanP.month, meanT.plot, ncol=1)
grid.text("A", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("B", x=unit(0, "npc"), y=unit(0.5, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
dev.off()

########################################################################################################################
#### figure 2 - precip, discharge changes
# timeseries of discharge, precip, groundwater, DOC arranged on top of one another
# Table of max, min sum of climatic and discharge variables
# As per figure 1 in Strohmeier, S, K H Knorr, M Reichert, and S Frei. 2013. 
# “Concentrations and Fluxes of Dissolved Organic Carbon in Runoff From a Forested Catchment: Insights From High Frequency Measurements.” …. doi:10.5194/bg-10-905-2013.

time.disc <- ggplot(discharge, aes(date, Q.m3.s, colour = logstatus)) + 
  geom_point(size = 0.4) +
  labs(x = expression("Date"),  
       y = expression(Q~(m^{3}~s^{-1}))) + 
  scale_color_manual(values=cbPalette[1:2]) + # colours
  theme(legend.position="") + 
  ggtitle("Discharge")

# do daily sum of precip for bar graph. Too complicated otherwise
dailysum.precip<- ddply(climate,.(format(climate$date, format = "%Y-%m-%d")),
                          summarise, sum.P = sum(as.numeric(Precip), na.rm = TRUE))
colnames(dailysum.precip) <- c("date", "dailyP.mm")
dailysum.precip$date <- as.POSIXct(strptime(dailysum.precip$date, format = "%Y-%m-%d"))
dailysum.precip <- wetdry.f(dailysum.precip)

# do plots - timeseries
time.precip <- ggplot(dailysum.precip, aes(as.Date(date), dailyP.mm, fill = factor(hydro))) + 
  geom_bar(stat="identity") +
  xlab("Date") + ylab("Precipitation (mm)") + 
  theme(legend.position="top") +
  scale_x_date(date_breaks = "year", 
               labels=date_format("%Y"))

time.soilT <- ggplot(climate, aes(date, SoilT)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("Soil Temperature (degC)") + theme()

time.airT <- ggplot(climate, aes(date, Tair)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("Air Temperature (degC)") + theme()

time.groundwater <- ggplot(groundwater, aes(date, water.height.average..mm.0904203)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("Groundwater Height (mm)") + theme()

time.DOC <- ggplot(subset(spectro.all, spectro.all$DOCcorr >=1), aes(date, DOCcorr, colour = logstatus)) + 
  geom_point(size = 0.4) +
  xlab("Date") + ylab("[DOC] (mg/L)") + 
  scale_color_manual(values=cbPalette[1:2]) + # colours
  theme(legend.position="top")

# Save as stacked figure
pdf(file=paste0(fig.dir,"/CRFigures_timeseries.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(time.disc, time.precip, time.DOC, ncol=1)
#grid.text("A", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
#grid.text("B", x=unit(0.5, "npc"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
dev.off()

########################            
## Boxplot of mean daily Q (L/s) - presented per month
# calculate the mean daily Q
meandaily.Q <- ddply(discharge,.(format(discharge$date, format='%Y-%m-%d')),
                     summarise, mean.Q=mean(Q.L.s, na.rm = TRUE))
colnames(meandaily.Q)[1] <- "date"
meandaily.Q$date <- as.POSIXct(strptime(meandaily.Q$date, format = "%Y-%m-%d"))
meandaily.Q$year_month <- paste(month(as.Date(meandaily.Q$date)), year(as.Date(meandaily.Q$date)))
meandaily.Q <- logstatus.f(meandaily.Q)

# boxplot the mean daily Q by month
boxplot.discharge <- ggplot(meandaily.Q, aes(x=date, y=mean.Q, group = year_month, fill=logstatus)) + 
  geom_boxplot() + 
  scale_fill_manual(breaks = c('pre', 'post'),
                    values = c(cbPalette[1], cbPalette[2]),
                    name="Logging\nStatus") +
  xlab("Date") + 
  #scale_x_date(breaks = "1 year", minor_breaks = "1 month") +
  ylab("Mean Daily Q (L/s)") + 
  #geom_line() + # line to show where 
  ggtitle("Mean Daily Q - By Month") +
  theme() +
  theme(legend.position=c(0.85, 0.8))

# save 
#pdf(file=paste0(fig.dir,"/CRFigures_boxplotQ.pdf"), width = 8.5, height = 11) #save figure
#grid.arrange(boxplot.discharge, ncol=1)
#dev.off() put with other runoff figures

# do monthly mean flow  
meanmonthly.Q <- ddply(discharge,.(format(discharge$date, format='%Y-%m')),
                       summarise, mean.Q=mean(Q.L.s, na.rm = TRUE))
colnames(meanmonthly.Q)[1] <- "date"
meanmonthly.Q$date1 <- as.Date(as.character(meanmonthly.Q$date))
#meanmonthly.Q$date1 <- as.POSIXct(strptime(meanmonthly.Q$date, format = '%Y-%m')) # not working!

meanmonth.Q <- ggplot(data=meanmonthly.Q, aes(x=date, y=mean.Q, group=1)) +
  geom_line() +
  geom_point() +
  #scale_x_date(breaks = "1 year", minor_breaks = "1 month") +
  expand_limits(y=0) +
  xlab("Date") + ylab("Monthly Mean Q (L/s)") +
  ggtitle("Monthly Mean Q") + 
  theme()

# save
pdf(file=paste0(fig.dir,"/CRFigures_monthlymeanQ.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(meanmonth.Q, ncol=1)
dev.off()

################################################
## Percent change in discharge table  - how does harvest affect runoff generation?
# Ref: Sørensen, Rasmus, Eva Ring, Markus Meili, Lars Högbom, Jan Seibert, Thomas Grabs, Hjalmar Laudon, and Kevin Bishop. 2009. “Forest Harvest Increases Runoff Most During Low Flows in Two Boreal Streams.” AMBIO: a Journal of the Human Environment 38 (7): 357–63. doi:10.1579/0044-7447-38.7.357.
# Laudon, Hjalmar, Johannes Hedtjärn, Jakob Schelker, Kevin Bishop, Rasmus Sørensen, and Anneli Ågren. 2009. “Response of Dissolved Organic Carbon Following Forest Harvesting in a Boreal Forest.” AMBIO: a Journal of the Human Environment 38 (7): 381–86. doi:10.1579/0044-7447-38.7.381.
# calculate mean and sd Q preharvest - by month
pre.monthlyQ <- ddply(subset(discharge, discharge$logstatus =="pre" ),.(format(discharge$date, format='%m')),
                      summarise, mean.Q = mean(Q.L.s, na.rm = TRUE), sd.Q = sd(Q.L.s, na.rm = TRUE))
# caclulate mean and sd Q postharvest - by month
post.monthlyQ <- ddply(subset(discharge, discharge$logstatus =="post" ),.(format(discharge$date, format='%m')),
                       summarise, mean.Q = mean(Q.L.s, na.rm = TRUE), sd.Q = sd(Q.L.s, na.rm = TRUE))
# calculate percent difference and absolute difference - by month
per.Q <- (pre.monthlyQ[,2] - post.monthlyQ[,2])/pre.monthlyQ[,2] *100
abs.Q <- pre.monthlyQ[,2] - post.monthlyQ[,2]

## Calculate the runoff from discharge - how much water exported by month and year (mm/month, mm/year)
# Calculate specific discharge-mm/s. 1 hec = 10000 m2
discharge$Q.m3.s <- discharge$Q.L.s/1000 # convert from L/s to m3/s
discharge$Q.m.min <- discharge$Q.m3.s/area.hect * 1/10000 * 60

# calculate the mm/day - calculate time interval in min between measurements 30 min measurements
# do by day - get daily mean discharge
daily.Q <- ddply(discharge,.(format(discharge$date, format='%Y-%m-%d')),
                 summarise, meandaily.Q = mean(Q.m.min, na.rm = TRUE))
daily.Q$Q.mm.day <- daily.Q$meandaily.Q * (60*24) *1000
colnames(daily.Q)[1:2] <- c("date", "meandailyQm.min")
daily.Q$date <- as.POSIXct(strptime(daily.Q$date, format = "%Y-%m-%d"))

#  interpolate disharge for missing data - 1 day intervals. Make daily timeseries, and merge
#daily.Q.zoo <- zoo(as.POSIXct(unique(daily.Q$date)), unique(daily.Q$Q.m.day))
date <- as.data.frame(as.POSIXct(strptime(seq(min(as.Date(daily.Q$date)), max(as.Date(daily.Q$date)), by = "day"), 
                                          format = "%Y-%m-%d")))
colnames(date)[1] <- "date"
daily.Q.all <- merge(daily.Q, date, by = "date", all = TRUE) # merge together

#  interpolate disharge for missing data - 1 day intervals
daily.Q.all$splineQ.mm.day <- na.spline((daily.Q.all$Q.mm.day)) #interpolation by spline (polynmial)
daily.Q.all$approxQ.mm.day <- na.approx((daily.Q.all$Q.mm.day)) #interpolation by linear approximation

plot(daily.Q.all$splineQ.mm.day) # check both interpolations
plot(daily.Q.all$approxQ.mm.day) # check both interpolations

# calculate mm/month by summing the Q (mm/day) by month
monthly.Q <- ddply(daily.Q.all,.(format(daily.Q.all$date, format='%Y-%m')),
                   summarise, summonthly.Q.mmmonth = sum(approxQ.mm.day, na.rm = TRUE))
monthly.Q$date <- as.POSIXct(strptime(paste(monthly.Q[,1], "-01", sep = ""), format = "%Y-%m-%d"))
# Take the mean Q (m/month) in the pre and post period
# add column for pre/post condition
monthly.Q <- logstatus.f(monthly.Q)

# do monthly mean by pre and most
# calculate mean and sd Q preharvest - by month
pre.monthlyQ.mmmonth <- ddply(subset(monthly.Q, monthly.Q$logstatus =="pre" ),.(format(monthly.Q$date, format='%m')),
                              summarise, mean.Q.pre = mean(summonthly.Q.mmmonth, na.rm = TRUE), 
                              sd.Q.pre = sd(summonthly.Q.mmmonth, na.rm = TRUE))

# caclulate mean and sd Q postharvest - by month
post.monthlyQ.mmmonth <- ddply(subset(monthly.Q, monthly.Q$logstatus =="post" ),.(format(monthly.Q$date, format='%m')),
                               summarise, mean.Q.post = mean(summonthly.Q.mmmonth, na.rm = TRUE), 
                               sd.Q.post = sd(summonthly.Q.mmmonth, na.rm = TRUE))
# calculate percent difference and absolute difference - by month
per.Q.mmday <- (pre.monthlyQ.mmmonth[,2] - post.monthlyQ.mmmonth[,2])/pre.monthlyQ.mmmonth[,2] *100
abs.Q.mmday <- pre.monthlyQ.mmmonth[,2] - post.monthlyQ.mmmonth[,2]
Qmday <- cbind(pre.monthlyQ.mmmonth, post.monthlyQ.mmmonth,per.Q.mmday, abs.Q.mmday)
write.csv(Qmday, file = paste0(fig.dir, "/Qmdayper.csv")) # write file

# calculate the percentage missing values

# calculate the yearly runoff - 2008-2013
yearly.runoff <- ddply(monthly.Q,.(format(monthly.Q$date, format='%Y')),
                       summarise, sum.Q.yearly = sum(summonthly.Q.mmmonth, na.rm = TRUE))
yearly.precip <- ddply(climate,.(format(climate$date, format='%Y')),
                       summarise, sum.P.yearly = sum(as.numeric(Precip), na.rm = TRUE))

# calculate the mean percent that each month contributes to the yearly precip
all.monthly.Q<- ddply(monthly.Q,.(format(monthly.Q$date, format='%m')),
                      summarise, mean.Q.all = mean(summonthly.Q.mmmonth, na.rm = TRUE), 
                      sd.Q.all = sd(summonthly.Q.mmmonth, na.rm = TRUE))
all.yearly.Q <- mean(yearly.runoff[2:6,2])
all.monthly.Q$percenttotal <- all.monthly.Q[,2]/all.yearly.Q*100

# calcule the mean monthly precipitation (mm/month)

# culate the mean daily T (by month)

##
# Wetting ratio Precip/Q - changes over time
# calculate the mm/day of precipitation
# 
date.seq <- as.data.frame(seq(as.POSIXct("2009-01-01 00:00:00"), as.POSIXct("2014-12-31 16:00:00"), by = "hour"), 
                          format = "%Y-%m-%d %H:%M:%S")
colnames(date.seq) <- "date"
# calculate the hourly precipitation
hour.precip <- ddply(climate,.(format(climate$date, format='%Y-%m-%d %H')),
                     summarise, Precip.hourly = mean(as.numeric(Precip), na.rm = TRUE))
hour.precip$date <- as.POSIXct(strptime(paste(hour.precip[,1], ":00", sep = ""), format = "%Y-%m-%d %H:%M"))
# Merge the mean hourly with the date 
all.precip <- merge(hour.precip, date.seq, by = 'date', all = TRUE)
# interpolate
#  interpolate disharge for missing data - 1 day intervals
all.precip$splineP.mm <- na.spline(all.precip$Precip.hourly) #interpolation by spline (polynmial)
all.precip$approxP.mm <- na.approx(all.precip$Precip.hourly) #interpolation by linear approximation

# sum per month
daily.P <- ddply(all.precip,.(format(all.precip$date, format='%Y-%m-%d')),
                 summarise, sumdaily.P = sum(Precip.hourly, na.rm = TRUE))
monthly.P <- ddply(all.precip,.(format(all.precip$date, format='%Y-%m')),
                   summarise, summonthly.P = sum(Precip.hourly, na.rm = TRUE))
colnames(monthly.P) <- c('date', "monthlyP.mm")
monthly.P$date <- as.POSIXct(strptime(paste(monthly.P$date, "-01", sep = ""), format = "%Y-%m-%d"))

# merge to monthly Q (mm/month)
monthly.wetting <- merge(monthly.P, monthly.Q, by = 'date', all = FALSE)
monthly.wetting$wetratio <- monthly.wetting[,4]/monthly.wetting[,2] #calculate Q/P (mm/month)
monthly.wetting$date <- as.POSIXct(strptime(paste(monthly.wetting$date, "-01", sep = ""), format = "%Y-%m-%d"))

# calculate mean and sd Q preharvest - by month
pre.wetratio <- ddply(subset(monthly.wetting, monthly.wetting$logstatus =="pre" ),.(format(as.Date(monthly.wetting$date), format='%m')),
                      summarise, mean.ratio.pre = mean(wetratio, na.rm = TRUE), sd.ratio.pre = sd(wetratio, na.rm = TRUE))

# caclulate mean and sd Q postharvest - by month
post.wetratio <- ddply(subset(monthly.wetting, monthly.wetting$logstatus =="post" ),.(format(as.Date(monthly.wetting$date), format='%m')),
                       summarise, meam.monwet = mean(wetratio, na.rm = TRUE), sd.monwet = sd(wetratio, na.rm = TRUE))

# calculate percent difference and absolute difference - by month
per.p.mmday <- (pre.wetratio[,2] - post.wetratio[,2])/pre.wetratio[,2] *100
abs.p.mmday <- pre.wetratio[,2] - post.wetratio[,2]
pmday <- cbind(pre.wetratio, post.wetratio,per.p.mmday, abs.p.mmday)
write.csv(pmday, file = paste0(fig.dir, "/wetratio.csv")) # write file

# Figure - mean monthly discharge (mm/month) timeseries with precipitation
Q.P.merged <- merge(monthly.P, monthly.Q, by = 'date', all = TRUE)

mmrunoff.plot <- ggplot(data=Q.P.merged, aes(x=date, y=summonthly.Q.mmmonth, group=1)) +
  geom_line() +
  geom_point() +
  #scale_x_date(breaks = "1 year", minor_breaks = "1 month") +
  expand_limits(y=0) +
  xlab("Date") + ylab("Monthly Runoff (mm/month)") +
  ggtitle("Monthly Mean Runoff") + 
  theme() 

mmp.month.plot <- ggplot(data=Q.P.merged, aes(x=date, y=monthlyP.mm)) +
  geom_bar(stat="identity") + 
  xlab("Date") + ylab("Monthly Precipitation (mm/month)") +
  ggtitle("Monthly  Precipitation") 
# save all figures - precip, daily mean discharge, and runoff
pdf(file=paste0(fig.dir,"/CRFigures_runoff.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(mmp.month.plot, boxplot.discharge, mmrunoff.plot, ncol = 1)
grid.text("A", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("B", x=unit(0, "npc")+ unit(2,"mm"), y=unit(2/3, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("C", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1/3, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
dev.off()

########## Do figure specifically for the section about change in discharge upon logging
# precip, Q (L/s), mean daily Q(L/s), runoff (mm/month)
# Do all by day? Better to compare?
mmp.month.plot <- ggplot(data=Q.P.merged, aes(x=date, y=monthlyP.mm, fill = factor(logstatus))) +
  geom_bar(stat="identity", colour = "black") + 
  scale_fill_manual(values = c(cbPalette[1], cbPalette[2]),
                    name="Log Status",
                    breaks=c("pre", "post"),
                    labels=c("Pre", "Post")) + # change colours to colour blind
  xlab("Date") + ylab("Monthly Precipitation (mm/month)") +
  ggtitle("Monthly  Precipitation") +
  theme(legend.position="top") #+ 
  #scale_x_date(limits=c(min = as.Date(Q.P.merged$date[1]), max=as.Date("2014-12-01")))

time.disc
boxplot.discharge <- boxplot.discharge + theme(legend.position="none")

# Monthly runoff barplot
mmrunoff.plot <- ggplot(data=Q.P.merged, aes(x=date, y=summonthly.Q.mmmonth, fill = factor(logstatus))) +
  geom_bar(stat="identity", colour = "black") + #, fill = c(cbPalette[1], cbPalette[2])) + #colour doesn't work!!
  #scale_x_date(breaks = "1 year", minor_breaks = "1 month") +
  scale_fill_manual(values = c(cbPalette[1], cbPalette[2])) +
  expand_limits(y=0) +
  xlab("Date") + ylab("Monthly Runoff (mm/month)") +
  ggtitle("Monthly Runoff") + 
  theme(legend.position="none")

# save figures - precip, daily mean discharge, and runoff - FIGURE 2*******
pdf(file=paste0(fig.dir,"/CRFigures_Fig2ChangQ.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(mmp.month.plot, time.disc, boxplot.discharge, mmrunoff.plot, ncol = 1)
grid.text("A", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("B", x=unit(0, "npc")+ unit(2,"mm"), y=unit(3/4, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("C", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1/2, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("D", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1/4, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
dev.off()

pdf(file=paste0(fig.dir,"/CRFigures_Fig2ChangQ_nobox.pdf"), width = 8.5, height = 11) #save figure without boxplot
grid.arrange(mmp.month.plot, time.disc, mmrunoff.plot, ncol = 1)
grid.text("A", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("B", x=unit(0, "npc")+ unit(2,"mm"), y=unit(2/3, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("C", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1/3, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
dev.off()

#####################################################################################
## Effect of harvest on discharge - cumulative distribution; and flood frequency analysis
# Calculate specific discharge - Q/watershed area 
# 8.64 is the conversion factor from L s 1 ha 1 to mm d 1.From: Sørensen, Rasmus, Eva Ring, Markus Meili, Lars Högbom, Jan Seibert, Thomas Grabs, Hjalmar Laudon, and Kevin Bishop. 2009. “Forest Harvest Increases Runoff Most During Low Flows in Two Boreal Streams.” AMBIO: a Journal of the Human Environment 38 (7): 357–63. doi:10.1579/0044-7447-38.7.357.
#discharge$Q.mm.day <- discharge$Q.L.s/area.hect * 8.64

## Cumulative distribution function of discharge - by year
CDF.runoff.year <- ggplot(discharge, aes(x = Q.mm.d)) + 
  stat_ecdf(aes(group = format(discharge$date, format='%Y'), colour = format(discharge$date, format='%Y'))) +
  scale_fill_manual(values=c(cbPalette[1:9]),
                    name="Year") +
  #theme(legend.title=element_text("Year")) +
  scale_color_manual(breaks=c("2007","2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015"),
                     values=c(cbPalette[1:9])) +
  theme(legend.position=c(.8, .3)) +
  ggtitle("Discharge by Year") +
  scale_y_continuous(name="Cumulative Distribution") +
  scale_x_continuous(name="Q (mm/day)", limits = c(0, 40))

## Cumulative distribution function of discharge - by year
CDF.runoff.log <- ggplot(discharge, aes(x = Q.mm.d, color = logstatus)) + 
  stat_ecdf(aes(group = logstatus)) +
  scale_fill_manual(values=c(cbPalette[1:2]),
                    name="Logging\nStatus",
                    labels=c("Pre-Harvest", "Post-Harvest")) +
  #theme(legend.title=element_text("Year")) +
  scale_color_manual(breaks=c("post","pre"),
                     values=c(cbPalette[1],cbPalette[2])) +
  theme(legend.position=c(.8, .3)) +
  ggtitle("Discharge by Logging Status") +
  scale_y_continuous(name="Cumulative Distribution") +
  scale_x_continuous(name="Q (mm/day)", limits = c(0, 40))

## Frequency Distribution - by year
frequ.year <- ggplot(discharge, aes(Q.mm.d, fill = format(discharge$date, format='%Y'))) + 
  geom_density(alpha = 0.2) +
  scale_fill_manual(values=c(cbPalette[1:9]),
                    name="Year") + 
  theme(legend.position=c(.8, .6)) +
  ggtitle("Discharge by Year") +
  scale_y_continuous(name="Density") +
  scale_x_continuous(name="Q (mm/day)", limits = c(0, 40))

# Frequency distribution - preharvest verusus post
frequ.log <- ggplot(discharge, aes(Q.mm.d, fill = logstatus)) + geom_density(alpha = 0.2) +
  scale_fill_manual(values=c(cbPalette[1:2]),
                    name="Logging\nStatus",
                    labels=c("Post-Harvest", "Pre-Harvest")) + 
  theme(legend.position=c(.8, .8)) +
  ggtitle("Discharge by Pre-Post Status") +
  scale_y_continuous(name="Density") +
  scale_x_continuous(name="Q (mm/day)", limits = c(0, 40)) + 
  theme()

# save distribution functions - all four in one document
pdf(file=paste0(fig.dir,"/CRFigures_distributionrunoff.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(CDF.runoff.year, CDF.runoff.log, frequ.year, frequ.log, ncol=2)
dev.off()  
##########
## Flood frequency Analysis
# http://www.headwateranalytics.com/blog/flood-frequency-analysis-in-r
# https://github.com/eheisman/hydroutils
# Make a pre and post harvest flodd exheedance curve
# Also, one plot with both with the pre and post discharge in different colours?

# extract instantaneous maxima flow series
#input_data.pre <- subset(discharge$Q.L.s, discharge$logstatus == "pre")

# choose a distribution 
# (see ?dist.list for available distributions)
# log Pearson 3 is not one of them but is a log-transformed equivalent of pe3
# note: this script recognizes 'lp3' to stand for log Pearson Type 3
#dist <- "gev" doesn't work for our data! different distribution
#dist <- 'lp3'

# fit frequency distribution
#source("/Users/user/SpecScripts/frequ_function.R")
#fa.pre <- FrequencyAnalysis(series=input_data.pre, distribution=dist)

# estimate 90% confidence intervals
#source("/Users/user/SpecScripts/bootstrap_function.R")
#ci.pre <- BootstrapCI(series=input_data.pre,   # flow data
#                  distribution=dist,   # distribution
#                  n.resamples = 2E3, # number of re-samples to conduct
#                  ci = 0.90)           # confidence interval level
# generate frequency plot
#source("/Users/user/SpecScripts/frequplot_function.R")
#frequ.plot.pre <- frequencyPlot(series=input_data.pre, ci.pre$ci)

# postharvest
#input_data.post <- subset(discharge$Q.L.s, discharge$logstatus == "post")
# choose a distribution 
# (see ?dist.list for available distributions)
# log Pearson 3 is not one of them but is a log-transformed equivalent of pe3
# note: this script recognizes 'lp3' to stand for log Pearson Type 3
#dist <- "gev" doesn't work for our data! different distribution
#dist <- 'lp3'

# fit frequency distribution
#source("/Users/user/SpecScripts/frequ_function.R")
#fa.post <- FrequencyAnalysis(series=input_data.post, distribution=dist)

# estimate 90% confidence intervals
#source("/Users/user/SpecScripts/bootstrap_function.R")
#ci.post <- BootstrapCI(series=input_data.post,   # flow data
 #                     distribution=dist,   # distribution
#                      n.resamples = 2E3, # number of re-samples to conduct
#                      ci = 0.90)           # confidence interval level
# generate frequency plot
#source("/Users/user/SpecScripts/frequplot_function.R")
#frequ.plot.post <- frequencyPlot(series=input_data.post, ci.post$ci)

#save the pre and post distributions
#pdf(file=paste0(fig.dir,"/CRFigures_floodfrequ.pdf"), width = 8.5, height = 11) #save figure
#grid.arrange(frequ.plot.pre, frequ.plot.post, ncol=2)
#grid.text("Pre-Harvest", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
#grid.text("Post-Harvest", x=unit(00.5, "npc"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
#dev.off()

#### Flow duration curves
# On daily mean Q(m3/s)

# Get daily mean Q (m3/s)
daily.Q$Q.m3s <- ddply(discharge,.(format(discharge$date, format='%Y-%m-%d')),
                 summarise, meandaily.Q = mean(Q.m3.s, na.rm = TRUE))[2]
# Convert to m3/day
daily.Q$Q.m3day <- daily.Q[4] *60*60*24
# For pre and post harvest period
daily.Q <- logstatus.f(daily.Q)
prelog.dailyQ <- subset(daily.Q, daily.Q$logstatus == "pre")
postlog.dailQ <- subset(daily.Q, daily.Q$logstatus == "post")

# Flow duration curves by pre/post
pre.post <- cbind.fill(prelog.dailyQ[,5], postlog.dailQ[,5], fill = NaN)
colnames(pre.post) <- c("Qm3day.pre", "Qm3day.post")
fdc.prepost <- fdc(pre.post, lQ.thr=0.7, hQ.thr=0.2, plot=TRUE, log="y",
    main= "Flow Duration Curve", xlab="% Time flow equalled or exceeded",
    ylab="Q, [m3/day]") 

# FDC by year
y2008 <- subset(daily.Q, format(daily.Q$date, format = "%Y") == "2008")
y2009 <- subset(daily.Q, format(daily.Q$date, format = "%Y") == "2009")
y2010 <- subset(daily.Q, format(daily.Q$date, format = "%Y") == "2010")
y2011 <- subset(daily.Q, format(daily.Q$date, format = "%Y") == "2011")
y2012 <- subset(daily.Q, format(daily.Q$date, format = "%Y") == "2012")
y2013 <- subset(daily.Q, format(daily.Q$date, format = "%Y") == "2013")

# bind as one
FDC.year <- cbind.fill(y2008[,5], y2009[,5], y2010[,5], y2011[,5], y2012[,5], y2013[,5], fill = NaN)
colnames(FDC.year) <- c("2008", "2009", "2010", "2011", "2012", "2013")
fdc.year <- fdc(FDC.year, lQ.thr=0.7, hQ.thr=0.2, plot=TRUE, log="y",
                   main= "Flow Duration Curve", xlab="% Time flow equalled or exceeded",
                   ylab="Q, [m3/day]") 

# wet/dry for pre and post
daily.Q <- wetdry.f(daily.Q) # add in column for wet/dry period
pre.wet <- subset(daily.Q, daily.Q$logstatus == "pre" & daily.Q$hydro == "wet")
pre.dry <- subset(daily.Q, daily.Q$logstatus == "pre" & daily.Q$hydro == "dry")
post.wet <- subset(daily.Q, daily.Q$logstatus == "post" & daily.Q$hydro == "wet")
post.dry <- subset(daily.Q, daily.Q$logstatus == "post" & daily.Q$hydro == "dry")
# bind
hydro.prepost <- cbind.fill(pre.wet[,5],pre.dry[,5], post.wet[,5], post.dry[,5], fill = NaN)
colnames(hydro.prepost) <- c("pre/wet", "pre/dry", "post/wet", "post/dry")
fdc.hydroprepost <- fdc(hydro.prepost, lQ.thr=0.7, hQ.thr=0.2, plot=TRUE, log="y",
                main= "Flow Duration Curve", xlab="% Time flow equalled or exceeded",
                ylab="Q, [m3/day]") 

# Save all as pdf
pdf(file=paste0(fig.dir,"/CRFigures_FDC.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(fdc(pre.post, lQ.thr=0.7, hQ.thr=0.2, plot=TRUE, log="y",
                                main= "Flow Duration Curve", xlab="% Time flow equalled or exceeded",
                                ylab="Q, [m3/day]"), 
             fdc(FDC.year, lQ.thr=0.7, hQ.thr=0.2, plot=TRUE, log="y",
                                          main= "Flow Duration Curve", xlab="% Time flow equalled or exceeded",
                                          ylab="Q, [m3/day]"), 
             fdc(hydro.prepost, lQ.thr=0.7, hQ.thr=0.2, plot=TRUE, log="y",
                                                            main= "Flow Duration Curve", xlab="% Time flow equalled or exceeded",
                 col=cbPalette[4:7],ylab=expression(Q~(m^{3}~day^{-1})), cex.lab=2), 
             ncol = 1)
grid.text("A", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("B", x=unit(0, "npc")+ unit(2,"mm"), y=unit(2/3, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("C", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1/3, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
dev.off()

pdf(file=paste0(fig.dir,"/CRFigures_FDC_preposthydro.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(fdc(hydro.prepost, lQ.thr=0.7, hQ.thr=0.2, plot=TRUE, log="y",
    main= "Flow Duration Curve", xlab="% Time flow equalled or exceeded",
    col=cbPalette[4:7], ylab=expression(Q~(m^{3}~day^{-1})), cex.lab=2,new = FALSE, leg.pos="bottomleft"))
dev.off()

#############################################################################

#############################################################################
## How does harvest alter DOC concentration?
## Figure *** -P, Q, DOC, DOC box plot

time.disc <- ggplot(discharge, aes(date, Q.m3.s, colour = logstatus)) + 
  geom_point(size = 0.4) +
  labs(x = expression("Date"),  
       y = expression(Q~(m^{3}~s^{-1}))) +
  scale_color_manual(values=cbPalette[1:2],
                     name="Log Status",
                     breaks=c("pre", "post"),
                     labels=c("Pre", "Post")) + # colours
  theme(legend.position="top")


# do daily sum of precip for bar graph. Too complicated otherwise
dailysum.precip<- ddply(climate,.(format(climate$date, format = "%Y-%m-%d")),
                        summarise, sum.P = sum(as.numeric(Precip), na.rm = TRUE))
colnames(dailysum.precip) <- c("date", "dailyP.mm")
dailysum.precip$date <- as.POSIXct(strptime(dailysum.precip$date, format = "%Y-%m-%d"))
dailysum.precip <- wetdry.f(dailysum.precip)

# do plots - timeseries
time.precip <- ggplot(dailysum.precip, aes(as.Date(date), dailyP.mm, fill = factor(hydro))) + 
  geom_bar(stat="identity") +
  xlab("Date") + ylab("Daily Precipitation (mm)") + 
  theme(legend.position="top") +
  scale_x_date(date_breaks = "year", 
               labels=date_format("%Y"))

time.DOC <- ggplot(subset(spectro.all, spectro.all$DOCcorr >=1), aes(date, DOCcorr, colour = logstatus)) + 
  geom_point(size = 0.4) +
  xlab("Date") + ylab("[DOC] (mg/L)") + 
  scale_color_manual(values=cbPalette[1:2],
                     name="Log Status",
                     breaks=c("pre", "post"),
                     labels=c("Pre", "Post")) + # colours
  theme(legend.position="top")

# boxplot of 30 minute DOC measurements, compiled by month. Deliniate the pre and post harvest period
DOC.box <- subset(spectro.all, spectro.all$DOCcorr >=1)
boxplot.DOC30 <- ggplot(DOC.box, aes(x=date, y=DOCcorr, group = format(DOC.box$date, format="%Y-%m"), fill=logstatus)) + 
  geom_boxplot(outlier.shape = NA) + # don't show outliers
  scale_fill_manual(breaks = c('pre', 'post'),
                    values = c(cbPalette[1], cbPalette[2]),
                    name="Logging\nStatus") +
  xlab("Date") + 
  ylab("[DOC] (mg/L)") + 
  #geom_line() + # line to show where logging occureed
  ggtitle("DOC Measurements") +
  theme() +
  theme(legend.position=c(0.93, 0.87))

### Save figure
pdf(file=paste0(fig.dir,"/CRFigures_DOCtimeseries.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(time.disc,(time.DOC+theme(legend.position= "NONE")), (boxplot.DOC30+theme(legend.position= "NONE")),  ncol = 1)
grid.text("A", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("B", x=unit(0, "npc")+ unit(2,"mm"), y=unit(2/3, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("C", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1/3, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
#grid.text("D", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1/4, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
#excude precip for now
dev.off()

#########
# Cumulative distribution function - DOC concentrations in the pre/post wet dry period (as per discharge)
## Cumulative distribution function of DOC - partitioned according to pre/post, wet/dry (30 minute measurements)
# create partitions in 30 minute data
DOC.partition <- logstatus.f(spectro.all)
DOC.partition <- wetdry.f(DOC.partition)
DOC.partition$hydro.log <- paste(DOC.partition$hydro, DOC.partition$logstatus, sep = "")

CDF.DOC.part <- ggplot(DOC.partition, aes(x = DOCcorr)) + 
  stat_ecdf(aes(group = hydro.log, colour = hydro.log)) +
  scale_fill_manual(values=c(cbPalette[1:4]),
                    name="Status") +
  #theme(legend.title=element_text("Year")) +
  scale_color_manual(breaks=c("drypost","drypre", "wetpost", "wetpre"),
                     values=c(cbPalette[1:4]), 
                     name = "Status") +
  theme(legend.position=c(0.2, 0.8)) +
  ggtitle("") +
  scale_y_continuous(name="Cumulative Distribution") +
  scale_x_continuous(name="[DOC] (mg/L)") + 
  coord_flip() + #flip x and y axis + 
  theme()

## save both together as figure
pdf(file=paste0(fig.dir,"/CRFigures_deltaDOCconc.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(boxplot.DOC30, CDF.DOC.part,  ncol = 1)
grid.text("A", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("B", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1/2, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
dev.off()
############################
# Test whether mean DOC is different in the pre/post period by month
#"nonparametric Kruskal-Wallis one way analysis of variance on ranks (ANOVA-R), which was combined with Dunn’s test to investigate significant differences in median values."
# From Schelker, J, K Eklöf, and K Bishop. 2012. “Effects of Forestry Operations on Dissolved Organic Carbon Concentrations and Export in Boreal First‐Order Streams.” Journal of Geophysical …. doi:10.1029/2011JG001827.
# compare the monthly mean DOC in the post period to that of the pre harvest period

# First, test whether data is normal or not.
## Have a look at the densities
plot(density(spectro.all$DOCcorr, na.rm = TRUE))

## Perform the test
shapiro.test(spectro.all$DOCcorr)

## Plot using a qqplot
qqnorm(spectro.all$DOCcorr);qqline(spectro.all$DOCcorr, col = 2)
# non-normal!

## nonparametric Kruskal-Wallis one way analysis of variance on ranks (ANOVA-R), which was combined with Dunn’s test to investigate significant differences in median values. -- Highlighted Jul 14, 2016
# http://www.r-bloggers.com/kruskal-wallis-one-way-analysis-of-variance-2/
# http://rcompanion.org/rcompanion/d_06.html

# first compile the months together (all jan, all feb, all march, etc.) into 
jan.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "01")
feb.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "02")
march.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "03")
april.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "04")
may.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "05")
june.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "06")
july.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "07")
aug.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "08")
sept.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "09")
oct.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "10")
nov.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "11")
dec.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "12")

# 
years.DOC <- function(data){
  test <- data[,c(as.numeric(match("date",names(data))), as.numeric(match("DOCcorr",names(data))))]
  sub <- cbind.fill(subset(test, format(test$date, format='%Y') == "2009"), 
                            subset(test, format(test$date, format='%Y') == "2010"),
                            subset(test, format(test$date, format='%Y') == "2011"),
                            subset(test, format(test$date, format='%Y') == "2012"),
                            subset(test, format(test$date, format='%Y') == "2013"),
                            subset(test, format(test$date, format='%Y') == "2014"), fill = NaN)
  return(sub)
}

jan.DOC.sub <- years.DOC(data = jan.DOC)
feb.DOC.sub <- years.DOC(feb.DOC)
march.DOC.sub <- years.DOC(march.DOC)
april.DOC.sub <- years.DOC(april.DOC)
may.DOC.sub <- years.DOC(may.DOC)
june.DOC.sub <- years.DOC(june.DOC)
july.DOC.sub <- years.DOC(july.DOC)
aug.DOC.sub <- years.DOC(aug.DOC)
sept.DOC.sub <- years.DOC(sept.DOC)
oct.DOC.sub <- years.DOC(oct.DOC)
nov.DOC.sub <- years.DOC(nov.DOC)
dec.DOC.sub <- years.DOC(dec.DOC)

#
jan.DOC$year <- format(jan.DOC$date, format='%Y')
jan.DOCselect <- jan.DOC[,c(1, as.numeric(match("DOCcorr",names(jan.DOC))), 
                            as.numeric(match("year",names(jan.DOC))))]
kruskal.test(jan.DOC.sub[,c(8,10,12)]) 
library(dunn.test)
dunn.test(jan.DOC.sub[,c(6,8,10,12)], g=DOCcorr, kw=TRUE)
kruskal.test(feb.DOC.sub[,c(6,8,10,12)]) 
kruskal.test(jan.DOC.sub[,c(4,6,8,10,12)]) 
kruskalmc(jan.DOC.sub$date, jan.DOC.sub[,c(4,6,8,10,12)]) 

#library(FSA)
#jan.DOC.test <- dunnTest(DOCcorr ~ year,
#         data=jan.DOCselect,
#         method="none")    # Can adjust p-values; 
# See ?p.adjust for options

#library(DescTools)
#DunnTest(x = jan.DOCselect$year,
#         g = jan.DOCselect$DOCcorr,
#         method="none")

# really slow! Try on daily means rather than 30 minute data
#First, aggregate to the daily mean, then look at whether daily means are different
DOC.daily <- ddply(subset(spectro.all, spectro.all$DOCcorr >=1),.(format(spectro.all$date, format='%Y-%m-%d')),
                     summarise, meandaily.DOC = mean(DOCcorr, na.rm = TRUE)) 
colnames(DOC.daily)[1] <- "date"
DOC.daily$date <- as.POSIXct(strptime(DOC.daily$date, format = "%Y-%m-%d"))
jan.DOC <- subset(DOC.daily, format(DOC.daily$date, format='%m') == "01")
feb.DOC <- subset(DOC.daily, format(DOC.daily$date, format='%m') == "02")
march.DOC <- subset(DOC.daily, format(DOC.daily$date, format='%m') == "03")
april.DOC <- subset(DOC.daily, format(DOC.daily$date, format='%m') == "04")
may.DOC <- subset(DOC.daily, format(DOC.daily$date, format='%m') == "05")
june.DOC <- subset(DOC.daily, format(DOC.daily$date, format='%m') == "06")
july.DOC <- subset(DOC.daily, format(DOC.daily$date, format='%m') == "07")
aug.DOC <- subset(DOC.daily, format(DOC.daily$date, format='%m') == "08")
sept.DOC <- subset(DOC.daily, format(DOC.daily$date, format='%m') == "09")
oct.DOC <- subset(DOC.daily, format(DOC.daily$date, format='%m') == "10")
nov.DOC <- subset(DOC.daily, format(DOC.daily$date, format='%m') == "11")
dec.DOC <- subset(DOC.daily, format(DOC.daily$date, format='%m') == "12")

# compile nov and december in the pre harvest period
nov.DOC$date <- as.Date(gsub("2009", "2010", nov.DOC$date), format = "%Y-%m-%d")
dec.DOC$date <- as.Date(gsub("2009", "2010", dec.DOC$date), format = "%Y-%m-%d")

# Function for doing dunn test, specifying method
Dunns.test <- function(data, method){
  data$year <- format(data$date, format='%Y')
  Dunn <- dunnTest(data$meandaily.DOC ~ as.factor(data$year),
                           data=data,
                           method=method)
  return(Dunn)
}

jan.DOC.test <- Dunns.test(jan.DOC, "by")
feb.DOC.test <- Dunns.test(feb.DOC, "by")
march.DOC.test <- Dunns.test(march.DOC, "by")
april.DOC.test <- Dunns.test(april.DOC, "by")
may.DOC.test <- Dunns.test(may.DOC, "by")
june.DOC.test <- Dunns.test(june.DOC, "by")
july.DOC.test <- Dunns.test(july.DOC, "by")
aug.DOC.test <- Dunns.test(aug.DOC, "by")
sept.DOC.test <- Dunns.test(sept.DOC, "by")
oct.DOC.test <- Dunns.test(oct.DOC, "by")
nov.DOC.test <- Dunns.test(nov.DOC, "by")
dec.DOC.test <- Dunns.test(dec.DOC, "by")

# take the mean of the pre-logging period by month to compare to months in the post-logging period
DOC.monthly <- ddply(spectro.all,.(format(spectro.all$date, format='%Y-%m')),
                   summarise, meanmonthly.DOC = mean(DOCcorr, na.rm = TRUE),
                   sdmonthly.DOC = sd(DOCcorr, na.rm = TRUE)) 
write.csv(DOC.monthly, file=paste0(fig.dir,"/CRmonthlymeanDOCconc.csv")) # write mean DOC

# plot the mean daily DOC as a boxplot
DOC.daily <- logstatus.f(DOC.daily)

boxplot.DOCdaily <- ggplot(DOC.daily, aes(x=date, y=meandaily.DOC, group = format(DOC.daily$date, format="%Y-%m"), fill=logstatus)) + 
  geom_boxplot(outlier.shape = NA) + # don't show outliers
  scale_fill_manual(breaks = c('pre', 'post'),
                    values = c(cbPalette[1], cbPalette[2]),
                    name="Logging\nStatus") +
  xlab("Date") + 
  ylab("Daily Mean DOC (mg/L)") + 
  #geom_line() + # line to show where logging occureed
  ggtitle("DOC Measurements") +
  theme() +
  theme(legend.position=c(0.93, 0.87))

### Save figure
pdf(file=paste0(fig.dir,"/CRFigures_DOCtimeseries_daily.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(time.disc,time.DOC+theme(legend.position= "NONE"), boxplot.DOCdaily+theme(legend.position= "NONE"),  ncol = 1)
grid.text("A", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("B", x=unit(0, "npc")+ unit(2,"mm"), y=unit(2/3, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("C", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1/3, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
#grid.text("D", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1/4, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
#excude precip for now
dev.off()

######## REDO THIS!!! Use this instead to look at significance of changes in mean monthly DOC.
#https://cran.r-project.org/web/packages/PMCMR/vignettes/PMCMR.pdf

######## Do table of mean DOC by month within the pre/post period
# Calculate mean DOC for pe/post by month
spectro.all$month <- format(spectro.all$date,"%m")
stats.logging.DOCconc <- ddply(spectro.all, c("month", "logstatus"), summarise,
                       mean.DOC = mean(DOCcorr, na.rm = TRUE), sd.DOC=sd(DOCcorr, na.rm = TRUE), 
                       max.DOC = max(DOCcorr, na.rm=TRUE), min.DOC = min(DOCcorr, na.rm=TRUE))
write.csv(stats.logging.DOCconc, file = paste0(fig.dir, "/DOCconc_prepoststats.csv"))
# do test for significance
jan.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "01")
jan.kwt <- kruskal.test(DOCcorr ~ as.factor(logstatus), data = jan.DOC[!is.na(jan.DOC$logstatus),])
feb.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "02")
feb.kwt <- kruskal.test(DOCcorr ~ as.factor(logstatus), data = feb.DOC[!is.na(feb.DOC$logstatus),])
march.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "03")
march.kwt <- kruskal.test(DOCcorr ~ as.factor(logstatus), data = march.DOC[!is.na(march.DOC$logstatus),])
april.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "04")
april.kwt <- kruskal.test(DOCcorr ~ as.factor(logstatus), data = april.DOC[!is.na(april.DOC$logstatus),])
may.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "05")
may.kwt <- kruskal.test(DOCcorr ~ as.factor(logstatus), data = may.DOC[!is.na(may.DOC$logstatus),])
june.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "06")
june.kwt <- kruskal.test(DOCcorr ~ as.factor(logstatus), data = june.DOC[!is.na(june.DOC$logstatus),])
july.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "07")
july.kwt <- kruskal.test(DOCcorr ~ as.factor(logstatus), data = july.DOC[!is.na(july.DOC$logstatus),])
aug.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "08")
aug.kwt <- kruskal.test(DOCcorr ~ as.factor(logstatus), data = aug.DOC[!is.na(aug.DOC$logstatus),])
sept.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "09")
sept.kwt <- kruskal.test(DOCcorr ~ as.factor(logstatus), data = sept.DOC[!is.na(sept.DOC$logstatus),])
oct.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "10")
oct.kwt <- kruskal.test(DOCcorr ~ as.factor(logstatus), data = oct.DOC[!is.na(oct.DOC$logstatus),])
nov.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "11")
nov.kwt <- kruskal.test(DOCcorr ~ as.factor(logstatus), data = nov.DOC[!is.na(nov.DOC$logstatus),])
dec.DOC <- subset(spectro.all, format(spectro.all$date, format='%m') == "12")
dec.kwt <- kruskal.test(DOCcorr ~ as.factor(logstatus), data = dec.DOC[!is.na(dec.DOC$logstatus),])
DOCconc.kwt <- rbind(jan.kwt[1:3], feb.kwt[1:3], march.kwt[1:3], april.kwt[1:3], may.kwt[1:3], june.kwt[1:3], 
                     july.kwt[1:3], aug.kwt[1:3], sept.kwt[1:3], nov.kwt[1:3], dec.kwt[1:3]) 
write.csv(DOCconc.kwt, file = paste0(fig.dir, "/DOCconc_KWT.csv"))
# calculate the overall difference/pre post

# calculate the percent change - pre and post by month
stats.DOCconc.all <- ddply(spectro.all, c("logstatus"), summarise,
                               mean.DOC = mean(DOCcorr, na.rm = TRUE), sd.DOC=sd(DOCcorr, na.rm = TRUE), 
                               max.DOC = max(DOCcorr, na.rm=TRUE), min.DOC = min(DOCcorr, na.rm=TRUE))
###################################################################################################################
# How does harvest affect DOC flux?
## Calculation of DOC Flux
# Method 1: daily average C = daily average DOC x average Q x 24 hours
#calculation of DOC flux

## calculate carbon flux - every 30 minutes
# merge Q (L/s) into the spectro.all to get a discharge column in this dataframe
# First, make sure discharge is in same timestamp as DOC using openair

DOC.Q <- merge(spectro.all, discharge, by = 'date', all = FALSE)
DOC.Q$cflux.mgs = (DOC.Q$DOCcorr * DOC.Q$Q.L.s )
DOC.Q$cflux.mgshec <- DOC.Q$cflux.mgs / area.hect # get the cflux normalized for watershed area (in hectares)
DOC.Q$cflux.gshec <- DOC.Q$cflux.mgshec/1000
c.flux= DOC.Q[complete.cases(DOC.Q[,as.numeric(match("DOCcorr",names(DOC.Q)))]),] # get the dataframe with complete cases of DOCcorr

# calculating cflux by first calculating the time intervals.. no interpolation
#n <- dim(c.flux)[1]
#for (i in 2:n) {
#  timeinterval <- c.flux$date[i]-c.flux$date[(i-1)]
#  c.flux$timeinterval <- timeinterval
#}

#for (i in 1:n){
##  cflux.mg = (c.flux$DOCcorr[i] * c.flux$Q.L.s[i] * c.flux$timeinterval[i]*60)
#  c.flux$cflux.mg[i] <- cflux.mg
#}

## Cflux - by 30 minute intervals interpolation
# interpolate for missing measurements 
date <- as.data.frame(as.POSIXct(strptime(seq(min(c.flux$date), max(c.flux$date), by = "30 mins"), 
                                          format = "%Y-%m-%d %H:%M")))
colnames(date)[1] <- "date"

cflux.all <- merge(c.flux, date, by = "date", all = TRUE) # merge together, creating NAs for missing data

#  interpolate disharge for missing data - 1 day intervals
cflux.all$splinecflux.mgs <- na.spline(cflux.all$cflux.mgs) #interpolation by spline (polynmial)
cflux.all$approxcflux.mgs <- na.approx(cflux.all$cflux.mgs) #interpolation by linear approximation

#get the percentage of missing data
permissing.30min <- length(which(is.na(cflux.all$cflux.mgs)))/length(cflux.all$cflux.mgs)*100 #30% missing 30 minute data!

## Cflux - by  hourly intervals (interpolation)
#do hourly to prevent the effect of missing data
hourly.DOC <- as.data.frame(ddply(c.flux,.(format(c.flux$date, format='%Y-%m-%d %H')),
                 summarise, hourly.DOCcorr = mean(DOCcorr, na.rm = TRUE)))
# interpolate for missing measurements 
date <- as.data.frame(as.POSIXct(strptime(seq(min(c.flux$date), max(c.flux$date), by = "hour"), 
                                          format = "%Y-%m-%d %H:%M")))
colnames(date)[1] <- "date"

cflux.all <- merge(test, date, by = "date", all = TRUE) # merge together, creating NAs for missing data

#  interpolate disharge for missing data - 1 day intervals
cflux.all$splinecflux.mgs <- na.spline(cflux.all$cflux.mgs) #interpolation by spline (polynmial)
cflux.all$approxcflux.mgs <- na.approx(cflux.all$cflux.mgs) #interpolation by linear approximation

#get the percentage of missing data
permissing.30min <- length(which(is.na(cflux.all$cflux.mgs)))/length(cflux.all$cflux.mgs)*100 #29% missing 30 minute data!

# do by day - to compensate for missing data
daily.DOC <- ddply(c.flux,.(format(c.flux$date, format='%Y-%m-%d')),
                    summarise, daily.cflux.mgs = mean(cflux.mgs, na.rm = TRUE)) # get the daily mean cflux
colnames(daily.DOC)[1] <- "date"
daily.DOC$date <- as.POSIXct(strptime(daily.DOC$date, format = "%Y-%m-%d"))
# get a timeseries of dates for interpolation
date.seq <- as.data.frame(as.POSIXct(strptime(seq(min(daily.DOC$date), max(daily.DOC$date), by = "day"), 
                                          format = "%Y-%m-%d")))
colnames(date.seq)[1] <- "date"
# merge into daily DOC to creat NAs for interpolation
cflux.day <- merge(daily.DOC, date.seq, by = "date", all = TRUE) # merge together, creating NAs for missing data
# convert mg/s to mg/day
cflux.day$daily.cflux.mgday <- cflux.day$daily.cflux.mgs*60*60*24

#  interpolate disharge for missing data - 1 day intervals
# cut missing beginning
cflux.day <- cflux.day[-c(1:22),]
cflux.day$daily.cflux.mgday.spline <- na.spline(cflux.day$daily.cflux.mgday) #interpolation by spline (polynmial)
cflux.day$daily.cflux.mgday.linear <- na.approx(cflux.day$daily.cflux.mgday) #interpolation by linear approximation

# percent of missing data - daily average DOC flux (mg/day)
permissing.day <- length(which(is.na(cflux.day$daily.cflux.mgday)))/length(cflux.day$daily.cflux.mgday)*100 #30% missing 30 minute data!

# sum linear approximation to get monthly DOC export
monthly.DOCflux <- ddply(cflux.day,.(format(cflux.day$date, format='%Y-%m')),
                                      summarise, monthlysumCs = sum(daily.cflux.mgday.linear, na.rm = TRUE)) # get the daily mean cflux
colnames(monthly.DOCflux)[1] <- 'date'

monthly.DOCflux$date <- as.POSIXct(strptime(paste(monthly.DOCflux$date, "-01", sep = ""), format = "%Y-%m-%d"))
monthly.DOCflux$cflux.gmonth = monthly.DOCflux$monthlysumCs/1000
monthly.DOCflux$cflux.gmonthhect <- monthly.DOCflux$monthlysumCs/area.hect
monthly.DOCflux$cflux.kgmonth = monthly.DOCflux$cflux.gmonth/1000
monthly.DOCflux$cflux.kgmonthhect <- monthly.DOCflux$cflux.kgmonth/area.hect

# add in pre/post harvest
monthly.DOCflux <- logstatus.f(monthly.DOCflux)
monthly.DOCflux$logstatus

# do boxplot for c-export by month. Pre and post harvest coloured.
boxplot.monthcexport <- ggplot(monthly.DOCflux, aes(as.Date(date), cflux.kgmonthhect)) + 
  geom_bar(aes(fill=logstatus),   # fill depends on cond2
           stat="identity",
           colour="black") +    # Black outline for all
  scale_fill_manual(breaks = c('pre', 'post'),
                    values = c(cbPalette[1], cbPalette[2]),
                    name="Logging\nStatus") +
  labs(x = expression("Date"),  
       y = expression(DOCFlux~(kg~month^{-1}~hec^{-1})))
  #xlab("Date") + ylab("DOC Flux (kg/monthhect)") + 
  theme() +
  scale_x_date(date_breaks = "year", 
               labels=date_format("%Y")) +
  theme(legend.position=c(0.8, 0.8))

# save monthly doc flux barplot figures
pdf(file=paste0(fig.dir,"/CRFigures_DOCboxplot_CFluxmonth.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(boxplot.monthcexport, ncol = 1)
dev.off()

## do boxplot of 30 minute DOC flux c - aggregared by month (as per DOC concentration)
DOC.Q <- logstatus.f(DOC.Q)
boxplot.DOCflux30 <- ggplot(DOC.Q, aes(x=date, y=cflux.gshec, group = format(spectro.all$date, format="%Y-%m"), fill=logstatus)) + 
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.3),  width = 0.3) + # don't show outliers
  scale_fill_manual(breaks = c('pre', 'post'),
                    values = c(cbPalette[1], cbPalette[2]),
                    name="Logging\nStatus") +
  xlab("Date") + 
  ylab("30 Minute DOC flux (g/(s hectares)") + 
  #geom_line() + # line to show where logging occureed
  ggtitle("DOC Flux (30 minute measurements") +
  theme() +
  theme(legend.position=c(0.93, 0.87))

## do cumulative distribution function - Cflux pre/post, and wet/dry
DOC.Q <- wetdry.f(DOC.Q)
DOC.Q$hydro.log <- paste(DOC.Q$hydro, DOC.Q$logstatus, sep = "")

CDF.DOC.flux <- ggplot(DOC.Q, aes(x = cflux.gshec)) + 
  stat_ecdf(aes(group = hydro.log, colour = hydro.log)) +
  scale_fill_manual(values=c(cbPalette[1:4]),
                    name="Status") +
  #theme(legend.title=element_text("Year")) +
  scale_color_manual(breaks=c("drypost","drypre", "wetpost", "wetpre"),
                     values=c(cbPalette[1:4]), 
                     name = "Status") +
  theme(legend.position=c(0.2, 0.8)) +
  ggtitle("") +
  scale_y_continuous(name="Cumulative Distribution") +
  scale_x_continuous(name="DOC flux (g/(s*hectare))") + 
  coord_flip() + #flip x and y axis + 
  theme()

# save as figure
pdf(file=paste0(fig.dir,"/CRFigures_CFlux.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(boxplot.monthcexport, boxplot.DOCflux30, CDF.DOC.flux, ncol = 1)
grid.text("A", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("B", x=unit(0, "npc")+ unit(2,"mm"), y=unit(2/3, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("C", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1/3, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
dev.off()

## save figure - DOC concentration and flux on one
maxWidth <- grid::unit.pmax(boxplot.DOC30$widths[2:5], boxplot.monthcexport$widths[2:5])
time.disc$widths[2:5] <- as.list(maxWidth)
boxplot.DOC30$widths[2:5] <- as.list(maxWidth)
time.DOC$widths[2:5] <- as.list(maxWidth)
boxplot.monthcexport$widths[2:5] <- as.list(maxWidth)
# make sure that plot y axis will line up.
p1 <- ggplot_gtable(ggplot_build(time.disc))
p2 <- ggplot_gtable(ggplot_build(time.DOC+theme(legend.position= "NONE")))
p3 <- ggplot_gtable(ggplot_build(boxplot.DOC30+theme(legend.position= "NONE")))
p4 <- ggplot_gtable(ggplot_build(boxplot.monthcexport+theme(legend.position= "NONE")))

maxWidth = unit.pmax(p1$widths[2:3], p2$widths[2:3], p3$widths[2:3], p4$widths[2:3])
p1$widths[2:3] <- maxWidth
p2$widths[2:3] <- maxWidth
p3$widths[2:3] <- maxWidth
p4$widths[2:3] <- maxWidth

pdf(file=paste0(fig.dir,"/CRFigures_DOCFlux.pdf"), width = 8.5, height = 11) #save figure
grid.arrange(p1,p2, p3, p4,  ncol = 1)
grid.text("A", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("B", x=unit(0, "npc")+ unit(2,"mm"), y=unit(3/4, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("C", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1/2, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("D", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1/3.42, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
dev.off()



# sum daily flux by year to get yearly fluxes
yearly.DOCflux <- ddply(cflux.day,.(format(cflux.day$date, format='%Y')),
                         summarise, yearlymgyear = sum(daily.cflux.mgday.linear, na.rm = TRUE)) # get the daily mean cflux
# get in kg/yearhect
yearly.DOCflux$yearlygyear <- yearly.DOCflux$yearlymgyear/1000
yearly.DOCflux$yearlykgyear <- yearly.DOCflux$yearlygyear/1000
yearly.DOCflux$yearlykgyearhect <- yearly.DOCflux$yearlykgyear / area.hect

# what percentage of flux is explained by DOC conc versus discharge?

fit <- lm(DOC.Q$cflux.gshec ~ DOC.Q$Q.L.s +DOC.Q$DOCcorr)

af <- anova(fit)
afss <- af$"Sum Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))

############# testing for significance
## test 30 minute data for significant differences
plot(density(DOC.Q$cflux.gshec, na.rm = TRUE))

## Perform the test
shapiro.test(DOC.Q$cflux.gshec)

## Plot using a qqplot
qqnorm(DOC.Q$cflux.gshec);qqline(DOC.Q$cflux.gshec, col = 2)
# non-normal again

# test daily flux 
# really slow! Try on daily means rather than 30 minute data
#First, aggregate to the daily mean, then look at whether daily means are different
DOC.daily.flux <- ddply(DOC.Q,.(format(DOC.Q$date, format='%Y-%m-%d')),
                   summarise, meandaily.DOC.gshe = mean(cflux.gshec, na.rm = TRUE)) 

colnames(DOC.daily.flux)[1] <- "date"
DOC.daily.flux$date <- as.POSIXct(strptime(DOC.daily.flux$date, format = "%Y-%m-%d"))

jan.DOC <- subset(DOC.daily.flux, format(DOC.daily.flux$date, format='%m') == "01")
feb.DOC <- subset(DOC.daily.flux, format(DOC.daily.flux$date, format='%m') == "02")
march.DOC <- subset(DOC.daily.flux, format(DOC.daily.flux$date, format='%m') == "03")
april.DOC <- subset(DOC.daily.flux, format(DOC.daily.flux$date, format='%m') == "04")
may.DOC <- subset(DOC.daily.flux, format(DOC.daily.flux$date, format='%m') == "05")
june.DOC <- subset(DOC.daily.flux, format(DOC.daily.flux$date, format='%m') == "06")
july.DOC <- subset(DOC.daily.flux, format(DOC.daily.flux$date, format='%m') == "07")
aug.DOC <- subset(DOC.daily.flux, format(DOC.daily.flux$date, format='%m') == "08")
sept.DOC <- subset(DOC.daily.flux, format(DOC.daily.flux$date, format='%m') == "09")
oct.DOC <- subset(DOC.daily.flux, format(DOC.daily.flux$date, format='%m') == "10")
nov.DOC <- subset(DOC.daily.flux, format(DOC.daily.flux$date, format='%m') == "11")
dec.DOC <- subset(DOC.daily.flux, format(DOC.daily.flux$date, format='%m') == "12")

# compile nov and december in the pre harvest period
nov.DOC$date <- as.Date(gsub("2009", "2010", nov.DOC$date), format = "%Y-%m-%d")
dec.DOC$date <- as.Date(gsub("2009", "2010", dec.DOC$date), format = "%Y-%m-%d")

jan.DOC.test <- Dunns.test(jan.DOC, "by")
feb.DOC.test <- Dunns.test(feb.DOC, "by")
march.DOC.test <- Dunns.test(march.DOC, "by")
april.DOC.test <- Dunns.test(april.DOC, "by")
may.DOC.test <- Dunns.test(may.DOC, "by")
june.DOC.test <- Dunns.test(june.DOC, "by")
july.DOC.test <- Dunns.test(july.DOC, "by")
aug.DOC.test <- Dunns.test(aug.DOC, "by")
sept.DOC.test <- Dunns.test(sept.DOC, "by")
oct.DOC.test <- Dunns.test(oct.DOC, "by")
nov.DOC.test <- Dunns.test(nov.DOC, "by")
dec.DOC.test <- Dunns.test(dec.DOC, "by")

# take the mean of the pre-logging period by month to compare to months in the post-logging period
DOC.monthly <- ddply(DOC.Q,.(format(DOC.Q$date, format='%Y-%m')),
                     summarise, meanmonthly.DOCflux = mean(cflux.gshec, na.rm = TRUE),
                     sdmonthly.DOCflux = sd(cflux.gshec, na.rm = TRUE)) 
write.csv(DOC.monthly, file=paste0(fig.dir,"/CRmonthlymeanDOCflux.csv")) # write mean DOC

######## Do table of mean daily DOC flux by month within the pre/post period + test for significance
# Calculate mean DOC for pe/post by month
DOC.daily.flux$month <- format(DOC.daily.flux$date,"%m")
DOC.daily.flux <- logstatus.f(DOC.daily.flux)
DOC.daily.flux <- wetdry.f(DOC.daily.flux)

stats.logging.DOCflux <- ddply(DOC.daily.flux, c("month", "logstatus"), summarise,
                               mean.DOC = mean(meandaily.DOC.gshe, na.rm = TRUE), sd.DOC=sd(meandaily.DOC.gshe, na.rm = TRUE), 
                               max.DOC = max(meandaily.DOC.gshe, na.rm=TRUE), min.DOC = min(meandaily.DOC.gshe, na.rm=TRUE))

write.csv(stats.logging.DOCflux, file = paste0(fig.dir, "/DOCflux_prepoststats.csv"))
# do test for significance kruskal wallis test 
jan.DOC <- subset(DOC.daily.flux, DOC.daily.flux$month == "01")
jan.kwt <- kruskal.test(meandaily.DOC.gshe ~ as.factor(logstatus), data = jan.DOC)

feb.DOC <- subset(DOC.daily.flux, DOC.daily.flux$month == "02")
feb.kwt <- kruskal.test(meandaily.DOC.gshe ~ as.factor(logstatus), data = feb.DOC[!is.na(feb.DOC$logstatus),])
march.DOC <- subset(DOC.daily.flux, DOC.daily.flux$month == "03")
march.kwt <- kruskal.test(meandaily.DOC.gshe ~ as.factor(logstatus), data = march.DOC[!is.na(march.DOC$logstatus),])
april.DOC <- subset(DOC.daily.flux, DOC.daily.flux$month == "04")
april.kwt <- kruskal.test(meandaily.DOC.gshe ~ as.factor(logstatus), data = april.DOC[!is.na(april.DOC$logstatus),])
may.DOC <- subset(DOC.daily.flux, DOC.daily.flux$month == "05")
may.kwt <- kruskal.test(meandaily.DOC.gshe ~ as.factor(logstatus), data = may.DOC[!is.na(may.DOC$logstatus),])
june.DOC <- subset(DOC.daily.flux, DOC.daily.flux$month == "06")
june.kwt <- kruskal.test(meandaily.DOC.gshe ~ as.factor(logstatus), data = june.DOC[!is.na(june.DOC$logstatus),])
july.DOC <- subset(DOC.daily.flux, DOC.daily.flux$month == "07")
july.kwt <- kruskal.test(meandaily.DOC.gshe ~ as.factor(logstatus), data = july.DOC[!is.na(july.DOC$logstatus),])
aug.DOC <- subset(DOC.daily.flux, DOC.daily.flux$month == "08")
aug.kwt <- kruskal.test(meandaily.DOC.gshe ~ as.factor(logstatus), data = aug.DOC[!is.na(aug.DOC$logstatus),])
sept.DOC <- subset(DOC.daily.flux, DOC.daily.flux$month == "09")
sept.kwt <- kruskal.test(meandaily.DOC.gshe ~ as.factor(logstatus), data = sept.DOC[!is.na(sept.DOC$logstatus),])
oct.DOC <- subset(DOC.daily.flux, DOC.daily.flux$month == "10")
oct.kwt <- kruskal.test(meandaily.DOC.gshe ~ as.factor(logstatus), data = oct.DOC[!is.na(oct.DOC$logstatus),])
nov.DOC <- subset(DOC.daily.flux, DOC.daily.flux$month == "11")
nov.kwt <- kruskal.test(meandaily.DOC.gshe ~ as.factor(logstatus), data = nov.DOC[!is.na(nov.DOC$logstatus),])
dec.DOC <- subset(DOC.daily.flux, DOC.daily.flux$month == "12")
dec.kwt <- kruskal.test(meandaily.DOC.gshe ~ as.factor(logstatus), data = dec.DOC[!is.na(dec.DOC$logstatus),])
# bind together
DOCflux.kwt <- rbind(jan.kwt[1:3], feb.kwt[1:3], march.kwt[1:3], april.kwt[1:3], may.kwt[1:3], june.kwt[1:3], 
                     july.kwt[1:3], aug.kwt[1:3], sept.kwt[1:3], dec.kwt[1:3]) 
write.csv(DOCflux.kwt, file = paste0(fig.dir, "/DOCflux_KWT.csv"))
# calculate the overall difference/pre post

# calculate the percent change - pre and post by month
stats.DOCflux.all <- ddply(DOC.daily.flux, c("logstatus"), summarise,
                           mean.DOC = mean(meandaily.DOC.gshe, na.rm = TRUE), sd.DOC=sd(meandaily.DOC.gshe, na.rm = TRUE), 
                           max.DOC = max(meandaily.DOC.gshe, na.rm=TRUE), min.DOC = min(meandaily.DOC.gshe, na.rm=TRUE))


#### Figure 4
# Hysterisis loops for event pre/post harvest. 
# as per Butturini, A, M Alvarez, and S Bernal. 2008. “Diversity and Temporal Sequences of Forms of DOC and NO3‐Discharge Responses in an Intermittent Stream: Predictable or Random Succession?.” Journal of …. doi:10.1029/2008JG000721/full.
# identification of rainfall events
# relationship between discharge and precipitation
# merge in the rainfall from climate
precip.discharge <- as.data.frame(merge(climate[,c(1,3)], discharge.ave[,c(1,2)], by = "date", all = FALSE))
linearQ.precip <- lm(as.numeric(precip.discharge$Precip) ~ precip.discharge$Q.L.s, na.rm = TRUE)
summary(linearQ.precip)
plot(as.numeric(precip.discharge$Precip) ~ precip.discharge$Q.L.s)
# Really low R2!!
plot(as.numeric(precip.discharge$Precip))
plot(precip.discharge$Q.L.s)

# plot the discharge by year in order to more closely 
discharge2009 <- ggplot(subset(discharge, format(discharge$date, format='%Y') == "2009"), 
                    aes(date, Q.L.s)) + geom_point(size = 0.4) +
  xlab("Date") + ylab("Q (L/s)") + theme()

discharge2010 <- ggplot(subset(discharge, format(discharge$date, format='%Y') == "2010"), 
                    aes(date, Q.L.s)) + 
  geom_point(size = 0.4) +
  xlab("Date") + ylab("Q (L/s)") + theme()

discharge2011 <- ggplot(subset(discharge, format(discharge$date, format='%Y') == "2011"), 
                       aes(date, Q.L.s)) + 
  geom_point(size = 0.4) +
  xlab("Date") + ylab("Q (L/s)") + theme()

discharge2012 <- ggplot(subset(discharge, format(discharge$date, format='%Y') == "2012"), 
                        aes(date, Q.L.s)) + 
  geom_point(size = 0.4) +
  xlab("Date") + ylab("Q (L/s)") + theme()

discharge2013 <- ggplot(subset(discharge, format(discharge$date, format='%Y') == "2013"), 
                        aes(date, Q.L.s)) + 
  geom_point(size = 0.4) +
  xlab("Date") + ylab("Q (L/s)") + theme()

discharge2014 <- ggplot(subset(discharge, format(discharge$date, format='%Y') == "2014"), 
                        aes(date, Q.L.s)) + 
  geom_point(size = 0.4) +
  xlab("Date") + ylab("Q (L/s)") + theme()

# breakpoint analysis to indicate when a hydrologic event occurs
#https://rpubs.com/MarkusLoew/12164
disch2014 <- ts(subset(discharge, format(discharge$date, format='%Y') == "2014")[,c(1,2)])

my.lm <- lm(Q.L.s ~ date, data = disch2014)
summary(my.lm)
library(segmented)
my.seg <- segmented(my.lm, 
                    seg.Z = ~ date, 
                    psi = NA, 
                    control = seg.control(K=50))

# try finding peaks within Q timeseries? peak Q to identify hydrologically significant e=vents within timeseries.
#http://www.inside-r.org/packages/cran/pracma/docs/findpeaks
require("pracma")

test <- findpeaks(disch2014[,2], minpeakheight = 30, thresh = 1)

# finding average dicharge for a event versus baseflow:
# http://www.lancs.ac.uk/~killick/Pub/KillickEckley2011.pdf - use to find average in baseflow conditions, as well as for 
require(changepoint)
pelt.disch2014 = cpt.mean(disch2014,method='PELT', test.stat='Normal')
plot(pelt.disch2014,type='l',cpt.col='blue',xlab='Index',cpt.width=4)
coef(pelt.disch2014)
# Threshold event for rainfall events. Get threshold for precipitation that causes a increase in discharge...

################
## Percent of DOC from top precipitation events pre and post harvest. 
################
# Show the amount of discharge generated and the DOC flux for an event pre and post harvest. Is this relationship changed?

###############
# Timing between peak precip, peak discharge and peak DOC as identified for hydrologically significant events identified previously.


#### Figure 5
# Timing of DOC and precipitation/groundwater levels: is DOC being mobilized into the stream faster after harvest?
# Figure 5 How to show this? Average time between peak precip to peak discharge, DOC?
# Cross correlation function
# http://www.ltrr.arizona.edu/~dmeko/notes_10.pdf
# https://www.researchgate.net/post/Can_anyone_help_to_determine_the_time_lag_between_flood_gauges_from_upstream_to_downstream

#### Figure 6
# What are the drivers of DOC concentration? Model of DOC and climatic/discharge/stream chemistry variables.
# Figure 6: Relationship between [DOC] and variables. Outline model of DOC concentration within the pre and post harvest period. Correlation table? Also, include soil temp?
# linear regression between DOC and different variables
# compile dataframe for analysis - DOC, precip, groundwater, Q, Q(stormflow), Q(baseflow), wet/dry, pre/post,  soilT, solarR
DOCmerged <- Reduce(function(x, y) merge(x, y, by = 'date', all=TRUE), 
                    list(spectro.all[,c(as.numeric(match("date", names(spectro.all))), as.numeric(match("DOCcorr", names(spectro.all))))], 
                         discharge[,c(as.numeric(match("date", names(discharge))), as.numeric(match("Q.L.s", names(discharge))), 
                                      as.numeric(match("baseflow", names(discharge))), as.numeric(match("quickflow", names(discharge))), 
                                      as.numeric(match("EC", names(discharge))),as.numeric(match("pH", names(discharge))),
                                      as.numeric(match("DO", names(discharge))),as.numeric(match("ORP", names(discharge))),
                                      as.numeric(match("Tw", names(discharge))))], 
                         climate[,c(1:7)],
                         groundwater.ave[,c(as.numeric(match("date", names(groundwater.ave))), 6)]))
DOCmerged <- logstatus.f(DOCmerged)
DOCmerged <- wetdry.f(DOCmerged)

# check normality of DOC
qqnorm(DOCmerged$DOCcorr);qqline(DOCmerged$DOCcorr, col = 2)
# check the log transformed DOC
qqnorm(log(DOCmerged$DOCcorr));qqline(log(DOCmerged$DOCcorr), col = 2) # still looks pretty funky. Log transformed data still not normal

# Generalized least square model - DOC and variables
#http://www.r-bloggers.com/in-depth-introduction-to-machine-learning-in-15-hours-of-expert-videos/
DOCmerged <- DOCmerged[complete.cases(DOCmerged[,as.numeric(match("DOCcorr",names(DOCmerged)))]),] # get only the cases where this is DOC data
DOCmerged$Precip <- as.numeric(DOCmerged$Precip)
DOCmerged$Tair<- as.numeric(DOCmerged$Tair)
# add in column for wet/dry. pre/post, month, time of day
mod1 <- gls(DOCcorr ~ Q.L.s * baseflow * Precip + SoilT * Tair * Solar , na.action = na.omit, data = DOCmerged) 
mod2 <- update(mod1, correlation=corARMA(q = 1)) # not working!! hmm...
   
lmmod1 <- lm(DOCcorr ~ Q.L.s + Precip, data = DOCmerged)
summary(lmmod1)
lmmod2 <- lm(DOCcorr ~ Precip, data = DOCmerged)
summary(lmmod2)
lmmod3 <- lm(DOCcorr ~ SoilT, data = DOCmerged)
summary(lmmod3)
lmmod4 <- lm(DOCcorr ~ bt, data = DOCmerged)
summary(lmmod4)

# Use ANOVA to look at models
anova(lmmod1,lmmod4)
summary(mod1)
# Idea: Correlation table between DOC and different climatic variables....
## logistic regression for pre/post; wet/dry. Can use data to predict pre/post; wet'dry?
#http://www.r-bloggers.com/how-to-perform-a-logistic-regression-in-r/

fit <- lm(DOCmerged$DOCcorr ~ DOCmerged$Q.L.s + DOCmerged$baseflow + DOCmerged$quickflow +
          DOCmerged$Tair + DOCmerged$Precip + DOCmerged$SoilT + 
          DOCmerged$Tw)
af <- anova(fit)
afss <- af$"Sum Sq"
all.fit <- cbind(af,PctExp=afss/sum(afss)*100)

## Do for pre/post
DOCpre <- subset(DOCmerged, DOCmerged$logstatus == "pre")
fit.pre <- lm(DOCpre$DOCcorr ~ DOCpre$Q.L.s + DOCpre$baseflow + DOCpre$quickflow +
            DOCpre$Tair + DOCpre$Precip + DOCpre$SoilT 
          + DOCpre$Tw)
af.pre <- anova(fit.pre)
afss.pre <- af.pre$"Sum Sq"
pre.fit <- cbind(af.pre,PctExp=afss.pre/sum(afss.pre)*100)

# post
DOCpost <- subset(DOCmerged, DOCmerged$logstatus == "post")
fit.post <- lm(DOCpost$DOCcorr ~ DOCpost$Q.L.s + DOCpost$baseflow + DOCpost$quickflow +
                 DOCpost$Tair + DOCpost$Precip + DOCpost$SoilT + DOCpost$Tw)
af.post <- anova(fit.post)
afss.post <- af.post$"Sum Sq"
post.fit <- cbind(af.post,PctExp=afss.post/sum(afss.post)*100)

## do for wet/pre
DOCwet <- subset(DOCmerged, DOCmerged$hydro == "wet")
fit.wet <- lm(DOCwet$DOCcorr ~ DOCwet$Q.L.s + DOCwet$baseflow + DOCwet$quickflow +
                DOCwet$Tair + DOCwet$Precip + DOCwet$SoilT 
                + DOCwet$Tw)
af.wet <- anova(fit.wet)
afss.wet <- af.wet$"Sum Sq"
wet.fit <- cbind(af.wet,PctExp=afss.wet/sum(afss.wet)*100)

DOCdry <- subset(DOCmerged, DOCmerged$hydro == "dry")
fit.dry <- lm(DOCdry$DOCcorr ~ DOCdry$Q.L.s + DOCdry$baseflow + DOCdry$quickflow +
                DOCdry$Tair + DOCdry$Precip + DOCdry$SoilT 
              + DOCdry$Tw)
af.dry <- anova(fit.dry)
afss.dry <- af.dry$"Sum Sq"
dry.fit <- cbind(af.dry,PctExp=afss.dry/sum(afss.dry)*100)

# do for wet/pre/post, dry/pre/post
DOCdrypre <- subset(DOCmerged, DOCmerged$hydro == "dry" & DOCmerged$logstatus == "pre")
fit.drypre <- lm(DOCdrypre$DOCcorr ~ DOCdrypre$Q.L.s + DOCdrypre$baseflow + DOCdrypre$quickflow +
                DOCdrypre$Tair + DOCdrypre$Precip + DOCdrypre$SoilT)
af.drypre <- anova(fit.drypre)
afss.drypre <- af.drypre$"Sum Sq"
drypre.fit <- cbind(af.drypre,PctExp=afss.drypre/sum(afss.drypre)*100)

##drypost
DOCdrypost <- subset(DOCmerged, DOCmerged$hydro == "dry" & DOCmerged$logstatus == "post")
fit.drypost <- lm(DOCdrypost$DOCcorr ~ DOCdrypost$Q.L.s + DOCdrypost$baseflow + DOCdrypost$quickflow +
                   DOCdrypost$Tair + DOCdrypost$Precip + DOCdrypost$SoilT)
af.drypost <- anova(fit.drypost)
afss.drypost <- af.drypost$"Sum Sq"
drypost.fit <- cbind(af.drypost,PctExp=afss.drypost/sum(afss.drypost)*100)
#wetpre
DOCwetpre <- subset(DOCmerged, DOCmerged$hydro == "wet" & DOCmerged$logstatus == "pre")
fit.wetpre <- lm(DOCwetpre$DOCcorr ~ DOCwetpre$Q.L.s + DOCwetpre$baseflow + DOCwetpre$quickflow +
                   DOCwetpre$Tair + DOCwetpre$Precip + DOCwetpre$SoilT)
af.wetpre <- anova(fit.wetpre)
afss.wetpre <- af.wetpre$"Sum Sq"
wetpre.fit <- cbind(af.wetpre,PctExp=afss.wetpre/sum(afss.wetpre)*100)
# wetpost
DOCwetpost <- subset(DOCmerged, DOCmerged$hydro == "wet" & DOCmerged$logstatus == "post")
fit.wetpost <- lm(DOCwetpost$DOCcorr ~ DOCwetpost$Q.L.s + DOCwetpost$baseflow + DOCwetpost$quickflow +
                    DOCwetpost$Tair + DOCwetpost$Precip + DOCwetpost$SoilT)
af.wetpost <- anova(fit.wetpost)
afss.wetpost <- af.wetpost$"Sum Sq"
wetpost.fit <- cbind(af.wetpost,PctExp=afss.wetpost/sum(afss.wetpost)*100)

# write table

write.csv(rbind(all.fit, pre.fit, post.fit, wet.fit, dry.fit, drypre.fit, drypost.fit, wetpre.fit, wetpost.fit),
  file = paste0(fig.dir, "/CRDOCconc_pertable.csv"))

########### Correlation table - between DOC concentration and environmental variables
# Correlation between
#cor.DOC <- na.omit(DOCmerged)
colnames(DOCmerged)[17] <- 'Groundwater'
corr.matrix <- cor(as.data.frame(DOCmerged[,c(2:12, 14:17)]), use="pairwise.complete.obs") # correlation matric between variables
# plot the correlation matrix of the spectral parameters
pdf(file=paste0(fig.dir,"/CR_DOC_corrplot.pdf"), width = 11, height = 8.5)
corrplot(corr.matrix, method = "circle") #plot matrix
dev.off()
write.csv(corr.matrix, paste0(fig.dir, "/CRDOC_corrplot.csv"))
library("Hmisc")
corr.matrix.2 <- rcorr(as.matrix(DOCmerged[,c(2:12, 14:17)]), type = c("pearson","spearman"))

##################
# Show diurnal cycles
# plot July and August for the time
julyaug <- subset(cflux.all, format.Date(cflux.all$date, "%m") == "07" | format.Date(cflux.all$date, "%m")=="08")
# add in air temp from biomet
julyaug <- merge(julyaug, climate, by = "date", all = FALSE)
julyaug$day <-  julyaug$date
year(julyaug$day) <- 0
julyaug.DOC.plot <- ggplot(subset(julyaug, julyaug$DOCcorr >=1), aes(day, DOCcorr)) + 
  geom_point(size = 0.4) +
  xlab("Date") + ylab("[DOC] (mg/L)") + 
  #scale_color_manual(values=cbPalette[1:2]) + # colours
  theme(legend.position="top") + theme() + ggtitle("Diurnal Cycles: July - August [DOC] ")

#' Create the multiple plots - DOC, Q and airT
plot1 <- ggplot(subset(julyaug, julyaug$DOCcorr >=1), aes(day, DOCcorr)) + 
  geom_point(size = 0.4) +
  #geom_point(data = subset(julyaug, julyaug$DOCcorr >=1), aes(x = day, y = Q.m3.s, color = "blue"))  + # add in discharge
  labs(x = expression(""),  
       y = expression('[DOC]'~(mg~L^{-1}))) +
  #scale_color_manual(values=cbPalette[1:2]) + # colours
  theme(legend.position="top") + theme() + ggtitle("Diurnal Cycles")
plottest <- subset(julyaug, julyaug$DOCcorr >=1)
# by hour - DOC
jan.DOC.hour <- ggplot(subset(julyaug, julyaug$DOCcorr >=1), aes(format(plottest$date, "%H:%M"), DOCcorr)) + 
  geom_point(size = 0.4) +
  labs(x = expression(""),  
       y = expression('[DOC]'~(mg~L^{-1}))) +
  #scale_color_manual(values=cbPalette[1:2]) + # colours
  #scale_x_datetime(limits = as.Date(jan.post$date), breaks=date_breaks("4 hour"), labels=date_format("%H:%M")) + 
  theme(legend.position="top") + theme(legend.position ="") + ggtitle("") +
  scale_y_continuous(limits = c(0.8, 3)) 

plot2 <- ggplot(subset(julyaug, julyaug$DOCcorr >=1), aes(day, Q.m3.s)) + 
  geom_point(size = 0.4) +
  labs(x = expression(""),  
       y = expression(Q~(m^{3}~s^{-1})))

plot3 <- ggplot(subset(julyaug, julyaug$DOCcorr >=1), aes(day, Tair)) + 
  geom_line() +
  xlab("Date") + ylab(expression("Mean Air Temperature ("~degree~"C)")) 

#AWESOME diurnal cycles
# save figures
pdf(file=paste0(fig.dir,"/CR_DOC_diurnaljulyaug.pdf"), width = 11, height = 8.5)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), ggplotGrob(plot3), size = "max")) # align the plots by date
dev.off()

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
julDOC.plot.DOC <- ggplot(subset(julDOC, julDOC$DOCcorr >=1), aes(hour, DOCcorr)) + 
  geom_point(size = 0.4) +
  xlab("Date") + ylab("[DOC] (mg/L)") + 
  #scale_color_manual(values=cbPalette[1:2]) + # colours
  theme(legend.position="top") + theme() + ggtitle("Diurnal Cycles: July 25 - 31 [DOC] ")

plot1 <- ggplot(subset(julDOC, julDOC$DOCcorr >=1), aes(hour, DOCcorr)) + 
  geom_point(size = 0.4) +
  #geom_point(data = subset(julyaug, julyaug$DOCcorr >=1), aes(x = day, y = Q.m3.s, color = "blue"))  + # add in discharge
  xlab("") + ylab("[DOC] (mg/L)") + 
  #scale_color_manual(values=cbPalette[1:2]) + # colours
  theme(legend.position="top") + theme() + ggtitle("Diurnal Cycles: July - August [DOC] ")

plot2 <- ggplot(subset(julDOC, julDOC$DOCcorr >=1), aes(hour, Q.m3.s)) + 
  geom_point(size = 0.4) +
  labs(x = expression(""),  
       y = expression(Q~(m^{3}~s^{-1}))) 
  
plot3 <- ggplot(subset(julDOC, julDOC$DOCcorr >=1), aes(hour, Tair)) + 
  geom_point(size = 0.4) +
  xlab("Date") + ylab("Air Temperature (degC)")

pdf(file=paste0(fig.dir,"/CR_DOC_diurnaljuly.pdf"), width = 11, height = 8.5)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), ggplotGrob(plot3), size = "last")) # align the plots by date
dev.off()


#### Q/groundwater/DOC/
# not enough groundwater??

#####################################################################
# Supplemental
########## timeseries analysis of discharge/DOC concentration
# https://a-little-book-of-r-for-time-series.readthedocs.io/en/latest/src/timeseries.html
ts.DOC.30min <- ts(spectro.all$DOCcorr, frequency=17520, start=c(2009,11,18,06,00,00))
plot.ts(ts.DOC.30min)
plot.ts(ts.DOC.30min)
# autocorrleaction behavior
acf(ts.DOC.30min, main ="30 min data", na.action = na.pass)

# go to daily mean DOC interpolation
# daily DOC concentartiol
dailyDOC.conc <- ddply(DOCmerged,.(format(DOCmerged$date, format='%Y-%m-%d')),
                       summarise, meandailyDOC= mean(DOCcorr, na.rm = TRUE))
colnames(dailyDOC.conc)[1] <- "date"
dailyDOC.conc$date <-as.POSIXct(strptime(dailyDOC.conc$date, format = "%Y-%m-%d")) 

# interpolate for missing measurements 
date <- as.data.frame(as.POSIXct(strptime(seq(min(dailyDOC.conc$date), max(dailyDOC.conc$date), by = "day"), 
                                          format = "%Y-%m-%d")))
colnames(date)[1] <- "date"

DOC.daily.all <- merge(dailyDOC.conc, date, by = "date", all = TRUE) # merge together, creating NAs for missing data
# percent missing
permissing.day <- length(which(is.na(DOC.daily.all$meandailyDOC)))/length(DOC.daily.all$meandailyDOC)*100 #30% missing 30 minute data!
#  interpolate disharge for missing data - 1 day intervals
DOC.daily.all$splineDOC.mgL <- na.spline(DOC.daily.all$meandailyDOC) #interpolation by spline (polynmial)
DOC.daily.all$approxDOC.mgL <- na.approx(DOC.daily.all$meandailyDOC) #interpolation by linear approximation

# create time series data with the interpolated daily DOC
daily.ts <- ts(DOC.daily.all$approxDOC.mgL, start=c(2009,11,18,06,00,00), frequency=365)
acf(daily.ts , main ="daily data", na.action = na.pass)

DOC.daily.comp <- stl(daily.ts, "periodic") #decompose into seasonal, trend and irregular components
DOC.daily.com <- as.matrix(DOC.daily.comp)
DOC.daily.com <- data.frame(DOC.daily.comp)
seasonal.daily   <- data.frame(DOC.daily.comp$components[,1])
trend.daily      <- DOC.daily.comp$approxDOC.mgL[,2]
remainder.daily  <- DOC.daily.comp$approxDOC.mgL[,3]

decomposed <- stl(daily.ts, s.window="periodic")
seasonal   <- decomposed$daily.ts[,"trend"]
trend      <- decomposed$time.series[,2]
remainder  <- decomposed$time.series[,3]
plot(DOC.daily.comp) # plot components
plot(daily.ts)
lines(trend.daily)
# - no clear trend or seasonality

## Try same on daily cflux
# daily cflux
# create time series data with the interpolated daily DOC flux
daily.ts <- ts(cflux.day$daily.cflux.mgday.linear, start=c(2009,11,18,06,00,00), frequency=365)
acf(daily.ts , main ="daily data", na.action = na.pass)

DOC.daily.comp <- stl(daily.ts, "periodic") #decompose into seasonal, trend and irregular components
DOC.daily.com <- as.matrix(DOC.daily.comp)
DOC.daily.com <- data.frame(DOC.daily.comp)
seasonal.daily   <- data.frame(DOC.daily.comp$components[,1])
trend.daily      <- DOC.daily.comp$approxDOC.mgL[,2]
remainder.daily  <- DOC.daily.comp$approxDOC.mgL[,3]

decomposed <- stl(daily.ts, s.window="periodic")
seasonal   <- decomposed$daily.ts[,"trend"]
trend      <- decomposed$time.series[,2]
remainder  <- decomposed$time.series[,3]
plot(DOC.daily.comp) # plot components
plot(daily.ts)
lines(trend.daily)

# try on weekly cflux
# try on monthly cflux
monthlyDOC.flux <- ddply(cflux.day,.(format(cflux.day$date, format='%Y-%m')),
                       summarise, monthlyDOCflux= sum(daily.cflux.mgday.linear, na.rm = TRUE))

monthly.ts <- ts(monthlyDOC.flux$monthlyDOCflux, start=c(2009,12,01,00,00,00), frequency=12)
acf(monthly.ts , main ="monthly data", na.action = na.pass)

DOC.monthly.comp <- stl(monthly.ts, "periodic") #decompose into seasonal, trend and irregular components
DOC.monthly.com <- as.matrix(DOC.monthly.comp)
DOC.monthly.com <- data.frame(DOC.monthly.comp)
seasonal.monthly   <- data.frame(DOC.monthly.comp$components[,1])
trend.monthly      <- DOC.monthly.comp$approxDOC.mgL[,2]
remainder.monthly  <- DOC.monthly.comp$approxDOC.mgL[,3]

decomposed <- stl(monthly.ts, s.window="periodic")
seasonal   <- decomposed$monthly.ts[,"trend"]
trend      <- decomposed$time.series[,2]
remainder  <- decomposed$time.series[,3]
plot(DOC.monthly.comp) # plot components
plot(monthly.ts)
lines(trend.monthly)

#### discharge - daily Q
daily.ts <- ts(daily.Q.all$approxQ.mm.day, start=c(2007,11,15,00,00,00), frequency=365)
acf(daily.ts , main ="daily data", na.action = na.pass)

Q.daily.comp <- stl(daily.ts, "periodic") #decompose into seasonal, trend and irregular components
Q.daily.com <- as.matrix(Q.daily.comp)
Q.daily.com <- data.frame(Q.daily.comp)
seasonal.daily   <- data.frame(Q.daily.comp$components[,1])
trend.daily      <- Q.daily.comp$approxQ.mgL[,2]
remainder.daily  <- Q.daily.comp$approxQ.mgL[,3]

decomposed <- stl(daily.ts, s.window="periodic")
seasonal   <- decomposed$daily.ts[,"trend"]
trend      <- decomposed$time.series[,2]
remainder  <- decomposed$time.series[,3]
plot(Q.daily.comp) # plot components
plot(daily.ts)
lines(trend.daily)

# monthly Q
monthlyQ <- ddply(daily.Q.all,.(format(daily.Q.all$date, format='%Y-%m')),
                         summarise, monthlyQflux= sum(approxQ.mm.day, na.rm = TRUE))

monthly.ts <- ts(monthlyQ$monthlyQflux, start=c(2007,11,00,00,00,00), frequency=12)
acf(monthly.ts , main ="monthly data", na.action = na.pass)

Q.monthly.comp <- stl(monthly.ts, "periodic") #decompose into seasonal, trend and irregular components
Q.monthly.com <- as.matrix(Q.monthly.comp)
Q.monthly.com <- data.frame(Q.monthly.comp)
seasonal.monthly   <- data.frame(Q.monthly.comp$components[,1])
trend.monthly      <- Q.monthly.comp$approxQ.mgL[,2]
remainder.monthly  <- Q.monthly.comp$approxQ.mgL[,3]

decomposed <- stl(monthly.ts, s.window="periodic")
seasonal   <- decomposed$monthly.ts[,"trend"]
trend      <- decomposed$time.series[,2]
remainder  <- decomposed$time.series[,3]
plot(Q.monthly.comp) # plot components
plot(monthly.ts)
lines(trend.monthly)
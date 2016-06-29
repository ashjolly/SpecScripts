#
#
#Climate files compilation - taken from ParFp analysis file 28Nov2014
#Done to compile climate files from UBC Biomet
# Done to decrease the size and complexity of the ParFP file R file!
#############


###################
rm(list = ls())
ls()

# set path
setwd("/Users/user/Documents/CR_Data/CR_Climate") 
# confirm path
getwd()


# Automatic way - read all of the .csv files corresponding to climate
library(reshape)
library(plyr)

filelist_climate <- list.files(pattern = ".csv$")
climatenames = c("date","Tair","RH","Solar","SoilT","wind", "Precip")
climate <- ldply(filelist_climate, read.csv, header=FALSE, skip=1, col.names=climatenames, na.strings = "NaN")

# change date and time
climate$date <- as.POSIXct(strptime(climate$date, format = "%Y-%m-%d %H:%M", tz = "America/Los_Angeles"))
#sort according to date
climate<- climate[order((climate$date)),] 
climate <- unique(climate)

toNumerics <- function(Date) {
  stopifnot(inherits(Date, c("Date", "POSIXt")))
  day <- as.numeric(strftime(Date, format = "%d"))
  week <- as.numeric(strftime(Date, format = "%U"))
  month <- as.numeric(strftime(Date, format = "%m"))
  year <- as.numeric(strftime(Date, format = "%Y"))
  list(year = year, month = month, week = week, day = day)
}

myDatTime<-toNumerics(as.Date(climate$date))

climate$year<-myDatTime$year

myDat<-within(climate, {
  Precip.cum.mm <- ave(Precip, year, FUN = cumsum,na.rm=FALSE)
})


Precip.cum <- ddply(climate, c("year"), summarise,Precip.cum = sum(as.numeric(Precip),na.rm = TRUE))

xyy.plot <- function(x, y1, y2, y1.par = NULL, y2.par = NULL, first = NULL, y1.last = NULL, y2.last = NULL, ...)
{
  options(warn = -1)
  if (is.null(y1.par)) y1.par <- list()
  if (is.null(y2.par)) y2.par <- list()   
  par(mar = c(5, 4, 4, 5))   
  if (!is.null(first)) eval(first)
  do.call(plot, modifyList(list(x = x, y = y1, xlab = deparse(substitute(x)), 
                                ylab = deparse(substitute(y1)), col = 1, ...), y1.par))
  if (!is.null(y1.last)) eval(y1.last)
  par(new = TRUE)
  do.call(plot, modifyList(list(x = x, y = y2, axes = FALSE, xlab = "", 
                                ylab = "", col = 2), y2.par)) 
  do.call(axis, modifyList(list(side = 4, col = 2, col.ticks = 2, col.axis = 2), y2.par))
  do.call(mtext, modifyList(list(text = deparse(substitute(y2)), 
                                 side = 4, line = 3, col = 2), y2.par))
  if (!is.null(y2.last)) eval(y2.last) 
}

#try xy y plotting
with(climate,xyy.plot(date,Precip,Precip.cum.mm,type = 'h'))


write.table(Precip.cum, file = "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Stream/cumul_precip.csv",
            sep = ",", col.names = NA, qmethod = "double")

write.table(climate, file = "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Stream/climate_test.csv",
            sep = ",", col.names = NA, qmethod = "double")
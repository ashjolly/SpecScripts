# Script for compiling climate data
# CR project - data from biomet
# 28 June0216
##########

climate.comp <- function(climate.dir) {
  # Automatic way - read all of the .csv files corresponding to climate
  setwd(climate.dir)
  filelist_climate <- list.files(path = climate.dir, pattern = ".csv$")
  climatenames = c("date","Tair","RH","Solar","SoilT","wind", "Precip")
  climate <- ldply(filelist_climate, read.csv, header=FALSE, skip=1, col.names=climatenames, na.strings = "NaN", 
                   stringsAsFactors = FALSE)
  # change date and time
  climate$date <- as.POSIXct(strptime(climate$date, format = "%Y-%m-%d %H:%M", tz = "America/Los_Angeles"))
  #sort according to date
  climate<- climate[order((climate$date)),] 
  climate <- unique(climate) # get rid of duplicates
  
  # add in the data from Iain - snow, and the 2014 data
  filelist_climateIain <- list.files(path = paste(climate.dir,"HDF11_Climate", sep= "/"), pattern = ".csv$")
  setwd(paste(climate.dir,"HDF11_Climate", sep= "/"))
  cnames <- c("first", "date","Tair","pressure","ave","Precip","rainsnow", "date2")
  climate.I <- ldply(filelist_climateIain, read.csv, header=FALSE, skip=1, col.names=cnames, na.strings = "NaN", 
                   stringsAsFactors = FALSE)
  # convert date to date
  climate.I$date <- as.POSIXct(strptime(climate.I$date, format = "%Y-%m-%d %H:%M", tz = "GMT"))
  # convert time zone from GMT to pacific standard
  attributes(climate.I$date)$tzone <- "America/Los_Angeles" 
  #sort according to date
  climate.I<- climate.I[order((climate.I$date)),] 
  climate.I <- unique(climate.I) # get rid of duplicates
  
  # merge the snow/precip into the climate variable for all of the time period
  climate.test <- merge(climate, climate.I[,c(2,7)], by = "date", all = TRUE)
  # merge in the 2014 data missing in the data from Zoran - missing from 26August2014 onwards
  climate.missing <- subset(climate.I, climate.I$date >= "2014-08-26")

  climate.test2 <- merge(climate.test, climate.missing[,c(2,3,6,7)], by = c("date", "Tair", "Precip", "rainsnow"), all = TRUE)
  #sort according to date
  climate.test2<- climate.test2[order((climate.test2$date)),] 
  climate.test2 <- subset(climate.test2,!duplicated(climate.test2$date)) # get rid of duplicate dates
  return(climate.2)
  }
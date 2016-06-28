# Script for compiling climate data
# CR project - data from biomet
# 28 June0216
##########

climate.comp <- function(climate.dir) {
  # Automatic way - read all of the .csv files corresponding to climate
  setwd(climate.dir)
  filelist_climate <- list.files(path = climate.dir, pattern = ".csv$")
  climatenames = c("date","Tair","RH","Solar","SoilT","wind", "Precip")
  climate <- ldply(filelist_climate, read.csv, header=FALSE, skip=1, col.names=climatenames, na.strings = "NaN")

  # change date and time
  climate$date <- as.POSIXct(strptime(climate$date, format = "%Y-%m-%d %H:%M", tz = "America/Los_Angeles"))
  #sort according to date
  climate<- climate[order((climate$date)),] 
  climate <- unique(climate) # get rid of duplicates
  return(climate)
}
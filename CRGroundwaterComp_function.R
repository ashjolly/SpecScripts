# 
# Function for compiling groundwater level changes from the tru track
# Riparian groundwater level changes - trutrack closest to the weir
# 28 June 2016
# Ashlee Jollymore's PhD
##############

groundwater.comp <- function(ground.dir){
  ####### groundwater changes
  # groundwater tru track beside weir
  groundwater <- read.csv(file = paste(ground.dir, "/SN_0904203.csv", sep = ""),
                        head=TRUE,sep=",") 
  groundwater$date <- as.POSIXct(strptime(groundwater$Date, format = "%Y-%m-%d  %H:%M:%S", tz="America/Los_Angeles"))
  # segregate groundwater into quick and baseflow according to ecohydro as per discharge
  ground.complete= groundwater[complete.cases(groundwater[,5]),]
  #ground.b <-BaseflowSeparation(ground.complete$waterheightpoint.mm0904203, filter_parameter = 0.925, passes = 3)
  #ground.b$groundbase <- ground.b$bt
  #ground.b$groundquick <- ground.b$qft
  #ground.complete <- cbind(ground.complete, ground.b)

  #Apply openair to create 30 min dataset
  ground.ave <- timeAverage(ground.complete, avg.time = "30 min", start.date = "2010-09-22 18:00:00", fill = TRUE)
 return(ground.ave)
  # choose only the groundwater height, quickflow and baseflow for merging with spectro data
  #ground <- ground.ave[,c(1,6,12,13)]
}
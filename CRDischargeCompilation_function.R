# function for compiling Q from multiple files
# made on 25nov2014 aj
# revised 28 june 2016. Converted to function
# aim:
# 1. Compile Q from multiple files to have one output file
##############

discharge.comp <- function(discharge.dir){
  top <- c('date', 'Q.L.s')
  ########
  discharge.pre2009 = read.csv(file=paste(discharge.dir,"discharge_pre2009_2.csv", sep = "/"), 
                      head=FALSE,skip = 1, sep=",", col.names = top)
  discharge.pre2009$date <- as.POSIXct(strptime(discharge.pre2009$date, format = "%y-%m-%d %H:%M", tz="GMT"))
  discharge1 = read.csv(file = paste(discharge.dir,"CR.Q.Feb2010_Oct2013_2.csv", sep = "/"), 
                    head=FALSE,skip = 1, sep=",", col.names = top)
  discharge1$date <- as.POSIXct(strptime(discharge1$date, format = "%m/%d/%y %H:%M", tz="GMT"))
  discharge2 = read.csv(file = paste(discharge.dir,"/CRDischarge_May2012.csv", sep = "/"), 
                    head=FALSE,skip = 1, sep=",", col.names = top)
  discharge2$date <- as.POSIXct(strptime(discharge2$date, format = "%m/%d/%Y %H:%M", tz="GMT"))
  discharge3 = read.csv(file = paste(discharge.dir,"CR.Q.Feb2010_Feb2014.csv", sep = "/"),  
                    head=FALSE,skip = 1, sep=",", col.names = top)
  discharge3$date <- as.POSIXct(strptime(discharge3$date, format = "%m/%d/%y %H:%M", tz="GMT"))
  discharge4 = read.csv(file = paste(discharge.dir,"TT_0506623,compilation,2010_june2014_30min.csv", sep = "/"),  
                    head=FALSE,skip = 1, sep=",", col.names = top)
  discharge4$date <- as.POSIXct(strptime(discharge4$date, format = "%y-%m-%d %H:%M", tz="GMT"))

  #merge the discharge data from the three sources, and ensure that only unique files are there
  discharge = rbind(discharge.pre2009,discharge1,discharge2,discharge3, discharge4)
  discharge = unique(discharge)

  # sort according to date
  discharge$date <- as.POSIXct(strptime(discharge$date, format = "%Y-%m-%d %H:%M:%S", tz="GMT"))
  discharge <- discharge[order((discharge$date)),]  

  ############################ 
  # calculate baseflow versus stormflow over time
  # uses baeflow calculator in R using BaseflowSeparation in EcoHydrology package
  #remove NA's prior to using baseflow
  discharge.complete= discharge[complete.cases(discharge[,2]),]
  baseflow <-BaseflowSeparation(discharge.complete$Q.L.s, filter_parameter = 0.925, passes = 3)
  #merge baseflow with streamflow
  discharge <- cbind(discharge.complete, baseflow)
  # percent of streamflow due to baseflow and quickflow
  discharge$bfper <- (discharge$bt)/discharge$Q.L.s *100
  discharge$qfper <- (discharge$qft)/discharge$Q.L.s *100

  # convert Q in L.s to m3.s
  m3.L = 0.001
  discharge$Q.m3.s <- discharge$Q.L.s *m3.L
  discharge$logQ.L.s <- log10(discharge$Q.L.s)

  # use open air to convert to 30 minute dataset
  #discharge.1 <- timeAverage(discharge, avg.time = "30 min", start.date = "2007-11-15 12:30:00", fill = FALSE)
  discharge.1 <- timeAverage(discharge, avg.time = "30 min", fill = TRUE)
  
  #convert from GST to PST
  attributes(discharge$date)$tzone <- "America/Los_Angeles"  
  attributes(discharge.1$date)$tzone <- "America/Los_Angeles"  
  return(discharge)
}

# function for trimming to beginning of the eem
# so that EEM begins and ends at same point
# 25nov2015
# Ashlee J's PhD
#########

eemtrimstart <- function(eem, start){
  # ensure that the ex ranges are the same for all of the data - 240 to 200 nm in 2 nm increments
  ex.temp <- colnames(eem)
  
  if(ex.temp[1] != start) {
    # if first value in ex.temp is not 240, trim 
    ex.length <- length(ex.temp)
    # find column where the exitation wavelength is 240 to cut from
    x240 = as.numeric(match(start,names(temp.EEMS)))
    temp.EEMS <- temp.EEMS[,c(x240:ex.length)]
  }
  return(temp.EEMS)
  }
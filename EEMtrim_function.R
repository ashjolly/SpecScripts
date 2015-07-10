# 
# eem trim function

# ashlee jollymore's phd project
# 10july2015
#############

EEMtrim <- function(eem, minex) {

  #trim so that exitation and emission goes from the same
  ex.temp <- colnames(temp.EEMS)
  
  # trim for min excitation - useful if you have wavelengths that start at different points
  if(ex.temp[1] != minex) {
    # if first value in ex.temp is not the min ex you specift, trim to this
    ex.length <- length(ex.temp)
    
    # find column where the exitation wavelength is new min ex to cut from
    xmin = as.numeric(match(minex,names(temp.EEMS)))
    trim.EEMS <- temp.EEMS[,c(xmin:ex.length)]
  } 
  
  return(trim.EEMS)
}


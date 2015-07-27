# 
# abs trim function

# ashlee jollymore's phd project
# 10july2015
#############

abstrim <- function(abs, minex) {
  
  #trim so that exitation and emission goes from the same
  ex.temp <- colnames(abs)
  
  x = dim(abs)[2]
  
  # trim for min excitation - useful if you have wavelengths that start at different points
  if(ex.temp[x] != minex) {
    # if last value in ex.temp is not the min ex you specift, trim to this
    
    # find column where the exitation wavelength is new min ex to cut from
    xmin = match(minex,colnames(abs))
    trim.abs <- as.data.frame(abs[c(1:xmin)])
  } 
  
  if (ex.temp[x] == minex){
    trim.abs <- abs
  }
  
  return(trim.abs)
}
###
# function for Raman correcting file
# Adapted from matlab script
# March 11 2015
############

#get correction factor
Ramancor <- function(blank) {
  Raman = blank[em375:em430, ex350]
  y = Raman
  x = em[em375:em430]
  n = length(x)-1
  summation = 0 
  iteration = 1
  
  for (i in 1:n) { 
    y0 = y[i] 
    y1 = y[i + 1] 
    dx = x[i+1] - x[i]
    summation = summation + dx * (y0 + y1)/2
    iteration = iteration+1
  }
  
  BaseRect = (y[1]+y[n])/2*(x[n]-x[1])
  RamanArea = summation - BaseRect
  
  return(RamanArea)
} 
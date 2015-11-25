###
# function for correcting EEMS according to area under Raman peak of water blank
# Adapted from matlab script
# March 11 2015#
# Ashlee Jollymore's PhD project
############

#get correction factor
Ramancor <- function(blank) {
  Raman = blank[em375:em430, ex350]
  y = Raman
  x = em[em375:em430]
  n = length(x)-1
  
  # first method - slightly underestimates the raman area
  #summation = 0 
  #iteration = 1
  
  #for (i in 1:n) { 
  #  y0 = y[i] 
  #  y1 = y[i + 1] 
  #  dx = x[i+1] - x[i]
  #  summation = summation + dx * (y0 + y1)/2
  #  iteration = iteration+1
  #}
  
  #BaseRect = (y[1]+y[n])/2*(x[n]-x[1])
  #RamanArea = summation - BaseRect
  
  ### Second method - more accurate when compared to eemR results 
  require(pracma)
  RamanArea = trapz(x,y)
  
  return(RamanArea)
} 
#
#
#
# Function for removing frst and second order Raleigh scatter from Aqualog EEMS
# Replaces with Nas
# 30June2015
# Ashlee Jollymore
#####################

raleigh <- function(eem, slitwidth){
  
  #Cut out the Rayleigh scattering
  #slitwidth or bandwidth (nmm) of the scan should be either 12 or 10
  
  Acut = eem
  ex = as.numeric(colnames(Acut))
  em = as.numeric(rownames(Acut))
  
  # identify first order Raleigh and replaced with 0's
  
  for (f in 1:length(ex)) {
    # identify first order, where em = ex +/- slitwidth
    tempem = (em>(ex[f]+slitwidth))
    Acut[tempem,f]= 0
  }
  
  #second order
  for (j in 1:length(ex)) {
  #temp = find(em>(ex(f)*2-slitwidth));
    # upper limit = ex*2 + slitwidth
    # lower limit = ex*2 - slitwidth
    upper = (ex[j]*2)+slitwidth
    lower = (ex[j]*2)-slitwidth
    tempem = as.character(subset(em, em >=lower & em <= upper))
    Acut[c(tempem),j]= 0
  }

  return(Acut)  

}
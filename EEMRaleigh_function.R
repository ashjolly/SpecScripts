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
  
  for (f in 1:length(ex)) {
    tempem = (em<(ex[f]+slitwidth))
    Acut[tempem,f]= 0
  }
  
  #second order
  for (j in 1:length(ex)) {
  #temp = find(em>(ex(f)*2-slitwidth));
    temp = (em>(ex[j]*2-slitwidth))
    Acut[temp,j]= 0
  }

  return(Acut)  

}
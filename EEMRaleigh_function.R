#
#
#
# Function for removing frst and second order Raleigh scatter from Aqualog EEMS
# Replaces with Nas
# 30June2015
# Ashlee Jollymore
#####################

raleigh <- function(eem, slitwidth, ex, em){
  
  #Cut out the Rayleigh scattering
  #slitwidth or bandwidth (nmm) of the scan should be either 12 or 10
  
  Acut = eem
  
  for (f in 1:length(ex)){
  #temp = find(em<(ex(f)+slitwidth));
  temp =  match(em %in% (em>(ex(f)+slitwidth)))
  Acut(temp,f)=NaN
  }
  
  #second order
  for (f in 1:length(ex)) {
  #temp = find(em>(ex(f)*2-slitwidth));
  temp = match(em %in% (em>(ex(f)*2-slitwidth)))
  Acut(temp,f)=NaN;
  }

  return(Acut)  
  
}
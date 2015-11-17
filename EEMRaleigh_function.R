# Function for removing frst and second order Raleigh scatter from Aqualog EEMS
# Replaces with Nas
# 30June2015
# Revised 6Aug2015 to add interpolation function, so that Raleigh lines are not zeros, but are interpolated
# Ashlee Jollymore
#####################

raleigh <- function(eem, slitwidth1, slitwidth2, R1){
  
  #Cut out the Rayleigh scattering
  #slitwidth or bandwidth (nmm) of the scan should be either 12 or 10. May have to change dending on whether the script cuts out enough of the band.
  
  Acut = eem
  ex = seq(240, 800, by = 2)
  em = as.numeric(rownames(Acut))
  
  # identify first order Raleigh and replaced with 0's
  
  for (f in 1:length(ex)) {
    # identify first order, where em = ex +/- slitwidth
    
    upper = (ex[f])+slitwidth1
    lower = (ex[f])-slitwidth1
    tempem = as.character(subset(em, em >=lower & em <= upper))
    Acut[tempem,f]= 0
  }
  
  #second order
  for (j in 1:length(ex)) {
  #temp = find(em>(ex(f)*2-slitwidth));
    # upper limit = ex*2 + slitwidth
    # lower limit = ex*2 - slitwidth
    upper = (ex[j]*2)+slitwidth2
    lower = (ex[j]*2)-slitwidth2
    tempem = as.character(subset(em, em >=lower & em <= upper))
    Acut[c(tempem),j]= 0
  }
  
  ################# Raman
  # identify Raman lines and replace with zeros or interpolate
  # this is a single line with ex = 350, emission from 
  # Note that the majority of this should be taken care of by subtracting the blank, but this region can be excised anyways
  
  
  
  ################# Interpolation
  # as per the eemscat.m script by Rasmus Bro et al.
  # interpolation code from this script:
  # mm=interp1(ax2cut,Eendcut,ax2,'pchip','extrap'); %%% interpolation using cubic option
  
  # Interpolate first order Raleigh lines if R1 = yes 
  
  if (R1 == "yes"){
    for (j in 1:length(ex)){
    
      upper = (ex[j])+slitwidth1
      lower = (ex[j])-slitwidth1
      tempem = as.character(subset(em, em >=lower & em <= upper))
    
     # Cut out the zeros, replace with NaN
      Acut[c(tempem),j]= NaN
    
      # interpolate 
      # gap fill using na.spline function n zoo package. Uses polynomical gap filling
      Acut[,j]= na.spline(Acut[,j])
    }
  }

  
  # Interpolate second order Raleigh lines
  for (j in 1:length(ex)){
    
    upper = (ex[j]*2)+slitwidth2
    lower = (ex[j]*2)-slitwidth2
    tempem = as.character(subset(em, em >=lower & em <= upper))
    
    # Cut out the zeros, replace with NaN
    Acut[c(tempem),j]= NaN

    # interpolate 
    # gap fill using na.spline function n zoo package. Uses polynomical gap filling
    Acut[,j]= na.spline(Acut[,j])
  }

  return(Acut)  

}
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
    tempem = (em<(ex[f]+slitwidth))
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
  
  # identify Raman lines and replace with zeros
  # this is a single line with ex = 350, emission from 
  # Note that the majority of this should be taken care of by subtracting the blank, but this region can be excised anyways
  
  
  ################# Interpolation
  # as per the eemscat.m script by Rasmus Bro et al.
  # interpolation code from this script:
  # mm=interp1(ax2cut,Eendcut,ax2,'pchip','extrap'); %%% interpolation using cubic option

  # Interpolate second order Raleigh lines
  library('pracma')
  for (j in 1:length(ex)){
    
    upper = (ex[j]*2)+slitwidth
    lower = (ex[j]*2)-slitwidth
    tempem = as.character(subset(em, em >=lower & em <= upper))
    
    # Cut out the zeros, replace with NaN
    Acut[c(tempem),j]= 0
    temp.data = cbind(em, Acut[,j])
    rownames(temp.data) = em
    
    #interpolate 
    #interp1(x, y, xi = x,
    #        method = c("constant", "linear", "nearest", "spline", "cubic"))
    test = interp1(temp.data[,1], temp.data[,2], temp.data[(tempem),2], method = c("cubic"), extrap = 'extrap')
    
    Acut[tempem,j]= interp1(em, Acut[,j], as.numeric(tempem), 
                            method = c("cubic"), 'extrap')
    
    
  }

  return(Acut)  

}
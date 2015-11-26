# Function for removing frst and second order Raleigh scatter from Aqualog EEMS
# Replaces with Nas
# 30June2015
# Revised 6Aug2015 to add interpolation function, so that Raleigh lines are not zeros, but are interpolated
# Ashlee Jollymore
#####################

raleigh <- function(eem, slitwidth1, slitwidth2, ramanwidth, R1){
  
  Acut = eem
  ex = seq(240, 800, by = 2)
  em = as.numeric(rownames(Acut)) 
  colnames(Acut) <- ex
  rownames(Acut) <- em
  
  ################# Raman
  # identify Raman lines and replace with zeros or interpolate
  # this is a single line with ex = 350 +/5 nm, em
  # Note that the majority of this should be taken care of by subtracting the blank, but this region can be excised anyways
  
  # from eemR package
  # For water, the Raman peak appears at a wavenumber 3600 cm lower than the
  # incident wavenumber. For excitation at 280 nm, the Raman peak from water
  # occurs at 311 nm. Source : Principles of Fluorescence Spectroscopy (2006) -
  # Third Edition.pdf
  
  ## Convert to wavenumber
  ex_wave_number = 1 / ex * 10000000
  ## For water
  raman_peaks = ex_wave_number - 3600     # I think Stedmon use 3400 TODO
  raman_peaks = 10000000 / raman_peaks    # conversion back to wavelengths
  raman_peaks2 = raman_peaks*2            #second order raman
  
  for (f in 1:length(raman_peaks)) {
    # identify first order, where em = ex +/- slitwidth
    upper = (raman_peaks[f])+ramanwidth
    lower = (raman_peaks[f])-ramanwidth
    tempem = as.character(subset(em, em >=lower & em <= upper))
    Acut[c(tempem),f]= NaN
    
    # second order
    upper = (raman_peaks[f])*2+ramanwidth
    lower = (raman_peaks[f])*2-ramanwidth
    tempem = as.character(subset(em, em >=lower & em <= upper))
    Acut[tempem,f]= NaN
  }
  
  # use na.spline to interpolate Raman bands
  Acut <- as.data.frame(na.spline(Acut))
  rownames(Acut) <- em
  
  ############## Cut out the Rayleigh scattering
  #slitwidth or bandwidth (nmm) of the scan should be either 12 or 10. May have to change dending on whether the script cuts out enough of the band.
  ################# Interpolation of Raleigh
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
  
  # if no interpolation, replace first order Raleigh with 0s
  if (R1 =="no"){
    for (j in 1:length(ex)){
      upper = (ex[j])+slitwidth1
      lower = (ex[j])-slitwidth1
      tempem = as.character(subset(em, em >=lower & em <= upper))
      
      # Cut out the zeros, replace with NaN
      Acut[c(tempem),j]= 0
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
  
  ################# Replace instances where the Em <= Ex with 0s
  # find where em is equal to ex.

  for (j in 1:dim(Acut)[1]){
    temp <- rownames(Acut)[j]
    test <- Acut[j,]
    Acut[j,temp <= colnames(test)] <- 0
    }
  

  # r is putting in very small numbers when it should be zero. Replace very small numbers with 0's
  Acut[ Acut < 1e-10 ] <- 0
  
  return(Acut)  
}
###
# function for calculating the IFE and correcting Raman masked file
# 
# Adapted from CM code, as well as from following references:
# 1. McKnight et al, 2002 Lim Ocean, Vol46(1), DOI: 10.4319/lo.2001.46.1.0038
# 2. 
# 30June2015
# Note that this gives the IFE correction factor, with the actual correction calculated in the loop
############

# EEM.IFC = EEM*10.^(0.5*IFC) #perform inner filter calculation

innerfilter <- function(eem, abs, emfluor, exfluor) {
  
  # em fluor = emission wavelength at which fluorescence EEM was collected
  # ex fluor = excitation wavelength at which fluorescence EEM was collected
  
  # Need to interpolate abs wavelength to reflect wavelength at which fluorescence was collected
  
  waves = as.numeric(gsub("X", "", colnames(abs)))  # wavelengths at which abs collected
 
  data = abs # column containing absorbance data
  
  # Interpolation absorbance for emission wavelengths
  emforIFC = em # emission wavalengths from fluorescence
  
  # note that CM only takes from 1:125 of the emission wavelenghts.. apply to all?
  em_abs = as.data.frame(approx(x = waves, y = data, xout = emforIFC, method = "linear")) #em_abs = interp1(waves,data,emforIFC)'; - initial CM line

  #note that CM file also inserts 0.175 as first em_abs number? Code below
  #em_abs(2:125) = em_abs(1:124);
  #em_abs(1) = 0.1705;
  
  # Interpolation absorbance for excitation wavelengths
  exforIFC = ex
  ex_abs = as.data.frame(approx(x = waves, y = data, xout = exforIFC, method = "linear"))
  IFC <- matrix(ncol = (dim(ex_abs)[1]), nrow = (dim(em_abs)[1]))
  colnames(IFC) <- ex_abs[,1]
  rownames(IFC) <- em_abs[,1]
  
  #loop adapted from CM file
  for (f in 1:(dim(em_abs)[1])) {
     for (b in 1:(dim(ex_abs)[1])) {
         IFC[f,b]=ex_abs[f,2]+em_abs[b,2]
     }
  }  
  
  # Apply IFE (ODex + ODem) to EEMs
  # EEM.IFE = EEM *10 ^(0.5 *IFC) Need to have something to ensure that if IFC = Na, 
  
  IFC <- as.data.frame(IFC)
  IFC.2 <- IFC *0.5
  IFC.2 <- as.data.frame(10^IFC.2)
  
  # replace Nas in IFC with 1
  IFC.2[is.na(IFC.2)] <- 1
  
  # get corrected fluorescence
  EEM.IFC <- eem * IFC.2
  
  return(EEM.IFC)
  
  # will need to do something fro Nas arising because of differences in wavelengths

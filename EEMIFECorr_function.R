###
# function for calculating the IFE and correcting Raman masked file
# 
# Adapted from CM code, as well as from following references:
# 1. McKnight et al, 2002 Lim Ocean, Vol46(1), DOI: 10.4319/lo.2001.46.1.0038
# 2. 
# 30June2015
# Note that this gives the IFE correction factor, with the actual correction calculated in the loop
############


innerfilter <- function(eem, abs, em) {
  
  waves = colnames(abs)  # wavelengths at which abs collected
  data = abs # column containing absorbance data
  
  ex_abs = data # wavelengths at which the absorbance emission was collected
  
  # em data
  emforIFC = em # note that CM only takes from 1:125 of the emission wavelenghts.. apply to all?
  em_abs = approxfun(x = waves, y = data, xout = emforIFC) #em_abs = interp1(waves,data,emforIFC)'; - initial CM line

  #note that CM file also inserts 0.175 as first em_abs number? Code below
  #em_abs(2:125) = em_abs(1:124);
  #em_abs(1) = 0.1705;

  #loop adapted from CM file
  for (f in 1:length(em_abs)) {
     for (b in 1:length(ex_abs)) {
         IFC(f,b)=ex_abs(f)+em_abs(b)
     }
  }  
  
  return(IFC)
}
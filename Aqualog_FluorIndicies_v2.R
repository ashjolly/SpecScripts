###
# function for calculating fluorescence indicies from EEMS file
# 
# March 11 2015
############

Fluor <- function(eem) {
  
  #############
  #calculate the fluorescence index
  #eem = EEMcorr
  FI <- eem[em470,ex370]/eem[em520, ex370]; #Calculates the fluorescence index
  
  #############
  #calculate the humification index
  # Find boundaries for red humic peaks - em = 435; em = 480 
  # % find column where the excitation = 254 nm
  line254 = data.frame(eem[,ex254])
  RedHum = line254[(em435:em480),]
  BlueHum = line254[(em300:em345),]
  RedHum.x <- em[em435:em480]
  BlueHum.x <- em[em300:em345]
  
  #Finds the area under the peaks
  require(pracma)
  
  RedA_area = trapz(RedHum.x, RedHum) ## CHANGE THIS w * sum of heights 
  BlueA_area = trapz(BlueHum.x,BlueHum) ###
 
  RedA_sum = as.numeric(sum(RedHum))
  BlueA_sum = as.numeric(sum(BlueHum))
  
  #This first code calculates HIX with the formula from Ohno(2002)
  HIX_ohno_area = RedA_area/(RedA_area+BlueA_area)
  
  #This calculates HIX using the original formula from Zsolnay (1999).  Ensure
  #you have low concentrations of DOC to not have inner filter or
  #concentration affects
  #Using area under the curve
  HIX_Zsonlay_area = RedA_area/BlueA_area
  
  # HIX from Ohno, according to sums of intensities over the range rather than area under the curve
  HIX_Ohno_sum <- RedA_sum/(RedA_sum+BlueA_sum)
  
  # Hix from Zsolnay (1999), according to sum on intensities
  HIX_Zsonlay_sum <- RedA_sum/(BlueA_sum)
  
  ############
  #calculate the freshness index
  
  FrI = eem[em380, ex310]/max(eem[(em420:em436), ex310])
  
  ############
  fl.out <- cbind(FI, HIX_ohno_area, HIX_Zsonlay_area, HIX_Ohno_sum, HIX_Zsonlay_sum, FrI)
  return(fl.out)
}
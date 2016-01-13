###
# function for calculating fluorescence indicies from EEMS file
# 
# March 11 2015
############

Fluor <- function(eem) {
  
  em <- as.numeric(rownames(eem))
  #############
  #calculate the fluorescence index
  #eem = EEMcorr
  # From McKnight et al (2001), Cory and McKNight (2005)
  # FI indicates if precursor material for DOM is more microbial (1.8) versus terrestrially derived
  # this is from the modified formula
  
  FI <- eem[em470,ex370]/eem[em520, ex370]; #Calculates the fluorescence index
  
  #############
  # Calculate the humification index
  # Ohno (2002) + Zsolnay (1999)
  # Find boundaries for red humic peaks - em = 435; em = 480 
  # % find column where the excitation = 254 nm
  line254 = data.frame(eem[,ex254])
  RedHum = line254[(em435:em480),]
  BlueHum = line254[(em300:em345),]
  RedHum.x <- em[em435:em480]
  BlueHum.x <- em[em300:em345]
  
  #Finds the area under the peaks at 254 nm excitation
  require(pracma)
  
  RedA_area = trapz(RedHum.x, RedHum) ## CHANGE THIS w * sum of heights 
  BlueA_area = trapz(BlueHum.x,BlueHum) ###
 
  RedA_sum = as.numeric(sum(RedHum))
  BlueA_sum = as.numeric(sum(BlueHum))
  
  #This first code calculates HIX with the formula from Ohno (2002)
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
  # calculate the freshness index
  # Also called BIX. Ratio of Beta/alpha peak.
  # Indicates proportion of recently produced DOm. beta peak represents
  # recently created (microbial) DOM, while alpha peak is older, mostly decomposed DOM
  # According to Parlanti et al (2000), Huguet et al (2009), Wilson and Xenopoulos (2009)
  
  FrI = eem[em380, ex310]/max(eem[(em420:em436), ex310])
  
  ############
  # Peak A (humic)
  # references: 
  # Coble, P. G. (1996). Characterization of marine and terrestrial DOM in seawater using excitation-emission matrix spectroscopy. Marine Chemistry. http://doi.org/10.1016/0304-4203(95)00062-3
  # Peak A = max intensity between ex 240-270 nm; em = 380-470 nm
  
  peakA = max(eem[(em380:em470), (ex240:ex270)])
  
  ########### 
  # Peak C (humic)
  # references:
  # Coble, P. G. (1996). Characterization of marine and terrestrial DOM in seawater using excitation-emission matrix spectroscopy. Marine Chemistry. http://doi.org/10.1016/0304-4203(95)00062-3
  # Beggs, K. M. H., & Summers, R. S. (2011). Character and Chlorine Reactivity of Dissolved Organic Matter from a Mountain Pine Beetle Impacted Watershed. Environmental Science & Technology, 45(13), 5717–5724. http://doi.org/10.1021/es1042436
  # Peak c = max intensity between ex 300-340 nm; em 400-450
  
  peakC = max(eem[(em400:em450), (ex300:ex340)])
  
  ############
  # Peak B (protein like - tyrosine-like)
  # references:
  # Coble, P. G. (1996). Characterization of marine and terrestrial DOM in seawater using excitation-emission matrix spectroscopy. Marine Chemistry. http://doi.org/10.1016/0304-4203(95)00062-3
  # Beggs, K. M. H., & Summers, R. S. (2011). Character and Chlorine Reactivity of Dissolved Organic Matter from a Mountain Pine Beetle Impacted Watershed. Environmental Science & Technology, 45(13), 5717–5724. http://doi.org/10.1021/es1042436
  # Peak B = max intensity between ex 260-290 nm; em 300-320 nm
  
  peakB = max(eem[(em300:em320), (ex260:ex290)])
  
  ############
  # Peak T (protein, tryptophan- like)
  # references:
  # Coble, P. G. (1996). Characterization of marine and terrestrial DOM in seawater using excitation-emission matrix spectroscopy. Marine Chemistry. http://doi.org/10.1016/0304-4203(95)00062-3
  # Beggs, K. M. H., & Summers, R. S. (2011). Character and Chlorine Reactivity of Dissolved Organic Matter from a Mountain Pine Beetle Impacted Watershed. Environmental Science & Technology, 45(13), 5717–5724. http://doi.org/10.1021/es1042436
  # Peak T = max intensity between ex 260-290nm; em 326-350 nm
  
  peakT = max(eem[(em326:em350), (ex260:ex290)])
  
  ############
  # Peak T/Peak C ratio
  # reference:
  # Baker, A. (2001). Fluorescence Excitation−Emission Matrix Characterization of Some Sewage-Impacted Rivers. Environmental Science & Technology, 35(5), 948–953. http://doi.org/10.1021/es000177t
  # Used to identify the impact of sewer effluent on a river. 
  # Indicates biochemical oxygen demand relative to DOC
  
  # must interpolate to get ex275 from ex 274 and ex276...
  em.ex275 = (eem[,ex274])+(eem[,ex276]-eem[,ex274])*((275-274)/(276-274))
  peakt.peakC <- em.ex275[em350]/max((eem[(em410:em430), (ex320:ex340)]))
  
  ############ 
  # bind all indicies together
  fl.out <- cbind(FI, HIX_ohno_area, HIX_Zsonlay_area, HIX_Ohno_sum, HIX_Zsonlay_sum, FrI,peakA, peakC, peakB, peakT, peakt.peakC)
  return(fl.out)
}
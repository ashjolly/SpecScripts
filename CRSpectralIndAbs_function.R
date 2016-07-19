###
# function for calculating spectral indicies from s::can spectrolyzer
#
# Mean to be called from other scripts which require calculation of spectral indicies from s::can measurements
# Created 28Jan2015
# Based on MJ's code from Aug 2013 file = spectro.processing(Aug2013)_v2.R
#
# notes for using:
# 1. note that x is data set containing absorbance readings taken by spectrolyzer. Must be merged with par files
# 2. note that dataset must be complete. If it contains NaN's, this function will not work! (open air consideration)
# 3. All compilation (fp + par) and cleaning must occur prior to using this function
# 4. Also note that the exported data does not contain identifiers, and must be merged with original dataset
#   that contains data identifiers (done for flexibility between date and ID identifiers)
##


Abs.ind <- function(spec, pathlength){
  # 0. For many spectrometers, need to transform absorbance units to absorbance coefficients as a = 2.303A/l
  l <- pathlength / 1000 #35 mmm path length expressed in m
  ### Here, the spectrolyser .fp data is in Abs/m data, known as "alpha", equivalent to A/l or D/l (Abs/m)
  #which is related to the Napierian absorption coefficient as "alpha" * 2.303
  
  #1. SUVA
  alpha254 <- spec$X252.5 + (spec$X255 - spec$X252.5)*((254-252.5)/(255-252.5))
  SUVA = alpha254/spec$DOCcorr     
  
  # 2. CDOM absorption ratio at 250 to 365 (a250:a365) is already in the .fp *** Need to change alpha's to a's
  ## [see Spencer et al 2009, doi:10.1029/2008GL036831]
  ## this is also called E2:E3 - see Helms et al 2008, Limnology & Oceanography 53(3): 955-969
  # Spencer 2012, doi:10.1029/2011JG001928 talks about using Napierian, but it works out mathematically to be the same with alpha (Abs/m) or a [as alpha times log(10)] 
  
  e2e3 =  spec$X250/spec$X365
  
  # 3. indicator of humification, E4:E6; alternatively SUVA is used for this - see Helms 2008 pg 955
  e4e6 =  spec$X465/spec$X665
  
  # 4. Total CDOM absorption calculated as the integrated absorption from 250 to 450 nm [see Helms p959]
  x250 = as.numeric(match("X250",names(spec))) # column number where abs is 250 nm
  x450 = as.numeric(match("X450",names(spec))) # column number where abs is 450 nm
  
  CDOM.total =  rowSums(spec[,x250:x450])*200 # the 200 is for length of the integrated spectrum in nm
  
  
  # 5. Slope ratio as proxy for DOM Molecular Weight
  ## see Helms et al 2008, Limnology & Oceanography 53(3): 955-969
  ## and Spencer 2009, doi:10.1029/2008GL036831
  ## and Spencer 2012, doi:10.1029/2011JG001928
  slope295_275 = -log(((spec$X295/spec$X275)/(295 - 275)), base = exp(1))
  slope400_350 = -log(((spec$X400/spec$X350)/(400-350)), base = exp(1))
  slope_ratio =  slope295_275/slope400_350 #this isn't really the slope!
  
  
  #library(gnm)
  library(nlmrt)
  
  #From Spencer 2007 (DOI: 10.1002/hyp.6887): A steep spectral slope (e.g. closer to 0.02) indicates low molecular weight material or decreasing aromaticity while a shallower spectral slope (e.g. closer to 0.01) indicates humic-like or higher molecular weight material with a higher aromatic content (Blough and Del Vecchio, 2002).
  
  sloperatio.f <- function(data){
    # column numbers that are working with for slope calculations
    x275 = as.numeric(match("X275",names(data))) #column number at 275 nm
    x295 = as.numeric(match("X295",names(data))) #column number at 295 nm
    x290 = as.numeric(match("X290",names(data))) #column number at 290 nm
    x350 = as.numeric(match("X350",names(data))) #column number at 350 nm
    x400 = as.numeric(match("X400",names(data))) #column number at 400 nm
    
    #working with slope calculations on natural log-transformed data
    S1.alpha.ln <- log(data[(x275:x295)],base=exp(1)) # 275nm - 295nm
    S2.alpha.ln <- log(data[(x290:x350)],base=exp(1)) # 290nm - 350nm
    S3.alpha.ln <- log(data[(x350:x400)],base=exp(1)) # 350nm - 400nm
    
    S1.alpha <- data[(x275:x295)] # 275nm - 295nm
    S2.alpha <- data[(x290:x350)] # 290nm - 350nm
    S3.alpha <- data[(x350:x400)] # 350nm - 400nm
 
    S1.nm <- seq(275,295,2.5)
    S2.nm <- seq(290,350,2.5)
    S3.nm <- seq(350,400,2.5)
    
    #here, formulated using the regression against nm; code as against index1 if should use a one-dimensional regression, but seems nm is needed - see Fichot 2012 doi:10.4319/lo.2012.57.5.1453, eqn 1
    S1 <- coef(lm(as.numeric(S1.alpha.ln) ~ S1.nm))[2]   
    S2 <- coef(lm(as.numeric(S2.alpha.ln) ~ S2.nm))[2]
    S3 <- coef(lm(as.numeric(S3.alpha.ln) ~ S3.nm))[2]
    
    #from sigmaplot: Equation: Exponential Decay, Single, 2 Parameter
    #   f = a*exp(-b*x)
    
    # Best looks like to implement the nonlinear fit of an exponential function (Spencer 2012)
    # here - y is -- y <- as.numeric(S1.alpha[i,]) ; x is S1.nm <- seq(275,295,2.5)
    
    # S1 = nonlinear fit of an exponential function for S over n275 - n295
    y <- as.numeric(S1.alpha)
    x <- S1.nm
    temp <- data.frame(y=y,x=x)
    regmod <- "y ~ a * exp(-b * x)"
    ones <- c(a=1, b=1) # all ones start
    test.start <- c(a=100, b=0.01)
    anmrtx <- try(nlxb(regmod, start=test.start, trace=FALSE, data=temp))
    S1 <- as.numeric(anmrtx$coef[2])
    
    # S2 = nonlinear fit of an exponential function for S over n290 - n350
    y <- as.numeric(S2.alpha)
    x <- S2.nm
    temp <- data.frame(y=y,x=x)
    regmod <- "y ~ a * exp(-b * x)"
    ones <- c(a=1, b=1) # all ones start
    test.start <- c(a=100, b=0.01)
    anmrtx <- try(nlxb(regmod, start=test.start, trace=FALSE, data=temp))
    S2 <- as.numeric(anmrtx$coef[2])
    
    # S3 = nonlinear fit of an exponential function for S over n350 - n400
    y <- as.numeric(S3.alpha)
    x <- S3.nm
    temp <- data.frame(y=y,x=x)
    regmod <- "y ~ a * exp(-b * x)"
    ones <- c(a=1, b=1) # all ones start
    test.start <- c(a=100, b=0.01)
    anmrtx <- try(nlxb(regmod, start=test.start, trace=FALSE, data=temp))
    S3 <- as.numeric(anmrtx$coef[2])
    
    # Calcuate slope ratio
    SR <- S1/S3
    # Bind all of the indicies together
    slope.coeff <- cbind(S1, S2, S3, SR)
    return(slope.coeff)
  }
  
  data.spec <- data.frame(matrix(NA, nrow = 0, ncol = 4))
  for (i in 1:dim(spec)[1]){
    data.spec[i,] <- sloperatio.f(data = spec[i,])
  }

  spectral <- as.data.frame(cbind(SUVA, e2e3, e4e6, CDOM.total, slope_ratio,  data.spec))# bind together all calculated indicies by column
  colnames(spectral) <= c("SUVA", "e2e3", "e4e6", "CDOM.total", "slope_ratio", "S1", "S2", "S3", "SR")
  return(spectral)
}




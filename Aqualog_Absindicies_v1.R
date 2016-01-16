###
# function for calculating absorbance indicies from absorbance file
# 
# March 11 2015
############

Abs <- function(absorbance) {

  #first absorbance at 254 nm
  abs254 <- absorbance$X254
  
  # 2. CDOM absorption ratio at 250 to 365 (a250:a365) is already in the .fp *** Need to change alpha's to a's
  ## [see Spencer et al 2009, doi:10.1029/2008GL036831]
  ## this is also called E2:E3 - see Helms et al 2008, Limnology & Oceanography 53(3): 955-969
  # Spencer 2012, doi:10.1029/2011JG001928 talks about using Napierian, but it works out mathematically to be the same with alpha (Abs/m) or a [as alpha times log(10)] 
  
  #interpolate for Aqualog absorbance, which is every two nm
  absorbance$X365 <- ((absorbance$X366 - absorbance$X364)/(366-364))+absorbance$X364
  absorbance$X465 <- ((absorbance$X466 - absorbance$X464)/(366-364))+absorbance$X464
  absorbance$X665 <- ((absorbance$X666 - absorbance$X664)/(666-664))+absorbance$X664  
    
  e2e3 =  absorbance$X250/absorbance$X365
  
  # 3. indicator of humification, E4:E6; alternatively SUVA is used for this - see Helms 2008 pg 955
  e4e6 =  absorbance$X465/absorbance$X665
  
  # 4. Total CDOM absorption calculated as the integrated absorption from 250 to 450 nm 
  # Reference :
  # Helms et al 2008, Limnology & Oceanography 53(3): 955-969 [see p959]
  CDOM.total =  as.numeric(rowSums(absorbance[,grep("X250", colnames(absorbance)):grep("X450", colnames(absorbance))])*(450-250)) # the 200 is for length of the integrated absorbancetrum in nm
  
  # try second method of calculating area under CDOM curve... using trapz function from 
  require(pracma)
  CDOM.total.int = trapz(as.numeric(seq(250,450, by = 2)), 
                         as.numeric(absorbance[,grep("X250", colnames(absorbance)):grep("X450", colnames(absorbance))]))
  
  # 5. Slope ratio as proxy for DOM Molecular Weight
  ## see Helms et al 2008, Limnology & Oceanography 53(3): 955-969
  ## and Spencer 2009, doi:10.1029/2008GL036831
  ## and Spencer 2012, doi:10.1029/2011JG001928
  slope294_276 = -log(((absorbance$X294/absorbance$X276)/(294 - 276)), base = exp(1))
  slope400_350 = -log(((absorbance$X400/absorbance$X350)/(400-350)), base = exp(1))
  slope_ratio =  slope294_276/slope400_350 #this isn't really the slope!
  
  #library(gnm)
  library(nlmrt)
  
  # From Spencer 2007 (DOI: 10.1002/hyp.6887): A steep absorbancetral slope (e.g. closer to 0.02) 
  # indicates low molecular weight material or decreasing aromaticity while a shallower absorbancetral 
  # slope (e.g. closer to 0.01) indicates humic-like or higher molecular weight material with a higher 
  # aromatic content (Blough and Del Vecchio, 2002).
  
  # column numbers that are working with for slope calculations
  x274 = as.numeric(match("X274",names(absorbance))) #column number at 275 nm
  x296 = as.numeric(match("X296",names(absorbance))) #column number at 295 nm
  x290 = as.numeric(match("X290",names(absorbance))) #column number at 290 nm
  x350 = as.numeric(match("X350",names(absorbance))) #column number at 350 nm
  x400 = as.numeric(match("X400",names(absorbance))) #column number at 400 nm
  
  #working with slope calculations on natural log-transformed data
  S1.alpha.ln <- log(absorbance[,(x274:x296)],base=exp(1)) # 274nm - 296nm
  S2.alpha.ln <- log(absorbance[,(x290:x350)],base=exp(1)) # 290nm - 350nm
  S3.alpha.ln <- log(absorbance[,(x350:x400)],base=exp(1)) # 350nm - 400nm
  
  S1.alpha <- absorbance[,(x274:x296)] # 274nm - 296nm
  S2.alpha <- absorbance[,(x290:x350)] # 290nm - 350nm
  S3.alpha <- absorbance[,(x350:x400)] # 350nm - 400nm
  
  S1.nm <- seq(274,296,2)
  S2.nm <- seq(290,350,2)
  S3.nm <- seq(350,400,2)
  
  n <- dim(absorbance)[1]  # the number of rows of data on which to do the slope regressions
  n <- as.numeric(n)
  #slopes <- seq(1,n,1)
  #slopes <- data.frame(slopes)
  #x = data.frame(x)
  for(i in 1:n) {
    #here, formulated using the regression against nm; code as against index1 if should use a one-dimensional regression, but seems nm is needed - see Fichot 2012 doi:10.4319/lo.2012.57.5.1453, eqn 1
    S1 <- coef(lm(as.numeric(S1.alpha.ln[i,]) ~ S1.nm))[2]   
    S2 <- coef(lm(as.numeric(S2.alpha.ln[i,]) ~ S2.nm))[2]
    S3 <- coef(lm(as.numeric(S3.alpha.ln[i,]) ~ S3.nm))[2]
    
    #from sigmaplot: Equation: Exponential Decay, Single, 2 Parameter
    #   f = a*exp(-b*x)
    
    # Best looks like to implement the nonlinear fit of an exponential function (Spencer 2012)
    # here - y is -- y <- as.numeric(S1.alpha[i,]) ; x is S1.nm <- seq(275,295,2.5)
    
    # S1 = nonlinear fit of an exponential function for S over n275 - n295
    y <- as.numeric(S1.alpha[i,])
    x <- S1.nm
    temp <- data.frame(y=y,x=x)
    regmod <- "y ~ a * exp(-b * x)"
    ones <- c(a=1, b=1) # all ones start
    test.start <- c(a=100, b=0.01)
    anmrtx <- try(nlxb(regmod, start=test.start, trace=FALSE, data=temp))
    S1 <- as.numeric(anmrtx$coef[2])
    
    # S2 = nonlinear fit of an exponential function for S over n290 - n350
    y <- as.numeric(S2.alpha[i,])
    x <- S2.nm
    temp <- data.frame(y=y,x=x)
    regmod <- "y ~ a * exp(-b * x)"
    ones <- c(a=1, b=1) # all ones start
    test.start <- c(a=100, b=0.01)
    anmrtx <- try(nlxb(regmod, start=test.start, trace=FALSE, data=temp))
    S2 <- as.numeric(anmrtx$coef[2])
    
    # S3 = nonlinear fit of an exponential function for S over n350 - n400
    y <- as.numeric(S3.alpha[i,])
    x <- S3.nm
    temp <- data.frame(y=y,x=x)
    regmod <- "y ~ a * exp(-b * x)"
    ones <- c(a=1, b=1) # all ones start
    test.start <- c(a=100, b=0.01)
    anmrtx <- try(nlxb(regmod, start=test.start, trace=FALSE, data=temp))
    S3 <- as.numeric(anmrtx$coef[2])
    
    absorbance$S1[i] <- S1
    absorbance$S2[i] <- S2
    absorbance$S3[i] <- S3
    absorbance$SR[i] <- S1/S3
  }
  
  absorbancetral <- cbind(abs254, e2e3, e4e6, CDOM.total, CDOM.total.int,slope_ratio, absorbance$S1, absorbance$S2, absorbance$S3, absorbance$SR)# bind together all calculated indicies by column
  # set column names
  colnames(absorbancetral) <- c("abs254", "e2e3", "e4e6", "CDOM.total", "CDOM.total.int", "slope_ratio", "S1", "S2", "S3", "SR")
  
  return(absorbancetral)
}
# Script for running loop that will calculate indicies from corrected EEMS
# This function will basically run two functions for calculating absrbance and fluorescence indicies
# Will take corrected EEM files iwthin a folder, and apply functions to calculate abs and fluor indicies
# Avoids having the loop within the script - idea is to make corrections more adaptable between different projects and easier to apply
# Started 9 August 2015
# Ashlee Jollymore's PhD project
################

calc.indicies <- function(filelist_EEMScor, directoryCorrectedAbs, directorynoncorabs, directoryEEMs, ex.wavelengths, em.wavelengths){

  n = dim(filelist_EEMScor)[1]
  
  Spectral.Indicies = data.frame(matrix(vector(), n,35)) #creating an empty vector
  
  for (i in 1:n){
    
    ###########
    # Calculating absorbance indicies
    # load the Abs file
    setwd(directoryCorrectedAbs)
    abs.temp.corr <-as.data.frame(read.delim(as.character(filelist_EEMScor[i,3]), 
                                        header= TRUE, sep = ",", stringsAsFactors=FALSE))
    
    # load the uncorrected ABS file
    setwd("/Users/user/SpecScripts") 
    source("EEMAbsLoadTrim_function.R")
    abs.temp.uncorr <- ABStrim(graphheadings = filelist_EEMScor, samplewd = directorynoncorabs, loopnum = i, column = 4)
    
    # Calculate absoprtion coefficients from blank corrected data
    # References:
    # Helms, J. R., Stubbins, A., & Ritchie, J. D. (2008). Absorption spectral slopes and slope ratios as indicators of molecular weight, source, and photobleaching of chromophoric dissolved organic matter. Limnology and ….
    
    l <- 10 / 1000 #10 mm (1cm) path length expressed in m
    Naperian = 2.303 * abs.temp.uncorr/l # naperian
    Decadic <- abs.temp.uncorr/l # decadic absorbance. Note that this is what spetrolyzer gives. USE THIS for the 
    
    # call function to calculate absorbance indicies
    setwd("/Users/user/SpecScripts") 
    source("Aqualog_Absindicies_v1.R")
    
    #call the function to calculate indicies
    Abs.ind.dec <- Abs(absorbance = Decadic)
    Abs.ind.Nap <- Abs(absorbance = Naperian)
    
    ##########
    # Calculating fluorescence indicies
    setwd(directoryEEMs)
    EEMcorr <-as.data.frame(read.delim(as.character(filelist_EEMScor[i,2]), 
                                       header= TRUE, sep = ",", stringsAsFactors=FALSE))
    
    # wavelengths for Fluorescence indicies calculation
    # ex wavelengths
    ex370 <- as.numeric(grep(ex.wavelengths[c("ex.370"),], colnames(EEMcorr)))
    ex254 <- as.numeric(grep(ex.wavelengths[c("ex.254"),], colnames(EEMcorr)))
    ex310 <- as.numeric(grep(ex.wavelengths[c("ex.310"),], colnames(EEMcorr)))
    ex274 <- as.numeric(grep(ex.wavelengths[c("ex.274"),], colnames(EEMcorr)))
    ex276 <- as.numeric(grep(ex.wavelengths[c("ex.276"),], colnames(EEMcorr)))
    ex320 <- as.numeric(grep(ex.wavelengths[c("ex.320"),], colnames(EEMcorr)))
    ex340 <- as.numeric(grep(ex.wavelengths[c("ex.340"),], colnames(EEMcorr)))
    ex260 <- as.numeric(grep(ex.wavelengths[c("ex.260"),], colnames(EEMcorr)))
    ex290 <- as.numeric(grep(ex.wavelengths[c("ex.290"),], colnames(EEMcorr)))
    ex240 <- as.numeric(grep(ex.wavelengths[c("ex.240"),], colnames(EEMcorr)))
    ex270 <- as.numeric(grep(ex.wavelengths[c("ex.270"),], colnames(EEMcorr)))
    ex300 <- as.numeric(grep(ex.wavelengths[c("ex.300"),], colnames(EEMcorr)))
  
    # em wavelengths
    em470 <- as.numeric(grep(em.wavelengths[c("em.470"),], rownames(EEMcorr)))
    em520 <- as.numeric(grep(em.wavelengths[c("em.520"),], rownames(EEMcorr))) 
    em435 <- as.numeric(grep(em.wavelengths[c("em.435"),], rownames(EEMcorr))) 
    em480 <- as.numeric(grep(em.wavelengths[c("em.480"),], rownames(EEMcorr)))
    em300 <- as.numeric(grep(em.wavelengths[c("em.300"),], rownames(EEMcorr)))
    em345 <- as.numeric(grep(em.wavelengths[c("em.345"),], rownames(EEMcorr)))
    em380 <- as.numeric(grep(em.wavelengths[c("em.380"),], rownames(EEMcorr)))
    em420 <- as.numeric(grep(em.wavelengths[c("em.420"),], rownames(EEMcorr)))
    em436 <- as.numeric(grep(em.wavelengths[c("em.436"),], rownames(EEMcorr)))
    em350 <- as.numeric(grep(em.wavelengths[c("em.350"),], rownames(EEMcorr)))
    em410 <- as.numeric(grep(em.wavelengths[c("em.410"),], rownames(EEMcorr)))
    em430 <- as.numeric(grep(em.wavelengths[c("em.430"),], rownames(EEMcorr)))
    em320 <- as.numeric(grep(em.wavelengths[c("em.320"),], rownames(EEMcorr)))
    em326 <- as.numeric(grep(em.wavelengths[c("em.326"),], rownames(EEMcorr)))
    em430 <- as.numeric(grep(em.wavelengths[c("em.430"),], rownames(EEMcorr)))
    em400 <- as.numeric(grep(em.wavelengths[c("em.400"),], rownames(EEMcorr)))
    em450 <- as.numeric(grep(em.wavelengths[c("em.450"),], rownames(EEMcorr)))
    
    # call function
    setwd("/Users/user/SpecScripts") 
    source("Aqualog_Fluorindicies_v2.R")
    
    # call function that calculates fluorescent indicies
    Fluor.ind <- Fluor(eem = EEMcorr)
    
    ##########
    # bind fluor indicies with abs indicies as well as the sample id
    samplename <- as.character(filelist_EEMScor[i,1]) # column name where sample ID is 
    
    Spectral.Ind <- cbind(samplename, Abs.ind.dec, Abs.ind.Nap, Fluor.ind) 
    top <- colnames(Spectral.Ind)
    Spectral.Indicies[i,]  <- Spectral.Ind
    colnames(Spectral.Indicies) <- top
  }
  
  return(Spectral.Indicies)
  
}
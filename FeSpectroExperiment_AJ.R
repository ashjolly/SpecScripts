rm(list = ls())
ls()

# set path
iron.dir <- "/Users/user/Dropbox/PhD Work/PhD Data/Iron Experiment"
setwd("/Users/user/Dropbox/PhD Work/PhD Data/Iron Experiment")

#set working directory here
# confirm path
getwd()

# 1. Read in the new par file.  Note that .par files will have to be saved as a .csv file.
OM95=read.csv(paste(iron.dir, "avabs95.csv", sep = "/"), head=FALSE,sep=",", stringsAsFactors=FALSE)
OM1 = read.csv(paste(iron.dir, "avabs1.csv", sep = "/"), head=FALSE,sep=",", stringsAsFactors=FALSE)

# read in the original
raw.OM1 <- read.csv(paste(iron.dir, "lowC_dec2011raw.csv", sep = "/"), head=FALSE,sep=",", stringsAsFactors=FALSE)
colnames(raw.OM1) <- raw.OM1[2,]
colnames(raw.OM1)[4:length(colnames(raw.OM1))] <- paste("X", colnames(raw.OM1)[4:length(colnames(raw.OM1))], sep = "")
raw.OM1 <- raw.OM1[-1:-2,]

raw.OM9.5 <- read.csv(paste(iron.dir, "midC_dec2011raw.csv", sep = "/"), head=FALSE,sep=",", stringsAsFactors=FALSE)
colnames(raw.OM9.5) <- raw.OM9.5[2,]
colnames(raw.OM9.5)[4:length(colnames(raw.OM9.5))] <- paste("X", colnames(raw.OM9.5)[4:length(colnames(raw.OM9.5))], sep = "")
raw.OM9.5 <- raw.OM9.5[-1:-2,]

########### colour blind colour palettes
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

############

# subtract absorbances of iron additions from no iron
OM95 = t(OM95)
OM1 = t(OM1)
OM95 = OM95[-1,]
OM1 = OM1[-1,]

#plot averages
plot (as.numeric(OM95[,1]), as.numeric(OM95[,3]),
      ylab = expression("Abs (Abs/m)"),
	    xlab = expression ("Wavelength (nm)"),
	    main = expression("Absorbance - Mid C"),
	    ylim = c(0,85),  col = 'yellow',
      cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1
      )
      lines(OM95[,1], OM95[,4], col = "black")
      lines(OM95[,1], OM95[,5], col = "blue")
      lines(OM95[,1], OM95[,6], col = "red")
      lines(OM95[,1], OM95[,7], col = "green")
      lines(OM95[,1], OM95[,8], col = "orange")
      legend("topright", c("No Iron", "0.5","1", "1.5", "2.5", "3.5"),
      lty=c(1,1), # gives the legend appropriate symbols (lines)
      lwd=c(2.5,2.5),col=c("yellow","black","blue", "red", "green", "orange")) # gives the legend lines the correct color and width

  plot (as.numeric(OM1[,1]), as.numeric(OM1[,3]),
      ylab = expression("Abs (Abs/m)"),
	    xlab = expression ("Wavelength (nm)"),
	    main = expression("Absorbance - Low C"),
	    ylim = c(0,90),  col = 'yellow',
      cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1
      )
      lines(OM1[,1], OM1[,4], col = "black")
      lines(OM1[,1], OM1[,5], col = "blue")
      lines(OM1[,1], OM1[,6], col = "red")
      lines(OM1[,1], OM1[,7], col = "green")
      lines(OM1[,1], OM1[,8], col = "orange")
      legend("topright", c("No Iron", "0.5","1", "1.5", "2.5", "3.5"),
      lty=c(1,1), # gives the legend appropriate symbols (lines)
      lwd=c(2.5,2.5),col=c("yellow","black","blue", "red", "green", "orange")) # gives the legend lines the correct color and width

#OM9.5  differences
OM9.5.dif0.5 = as.numeric(OM95[,4]) - as.numeric(OM95[,3])
OM9.5.dif1 = as.numeric(OM95[,5]) - as.numeric(OM95[,3])
OM9.5.dif1.5 = as.numeric(OM95[,6]) - as.numeric(OM95[,3])
OM9.5.dif2.5 = as.numeric(OM95[,7]) - as.numeric(OM95[,3])
OM9.5.dif3.5 = as.numeric(OM95[,8]) - as.numeric(OM95[,3])

#OM1    differences
OM1.dif0.5 = as.numeric(OM1[,4]) - as.numeric(OM1[,3])
OM1.dif1 = as.numeric(OM1[,5]) - as.numeric(OM1[,3])
OM1.dif1.5 = as.numeric(OM1[,6]) - as.numeric(OM1[,3])
OM1.dif2.5 = as.numeric(OM1[,7]) - as.numeric(OM1[,3])
OM1.dif3.5 = as.numeric(OM1[,8]) - as.numeric(OM1[,3])

#plot differences
  plot (OM95[,1], OM9.5.dif0.5,
      ylab = expression("DIfference in Absorbance"),
	    xlab = expression ("Wavelength (nm)"),
	    main = expression("Change in absorbance - Mid C"),
	    ylim = c(0,50),
      cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1
      )
      lines(OM95[,1], OM9.5.dif1, col = "blue")
      lines(OM95[,1], OM9.5.dif1.5, col = "red")
      lines(OM95[,1], OM9.5.dif2.5, col = "green")
      lines(OM95[,1], OM9.5.dif3.5, col = "orange")
      legend("topright", c("0.5","1", "1.5", "2.5", "3.5"),
      lty=c(1,1), # gives the legend appropriate symbols (lines)
      lwd=c(2.5,2.5),col=c("black","blue", "red", "green", "orange")) # gives the legend lines the correct color and width

png(paste(iron.dir, "/differencegraph.png", sep = ""),    # create graphic for the         
          width = 5*300,        # 5 x 300 pixels
          height = 3*300,
          res = 300,            # 300 pixels per inch
          pointsize = 6)        # smaller font size
      
      plot (OM1[,1], OM1.dif0.5,
      ylab = expression("Difference in Absorbance (1/m)"),
	    xlab = expression ("Wavelength (nm)"),
	    main = expression("Change in Absorbance - Low C"),
	    ylim = c(0,50),
      col = "#D55E00", type = "l",
      cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1
      )
      lines(OM1[,1], OM1.dif1, col = "#E69F00")
      lines(OM1[,1], OM1.dif1.5, col = "#56B4E9")
      lines(OM1[,1], OM1.dif2.5, col = "#009E73")
      lines(OM1[,1], OM1.dif3.5, col = "#F0E442")
      abline()
      legend("topright", c("0.5 mg/L Fe(III)","1.0 mg/L Fe(III)", "1.5 mg/L Fe(III)", "2.5 mg/L Fe(III)", "3.5 mg/L Fe(III)"), cex = 1,pt.cex = 5,
      lty=c(1,1), # gives the legend appropriate symbols (lines)
      lwd=c(2.5,2.5),col=c("#D55E00","#E69F00", "#56B4E9", "#009E73", "#F0E442")) # gives the legend lines the correct color and width
dev.off()

      
# calculate spectral parameters
OM1 = cbind(OM1, OM1.dif0.5, OM1.dif1, OM1.dif1.5, OM1.dif2.5, OM1.dif3.5)
OM9.5 = cbind(OM95, OM9.5.dif0.5, OM9.5.dif1, OM9.5.dif1.5, OM9.5.dif2.5, OM9.5.dif3.5)

OM1 = t(OM1)
OM9.5 =t(OM9.5)
y = c("", "blank", 'no iron', '0.5','1','1.5','2.5','3.5','0.5diff','1diff','1.5diff','2.5diff','3.5diff' )
iron <- c(0,0.5,1,1.5,2.5,3.5)
rownames(OM1) = y
rownames(OM9.5) = y

lambda <- seq(200,750,2.5) 
lambda.dfr <- c(0,0,lambda)
lambda.a <- paste("alpha",lambda,sep="")
colnames(OM1) = as.character(lambda.a)
OM1 = OM1[-1,]
OM1 = as.data.frame(OM1)
colnames(OM9.5) = as.character(lambda.a)
OM9.5 = OM9.5[-1,]
OM9.5 = as.data.frame(OM9.5)
l <- 35 / 1000 #35 mmm path length expressed in m

# 1. SUVA -- needs a254, then divide a254 by DOC after merging .par and .fp files
OM9.5$alpha254 <- (254-252.5)*as.numeric(as.character(OM9.5$alpha255)) - as.numeric(as.character(OM9.5$alpha252.5))/2.5+as.numeric(as.character(OM9.5$alpha252.5))
OM1$alpha254 <- (254-252.5)*as.numeric(as.character(OM1$alpha255)) - as.numeric(as.character(OM1$alpha252.5))/2.5 + as.numeric(as.character(OM1$alpha252.5))
OM1$SUVA <- OM1$alpha254/1.0
OM9.5$SUVA <- OM9.5$alpha254/9.5
# on raw data
raw.OM1$alpha254 <-(254-252.5)*(raw.OM1$X255 - raw.OM1$X252.5)/2.5+(raw.OM1$X252.5)
raw.OM1$SUVA <- raw.OM1$alpha254

# 2. CDOM absorption ratio at 250 to 365 (a250:a365) is already in the .fp *** Need to change alpha's to a's
## [see Spencer et al 2009, doi:10.1029/2008GL036831]
## this is also called E2:E3 - see Helms et al 2008, Limnology & Oceanography 53(3): 955-969
OM9.5$E2E3 <- as.numeric(as.character(OM9.5$alpha250))/as.numeric(as.character(OM9.5$alpha365)) #aka a250:a365
OM1$E2E3 <- as.numeric(as.character(OM1$alpha250))/as.numeric(as.character(OM1$alpha365)) #aka a250:a365
raw.OM1$E2E3 <- raw.OM1$X250/raw.OM1$X365

# 3. indicator of humification, E4:E6; alternatively SUVA is used for this - see Helms 2008 pg 955
OM9.5$E4E6 <- as.numeric(as.character(OM9.5$alpha465))/as.numeric(as.character(OM9.5$alpha665))
OM1$E4E6 <- as.numeric(as.character(OM1$alpha465)) / as.numeric(as.character(OM1$alpha665))
raw.OM1$E4E6 <- raw.OM1$X465/raw.OM1$X665

# 4. Slope ratio as proxy for DOM Molecular Weight
## see Helms et al 2008, Limnology & Oceanography 53(3): 955-969
## and Spencer 2009, doi:10.1029/2008GL036831
## and Spencer 2012, doi:10.1029/2011JG001928
                                                                    
# Spectral slopes for the ranges 275–295; 290–350; 350-400
#alpha275 - V49
#alpha295 - V57
#alpha290 - V55
#alpha350 - V79
#alpha400 - V99

#library(gnm)
  install.packages("nlmrt")
library(nlmrt)
#From Spencer 2007 (DOI: 10.1002/hyp.6887): A steep spectral slope (e.g. closer to 0.02) indicates low molecular weight material or decreasing aromaticity while a shallower spectral slope (e.g. closer to 0.01) indicates humic-like or higher molecular weight material with a higher aromatic content (Blough and Del Vecchio, 2002).

#working with slope calcualtions on natural log-transformed data   OM1 only
S1.alpha.ln <- log(OM1[,(47:55)],base=exp(1)) # 275nm - 295nm; note that this is shifted over two for these files
S2.alpha.ln <- log(OM1[,(53:77)],base=exp(1)) # 290nm - 350nm ; same
S3.alpha.ln <- log(OM1[,(77:97)],base=exp(1)) # 350nm - 400nm  same shift

S1.alpha <- OM1[,(47:55)] # 275nm - 295nm
S2.alpha <- OM1[,(53:77)] # 290nm - 350nm
S3.alpha <- OM1[,(77:97)] # 350nm - 400nm

S1.nm <- seq(275,295,2.5)
S2.nm <- seq(290,350,2.5)
S3.nm <- seq(350,400,2.5)

n <- dim(OM1)[1]  # the number of rows of data on which to do the slope regressions
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
  
  OM1$S1[i] <- S1
	OM1$S2[i] <- S2
	OM1$S3[i] <- S3
	OM1$SR[i] <- S1/S3
	}



# Peak at 300ish nm in differential spectra find peak for the differences..
#OM1

  install.packages("zoo")
  library(zoo)

#find local max in OM1 samples  
 xz <- as.zoo(OM1.dif0.5)
 rollapply(xz, 3, function(OM1.dif0.5) which.max(OM1.dif0.5)==2)
 rollapply(xz, 3, function(OM1.dif0.5) which.min(OM1.dif0.5)==2) 
  rxz <- rollapply(xz, 3, function(OM1.dif0.5) which.max(OM1.dif0.5)==2)
a =  index(rxz)[coredata(rxz)]

 xz <- as.zoo(OM1.dif1)
 rollapply(xz, 3, function(OM1.dif1) which.max(OM1.dif1)==2)
 rollapply(xz, 3, function(OM1.dif1) which.min(OM1.dif1)==2) 
  rxz <- rollapply(xz, 3, function(OM1.dif1) which.max(OM1.dif1)==2)
b =  index(rxz)[coredata(rxz)]

 xz <- as.zoo(OM1.dif1.5)
 rollapply(xz, 3, function(OM1.dif1.5) which.max(OM1.dif1.5)==2)
 rollapply(xz, 3, function(OM1.dif1.5) which.min(OM1.dif1.5)==2) 
  rxz <- rollapply(xz, 3, function(OM1.dif1.5) which.max(OM1.dif1.5)==2)
c =  index(rxz)[coredata(rxz)]
 
 xz <- as.zoo(OM1.dif2.5)
 rollapply(xz, 3, function(OM1.dif2.5) which.max(OM1.dif2.5)==2)
 rollapply(xz, 3, function(OM1.dif2.5) which.min(OM1.dif2.5)==2) 
  rxz <- rollapply(xz, 3, function(OM1.dif2.5) which.max(OM1.dif2.5)==2)
d =  index(rxz)[coredata(rxz)]

 xz <- as.zoo(OM1.dif3.5)
 rollapply(xz, 3, function(OM1.dif3.5) which.max(OM1.dif3.5)==2)
 rollapply(xz, 3, function(OM1.dif3.5) which.min(OM1.dif3.5)==2) 
  rxz <- rollapply(xz, 3, function(OM1.dif3.5) which.max(OM1.dif3.5)==2)
e =  index(rxz)[coredata(rxz)]

#create string with all of the indicies where local max occur
maxOM1 = cbind(a,b,c,d,e)

#install.packages("stringr")
library(stringr)

maxOM1 = str_c(maxOM1, sep = "", collapse = NULL)
maxOM1 = unique(maxOM1)   #get rid of repeated  indicies
maxOM1 = as.numeric(maxOM1)

#create subset of OM1 data that contains wavelengths and absorbance for local max
maxOM1.subset = OM1[maxOM1,]

#plot max absorbance versus Fe concentration
Fe = c(0,0.5,1,1.5,2.5,3.5)
Fediff = c(0.5,1,1.5,2.5,3.5)

rownames(maxOM1.subset) = maxOM1.subset[,1]
maxOM1.subset = maxOM1.subset[,-1]
maxOM1.subset = maxOM1.subset[,-1]  #get rid of blank

#subset into absorbance and absorbance difference subsets
maxOM1.abs = maxOM1.subset[,1:6]
maxOM1.diff = maxOM1.subset[,7:11]

test1 = maxOM1.abs[2,]
test = coef(lm(as.numeric(test1) ~ Fe))[2]

#do linear regression on response of difference at each wavelength to see which wavlenegth has greatest response

n <- dim(maxOM1.abs)[1] 

for (i in 1:n) {
  mabs <-  coef(lm(as.numeric(maxOM1.abs[i,]) ~ as.numeric(Fe)))[2]
  mabs[i] <- mabs
  }

for (i in 1:n) {
  mdiff <-  coef(lm(as.numeric(maxOM1.diff[i,]) ~ as.numeric(Fediff)))[2]
  mdiff[i] <- mdiff
  }
 
 
 maxOM1.subset = t(maxOM1.subset)   
 
 # plot the e4e6 versus SUVA
 plot(OM1$SUVA ~ OM1$E4E6)
 plot(OM9.5$SUVA ~ OM9.5$E4E6)
 
 # plot variables versus Fe Concentrations
 OM1.sub <- cbind(OM1[2:7,], iron)
 OM9.5.sub <- cbind(OM9.5[2:7,], iron)
 
plot(OM1.sub$iron ~ OM1.sub$SUVA)
plot(OM9.5.sub$iron ~ OM9.5.sub$SUVA)

SUVA.lm <- lm(OM1.sub$SUVA~iron)
summary(SUVA.lm)
 
plot(iron ~ OM1.sub$E4E6)
E4E6.lm <- lm(OM1.sub$E4E6~iron)
summary(E4E6.lm)

####### sexy plotting
# calculate average for low carbon dataset
library(plyr)
colnames(raw.OM1)[1] <- 'Fe'
ave.lowC <- ddply(raw.OM1, 'Fe', summarise, 
                  mean.SUVA=mean(SUVA), sd.SUVA=sd(SUVA), N = length(SUVA), se.SUVA = sd.SUVA/N,
                  mean.E4E6=mean(E4E6), sd.E4E6=sd(E4E6), N = length(E4E6), se.E4E6 = sd.E4E6/N
                  )
saveRDS(ave.lowC, paste(iron.dir, "lowCaverages", sep = "")) 

library(ggpmisc)
library(ggplot2)
my.formula = y~x

## Iron versus SUVA - Low C
png(paste(iron.dir, "/SUVvsFe.png", sep = ""),    # create graphic for the         
    width = 5*300,        # 5 x 300 pixels
    height = 3*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

ggplot(ave.lowC[1:6,], aes(x=iron, y=mean.SUVA)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  ggtitle("1.0 mg/L C Standard") +
  labs(x="Iron Concentration (mg/L)",y="SUVA (L/mg*m)") +
  geom_errorbar(aes(ymin=mean.SUVA-sd.SUVA, ymax=mean.SUVA+sd.SUVA), width=.1) # add in error bars
dev.off()

# high C standard - SUVA
ggplot(OM9.5.sub, aes(x=iron, y=SUVA)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  ggtitle("9.5 mg/L C Standard") +
  labs(x="Iron Concentration (mg/L)",y="SUVA (L/mg*m)")

# low C standard - E4E6
png(paste(iron.dir, "/E4E6vsFe.png", sep = ""),    # create graphic for the         
    width = 5*300,        # 5 x 300 pixels
    height = 3*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

  ggplot(ave.lowC[1:6,], aes(x=iron, y=mean.E4E6)) +
    geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
    geom_point() +
    ggtitle("1.0 mg/L C Standard") +
    labs(x="Iron Concentration (mg/L)",y="E4E6 (1/m)") +
    geom_errorbar(aes(ymin=mean.E4E6-sd.E4E6, ymax=mean.E4E6+sd.E4E6), width=.1) # add in error bars
dev.off()
  
# high C standard - E4E6
ggplot(OM9.5.sub, aes(x=iron, y=E4E6)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  ggtitle("9.5 mg/L C Standard") +
  labs(x="Iron Concentration (mg/L)",y="E4E6 (1/m)")

## SUVA/e4e6 versus FE
raw.OM1$ratio <- raw.OM1$SUVA/raw.OM1$E4E6
OM9.5.sub$ratio <- OM9.5.sub$SUVA/OM9.5.sub$E4E6

ave.ratio <-  ddply(raw.OM1, 'Fe', summarise, 
           mean.ratio=mean(ratio), sd.ratio=sd(ratio), N = length(ratio), se.ratio = sd.ratio/N
  )

## plot ratio versus Fe
png(paste(iron.dir, "/ratiovsFe.png", sep = ""),    # create graphic for the         
    width = 5*300,        # 5 x 300 pixels
    height = 3*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

  ggplot(ave.ratio[1:6,], aes(x=iron, y=mean.ratio)) +
    geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
    geom_point() +
    ggtitle("1.0 mg/L C Standard") +
    labs(x="Iron Concentration (mg/L)",y="SUVA/E4E6 (L/mg)") +
    geom_errorbar(aes(ymin=mean.ratio-sd.ratio, ymax=mean.ratio+sd.ratio), width=.1) # add in error bars
dev.off()

# higher concentration
ggplot(OM9.5.sub, aes(x=iron, y=ratio)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  ggtitle("9.5 mg/L C Standard") +
  labs(x="log(Iron Concentration) (mg/L)",y="SUVA/e4e6")

# SUVA versu e4e6 - low conc
png(paste(iron.dir, "/SUVAvsE4E6.png", sep = ""),    # create graphic for the         
    width = 5*300,        # 5 x 300 pixels
    height = 3*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

  ggplot(ave.lowC[1:6,], aes(x=mean.E4E6, y=mean.SUVA)) +
    geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
    geom_point() +
    ggtitle("1.0 mg/L C Standard") +
    labs(x="E4E6 (1/m)",y="SUVA (L/mg*m)") +
    geom_errorbar(aes(ymin=mean.SUVA-sd.SUVA, ymax=mean.SUVA+sd.SUVA), width=.1) + # add in error bars
    geom_errorbarh(aes(xmin = mean.E4E6-se.E4E6,xmax = mean.E4E6+se.E4E6))
dev.off()

# SUVA versu e4e6 - higher conc
ggplot(OM9.5.sub, aes(x=E4E6, y=SUVA)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  ggtitle("9.5 mg/L C Standard") +
  labs(x="E4E6 (1/m)",y="SUVA (L/mg*m)")

 
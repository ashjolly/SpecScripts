## Overall aim of script
# to provide analysis for CR related data. This includes

# Timeseries analysis of DOC trends
# Multiple Linear Regression to model [DOC]
# God. Begun fuck knows when. Continued ad nauseum, but relaunched with aplomb and excitement 22sept2015
# Ashlee Jollymore's PhD project

##
# purpose: to analyze DOC accoridng to timeseries 
# references: 'The R Book" Michael J Crawley, WIly
# 

############################################################################################
# Initial stage
rm(list = ls())
ls()

###### Necessary toolboxes

###### Set paths for data
# set path
setwd("/Users/ashlee/Dropbox/par and fp compilation")

##### Read in data
# read in spectro.all data (all variables)
spectro.all<- read.csv(file = "/Users/ashlee/Dropbox/par and fp compilation/spectro.all.csv", head=TRUE,sep=",")
attach(spectro.all)
names(spectro.all)
spectro.all$date <- as.POSIXct(strptime(spectro.all$date, format = "%Y-%m-%d %H:%M:%S", tz = "America/Los_Angeles"))
data.frame(spectro.all)

############################################################################################
# Homemade functions

# function for three point moving average
tpma = function(test){
  y <- as.numeric(length(test)-2)
  for (i in 2:(as.numeric(length(test))-1)){
    y[i] <- (as.numeric(test[i-1])+as.numeric(test[i])+as.numeric(test[i+1]))/3
  }
  return (y)
}
############################################################################################
# Time series analysis - 
# see http://a-little-book-of-r-for-time-series.readthedocs.org/en/latest/src/timeseries.html
# http://cran.r-project.org/web/packages/deseasonalize/deseasonalize.pdf
# http://www.jstatsoft.org/v28/i02/paper

## Finding relationship between discharge, groundwater and precip'
# july 2014 note you have complete data for groundwater and discharge from 2012-2013 (2 year) segment
sp.20122013 <-  subset(spectro.all, as.Date(date) >= "2012-01-01 00:00:00" & as.Date(date) <= "2013-12-31 24:00:00")
ts.discharge <- ts(sp.20122013$Q.L.s, frequency=17520, start=c(2012,00,00,00,00,00))
plot(ts.discharge)
ts.ground <- ts(sp.20122013$waterheightpoint.mm0904203, frequency=17520, start=c(2012,00,00,00,00,00))
plot(ts.ground)
ts.precip <- ts(sp.20122013$Precip, frequency=17520, start=c(2012,00,00,00,00,00))
plot(ts.precip)

plot.ts(ts.ground)
par(new = TRUE)
plot.ts(ts.discharge, axes = F, col = "blue", xlab = "", ylab = "",)
axis(4)
mtext("groundwater height (mm)", side = 4, line =3)

## at seasonal decomposition of each
par(mfrow = c(1,1))
acf(ts.discharge, type='p')
acf(ts.ground, type = 'p')
acf(ts.precip, type = 'p')

test <- stl(ts.discharge, "periodic")
plot(test)

test.2 <- stl(ts.precip, "periodic")
plot(test)

## 3 point moving averages
ground.av <- ts(tpma(ts.ground), frequency=17520, start=c(2012,00,00,00,30,00))
plot(ground.av, ylim = c(0,1100))

#try on daily mean
daily.ground <- aggregate(spectro.all$waterheightpoint.mm0904203, list(format(spectro.all$date1)), FUN =mean, na.rm=TRUE)
ma.daily.ground <- ts(tpma(daily.ground$x), frequency=365, start=c(2010,09,22))
plot(ma.daily.ground)

# weekly mean
spectro.all$week <- format((spectro.all$date), "%Y-%U") #get column with weeks labelled
weekly.ground <- aggregate(spectro.all$waterheightpoint.mm0904203, list(format(spectro.all$week)), FUN =mean, na.rm=TRUE)
ma.weekly.ground <- ts(tpma(weekly.ground $x), frequency=52, start=c(2010,09))
plot(ma.weekly.ground)


############################################################################################
# Relationship between groundwater level and streamwater discharge
# plot discharge versus groundwater
install.packages('ggplot2')
require(ggplot2)
p <- ggplot(sp.20122013, aes(x = waterheightpoint.mm0904203, y = Q.L.s)) + geom_point()

##for fun: Polynom, third degree: ?poly
# how to use a polynom in a linear model
gw.complete= sp.20122013[complete.cases(sp.20122013[,'waterheightpoint.mm0904203']),]
Qground.lm3 <- lm(Q.L.s ~ poly(waterheightpoint.mm0904203, 3), data = gw.complete)
Qground.lm4 <- lm(Q.L.s ~ poly(waterheightpoint.mm0904203, 4), data = gw.complete)
Qground.lm5 <- lm(Q.L.s ~ poly(waterheightpoint.mm0904203, 5), data = gw.complete)
summary(Qground.lm3)
summary(Qground.lm4)
summary(Qground.lm5)

# plot
p <- p + geom_abline(intercept = ground.coef[1], 
                     slope = ground.coef[2], 
                     aes(colour = "overall"))
p <- p + geom_smooth(method = "lm",
                formula = y ~ poly(x, degree = 3), 
                se = FALSE, colour = "orange")
p<- p+geom_smooth(method = "lm",
                  formula = y ~ poly(x, degree = 4), 
                  se = FALSE, colour = "red")
p<- p+geom_smooth(method = "lm",
                  formula = y ~ poly(x, degree = 5), 
                  se = FALSE, colour = "blue")
coef(Qground.lm4)
#test quadratic model
height = seq(1,1000, by =10)
Q <- 1362.68741*height^4+2790.16892*height^3+3938.64732*height^2+4892.74652*height+38.39423
plot(height, Q)
#######
# Do breakpoint analysis to see hot to segment data
# as per http://rpubs.com/MarkusLoew/12164
#first on groundwater data
#ground <- spectro.all[,c('date','waterheightpoint.mm0904203')]
#plot(ground$date, ground$waterheightpoint.mm0904203)
# create a linear model 2012-2013 data
sp.20122013$date <- as.POSIXct(strptime(sp.20122013$date, format = "%Y-%m-%d %H:%M:%S", tz = "America/Los_Angeles"))
plot(sp.20122013$date, sp.20122013$waterheightpoint.mm0904203)
ground.lm <- lm(waterheightpoint.mm0904203 ~ date, data = sp.20122013)
summary(ground.lm)

# Extract  coefficients from the overall model
ground.coef <- coef(ground.lm)

# -------------------
# analyse breakpoints
# -------------------
# http://cran.r-project.org/doc/Rnews/Rnews_2008-1.pdf
#install.packages("segmented")
library(segmented)

# have to provide estimates for breakpoints.
# after looking a the data, 
ground.seg <- segmented(ground.lm, 
                    seg.Z = ~ date, 
                    psi = NA)
  #psi =list(date = c(4, 15)))
summary(ground.seg)
# get the breakpoints
ground.seg$psi
# get the slopes
slope(ground.seg)

# get the fitted data
ground.fitted <- fitted(ground.seg)
ground.model <- data.frame(date = ground$date, height = ground.fitted)

# plot the fitted model
ggplot(my.model, aes(x = Distance, y = Elevation)) + geom_line()

##############################
###DOC
#find start time for DOC data
#q = as.numeric(match("DOCcorr",names(spectro.all)))
#c.flux= spectro.all[complete.cases(spectro.all[,q]),]

# data analysis as timeseries - 30 min data
ts.DOC.30min <- ts(spectro.all$DOCcorr, frequency=17520, start=c(2009,11,18,06,00,00))
plot(ts.DOC.30min)





####
# three point moving average for DOC concentration - 30 min data (original)
ma.DOC <- ts(tpma(test = ts.DOC.30min), frequency=17520, start=c(2009,11,18,06,00,00))
plot(ts.DOC.30min, ylim = c(0,20))
lines(ma.DOC, col = "blue")

# three point moving average - daily averages
daily <- read.csv(file = "/Users/ashlee/Dropbox/par and fp compilation/dailymean.csv", head=TRUE,sep=",")
attach(daily)
names(daily)
daily.ts <- ts(daily$mean.DOC., start=c(2009,11,18,06,00,00), frequency=365)
plot.ts(daily.ts)

#
ma.DOC.daily <- ts(tpma(test = daily.ts),start=c(2009,11,18,06,00,00), frequency=365) 
plot(ts.DOC.30min, ylim = c(0,20), main = "3-Point Moving Average: Daily Mean [DOC]")
points(ma.DOC.daily, col = "blue", ylim = c(0,20), pch = 20)
legend ("topleft", c("30 minute data", "Three point average (daily mean)"),pch = c(20,20), col = c("black", "blue"), cex = 0.6, bty = "n")

##########
# three point moving average - weekly averages
#First calculate weekly means of DOC data
weekly <- aggregate(spectro.all$DOCcorr, list(format(spectro.all$week)), FUN =mean, na.rm=TRUE)
# convert to ts object
ts.weekly <- ts(weekly$x, frequency=52, start=c(2009,11,18,06,00,00))  # make weekly timeseries from mean
plot(ts.weekly) #plot just to see
#calculate three point moving average - weekly data
ma.DOC.weekly <- ts(tpma(test = ts.weekly), frequency=52, start=c(2009,11,18,06,00,00)) 
#Plot agains 30 minute data
plot(ts.DOC.30min, ylim = c(0,20), main = "3-Point Moving Average: Weekly Mean [DOC]")
points(ma.DOC.weekly, col = "blue", pch = 20)
legend ("topleft", c("30 minute data", "Three point average (weekly mean)"),pch = c(20,20), col = c("black", "blue"), cex = 0.6, bty = "n")

####
# three point moving average - monthly averages
#monthly <- aggregate(spectro.all$DOCcorr, list(format(spectro.all$date2)), FUN =mean, na.rm=TRUE)
monthly  <- aggregate(spectro.all$DOCcorr, list(format(spectro.all$date2)), FUN =mean, na.rm=TRUE)

ts.month<- ts(monthly$x, frequency = 12, start=c(2008,12))
plot(ts.month)
ma.monthly <- ts(tpma(test = ts.month),frequency=12, start=c(2008,12)) #calculate moving average time series
#plot
plot(ts.DOC.30min, ylim = c(0,20), main = "3-Point Moving Average: Monthly Mean [DOC]") 
points(ma.monthly, col = "blue", pch = 20)
legend ("topleft", c("30 minute data", "Three point average (monthly mean)"),pch = c(20,20), col = c("black", "blue"), cex = 0.6, bty = "n")

# make a nice time axis
#d.range <- range(y1$date)
#d.list <- seq(d.range[1], d.range[2], by='week')

#######
# Timeseries analysis: looking at seasonality, correlation etc
# correlation and autocorrlation

# autocorrleaction behavior
acf(ts.DOC.30min, main ="30 min data", na.action = na.pass)
acf(daily.ts , main ="daily data", na.action = na.pass)
acf(ts.weekly, main ="weekly data", na.action = na.pass)
acf(ts.month, main ="monthly mean", na.action = na.pass)

##########
#decompose into seasonal, trend
# references for looking at timeseries data with missing data
# http://arxiv.org/pdf/0811.0659.pdf

#install.packages("mgcv")
require(mgcv)
require(nlme)
mod <- gamm(zoo.all$DOCcorr ~ zoo.all$date, data = zoo.all,
            correlation = corCAR1(form = ~ zoo.all$date))

### Interpolate a vector time series and highlight the imputed data
# monthly data only
# note that interpTs only works with monthly data
# interpolate to get missing data
#install.packages('wq')
require(wq)
interp.month <- interpTs(ts.month, gap = 1)

# plot interp versus original
plot(ts.month, col = 'red')
points(interp.month, col = 'blue', pch = 19)

# try seasonal decomposition
interp.month= interp.month[complete.cases(interp.month)] #get rid of NAs at beginning
interp.month = ts(interp.month, frequency = 12, start=c(2009,11)) #convert back to time series

DOC.month.comp <- stl(interp.month, "periodic") #decompose into seasonal, trend and irregular components
DOC.monthly.com <- as.matrix(DOC.month.comp)
DOC.monthly.com <- data.frame(DOC.monthly.comp)
seasonal.monthly   <- data.frame(DOC.month.comp$components[,1])
trend.monthly      <- DOC.month.comp$interp.month[,2]
remainder.monthly  <- DOC.month.comp$interp.month[,3]

decomposed <- stl(interp.month, s.window="periodic")
seasonal   <- decomposed$interp.month[,"trend"]
trend      <- decomposed$time.series[,2]
remainder  <- decomposed$time.series[,3]
plot(DOC.month.comp) # plot components
plot(ts.month)
lines(trend.monthly)

#try interp on series median
interp2.month<- interpTs(ts.month, type = "series.median", gap = 1)
plot(ts.month, col = 'black')
points(interp2.month, col = 'blue', pch = 19)
points(test.month,col = 'red', pch = 19)
test2.month= test2.month[complete.cases(test2.month)]
test2.month = ts(test2.month, frequency = 12, start=c(2009,11))
DOC.month.comp <- stl(test2.month, "periodic") 
# same as before

#############
# from 
install.packages("gamair")
require(mgcv)
require(gamair)
data(ts.DOC.30min)
cairo2 <- within(cairo, Date <- as.Date(paste(year, month, day.of.month, 
                                              sep = "-")))
plot(temp ~ Date, data = cairo2, type = "l")


DOC.ts.comp$seasonal
plot(DOCtimeseriescomponent)
#install.packages("deseasonalize")
require(deseasonalize)
plot(DOCtimeseriescomponent)
DOCseasonallyadjusted <- DOCtimeseriescomponent - DOCtimeseriescomponent$seasonal


############
# data analysis as zoo object
# see http://cran.r-project.org/web/packages/zoo/vignettes/zoo.pdf
# http://blog.revolutionanalytics.com/2013/06/learning-time-series-with-r.html

library(zoo)
zoo.all <- zoo(spectro.all)

## checking of a strictly regular zoo series
class(zoo.all)
frequency(zoo.all) ## extraction of frequency attribute = 1... issue?
is.regular(zoo.all)

# plot timeseries
plot.ts(DOC)
plot(zoo.all$date, zoo.all$DOCcorr)

############################################################################################
# Compare DOC quality metrics to soil leechate qualities
# Aim is to quantify DOC quality metrics from lysimeter and soil extracts
# Determine what the trends are in DOC metrics with soil depth
# Used as apreliminary stage for quantifying flowpaths of DOC into stream - will compare trends to streamwater DOC quality


############################################################################################
# Multiple linear regression
# Aim is to model [DOC] according to climatic and hydrologic variables.
# [DOC] = f1*discharge + f2*solarradiation + f3*airtemp + f4*precip..etc
# First, use multiple linear regression approach to first find the variables that are most critical to model
# Second, model with multiple linear regression, using gradient descent to find the best fits for all of teh different factors.
# Premise - can we use climatic and environmental factors to predict [DOC]?

# Thoughts:
# May have to partition into pre and post logging periods. Could examine how factors change in terms of the effect of harvest
# Also, use this approach to predict the concentration of various DOM elements 
#        - question being can you model DOM fraction on the basis of the same variables as DOC, and what does this say about hydrologic flowpaths?
#        - compare to soil water DOC qualities to assign hydrologic flowpaths

# see http://www.r-bloggers.com/regression-via-gradient-descent-in-r/





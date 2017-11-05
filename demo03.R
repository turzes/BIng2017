#Install necessary package
#install.packages(TTR)
#install.packages(PerformanceAnalytics)
#install.packages(quantmod)
#install.packages(caTools)
#install.packages(caret)
#install.packages(reshape)
#install.packages(ggplot2)
#install.packages(pbapply)
#install.packages(lbfgs)


#########################################################
# PPPAENET: demo
# Bingbing Li
#########################################################
# Clear workspace
rm(list=ls()) 

library(TTR)
library(PerformanceAnalytics)
library(quantmod)
library(caTools)
library(caret)
library(reshape2)
library(ggplot2)
library(pbapply)
library(lbfgs)


# files location
filelocation <- "~/Documents/working directory/PPPrcodes/"
# Load C++ functions
#cfunfile <- paste(filelocation, 'cfuns.cpp', sep='')


# Load systematic investor toolbox
con = gzcon(file(paste(filelocation, 'sit.gz', sep=''), 'rb'))
source(con)
close(con)

source(paste(filelocation, 'utilityFunctions.R', sep=''))
source(paste(filelocation, 'signalFunctions.R', sep=''))
source(paste(filelocation, 'estimation.R', sep=''))
source(paste(filelocation, 'backtest.R', sep=''))

# Setup start/end dates for estimation
startdate <-'2014-01-01'
enddate <-'2017-05-11'
# Create dataset ia (a list contains prices,ret.vol)
ia <- test.create.historical.data(filelocation = filelocation, dat.from = "2008::")

# Compute Momentum signals
charlist.mom <- create.MOM.signals(ia)
# EMA1
charlist.ema1 <- create.EMA1.signals(ia, from = startdate)
# Donchian Channel
charlist.dc <- create.DC.signals(ia, from = startdate)
# Combine all signals into a list
charlist.tf <- c(charlist.mom,charlist.ema1,charlist.dc)
# Apply position sizing function: vp (volatility parity)
charlist.tf.vp <- volatilityWeighted(charlist.tf, ia$rets.vol)

# Load macro economic data
macrolist <- create.marco.signals(ia, filelocation = filelocation)
# Select macro economic data into a list
charlist.macro <- macrolist[c("IR","ICSA","AMTMNO")]
# Combine with trend signals
charlist <- c(charlist.tf.vp,charlist.macro)

# Bechmark weights: 1/N
BenchW <- ntop(ia$prices, ncol(ia$prices))
# Transaction cost: 5bps
cost <- 0.05/100

# Cross validation of trend and marco signals
ppp.cv.marco <- crossvalidate(charlist, BenchW, ia$prices, cost, gamma = 20, alpha = 0, 
                              startdate = startdate, enddate = enddate,trainDataRatio = c(0.8,0.2),
                              type = "CEG.out",nfolds = 10,lambda1 = NULL, minRatio = 1e-10, 
                              nLambda1 = 100, lambda2 = c(10^(-(0:7)),0), optcost = "NO_COST", inference.type = "boot",
                              name="Trend and Macro Vol Parity",show = TRUE)
# Summary of cross validation results
summary(ppp.cv.marco)
# Estimated coefficients
print(coef(ppp.cv.marco))
# Estimated standard errors of coefficients
print(coef.se(ppp.cv.marco))

# Cross validation of trend signals
ppp.cv.tf <- crossvalidate(charlist.tf.vp, BenchW, ia$prices, cost, gamma = 20, alpha = 0, 
                           startdate = startdate, enddate = enddate,trainDataRatio = c(0.8,0.2),
                           type = "CEG.out",nfolds = 10,lambda1 = NULL, minRatio = 1e-10, 
                           nLambda1 = 100, lambda2 = c(10^(-(0:7)),0), optcost = "NO_COST", inference.type = "boot",
                           name="Trend Vol Parity",show = TRUE)
# Summary of cross validation results
summary(ppp.cv.marco)
# Estimated coefficients
print(coef(ppp.cv.marco))
# Estimated standard errors of coefficients
print(coef.se(ppp.cv.marco))

# Backtest of trend and marco signals
ppp.bt.macro <- backtest(ppp.cv.marco,initialWindow = NULL, horizon = 1, skip = 0,fixedWindow = FALSE, oos = TRUE)
# Summary of backtest results
summary(ppp.bt.macro)
# Plot of estimated coefficients over time
plot.coef(ppp.bt.macro)
# Plot of backtest performance
X11()
plot(ppp.bt.macro)

ppp.bt.tf <- backtest(ppp.cv.tf,initialWindow = NULL, horizon = 1, skip = 0,fixedWindow = FALSE, oos = TRUE)
# Summary of backtest results
summary(ppp.bt.tf)
# Plot of estimated coefficients over time
plot.coef(ppp.bt.tf)
# Plot of backtest performance
X11()
plot(ppp.bt.tf)

# Compare of PPP models with different signals
models<- list()
models$trend <- ppp.bt.tf
models$trend.with.marco<- ppp.bt.macro

X11()
plotbtlist(models)
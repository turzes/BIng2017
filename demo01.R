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
#install.packages(Rcpp)
#install.packages(RcppArmadillo)

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

# Optional: Load C++ functions
# cfunfile <- paste(filelocation, 'cfuns.cpp', sep='')

# Load systematic investor toolbox
#con = gzcon(file(paste(filelocation, 'sit.gz', sep=''), 'rb'))
#source(con)
#close(con)

source(paste(filelocation, 'utilityFunctions.R', sep=''))
source(paste(filelocation, 'signalFunctions.R', sep=''))
source(paste(filelocation, 'estimation.R', sep=''))
source(paste(filelocation, 'backtest.R', sep=''))


# load price data
load(paste(filelocation, 'Qdata150517.RData', sep=''))

dataq[sapply(dataq, is.null)] <- NULL

dataq <- lapply(dataq, function(data){
  data1 <- xts(data[,-1], order.by = as.Date(as.character(data[,1]),"%Y-%m-%d"))
  outdata <- data1[,c("Open", "High","Low", "Settle","Volume","Prev. Day Open Interest")]   
  colnames(outdata) <- c("Open", "High","Low", "Close","Volume","Open.Interest")
  outdata
})

prices <- Reduce(cbind, lapply(dataq,function(x){x[,"Close"]}))

comxts <- prices[,1]
comxts[] <- NA

dataq <- lapply(dataq, function(data){
  outdata <- cbind(data,comxts)
  outdata <- na.locf(outdata)
  outdata[,-ncol(outdata)]
})

prices <- Reduce(cbind, lapply(dataq,function(x){x[,"Close"]}))
colnames(prices) <- names(dataq)

dat.from <- "2008::"

prices <- prices[dat.from]
idx <- apply(!is.na(prices), 2, all)

prices <- prices[,idx]
data <- setupListRange(dataq,dat.from, outtype = "xts")
data <- data[idx]

# Compute Returns
rets <- prices/mlag(prices) - 1
# Compute Returns volatiliy
rets.vol <- bt.apply.matrix(rets, runsd, k = 21)*sqrt(252)

# Setup start/end dates for estimation
startdate <-'2014-01-01'
enddate <-'2017-05-11'
# Setup start dates to construct trend signals (here, same as startdate)
from <- startdate

# Get index of startdate/enddate
dates <- index(prices)
if(is.character(startdate)) startdate <- which(dates ==  dates[dates>=startdate][1])
if(is.character(enddate)) enddate <- which(dates ==  dates[dates<=enddate][length(dates[dates<=enddate])]) 

# Compute Momentum signals
mom01 <- getMomentum(prices, nmonth = 1)
mom03 <- getMomentum(prices, nmonth = 3)
mom06 <- getMomentum(prices, nmonth = 6)
mom09 <- getMomentum(prices, nmonth = 9)
mom12 <- getMomentum(prices, nmonth = 12)
mom18 <- getMomentum(prices, nmonth = 18)
mom24 <- getMomentum(prices, nmonth = 24)
mom36 <- getMomentum(prices, nmonth = 36)
mom48 <- getMomentum(prices, nmonth = 48)

# EMA1
ema1.10<- getEMA1(prices,n = 10, from = from)
ema1.20<- getEMA1(prices,n = 20, from = from)
ema1.50 <- getEMA1(prices,n = 50, from = from)
ema1.100 <- getEMA1(prices,n = 100, from = from)
ema1.200 <- getEMA1(prices,n = 200, from = from)
ema1.300 <- getEMA1(prices,n = 300, from = from)
ema1.400 <- getEMA1(prices,n = 400, from = from)

# DC
dc.10 <- getDonchianChannel(data, n = 10, from = from)
dc.20 <- getDonchianChannel(data, n = 20, from = from)
dc.50 <- getDonchianChannel(data, n = 50, from = from)
dc.100 <- getDonchianChannel(data, n = 100, from = from)
dc.200 <- getDonchianChannel(data, n = 200, from = from)
dc.300 <- getDonchianChannel(data, n = 300, from = from)
dc.400 <- getDonchianChannel(data, n = 400, from = from)

charlist.mom <- list(mom01=mom01,mom03=mom03,mom06=mom06, 
                     mom09=mom09,mom12=mom12,mom24=mom24,
                     mom48=mom48)
charlist.ema1 <- list(ema1.10=ema1.10,ema1.20=ema1.20,ema1.50=ema1.50, 
                      ema1.100=ema1.100, ema1.200=ema1.200,
                      ema1.300=ema1.300, ema1.400=ema1.400)
charlist.dc <- list(dc.10=dc.10,dc.20=dc.20,dc.50=dc.50, 
                    dc.100=dc.100, dc.200=dc.200,dc.300=dc.300, 
                    dc.400=dc.400)

# Combine all signals into a list
charlist.tf <- c(charlist.mom,charlist.ema1,charlist.dc)
# Apply position sizing function: vp (volatility parity)
charlist.tf.vp <- volatilityWeighted(charlist.tf, rets.vol)

# Bechmark weights: 1/N
BenchW <- ntop(prices, ncol(prices))
# Transaction cost: 4bps
cost <- 0.05/100


# Compare estimation results of PPP without penalty and various PPP with AENET.
fit.lasso <- pppaenet(charlist.tf.vp, BenchW, prices, cost, gamma = 20, alpha = 0, 
                      startdate = startdate, enddate = enddate, trainDataRatio = c(0.8,0.2), 
                      lambda1 = 0.0000001,lambda2 = 0)

fit.ridge <- pppaenet(charlist.tf.vp, BenchW, prices, cost, gamma = 20, alpha = 0, 
                      startdate = startdate, enddate = enddate, trainDataRatio = c(0.8,0.2), 
                      lambda1 = 0,lambda2 = 0.0000001)

fit.enet <- pppaenet(charlist.tf.vp, BenchW, prices, cost, gamma = 20, alpha = 0, 
                     startdate = startdate, enddate = enddate, trainDataRatio = c(0.8,0.2), 
                     lambda1 = 0.0000001,lambda2 = 0.0000001)

fit.aenet <- pppaenet(charlist.tf.vp, BenchW, prices, cost, gamma = 20, alpha = 1, 
                      startdate = startdate, enddate = enddate, trainDataRatio = c(0.8,0.2), 
                      lambda1 = 0.0000001,lambda2 = 0.0000001)

# Get estimated coefficients by accessing the slotname
outMat <- cbind(t(fit.lasso@coefs),t(fit.ridge@coefs),t(fit.enet@coefs),t(fit.aenet@coefs))
colnames(outMat) <- c('Lasso','Ridge','ENET','AENET')
rownames(outMat) <- names(fit.lasso@data@signallist)

print(outMat)

# Estimate the path
ppp.lasso <- pppaenet(charlist.tf.vp, BenchW, prices, cost, gamma = 20, alpha = 0, 
                      startdate = startdate, enddate = enddate,trainDataRatio = c(0.8,0.2), 
                      lambda1 = NULL, minRatio = 1e-10, nLambda1 = 100, lambda2 = 0.0000001)

# Plot the path using plot.coef()
plot.coef(ppp.lasso)

# Clear workspace
rm(list=ls()) 

library(TTR)
library(PerformanceAnalytics)
library(quantmod)
library(caTools)
library(caret)
library(reshape)
library(ggplot2)
library(pbapply)
library(lbfgs)

# files location
filelocation <- "D:/PPPrcodes/"
# Load C++ functions
# cfunfile <- paste(filelocation, 'cfuns.cpp', sep='')


# Load systematic investor toolbox
con = gzcon(file(paste(filelocation, 'sit.gz', sep=''), 'rb'))
source(con)
close(con)

source(paste(filelocation, 'utilityFunctions.R', sep=''))
source(paste(filelocation, 'signalFunctions.R', sep=''))
source(paste(filelocation, 'estimation.R', sep=''))
source(paste(filelocation, 'backtest.R', sep=''))

# Setup start/end dates for estimation
startdate <-'1995-01-01'
enddate <-'2017-05-11'
startdate2<-startdate
enddate2<-enddate
# Create dataset ia (a list contains prices,ret.vol)
ia <- test.create.historical.data(filelocation = filelocation, dat.from = "1992::")

# Compute Trend signals

dates <- index(ia$prices)
if(is.character(startdate)) startdate <- which(dates ==  dates[dates>=startdate][1])
if(is.character(enddate)) enddate <- which(dates ==  dates[dates<=enddate][length(dates[dates<=enddate])]) 

charlist.tf <- create.timeseries.signals(ia)

CovList <- getRollingCov(ia$rets,lookback.len = 90, from = startdate)

# Bechmark weights: 1/N
BenchW <- ntop(ia$prices, ncol(ia$prices))
# Transaction cost: 5bps

frequency = 'weeks'
# find period ends, can be 'weeks', 'months', 'quarters', 'years'
period.ends = endpoints(ia$prices, frequency)
period.ends = period.ends[period.ends > 0]
prices <- ia$prices[period.ends,]
BenchW <- BenchW[period.ends,]
rets <- ia$rets[period.ends,]
rets.vol <- ia$rets.vol[period.ends,]
CovList <- CovList[period.ends]


dates2 <- index(prices)
if(is.character(startdate2)) startdate2 <- which(dates2 ==  dates2[dates2>=startdate2][1])
if(is.character(enddate2)) enddate2 <- which(dates2 ==  dates2[dates2<=enddate2][length(dates2[dates2<=enddate2])]) 


charlist.tf <- setupListRange(charlist.tf, period.ends, outtype = "xts")

charlist.tf.vp <- volatilityWeighted(charlist.tf, rets.vol)

charlist.tf <- c(charlist.tf, charlist.tf.vp)


BenchW[] <- 0

cost <- 0.10/100

lamdda2 <- c(1,1e-03,1e-05,1e-07,0)



# ppp.cv.tf <- crossvalidate(charlist.tf, BenchW, prices, cost, gamma = 50, alpha = 1, 
#                             startdate = startdate2, enddate = enddate2,trainDataRatio = c(0.8,0.2),
#                           type = "OMG.out",nfolds = 10,lambda1 = NULL, minRatio = 1e-10, 
#                         nLambda1 = 100, lambda2 = lamdda2, optcost = "EXACT_COST",name="Trend Vol Parity",show = TRUE)

# ppp.bt.tf.insample <- backtest(ppp.cv.tf,initialWindow = 265, horizon = 1, skip = 0,fixedWindow = FALSE, oos = FALSE, show = TRUE)

# ppp.bt.tf.outsample <- backtest(ppp.cv.tf,initialWindow = NULL, horizon = 1, skip = 0,fixedWindow = FALSE, oos = TRUE, show = TRUE)




# ppp.cv.tf.nocost <- crossvalidate(charlist.tf, BenchW, prices, cost, gamma = 50, alpha = 1, 
#                                    startdate = startdate2, enddate = enddate2,trainDataRatio = c(0.8,0.2),
#                                    type = "OMG.out",nfolds = 10,lambda1 = NULL, minRatio = 1e-10, 
#                                    nLambda1 = 100, lambda2 = lamdda2, optcost = "NO_COST",name="Trend Vol Parity",show = TRUE)
# 


# ppp.bt.tf.nocost.insample <- backtest(ppp.cv.tf.nocost,initialWindow = 265, horizon = 1, skip = 0,fixedWindow = FALSE, oos = FALSE, show = TRUE)

# ppp.bt.tf.nocost.outsample <- backtest(ppp.cv.tf.nocost,initialWindow = NULL, horizon = 1, skip = 0,fixedWindow = FALSE, oos = TRUE, show = TRUE)




# ppp.cv.tf.extracost <- crossvalidate(charlist.tf, BenchW, prices, cost, gamma = 50, alpha = 1, 
#                                      startdate = startdate2, enddate = enddate2,trainDataRatio = c(0.8,0.2),
#                                      type = "OMG.out",nfolds = 10,lambda1 = NULL, minRatio = 1e-10, 
#                                      nLambda1 = 100, lambda2 = lamdda2, optcost = "EXTRA_COST",name="Trend Vol Parity",show = TRUE)



# ppp.bt.tf.extracost.insample <- backtest(ppp.cv.tf.extracost,initialWindow = 265, horizon = 1, skip = 0,fixedWindow = FALSE, oos = FALSE, show = TRUE)

# ppp.bt.tf.extracost.outsample <- backtest(ppp.cv.tf.extracost,initialWindow = NULL, horizon = 1, skip = 0,fixedWindow = FALSE, oos = TRUE, show = TRUE)


load(file=paste(filelocation, 'myCrossValidation01.RData', sep=''))

ppp.tf.weights <- ppp.bt.tf.insample@ppp$orgweight
ppp.tf.nocost.weights <- ppp.bt.tf.nocost.insample@ppp$orgweight
ppp.tf.extracost.weights <- ppp.bt.tf.extracost.insample@ppp$orgweight



teststartdate <- ppp.bt.tf.insample@teststartdate
testenddate <- ppp.bt.tf.insample@testenddate

commission <- list(cps = 0.00, fixed =0, percentage =  cost)
data.bt <- preparedata(ia$prices[paste(ppp.bt.tf.insample@test.from,"::",ppp.bt.tf.insample@test.to,sep="")])


data.bt$weight[] = NA
data.bt$weight[index(ppp.tf.weights),] = ppp.tf.weights
ppp.tf.org.insample = bt.run.share(data.bt, commission = commission, silent=T)

data.bt$weight[] = NA
data.bt$weight[index(ppp.tf.nocost.weights),] = ppp.tf.nocost.weights
ppp.tf.nocost.org.insample = bt.run.share(data.bt, commission = commission, silent=T)

data.bt$weight[] = NA
data.bt$weight[index(ppp.tf.extracost.weights),] = ppp.tf.extracost.weights
ppp.tf.extracost.org.insample = bt.run.share(data.bt, commission = commission, silent=T)

models <- list()

models$TF.EXACT.COST <- ppp.tf.org.insample
models$TF.NO.COST <- ppp.tf.nocost.org.insample
models$TF.EXTRA.COST <- ppp.tf.extracost.org.insample
X11()
strategy.performance.snapshoot(models,one.page = T, title = "In-sample performance: PPPAENET Trend Following")



#######  out of sample

ppp.tf.weights <- ppp.bt.tf.outsample@ppp$orgweight

ppp.tf.nocost.weights <- ppp.bt.tf.nocost.outsample@ppp$orgweight
ppp.tf.extracost.weights <- ppp.bt.tf.extracost.outsample@ppp$orgweight



teststartdate <- ppp.bt.tf.outsample@teststartdate
testenddate <- ppp.bt.tf.outsample@testenddate

commission <- list(cps = 0.00, fixed =0, percentage =  cost)
data.bt <- preparedata(ia$prices[paste(ppp.bt.tf.outsample@test.from,"::",ppp.bt.tf.outsample@test.to,sep="")])

data.bt$weight[] = NA
data.bt$weight[index(ppp.tf.weights),] = ppp.tf.weights
ppp.tf.org.outsample = bt.run.share(data.bt, commission = commission, silent=T)

data.bt$weight[] = NA
data.bt$weight[index(ppp.tf.nocost.weights),] = ppp.tf.nocost.weights
ppp.tf.nocost.org.outsample = bt.run.share(data.bt, commission = commission, silent=T)

data.bt$weight[] = NA
data.bt$weight[index(ppp.tf.extracost.weights),] = ppp.tf.extracost.weights
ppp.tf.extracost.org.outsample = bt.run.share(data.bt, commission = commission, silent=T)

models <- list()
models$TF.EXACT.COST <- ppp.tf.org.outsample
models$TF.NO.COST <- ppp.tf.nocost.org.outsample
models$TF.EXTRA.COST <- ppp.tf.extracost.org.outsample
X11()
strategy.performance.snapshoot(models,one.page = T, title = "Out-of-sample performance: PPPAENET Trend Following")


weight.tf.tgt <- target.Weighted.strategy2(ppp.tf.weights, CovList[teststartdate:testenddate],targetvol = 10/100, max.portfolio.leverage = 6,from = 1)
data.bt$weight[] = NA
data.bt$weight[index(weight.tf.tgt),] = weight.tf.tgt
ppp.tf.outsample = bt.run.share(data.bt, commission = commission, silent=T)



weight.tf.nocost.tgt <- target.Weighted.strategy2(ppp.tf.nocost.weights, CovList[teststartdate:testenddate],targetvol = 10/100, max.portfolio.leverage = 6,from = 1)
data.bt$weight[] = NA
data.bt$weight[index(weight.tf.nocost.tgt),] = weight.tf.nocost.tgt
ppp.tf.nocost.outsample = bt.run.share(data.bt, commission = commission, silent=T)

weight.tf.extracost.tgt <- target.Weighted.strategy2(ppp.tf.extracost.weights, CovList[teststartdate:testenddate],targetvol = 10/100, max.portfolio.leverage = 6,from = 1)
data.bt$weight[] = NA
data.bt$weight[index(weight.tf.extracost.tgt),] = weight.tf.extracost.tgt
ppp.tf.extracost.outsample = bt.run.share(data.bt, commission = commission, silent=T)



bt.SGTrendIndex.outsample <- get.bt.benchmarkIndex(ia$SGTrendIndex.ret, index(ppp.tf.outsample$ret))
bt.SGCTAIndex.outsample <- get.bt.benchmarkIndex(ia$SGCTAIndex.ret, index(ppp.tf.outsample$ret))


models <- list()

#models$SG.CTA.Index <- bt.SGCTAIndex.outsample
models$SG.Trend.Index <- bt.SGTrendIndex.outsample

models$PPP.TF.NOCOST<-ppp.tf.nocost.outsample
models$PPP.TF.EXTRACOST<-ppp.tf.extracost.outsample
models$PPP.TF.EXACTCOST<-ppp.tf.outsample
x11()
plotbt2(models,tilttxt = "Out-of-sample performance: PPPAENET Trend Following v.s. CTA Index")



models <- list()
models$SG.Trend.Index <- bt.SGTrendIndex.outsample
models$PPP.TF.EXACTCOST <- ppp.tf.outsample
x11()

plotbt.custom.report.part1.2(models,trade.summary=T,
                             main = 'Out-of-sample: PPPAENET Trend Following v.s. CTA Index')

x11()
plotbt.transition.map(ppp.tf.outsample$weight,
                      name = 'PPPAENET Trend Following')

x11()
out <- plotbt.monthly.table(ppp.tf.outsample$equity,
                            smain = 'PPPAENET Trend Following Monthly Performance')


x11()

plotbt.holdings.time(ppp.tf.outsample$weight,
                     smain='PPPAENET Trend Following')


x11()
plot.coef(ppp.bt.tf.extracost.insample,tilttxt = "In-sample: Coefficient estimates EXTRA COST", coef.index = c(1:11,55:65))

x11()
plot.coef(ppp.bt.tf.nocost.insample,tilttxt = "In-sample: Coefficient estimates NO COST", coef.index = c(1:11,55:65))

x11()
plot.coef(ppp.bt.tf.insample,tilttxt = "In-sample: Coefficient estimates EXACT COST", coef.index = c(1:11,55:65))

x11()
plot.coef(ppp.bt.tf.extracost.outsample,tilttxt = "Out-of-sample: Coefficient estimates EXTRA COST", coef.index = c(1:11,55:65))

x11()
plot.coef(ppp.bt.tf.nocost.outsample,tilttxt = "Out-of-sample: Coefficient estimates NO COST", coef.index = c(1:11,55:65))

x11()
plot.coef(ppp.bt.tf.outsample,tilttxt = "Out-of-sample: Coefficient estimates EXACT COST", coef.index = c(1:11,55:65))

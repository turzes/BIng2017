###############################################################################
# Load Systematic Investor Toolbox (SIT)
# http://systematicinvestor.wordpress.com/systematic-investor-toolbox/
###############################################################################
#setInternet2(TRUE)
con = gzcon(url('http://www.systematicportfolio.com/sit.gz', 'rb'))
con = gzcon(file('c:/temp/sit.gz', 'rb'))
require(RCurl)
sit = getURLContent('https://github.com/systematicinvestor/SIT/raw/master/sit.gz', binary=TRUE, followlocation = TRUE, ssl.verifypeer = FALSE)
con = gzcon(rawConnection(sit, 'rb'))

source(con)
close(con)
#*****************************************************************
# Load historical data
#******************************************************************
load.packages('quantmod')

tickers = spl('SPY,TLT')

data <- new.env()
getSymbols(tickers, src = 'yahoo', from = '1980-01-01', env = data, auto.assign = T)
for(i in ls(data)) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align='remove.na', dates='2005::')


#*****************************************************************
# Setup
#******************************************************************
lookback.len = 60

prices = data$prices

models = list()

#*****************************************************************
# Simple Momentum
#******************************************************************
momentum = prices / mlag(prices, lookback.len)
data$weight[] = NA
data$weight$SPY[] = momentum$SPY > momentum$TLT
data$weight$TLT[] = momentum$SPY <= momentum$TLT
models$Simple  = bt.run.share(data, clean.signal=T)     

#*****************************************************************
# Probabilistic Momentum
#******************************************************************
confidence.level = 60/100
ret = prices / mlag(prices) - 1

ir = sqrt(lookback.len) * runMean(ret$SPY - ret$TLT, lookback.len) / runSD(ret$SPY - ret$TLT, lookback.len)
momentum.p = pt(ir, lookback.len - 1)

data$weight[] = NA
data$weight$SPY[] = iif(cross.up(momentum.p, confidence.level), 1, iif(cross.dn(momentum.p, (1 - confidence.level)), 0,NA))
data$weight$TLT[] = iif(cross.dn(momentum.p, (1 - confidence.level)), 1, iif(cross.up(momentum.p, confidence.level), 0,NA))
models$Probabilistic  = bt.run.share(data, clean.signal=T)  

#*****************************************************************
# Probabilistic Momentum + SPY Leverage
#******************************************************************
data$weight[] = NA
data$weight$SPY[] = iif(cross.up(momentum.p, confidence.level), 1, iif(cross.up(momentum.p, (1 - confidence.level)), 0,NA))
data$weight$TLT[] = iif(cross.dn(momentum.p, (1 - confidence.level)), 1, iif(cross.up(momentum.p, confidence.level), 0,NA))
models$Probabilistic.Leverage = bt.run.share(data, clean.signal=T) 

#*****************************************************************
# Create Report
#******************************************************************   
strategy.performance.snapshoot(models, T)

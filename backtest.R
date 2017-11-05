
Sys.setenv(TZ='Europe/London')

create.historical.data <- function(datalist, prices,rets, rets.vol){
  ia <- list()
  ia$datalist <- datalist
  ia$prices <- prices
  ia$rets <- rets
  ia$rets.vol <- rets.vol
  return(ia)
}

test.create.historical.data <- function(filelocation = NULL, dat.from = NULL){
  if(is.null(filelocation)) filelocation <- "D:/Rcode/Rcpp/"
  if(is.null(dat.from)) {
    dat.from <- "2008::"
    
  }
  
  # load data
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
  
  prices <- prices[dat.from]
  idx <- apply(!is.na(prices), 2, all)
  
  prices <- prices[,idx]
  data <- setupListRange(dataq,dat.from, outtype = "xts")
  data <- data[idx]
  
  # Returns
  rets <- prices/mlag(prices) - 1
  # Returns volatiliy
  # rets.vol <- bt.apply.matrix(rets, runsd, k = 21)*sqrt(252)
  rets.vol <- bt.apply.matrix(rets, runsd, k = 60)*sqrt(252)
  
  # vol.tmp <- sapply(data, function(x){
  #   return(volatility(x, calc = "parkinson", n = 60, N = 252, mean0 = TRUE))
  # })
  
  # rets.vol<- xts(vol.tmp, order.by = index(prices))
  
  # rets.vol[is.na(rets.vol)] = 0
  #n <- 90
  #rets.vol <- sqrt(252)*sqrt(bt.apply.matrix(rets^2, runSum, n - 1))/(n - 2)
  #rets.vol <- sqrt(252)*bt.apply.matrix(abs(rets), runMean, n - 1)

  ia <- create.historical.data(data, prices, rets, rets.vol)

  SGTrendIndex <- read.csv(paste(filelocation, 'Trend_Index_Historical.csv', sep=''))
  SGTrendIndex <- xts(SGTrendIndex[,-1], order.by = as.Date(SGTrendIndex[, 1],"%d/%m/%Y"))
  SGTrendIndex <- SGTrendIndex[,c(1,2)]
  colnames(SGTrendIndex) <- c("price","rets")
  
  SGTrendIndex <- SGTrendIndex[index(ia$prices),]
  SGTrendIndex.ret <- SGTrendIndex[,"rets"]
  
  SGCTAIndex <- read.csv(paste(filelocation, 'CTA_Historical.csv', sep=''))
  SGCTAIndex <- xts(SGCTAIndex[,-1], order.by = as.Date(SGCTAIndex[, 1],"%d/%m/%Y"))
  SGCTAIndex <- SGCTAIndex[,c(1,2)]
  colnames(SGCTAIndex) <- c("price","rets")
  
  SGCTAIndex <- SGCTAIndex[index(ia$prices),]
  SGCTAIndex.ret <- SGCTAIndex[,"rets"]
  
  ia$SGCTAIndex.ret <- SGCTAIndex.ret
  ia$SGTrendIndex.ret <- SGTrendIndex.ret
  
  return(ia)
}

create.marco.signals <- function(ia, filelocation = NULL,from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  macrofactors <- read.csv(paste(filelocation, 'DailyMacroFactors_zscore.csv', sep=''))
  macrofactors <- xts(macrofactors[,-1], order.by = as.Date(macrofactors[, 1],"%Y-%m-%d"))
  macrofactors[is.na(macrofactors)] <- 0
  #macrofactors <- read.csv(paste(filelocation, 'dailyMacroFactors_PrincipalComponents.csv', sep=''))
  #macrofactors <- xts(macrofactors[,-1], order.by = as.Date(macrofactors[, 1],"%Y-%m-%d"))
  
  charlist.macro <- vector("list",ncol(macrofactors))
  for(i in 1:ncol(macrofactors)){
    
    charlist.macro[[i]]<- getMacroFactors(macrofactors,index(prices),names(macrofactors)[i],ncol(prices))
  }
  names(charlist.macro) <- names(macrofactors)
  
  return(charlist.macro)
}



create.MOM.ts.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  # PPP: Momentum
  # Signal: 1 Months Returns
  
  momD05 <- getMomentum(prices, nmonth = 0.25, type = 'value', normalize = TRUE, normtype = 'ts')
  momD10 <- getMomentum(prices, nmonth = 0.5, type = 'value', normalize = TRUE, normtype = 'ts')
  
  
  mom01 <- getMomentum(prices, nmonth = 1, type = 'value', normalize = TRUE, normtype = 'ts')
  mom03 <- getMomentum(prices, nmonth = 3, type = 'value', normalize = TRUE, normtype = 'ts')
  mom06 <- getMomentum(prices, nmonth = 6, type = 'value', normalize = TRUE, normtype = 'ts')
  mom09 <- getMomentum(prices, nmonth = 9, type = 'value', normalize = TRUE, normtype = 'ts')
  mom12 <- getMomentum(prices, nmonth = 12, type = 'value', normalize = TRUE, normtype = 'ts')
  mom18 <- getMomentum(prices, nmonth = 18, type = 'value', normalize = TRUE, normtype = 'ts')
  mom24 <- getMomentum(prices, nmonth = 24, type = 'value', normalize = TRUE, normtype = 'ts')
  mom36 <- getMomentum(prices, nmonth = 36, type = 'value', normalize = TRUE, normtype = 'ts')
  mom48 <- getMomentum(prices, nmonth = 48, type = 'value', normalize = TRUE, normtype = 'ts')
  mom54 <- getMomentum(prices, nmonth = 54, type = 'value', normalize = TRUE, normtype = 'ts')
  mom60 <- getMomentum(prices, nmonth = 60, type = 'value', normalize = TRUE, normtype = 'ts')
  
  charlist.mom.ts <- list(momD05=momD05,momD10=momD10,mom01=mom01,mom03=mom03,mom06=mom06, mom09=mom09,mom12=mom12,mom24=mom24,mom48=mom48,mom54=mom54,mom60=mom60)
  
  
  
  return(charlist.mom.ts)
}

create.EMA1.ts.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  ema1.10<- getEMA1(prices,n = 10, type = 'value', normalize = T, normtype = 'ts')
  ema1.20<- getEMA1(prices,n = 20, type = 'value', normalize = T, normtype = 'ts')
  ema1.50 <- getEMA1(prices,n = 50, type = 'value', normalize = T, normtype = 'ts')
  ema1.100 <- getEMA1(prices,n = 100, type = 'value', normalize = T, normtype = 'ts')
  ema1.200 <- getEMA1(prices,n = 200, type = 'value', normalize = T, normtype = 'ts')
  ema1.300 <- getEMA1(prices,n = 300, type = 'value', normalize = T, normtype = 'ts')
  ema1.400 <- getEMA1(prices,n = 400, type = 'value', normalize = T, normtype = 'ts')
  ema1.600 <- getEMA1(prices,n = 600, type = 'value', normalize = T, normtype = 'ts')
  ema1.800 <- getEMA1(prices,n = 800, type = 'value', normalize = T, normtype = 'ts')
  
  
  charlist.ema1 <- list(ema1.10=ema1.10,ema1.20=ema1.20,ema1.50=ema1.50, ema1.100=ema1.100, ema1.200=ema1.200,ema1.300=ema1.300, ema1.400=ema1.400,ema1.600=ema1.600,ema1.800=ema1.800)
  
  return(charlist.ema1)
}

create.EMA2.ts.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  
  ema2.10.20 <- getEMA2(prices, n1 = 10, n2 = 20, type = 'value', normalize = T, normtype = 'ts')
  ema2.20.50 <- getEMA2(prices, n1 = 20, n2 = 50, type = 'value', normalize = T, normtype = 'ts')
  ema2.40.80 <- getEMA2(prices, n1 = 40, n2 = 80, type = 'value', normalize = T, normtype = 'ts')
  ema2.50.200 <- getEMA2(prices, n1 = 50, n2 = 200, type = 'value', normalize = T, normtype = 'ts')
  ema2.100.300 <- getEMA2(prices, n1 = 100, n2 = 300, type = 'value', normalize = T, normtype = 'ts')
  ema2.200.400 <- getEMA2(prices, n1 = 200, n2 = 400, type = 'value', normalize = T, normtype = 'ts')
  ema2.300.500<- getEMA2(prices, n1 = 300, n2 = 500, type = 'value', normalize = T, normtype = 'ts')
  ema2.400.800 <- getEMA2(prices, n1 = 400, n2 = 800, type = 'value', normalize = T, normtype = 'ts')
  
  
  charlist.ema2 <- list(ema2.10.20=ema2.10.20,ema2.20.50=ema2.20.50,ema2.40.80=ema2.40.80,ema2.50.200=ema2.50.200,ema2.100.300=ema2.100.300, ema2.200.400=ema2.200.400, ema2.300.500=ema2.300.500,ema2.400.800=ema2.400.800)
  
  return(charlist.ema2)
}


create.RSI.ts.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  
  rsi.10 <- getRSI(prices, n = 10, type = 'value', normalize = T, normtype = 'ts')
  rsi.20 <- getRSI(prices, n = 20, type = 'value', normalize = T, normtype = 'ts')
  rsi.50 <- getRSI(prices, n = 50, type = 'value', normalize = T, normtype = 'ts')
  rsi.100 <- getRSI(prices, n = 100, type = 'value', normalize = T, normtype = 'ts')
  rsi.200 <- getRSI(prices, n = 200, type = 'value', normalize = T, normtype = 'ts')
  rsi.300 <- getRSI(prices, n = 300, type = 'value', normalize = T, normtype = 'ts')
  rsi.400 <- getRSI(prices, n = 400, type = 'value', normalize = T, normtype = 'ts')
  rsi.600 <- getRSI(prices, n = 600, type = 'value', normalize = T, normtype = 'ts')
  rsi.800 <- getRSI(prices, n = 800, type = 'value', normalize = T, normtype = 'ts')
  
  charlist.rsi <- list(rsi.10=rsi.10,rsi.20=rsi.20,rsi.50=rsi.50, rsi.100=rsi.100, rsi.200=rsi.200,rsi.300=rsi.300, rsi.400=rsi.400,rsi.600=rsi.600,rsi.800=rsi.800)
  return(charlist.rsi)
}


create.TRIX.ts.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  
  
  
  trix.10 <- getTRIX(prices, n = 20, nSig = 10, type = 'value', normalize = T, normtype = 'ts')
  trix.20 <- getTRIX(prices, n = 20, nSig = 10, type = 'value', normalize = T, normtype = 'ts')
  trix.50 <- getTRIX(prices, n = 50, nSig = 10, type = 'value', normalize = T, normtype = 'ts')
  trix.100 <- getTRIX(prices, n = 100, nSig = 10, type = 'value', normalize = T, normtype = 'ts')
  trix.200 <- getTRIX(prices, n = 200, nSig = 10, type = 'value', normalize = T, normtype = 'ts')
  trix.300 <- getTRIX(prices, n = 300, nSig = 10, type = 'value', normalize = T, normtype = 'ts')
  trix.400 <- getTRIX(prices, n = 400, nSig = 10, type = 'value', normalize = T, normtype = 'ts')
  trix.600 <- getTRIX(prices, n = 600, nSig = 10, type = 'value', normalize = T, normtype = 'ts')
  trix.800 <- getTRIX(prices, n = 800, nSig = 10, type = 'value', normalize = T, normtype = 'ts')
  
  
  
  charlist.trix <- list(trix.10=trix.10,trix.20=trix.20,trix.50=trix.50, trix.100=trix.100, trix.200=trix.200,trix.300=trix.300, trix.400=trix.400,trix.600=trix.600,trix.800=trix.800)
  
  return(charlist.trix)
}

create.MACD.ts.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  
  
  # MACD STOCH
  macd.10.20 <- getMACD(prices, nFast = 10, nSlow = 20, nSig = 9, type = 'value', normalize = T, normtype = 'ts')
  macd.20.50 <- getMACD(prices, nFast = 20, nSlow = 50, nSig = 9, type = 'value', normalize = T, normtype = 'ts')
  macd.40.80 <- getMACD(prices, nFast = 40, nSlow = 80, nSig = 9, type = 'value', normalize = T, normtype = 'ts')
  macd.50.200 <- getMACD(prices, nFast = 50, nSlow = 200, nSig = 9, type = 'value', normalize = T, normtype = 'ts')
  macd.100.300 <- getMACD(prices, nFast = 100, nSlow = 300, nSig = 9, type = 'value', normalize = T, normtype = 'ts')
  macd.200.400 <- getMACD(prices, nFast = 200, nSlow = 400, nSig = 9, type = 'value', normalize = T, normtype = 'ts')
  macd.300.500 <- getMACD(prices, nFast = 300, nSlow = 500, nSig = 9, type = 'value', normalize = T, normtype = 'ts')
  macd.400.800 <- getMACD(prices, nFast = 400, nSlow = 800, nSig = 9, type = 'value', normalize = T, normtype = 'ts')
  
  
  
  charlist.macd <- list(macd.10.20=macd.10.20,macd.20.50=macd.20.50,macd.40.80=macd.40.80,macd.50.200=macd.50.200,macd.100.300=macd.100.300,macd.200.400=macd.200.400, macd.300.500=macd.300.500,macd.400.800=macd.400.800)
  
  return(charlist.macd)
}


create.MOM.xs.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  # PPP: Momentum
  # Signal: 1 Months Returns
  
  momD05 <- getMomentum(prices, nmonth = 0.25, type = 'value', normalize = TRUE, normtype = 'xs')
  momD10 <- getMomentum(prices, nmonth = 0.5, type = 'value', normalize = TRUE, normtype = 'xs')
  
  
  mom01 <- getMomentum(prices, nmonth = 1, type = 'value', normalize = TRUE, normtype = 'xs')
  mom03 <- getMomentum(prices, nmonth = 3, type = 'value', normalize = TRUE, normtype = 'xs')
  mom06 <- getMomentum(prices, nmonth = 6, type = 'value', normalize = TRUE, normtype = 'xs')
  mom09 <- getMomentum(prices, nmonth = 9, type = 'value', normalize = TRUE, normtype = 'xs')
  mom12 <- getMomentum(prices, nmonth = 12, type = 'value', normalize = TRUE, normtype = 'xs')
  mom18 <- getMomentum(prices, nmonth = 18, type = 'value', normalize = TRUE, normtype = 'xs')
  mom24 <- getMomentum(prices, nmonth = 24, type = 'value', normalize = TRUE, normtype = 'xs')
  mom36 <- getMomentum(prices, nmonth = 36, type = 'value', normalize = TRUE, normtype = 'xs')
  mom48 <- getMomentum(prices, nmonth = 48, type = 'value', normalize = TRUE, normtype = 'xs')
  mom54 <- getMomentum(prices, nmonth = 54, type = 'value', normalize = TRUE, normtype = 'xs')
  mom60 <- getMomentum(prices, nmonth = 60, type = 'value', normalize = TRUE, normtype = 'xs')
  
  charlist.mom.xs <- list(momD05=momD05,momD10=momD10,mom01=mom01,mom03=mom03,mom06=mom06, mom09=mom09,mom12=mom12,mom24=mom24,mom48=mom48,mom54=mom54,mom60=mom60)
  
  
  
  return(charlist.mom.xs)
}

create.EMA1.xs.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  ema1.10<- getEMA1(prices,n = 10, type = 'value', normalize = T, normtype = 'xs')
  ema1.20<- getEMA1(prices,n = 20, type = 'value', normalize = T, normtype = 'xs')
  ema1.50 <- getEMA1(prices,n = 50, type = 'value', normalize = T, normtype = 'xs')
  ema1.100 <- getEMA1(prices,n = 100, type = 'value', normalize = T, normtype = 'xs')
  ema1.200 <- getEMA1(prices,n = 200, type = 'value', normalize = T, normtype = 'xs')
  ema1.300 <- getEMA1(prices,n = 300, type = 'value', normalize = T, normtype = 'xs')
  ema1.400 <- getEMA1(prices,n = 400, type = 'value', normalize = T, normtype = 'xs')
  ema1.600 <- getEMA1(prices,n = 600, type = 'value', normalize = T, normtype = 'xs')
  ema1.800 <- getEMA1(prices,n = 800, type = 'value', normalize = T, normtype = 'xs')
  
  
  charlist.ema1 <- list(ema1.10=ema1.10,ema1.20=ema1.20,ema1.50=ema1.50, ema1.100=ema1.100, ema1.200=ema1.200,ema1.300=ema1.300, ema1.400=ema1.400,ema1.600=ema1.600,ema1.800=ema1.800)
  
  return(charlist.ema1)
}

create.EMA2.xs.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  
  ema2.10.20 <- getEMA2(prices, n1 = 10, n2 = 20, type = 'value', normalize = T, normtype = 'xs')
  ema2.20.50 <- getEMA2(prices, n1 = 20, n2 = 50, type = 'value', normalize = T, normtype = 'xs')
  ema2.40.80 <- getEMA2(prices, n1 = 40, n2 = 80, type = 'value', normalize = T, normtype = 'xs')
  ema2.50.200 <- getEMA2(prices, n1 = 50, n2 = 200, type = 'value', normalize = T, normtype = 'xs')
  ema2.100.300 <- getEMA2(prices, n1 = 100, n2 = 300, type = 'value', normalize = T, normtype = 'xs')
  ema2.200.400 <- getEMA2(prices, n1 = 200, n2 = 400, type = 'value', normalize = T, normtype = 'xs')
  ema2.300.500<- getEMA2(prices, n1 = 300, n2 = 500, type = 'value', normalize = T, normtype = 'xs')
  ema2.400.800 <- getEMA2(prices, n1 = 400, n2 = 800, type = 'value', normalize = T, normtype = 'xs')
  
  
  charlist.ema2 <- list(ema2.10.20=ema2.10.20,ema2.20.50=ema2.20.50,ema2.40.80=ema2.40.80,ema2.50.200=ema2.50.200,ema2.100.300=ema2.100.300, ema2.200.400=ema2.200.400, ema2.300.500=ema2.300.500,ema2.400.800=ema2.400.800)
  
  return(charlist.ema2)
}


create.RSI.xs.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  
  rsi.10 <- getRSI(prices, n = 10, type = 'value', normalize = T, normtype = 'xs')
  rsi.20 <- getRSI(prices, n = 20, type = 'value', normalize = T, normtype = 'xs')
  rsi.50 <- getRSI(prices, n = 50, type = 'value', normalize = T, normtype = 'xs')
  rsi.100 <- getRSI(prices, n = 100, type = 'value', normalize = T, normtype = 'xs')
  rsi.200 <- getRSI(prices, n = 200, type = 'value', normalize = T, normtype = 'xs')
  rsi.300 <- getRSI(prices, n = 300, type = 'value', normalize = T, normtype = 'xs')
  rsi.400 <- getRSI(prices, n = 400, type = 'value', normalize = T, normtype = 'xs')
  rsi.600 <- getRSI(prices, n = 600, type = 'value', normalize = T, normtype = 'xs')
  rsi.800 <- getRSI(prices, n = 800, type = 'value', normalize = T, normtype = 'xs')
  
  charlist.rsi <- list(rsi.10=rsi.10,rsi.20=rsi.20,rsi.50=rsi.50, rsi.100=rsi.100, rsi.200=rsi.200,rsi.300=rsi.300, rsi.400=rsi.400,rsi.600=rsi.600,rsi.800=rsi.800)
  return(charlist.rsi)
}


create.TRIX.xs.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  
  
  
  trix.10 <- getTRIX(prices, n = 20, nSig = 10, type = 'value', normalize = T, normtype = 'xs')
  trix.20 <- getTRIX(prices, n = 20, nSig = 10, type = 'value', normalize = T, normtype = 'xs')
  trix.50 <- getTRIX(prices, n = 50, nSig = 10, type = 'value', normalize = T, normtype = 'xs')
  trix.100 <- getTRIX(prices, n = 100, nSig = 10, type = 'value', normalize = T, normtype = 'xs')
  trix.200 <- getTRIX(prices, n = 200, nSig = 10, type = 'value', normalize = T, normtype = 'xs')
  trix.300 <- getTRIX(prices, n = 300, nSig = 10, type = 'value', normalize = T, normtype = 'xs')
  trix.400 <- getTRIX(prices, n = 400, nSig = 10, type = 'value', normalize = T, normtype = 'xs')
  trix.600 <- getTRIX(prices, n = 600, nSig = 10, type = 'value', normalize = T, normtype = 'xs')
  trix.800 <- getTRIX(prices, n = 800, nSig = 10, type = 'value', normalize = T, normtype = 'xs')
  
  
  
  charlist.trix <- list(trix.10=trix.10,trix.20=trix.20,trix.50=trix.50, trix.100=trix.100, trix.200=trix.200,trix.300=trix.300, trix.400=trix.400,trix.600=trix.600,trix.800=trix.800)
  
  return(charlist.trix)
}

create.MACD.xs.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  
  
  # MACD STOCH
  macd.10.20 <- getMACD(prices, nFast = 10, nSlow = 20, nSig = 9, type = 'value', normalize = T, normtype = 'xs')
  macd.20.50 <- getMACD(prices, nFast = 20, nSlow = 50, nSig = 9, type = 'value', normalize = T, normtype = 'xs')
  macd.40.80 <- getMACD(prices, nFast = 40, nSlow = 80, nSig = 9, type = 'value', normalize = T, normtype = 'xs')
  macd.50.200 <- getMACD(prices, nFast = 50, nSlow = 200, nSig = 9, type = 'value', normalize = T, normtype = 'xs')
  macd.100.300 <- getMACD(prices, nFast = 100, nSlow = 300, nSig = 9, type = 'value', normalize = T, normtype = 'xs')
  macd.200.400 <- getMACD(prices, nFast = 200, nSlow = 400, nSig = 9, type = 'value', normalize = T, normtype = 'xs')
  macd.300.500 <- getMACD(prices, nFast = 300, nSlow = 500, nSig = 9, type = 'value', normalize = T, normtype = 'xs')
  macd.400.800 <- getMACD(prices, nFast = 400, nSlow = 800, nSig = 9, type = 'value', normalize = T, normtype = 'xs')
  
  
  
  charlist.macd <- list(macd.10.20=macd.10.20,macd.20.50=macd.20.50,macd.40.80=macd.40.80,macd.50.200=macd.50.200,macd.100.300=macd.100.300,macd.200.400=macd.200.400, macd.300.500=macd.300.500,macd.400.800=macd.400.800)
  
  return(charlist.macd)
}

create.MOM.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  # PPP: Momentum
  # Signal: 1 Months Returns
  momD05 <- getMomentum(prices, nmonth = 0.25)
  momD10 <- getMomentum(prices, nmonth = 0.5)
  
  
  mom01 <- getMomentum(prices, nmonth = 1)
  mom03 <- getMomentum(prices, nmonth = 3)
  mom06 <- getMomentum(prices, nmonth = 6)
  mom09 <- getMomentum(prices, nmonth = 9)
  mom12 <- getMomentum(prices, nmonth = 12)
  mom18 <- getMomentum(prices, nmonth = 18)
  mom24 <- getMomentum(prices, nmonth = 24)
  mom36 <- getMomentum(prices, nmonth = 36)
  mom48 <- getMomentum(prices, nmonth = 48)
  mom54 <- getMomentum(prices, nmonth = 54)
  mom60 <- getMomentum(prices, nmonth = 60)
  
  
  
  charlist.mom <- list(momD05,momD10,mom01=mom01,mom03=mom03,mom06=mom06, mom09=mom09,mom12=mom12,mom24=mom24,mom48=mom48,mom54=mom54,mom60=mom60)
  
  return(charlist.mom)
}

create.EMA1.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  ema1.10<- getEMA1(prices,n = 10, from = from)
  ema1.20<- getEMA1(prices,n = 20, from = from)
  ema1.50 <- getEMA1(prices,n = 50, from = from)
  ema1.100 <- getEMA1(prices,n = 100, from = from)
  ema1.200 <- getEMA1(prices,n = 200, from = from)
  ema1.300 <- getEMA1(prices,n = 300, from = from)
  ema1.400 <- getEMA1(prices,n = 400, from = from)
  ema1.600 <- getEMA1(prices,n = 600, from = from)
  ema1.800 <- getEMA1(prices,n = 800, from = from)
  
  charlist.ema1 <- list(ema1.10=ema1.10,ema1.20=ema1.20,ema1.50=ema1.50, ema1.100=ema1.100, ema1.200=ema1.200,ema1.300=ema1.300, ema1.400=ema1.400,ema1.600=ema1.600,ema1.800=ema1.800)
  
  return(charlist.ema1)
}

create.EMA2.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  
  ema2.10.20 <- getEMA2(prices, n1 = 10, n2 = 20, from = from)
  ema2.20.50 <- getEMA2(prices, n1 = 20, n2 = 50, from = from)
  ema2.40.80 <- getEMA2(prices, n1 = 40, n2 = 80, from = from)
  ema2.50.200 <- getEMA2(prices, n1 = 50, n2 = 200, from = from)
  ema2.100.300 <- getEMA2(prices, n1 = 100, n2 = 300, from = from)
  ema2.200.400 <- getEMA2(prices, n1 = 200, n2 = 400, from = from)
  ema2.300.500 <- getEMA2(prices, n1 = 300, n2 = 500, from = from)
  ema2.400.800 <- getEMA2(prices, n1 = 400, n2 = 800, from = from)
  
  charlist.ema2 <- list(ema2.10.20,ema2.20.50,ema2.40.80,ema2.50.200,ema2.100.300, ema2.200.400, ema2.300.500,ema2.400.800=ema2.400.800)
  
  return(charlist.ema2)
}

create.RSI.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  
  rsi.10 <- getRSI(prices, n = 10, from = from)
  rsi.20 <- getRSI(prices, n = 20, from = from)
  rsi.50 <- getRSI(prices, n = 50, from = from)
  rsi.100 <- getRSI(prices, n = 100, from = from)
  rsi.200 <- getRSI(prices, n = 200, from = from)
  rsi.300 <- getRSI(prices, n = 300, from = from)
  rsi.400 <- getRSI(prices, n = 400, from = from)
  rsi.600 <- getRSI(prices, n = 600, from = from)
  rsi.800 <- getRSI(prices, n = 800, from = from)
  
  charlist.rsi <- list(rsi.10=rsi.10,rsi.20=rsi.20,rsi.50=rsi.50, rsi.100=rsi.100, rsi.200=rsi.200,rsi.300=rsi.300, rsi.400=rsi.400,rsi.600=rsi.600,rsi.800=rsi.800)
  return(charlist.rsi)
}

create.TRIX.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  
  
  trix.10 <- getTRIX(prices, n = 20, nSig = 10, from = from)
  trix.20 <- getTRIX(prices, n = 20, nSig = 10, from = from)
  trix.50 <- getTRIX(prices, n = 50, nSig = 10, from = from)
  trix.100 <- getTRIX(prices, n = 100, nSig = 10, from = from)
  trix.200 <- getTRIX(prices, n = 200, nSig = 10, from = from)
  trix.300 <- getTRIX(prices, n = 300, nSig = 10, from = from)
  trix.400 <- getTRIX(prices, n = 400, nSig = 10, from = from)
  trix.600 <- getTRIX(prices, n = 600, nSig = 10, from = from)
  trix.800 <- getTRIX(prices, n = 800, nSig = 10, from = from)
  
  charlist.trix <- list(trix.10=trix.10,trix.20=trix.20,trix.50=trix.50, trix.100=trix.100, trix.200=trix.200,trix.300=trix.300, trix.400=trix.400,trix.600=trix.600,trix.800=trix.800)
  
  return(charlist.trix)
}

create.DC.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  
  dc.10 <- getDonchianChannel(data, n = 10, from = from)
  dc.20 <- getDonchianChannel(data, n = 20, from = from)
  dc.50 <- getDonchianChannel(data, n = 50, from = from)
  dc.100 <- getDonchianChannel(data, n = 100, from = from)
  dc.200 <- getDonchianChannel(data, n = 200, from = from)
  dc.300 <- getDonchianChannel(data, n = 300, from = from)
  dc.400 <- getDonchianChannel(data, n = 400, from = from)
  dc.600 <- getDonchianChannel(data, n = 600, from = from)
  dc.800 <- getDonchianChannel(data, n = 800, from = from)
  
  charlist.dc <- list(dc.10=dc.10,dc.20=dc.20,dc.50=dc.50, dc.100=dc.100, dc.200=dc.200,dc.300=dc.300, dc.400=dc.400,dc.600=dc.600,dc.800=dc.800)
  
  return(charlist.dc)
}

create.BBANDS.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  bbands.10 <- getBBands(data, n = 10, from = from)
  bbands.20 <- getBBands(data, n = 20, from = from)
  bbands.50 <- getBBands(data, n = 50, from = from)
  bbands.100 <- getBBands(data, n = 100, from = from)
  bbands.200 <- getBBands(data, n = 200, from = from)
  bbands.300 <- getBBands(data, n = 300, from = from)
  bbands.400 <- getBBands(data, n = 400, from = from)
  bbands.600 <- getBBands(data, n = 600, from = from)
  bbands.800 <- getBBands(data, n = 800, from = from)
  
  charlist.bbands <- list(bbands.10=bbands.10,bbands.20=bbands.20,bbands.50=bbands.50, bbands.100=bbands.100, bbands.200=bbands.200,bbands.300=bbands.300, bbands.400=bbands.400,bbands.600=bbands.600,bbands.800=bbands.800)
  
  return(charlist.bbands)
}


create.OBV.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  
  obv.10 <- getOBV(data, n = 10, from = from)
  obv.20 <- getOBV(data, n = 20, from = from)
  obv.50 <- getOBV(data, n = 50, from = from)
  obv.100 <- getOBV(data, n = 100, from = from)
  obv.200 <- getOBV(data, n = 200, from = from)
  obv.300 <- getOBV(data, n = 300, from = from)
  obv.400 <- getOBV(data, n = 400, from = from)
  obv.600 <- getOBV(data, n = 600, from = from)
  obv.800 <- getOBV(data, n = 800, from = from)
  
  
  
  
  charlist.obv <- list(obv.10=obv.10,obv.20=obv.20,obv.50=obv.50, obv.100=obv.100, obv.200=obv.200,obv.300=obv.300, obv.400=obv.400,obv.600=obv.600,obv.800=obv.800)
  
  return(charlist.obv)
}

create.MACD.signals <- function(ia, from = NULL, to = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  prices <- ia$prices
  data <- ia$datalist
  
  
  
  # MACD STOCH
  macd.10.20 <- getMACD(prices, nFast = 10, nSlow = 20, nSig = 9, from = from)
  macd.20.50 <- getMACD(prices, nFast = 20, nSlow = 50, nSig = 9, from = from)
  macd.40.80 <- getMACD(prices, nFast = 40, nSlow = 80, nSig = 9, from = from)
  macd.50.200 <- getMACD(prices, nFast = 50, nSlow = 200, nSig = 9, from = from)
  macd.100.300 <- getMACD(prices, nFast = 100, nSlow = 300, nSig = 9, from = from)
  macd.200.400 <- getMACD(prices, nFast = 200, nSlow = 400, nSig = 9, from = from)
  macd.300.500 <- getMACD(prices, nFast = 300, nSlow = 500, nSig = 9, from = from)
  macd.400.800 <- getMACD(prices, nFast = 400, nSlow = 800, nSig = 9, from = from)
  
  charlist.macd <- list(macd.10.20=macd.10.20,macd.20.50=macd.20.50,macd.40.80=macd.40.80,macd.50.200=macd.50.200,macd.100.300=macd.100.300,macd.200.400=macd.200.400, macd.300.500=macd.300.500,macd.400.800=macd.400.800)

  return(charlist.macd)
}

create.trendfollowing.signals <- function(ia, from = NULL){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  charlist.mom <- create.MOM.signals(ia)
  charlist.ema1 <- create.EMA1.signals(ia, from = from)
  charlist.rsi <- create.RSI.signals(ia, from = from)
  charlist.trix <- create.TRIX.signals(ia, from = from)
  charlist.dc <- create.DC.signals(ia, from = from)
  charlist.bbands <- create.BBANDS.signals(ia, from = from)
  charlist.obv <- create.OBV.signals(ia, from = from)
  charlist.ema2 <- create.EMA2.signals(ia, from = from)
  charlist.macd <- create.MACD.signals(ia, from = from)
  
  
  
  charlist.tf <- c(charlist.mom,charlist.ema1,charlist.trix,charlist.dc,charlist.obv,charlist.rsi,charlist.bbands,charlist.ema2,charlist.macd)
  
  return(charlist.tf)
}

create.timeseries.signals <- function(ia){
  # Compute Trend signals
  charlist.mom <- create.MOM.ts.signals(ia)
  charlist.ema1 <- create.EMA1.ts.signals(ia)
  charlist.rsi <- create.RSI.ts.signals(ia)
  charlist.trix <- create.TRIX.ts.signals(ia)
  charlist.ema2 <- create.EMA2.ts.signals(ia)
  charlist.macd <- create.MACD.ts.signals(ia)
  
  # Timeseries trend signals (Momentum,EMA1,EMA2,MACD,TRIX,RSI)
  charlist.tf <- c(charlist.mom,charlist.ema1,charlist.ema2,
                   charlist.macd,charlist.rsi,charlist.trix)
  
  return(charlist.tf)
}

create.crosssectional.signals <- function(ia){
  # Compute Trend signals
  charlist.mom <- create.MOM.xs.signals(ia)
  charlist.ema1 <- create.EMA1.xs.signals(ia)
  charlist.rsi <- create.RSI.xs.signals(ia)
  charlist.trix <- create.TRIX.xs.signals(ia)
  charlist.ema2 <- create.EMA2.xs.signals(ia)
  charlist.macd <- create.MACD.xs.signals(ia)
  
  # Timeseries trend signals (Momentum,EMA1,EMA2,MACD,TRIX,RSI)
  charlist.tf <- c(charlist.mom,charlist.ema1,charlist.ema2,
                   charlist.macd,charlist.rsi,charlist.trix)
  
  return(charlist.tf)
}

test.create.signals <- function(ia, from = NULL, signaltype = 'ts'){
  # Setup start/end dates for estimation
  #startdate <-'2014-01-01'
  # from <- startdate
  
  #from = NULL
  
  if(signaltype == 'ts') charlist <- create.trend.signals(ia, from = from)
  if(signaltype == 'xs') charlist <- create.value.signals(ia, from = from)
  if(signaltype == 'macro') charlist <- create.marco.signals(ia, from = from)
  return(charlist)
}

apply.postion.sizing <- function(signallist, type = c("vp","vp.tgtvol","vp.corr","vp.corr.tgtvol"),...){
  
  type = match.arg(type[1], c("vp","vp.tgtvol","vp.corr","vp.corr.tgtvol"))
  
  xfun <- switch(type,
                 "vp" = riskWeighted,
                 "vp.tgtvol" = targetVolatility,
                 "vp.corr" = riskWeightedCorr,
                 "vp.corr.tgtvol" = targetVolatility3)
  
  msg = try( xfun(signallist,... ) , silent=TRUE)
  if (class(msg)[1] == 'try-error')
    warning(msg, '\n')
  return(msg)
}



test.postion.sizing2 <- function(startdate = '2014-01-01', enddate ='2017-05-11', signaltype = 'ts',gamma = 50, alphatest = 0,  CostData = NULL, optcost = FALSE, type = 'OMG.out',write.to.file = FALSE,filelocation = NULL){
  
  ia <- test.create.historical.data()
  
  charlist.tf <- test.create.signals(ia, from = startdate, signaltype = signaltype)
  
  BenchW <- ntop(ia$prices, ncol(ia$prices))
  
  if(CostData > 0){
    cost.test <- CostData
    CostData <- BenchW
    CostData[] <- NA
    CostData[] <- cost.test
    
  }
  
  
  charlist.vp <- apply.postion.sizing(charlist.tf,type = "vp", ia$rets.vol,ia$prices)
  
  charlist.vp.corr <- apply.postion.sizing(charlist.tf, type = "vp.corr", ia$rets, ia$rets.vol, from = startdate)
  
  charlist.all <- list(charlist.tf, charlist.vp, charlist.vp.corr)
  
  bt.models <- list()
  
  bt.names <- c("Org","vp", "vp.corr")
  est.out <- NULL
  
  for(i in 1:length(charlist.all)){
    
    ppp.cv <- crossvalidate(charlist.all[[i]], BenchW, ia$prices, CostData = CostData, gamma = gamma, alpha = alphatest, 
                            startdate = startdate, enddate = enddate,type = type,nfolds = 15, 
                            inference.type = "kfolds",name=bt.names[i])
    
    bt.models[[i]] <- backtest(ppp.cv)
    est.out <- c(est.out, capture.output(print(bt.models[[i]])))
    
  }
  names(bt.models) <- bt.names
  x11()
  plotbtlist(bt.models)
  
  if(write.to.file){
    cat(paste("Backtest of different position sizing functions",Sys.time()), est.out, file=paste(filelocation, 'Test_position_sizing.txt', sep=''), sep="\n", append=TRUE)
  }
}


test.postion.sizing <- function(startdate = '2014-01-01', enddate ='2017-05-11', signaltype = 'ts',gamma = 50, alphatest = 0,  CostData = NULL, optcost = FALSE, type = 'OMG.out',write.to.file = FALSE,filelocation = NULL){
  
  ia <- test.create.historical.data()
  
  charlist.tf <- test.create.signals(ia, from = startdate, signaltype = signaltype)
  
  charlist.macro <- test.create.signals(ia, from = startdate, signaltype = 'macro')
  
  BenchW <- ntop(ia$prices, ncol(ia$prices))
  
  if(CostData > 0){
    cost.test <- CostData
    CostData <- BenchW
    CostData[] <- NA
    CostData[] <- cost.test
    
  }
  
  
  charlist.vp <- apply.postion.sizing(charlist.tf,type = "vp", ia$rets.vol,ia$prices)
  charlist.vp.tgtvol <- apply.postion.sizing(charlist.tf, type = "vp.tgtvol",ia$prices, target = 10/100, 
                                             lookback.len = 21, max.portfolio.leverage = 300/100,
                                             risk.weighted = TRUE, vols = ia$rets.vol)
  charlist.vp.corr <- apply.postion.sizing(charlist.tf, type = "vp.corr", ia$rets, ia$rets.vol, from = startdate)
  charlist.vp.corr.tgtvol <- apply.postion.sizing(charlist.vp.corr, type = "vp.corr.tgtvol", ia$prices, target = 10/100, 
                                                  lookback.len = 21, max.portfolio.leverage = 300/100)
  
  
  charlist.tf <- c(charlist.tf, charlist.macro)
  charlist.vp <- c(charlist.vp, charlist.macro)
  charlist.vp.tgtvol <- c(charlist.vp.tgtvol, charlist.macro)
  charlist.vp.corr <- c(charlist.vp.corr, charlist.macro)
  charlist.vp.corr.tgtvol <- c(charlist.vp.corr.tgtvol, charlist.macro)
  
  charlist.all <- list(charlist.tf, charlist.vp, charlist.vp.tgtvol, charlist.vp.corr, charlist.vp.corr.tgtvol)
  
  bt.models <- list()
  
  bt.names <- c("Org","vp","vp.targtvol", "vp.corr", "vp.corr.tgtvol")
  est.out <- NULL
  
  for(i in 1:length(charlist.all)){
    
    ppp.cv <- crossvalidate(charlist.all[[i]], BenchW, ia$prices, CostData = CostData, gamma = gamma, alpha = alphatest, 
                            startdate = startdate, enddate = enddate,type = type,nfolds = 15, 
                            inference.type = "kfolds",name=bt.names[i])
    
    bt.models[[i]] <- backtest(ppp.cv)
    est.out <- c(est.out, capture.output(print(bt.models[[i]])))
    
  }
  names(bt.models) <- bt.names
  x11()
  plotbtlist(bt.models)
  
  if(write.to.file){
    cat(paste("Backtest of different position sizing functions",Sys.time()), est.out, file=paste(filelocation, 'Test_position_sizing.txt', sep=''), sep="\n", append=TRUE)
  }
}

test.risk.aversion <- function(gamma.test = c(5,10,15,50,100), startdate = '2014-01-01', enddate ='2017-05-11', signaltype = 'ts',pos.type = "vp", alphatest = 0,  CostData = NULL, type = 'OMG.out', write.to.file = FALSE,filelocation = NULL){
  
  ia <- test.create.historical.data()
  
  charlist.tf <- test.create.signals(ia, from = startdate, signaltype = signaltype)
  charlist.macro <- test.create.signals(ia, from = startdate, signaltype = 'macro')
  
  if(pos.type == "vp"){
    charlist.tf <- apply.postion.sizing(charlist.tf,type = "vp", ia$rets.vol,ia$prices)
  }
  
  if(pos.type == "vp.tgtvol"){
    
    charlist.tf <- apply.postion.sizing(charlist.tf, type = "vp.tgtvol",ia$prices, target = 10/100, 
                                        lookback.len = 21, max.portfolio.leverage = 300/100,
                                        risk.weighted = TRUE, vols = ia$rets.vol)
  }
  
  charlist.tf <- c(charlist.tf,charlist.macro)
  
  BenchW <- ntop(ia$prices, ncol(ia$prices))
  
  if(CostData > 0){
    cost.test <- CostData
    CostData <- BenchW
    CostData[] <- NA
    CostData[] <- cost.test
    
  }
  
  bt.models <- list()
  
  est.out <- NULL
  
  bt.names <- paste("gamma.",gamma.test,sep="")
  
  for(i in 1:length(gamma.test)){
    
    ppp.cv <- crossvalidate(charlist.tf, BenchW, ia$prices, CostData = CostData, gamma = gamma.test[i], alpha = alphatest, 
                            startdate = startdate, enddate = enddate,type = type,nfolds = 15, 
                            inference.type = "kfolds",name=bt.names[i])
    
    bt.models[[i]] <- backtest(ppp.cv)
    est.out <- c(est.out, capture.output(print(bt.models[[i]])))
    
  }
  names(bt.models) <- bt.names
  x11()
  plotbtlist(bt.models)
  
  if(write.to.file){
    cat(paste("Backtest of different risk aversions",Sys.time()), est.out, file=paste(filelocation, 'Test_risk_aversion.txt', sep=''), sep="\n", append=TRUE)
  }
}


test.anent.alpha <- function(alpha.test = c(0,1,2,4), startdate = '2014-01-01', enddate ='2017-05-11', signaltype = 'ts',pos.type = "vp", gamma = 50,  CostData = NULL, type = 'OMG.out', write.to.file = FALSE,filelocation = NULL){
  
  ia <- test.create.historical.data()
  
  charlist.tf <- test.create.signals(ia, from = startdate, signaltype = signaltype)
  charlist.macro <- test.create.signals(ia, from = startdate, signaltype = 'macro')
  
  if(pos.type == "vp"){
    charlist.tf <- apply.postion.sizing(charlist.tf,type = "vp", ia$rets.vol,ia$prices)
  }
  
  if(pos.type == "vp.tgtvol"){
    
    charlist.tf <- apply.postion.sizing(charlist.tf, type = "vp.tgtvol",ia$prices, target = 10/100, 
                                        lookback.len = 21, max.portfolio.leverage = 300/100,
                                        risk.weighted = TRUE, vols = ia$rets.vol)
  }
  
  charlist.tf <- c(charlist.tf,charlist.macro)
  
  
  BenchW <- ntop(ia$prices, ncol(ia$prices))
  
  if(CostData > 0){
    cost.test <- CostData
    CostData <- BenchW
    CostData[] <- NA
    CostData[] <- cost.test
    
  }
  
  bt.models <- list()
  
  est.out <- NULL
  
  bt.names <- paste("alpha.",alpha.test,sep="")
  
  for(i in 1:length(alpha.test)){
    
    ppp.cv <- crossvalidate(charlist.tf, BenchW, ia$prices, CostData = CostData, gamma = gamma, alpha = alpha.test[i], 
                            startdate = startdate, enddate = enddate,type = type,nfolds = 15, 
                            inference.type = "kfolds",name=bt.names[i])
    
    bt.models[[i]] <- backtest(ppp.cv)
    est.out <- c(est.out, capture.output(print(bt.models[[i]])))
    
  }
  names(bt.models) <- bt.names
  x11()
  plotbtlist(bt.models)
  
  if(write.to.file){
    cat(paste("Backtest of different alpha",Sys.time()), est.out, file=paste(filelocation, 'Test_aenet_alpha.txt', sep=''), sep="\n", append=TRUE)
  }
}


test.transaction.costs <- function(cost.test = c(0,0.03/100,0.05/100, 0.10/100), startdate = '2014-01-01', enddate ='2017-05-11', signaltype = 'ts',pos.type = "vp", gamma = 50, alphatest = 0, type = 'OMG.out', write.to.file = FALSE,filelocation = NULL){
  
  ia <- test.create.historical.data()
  
  charlist.tf <- test.create.signals(ia, from = startdate, signaltype = signaltype)
  charlist.macro <- test.create.signals(ia, from = startdate, signaltype = 'macro')
  
  if(pos.type == "vp"){
    charlist.tf <- apply.postion.sizing(charlist.tf,type = "vp", ia$rets.vol,ia$prices)
  }
  
  if(pos.type == "vp.tgtvol"){
    
    charlist.tf <- apply.postion.sizing(charlist.tf, type = "vp.tgtvol",ia$prices, target = 10/100, 
                                        lookback.len = 21, max.portfolio.leverage = 300/100,
                                        risk.weighted = TRUE, vols = ia$rets.vol)
  }
  
  charlist.tf <- c(charlist.tf, charlist.macro)
  
  BenchW <- ntop(ia$prices, ncol(ia$prices))
  
  
  bt.models <- list()
  
  est.out <- NULL
  
  bt.names <- paste("cost.",cost.test/(0.01/100),"bps",sep="")
  
  CostData <- BenchW
  
  for(i in 1:length(cost.test)){
    
    
    CostData[] <- NA
    CostData[] <- cost.test[i]
    
    ppp.cv <- crossvalidate(charlist.tf, BenchW, ia$prices, CostData = CostData, gamma = gamma, alpha = alphatest, 
                            startdate = startdate, enddate = enddate,type = type,nfolds = 15, 
                            inference.type = "kfolds",name=bt.names[i])
    
    bt.models[[i]] <- backtest(ppp.cv)
    est.out <- c(est.out, capture.output(print(bt.models[[i]])))
    
  }
  names(bt.models) <- bt.names
  x11()
  plotbtlist(bt.models)
  
  if(write.to.file){
    cat(paste("Backtest of different transaction costs",Sys.time()), est.out, file=paste(filelocation, 'Test_transaction_costs.txt', sep=''), sep="\n", append=TRUE)
  }
}



test.marco.signals <- function(startdate = '2014-01-01', enddate ='2017-05-11',gamma = 50, alphatest = 0,  CostData = NULL, optcost = FALSE, type = 'OMG.out',write.to.file = FALSE,filelocation = NULL){
  
  ia <- test.create.historical.data()
  
  charlist.tf <- test.create.signals(ia, from = startdate, signaltype = 'ts')
  
  charlist.marco <- test.create.signals(ia, from = startdate, signaltype = 'macro')
  
  BenchW <- ntop(ia$prices, ncol(ia$prices))
  
  if(CostData > 0){
    cost.test <- CostData
    CostData <- BenchW
    CostData[] <- NA
    CostData[] <- cost.test
    
  }
  
  
  charlist.vp <- apply.postion.sizing(charlist.tf,type = "vp", ia$rets.vol,ia$prices)
  
  charlist.vp.macro <- c(charlist.vp, charlist.marco)
  
  charlist.all <- list(charlist.vp, charlist.vp.macro)
  
  bt.models <- list()
  
  bt.names <- c("vp","vp.macro")
  est.out <- NULL
  
  for(i in 1:length(charlist.all)){
    
    ppp.cv <- crossvalidate(charlist.all[[i]], BenchW, ia$prices, CostData = CostData, gamma = gamma, alpha = alphatest, 
                            startdate = startdate, enddate = enddate,type = type,nfolds = 15, 
                            inference.type = "kfolds",name=bt.names[i])
    
    bt.models[[i]] <- backtest(ppp.cv,fixedWindow = F)
    est.out <- c(est.out, capture.output(print(bt.models[[i]])))
    
  }
  names(bt.models) <- bt.names
  x11()
  plotbtlist(bt.models)
  
  if(write.to.file){
    cat(paste("Backtest of macro",Sys.time()), est.out, file=paste(filelocation, 'Test_macro.txt', sep=''), sep="\n", append=TRUE)
  }
}

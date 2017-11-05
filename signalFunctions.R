
norm_vec <- function(x) sqrt(sum(x^2),na.rm = TRUE)
norm1_vec <- function(x) sum(abs(x),na.rm = TRUE)
norm2_vec <- function(x) sum((x),na.rm = TRUE)

normalizeVec <- function(out){
  
  outnorm <- norm1_vec(out)
  if(is.na(outnorm)) outnorm <- 0
  if(outnorm != 0) out <- out/outnorm
  out
}

normalizeSignal <- function(out, normtype){
  
  if(normtype == 'ts'){
    out <- (out)/apply(out,1,sd)
  }
  if(normtype == 'xs'){
    out <- (out - apply(out,1,mean,na.rm=TRUE))/apply(out,1,sd,na.rm=TRUE)
  }
  if(normtype == 'norm'){
    out <- out/apply(out,1,norm_vec)
  }
  out
}


getMomentum <- function(x, nmonth = 3, nday = 20, type = 'sign',normalize = F, normtype = 'ts', include.lag = TRUE){
  
  mom = x/mlag(x, nday*nmonth) - 1
  
  if(type == 'sign') {
    mom <- sign(mom)
    mom[is.na(mom)] = 0
  }
  
  if(type == 'value'){
    
    mom.vol = bt.apply.matrix(mom, runsd, k = 250)
    
    mom <- mom/mom.vol
  }

  if(normalize){
    
    mom <- normalizeSignal(mom, normtype)
    
  }
  mom[is.na(mom)] = 0
  
  if(include.lag) mom <- mlag(mom)
  
  mom
}


getEMA1 <- function(x, n = 50, percent = F, type = 'binary',normalize = FALSE, normtype = 'norm', from = NULL, include.lag = TRUE){
  x.ema <- bt.apply.matrix(x, EMA, n = n)
  
  
  if(type == "binary"){
    if(!is.null(from)){
      from <- paste("/",from,sep="")
      x[from,] <- NA
      x.ema[from,] <- NA
    }
    out <- iif(cross.up(x, x.ema), 1, iif(cross.dn(x, x.ema), -1, NA))
    out <- apply(out, 2, ifna.prev)
    out[is.na(out)] = 0
    out <- xts(out,order.by = index(x))
  }

    if(type == 'sign'){
      if(percent){
        out <- 100*(x/x.ema - 1)
      }else{
        out <- (x - x.ema)
      }
  
      out <- sign(out)
    }
  
  if(type == 'value'){
    
    if(percent){
      out <- 100*(x/x.ema - 1)
    }else{
      out <- (x - x.ema)
    }
    
    x.vol = bt.apply.matrix(out, runsd, k = 250)
    out <- out/x.vol
  }
  
  if(normalize){
    
    out <- normalizeSignal(out, normtype)
    
  }
  out[is.na(out)] = 0
  if(include.lag) out <- mlag(out)
  out
}



getEMA2 <- function(x, n1 = 50, n2 = 200, percent = F, type = 'binary',normalize = FALSE, normtype = 'norm', from = NULL, include.lag = TRUE){
  x.ema1 <- bt.apply.matrix(x, EMA, n = n1)
  x.ema2 <- bt.apply.matrix(x, EMA, n = n2)
  
  if(type == 'binary'){
    
    if(!is.null(from)){
      from <- paste("/",from,sep="")
      x.ema1[from,] <- NA
      x.ema2[from,] <- NA
    }
    
    out <- iif(cross.up(x.ema1, x.ema2), 1, iif(cross.dn(x.ema1, x.ema2), -1, NA))
    out <- apply(out, 2, ifna.prev)
    out[is.na(out)] = 0
    out <- xts(out,order.by = index(x))
  }
  
  if(type == 'sign'){
    
    if(!is.null(from)){
      from <- paste("/",from,sep="")
      x.ema1[from,] <- NA
      x.ema2[from,] <- NA
    }
    
    out <- iif(x.ema1 > x.ema2, 1, iif(x.ema1 < x.ema2, -1, NA))
    out <- apply(out, 2, ifna.prev)
    out[is.na(out)] = 0
    out <- xts(out,order.by = index(x))
  }
  
  if(type == 'value'){
    
    if(percent){
      out <- 100*(x.ema1/x.ema2 - 1)
    }else{
      out <- x.ema1 - x.ema2
    }
    
    x.vol = bt.apply.matrix(out, runsd, k = 250)
    
    out <- out /x.vol
  }
  
  if(normalize){
    
    out <- normalizeSignal(out, normtype)
    
  }
  out[is.na(out)] = 0
  if(include.lag) out <- mlag(out)
  out
}


computeTrix <- function(price,n = 20, nSig = 9, percent = TRUE, outtype = "TRIX"){
  
  TRIX(price, n, nSig)[,outtype]
}


getTRIX <- function(x, n = 20, nSig = 9, percent = TRUE, type = 'binary',normalize = FALSE, normtype = 'norm', from = NULL, include.lag = TRUE){
  x.trix <- bt.apply.matrix(x, computeTrix, n = n, nSig = nSig, percent = TRUE, outtype = "TRIX")
  x.signal<- bt.apply.matrix(x, computeTrix, n = n, nSig = nSig, percent = TRUE, outtype = "signal")

  
  if(type == 'binary'){
    
    if(!is.null(from)){
      from <- paste("/",from,sep="")
      x.trix[from,] <- NA
      x.signal[from,] <- NA
    }
    
    out <- iif(cross.up(x.trix, x.signal) & (x.trix > 0), 1, iif(cross.dn(x.trix, x.signal) & (x.trix < 0), -1, NA))
    out <- apply(out, 2, ifna.prev)
    out[is.na(out)] = 0
    out <- xts(out,order.by = index(x))
  }
  
  if(type == 'value'){
    
    out <- x.trix - x.signal
    x.vol = bt.apply.matrix(out, runsd, k = 250)
    out <- out/x.vol
  }
  
  if(normalize){
    
    out <- normalizeSignal(out, normtype)
    
  }
  out[is.na(out)] = 0
  if(include.lag) out <- mlag(out)
  
  out
}


getDonchianChannel <- function(datalist, n = 10, dc.lag = FALSE,normalize = FALSE, normtype = 'norm', from = NULL, include.lag = TRUE){
  
  datalist2 <- lapply(datalist, function(data){
    
    dc <- DonchianChannel(data[,c("High","Low")], n, dc.lag)
    
    if(!is.null(from)){
      from <- paste("/",from,sep="")
      data[from,] <- NA
      dc[from,] <- NA
    }
    
    out <- iif(data[,"Close"] == dc[,"high"] , 1, iif(data[,"Close"] == dc[,"low"] , -1, NA))
    
    out <- apply(out, 2, ifna.prev)
    out[is.na(out)] = 0
    out <- xts(out,order.by = index(data))
    out
  })
  
  out <- Reduce(cbind, datalist2)
  colnames(out) <- names(datalist)
  
  if(normalize){
    
    out <- normalizeSignal(out, normtype)
    
  }
  
  if(include.lag) out <- mlag(out)
  
  out
}


getBBands <- function(datalist, n = 20, maType="SMA", sd = 2, normalize = FALSE, normtype = 'norm', from = NULL, include.lag = TRUE){
  
  datalist2 <- lapply(datalist, function(data){
    
    bbands.HLC <- BBands(data[,c("High","Low","Close")], n = n, maType = maType, sd = sd)
    
    if(!is.null(from)){
      from <- paste("/",from,sep="")
      data[from,] <- NA
      bbands.HLC[from,] <- NA
    }
    
    #out <- iif(cross.dn(data[,"Close"], bbands.HLC[,"dn"]) , 1, iif(cross.up(data[,"Close"], bbands.HLC[,"up"]) , -1, iif(cross(data[,"Close"], bbands.HLC[,"mavg"]), 0, NA)))
    out <- iif(cross.dn(data[,"Close"], bbands.HLC[,"dn"]) , 1, iif(cross.up(data[,"Close"], bbands.HLC[,"up"]) , -1, NA))
    
    out <- apply(out, 2, ifna.prev)
    out[is.na(out)] = 0
    out <- xts(out,order.by = index(data))
    out
  })
  out <- Reduce(cbind, datalist2)
  colnames(out) <- names(datalist)
  
  if(normalize){
    
    out <- normalizeSignal(out, normtype)
    
  }
  
  if(include.lag) out <- mlag(out)
  
  out
}

computeMACD <- function(x, nFast = 12, nSlow = 26, nSig = 9, maType = "EMA", percent = TRUE, type = "binary", from = NULL){
  
  x.macd <- MACD(x, nFast = nFast, nSlow = nSlow, nSig = nSig, maType = maType,percent = percent)
  
  if(type == 'binary'){
    
    if(!is.null(from)){
      from <- paste("/",from,sep="")
      #x[from,] <- NA
      x.macd[from,] <- NA
    }
    
    #out <- iif(cross.up(x.macd[,'macd'],x.macd[,'signal']),1,iif(cross.dn(x.macd[,'macd'],x.macd[,'signal']) , -1, NA))
    
    out <- iif(cross.up(x.macd[,'macd'],x.macd[,'signal']) & (x.macd[,'macd'] > 0),1,iif(cross.dn(x.macd[,'macd'],x.macd[,'signal']) | (x.macd[,'macd'] < 0), -1, NA))
    
    out <- apply(out, 2, ifna.prev)
    out[is.na(out)] = 0
    out <- xts(out,order.by = index(x))
  }
  
  if(type == 'value'){
    
    out <- x.macd[,'macd'] - x.macd[,'signal']
    
    x.vol = bt.apply.matrix(out, runsd, k = 250)
    out <- out/x.vol
  }
  
  
  out
}

getMACD <- function(x, nFast = 12, nSlow = 26, nSig = 9, maType = "EMA", percent = TRUE, type = "binary", normalize = FALSE, normtype = 'norm', from = NULL, include.lag = TRUE){
  
  out <- bt.apply.matrix(x, computeMACD, usecoredata = FALSE, nFast = nFast, nSlow = nSlow, nSig = nSig, maType = maType,percent = percent, type = type, from = from)
  
  if(normalize){
    
    out <- normalizeSignal(out, normtype)
    
  }
  out[is.na(out)] = 0
  if(include.lag) out <- mlag(out)
  
  out
}


getRSI <- function(x, n = 14, maType="EMA", b = 0.3, type = "binary", normalize = FALSE, normtype = "norm", from = NULL, include.lag = TRUE){
  
  x.rsi <- bt.apply.matrix(x, RSI, n = n, maType = maType)
  x.rsi <- x.rsi/100
  
  if(type == 'binary'){
    
    if(!is.null(from)){
      from <- paste("/",from,sep="")
      #x[from,] <- NA
      x.rsi[from,] <- NA
    }
    
    
    out <- iif(cross.dn(x.rsi, b), 1, iif(cross.up(x.rsi, 1- b), -1, NA))
    
    #out <- iif(cross.up(x.rsi, b), 1, iif(cross.dn(x.rsi, 1- b), -1, NA))

    out <- apply(out, 2, ifna.prev)
    out[is.na(out)] = 0
    out <- xts(out,order.by = index(x))
    
  }
  
  if(type == 'value'){
    
    x.rsi <- x.rsi - 0.5
    rsi.vol = bt.apply.matrix(x.rsi, runsd, k = 250)
    out <- x.rsi/rsi.vol
    
  }
  
  
  if(normalize){
    
    out <- normalizeSignal(out, normtype)
    
  }
  out[is.na(out)] = 0
  if(include.lag) out <- mlag(out)
  
  out
  
}



getStoch <- function(datalist, nFastK = 14, nFastD = 3, nSlowD = 3, maType = "SMA", b = 0.3, bounded = TRUE, smooth = 1,normalize = FALSE, normtype = 'norm', from = NULL){
  
  datalist2 <- lapply(datalist, function(data){
    
    x.stoch <- stoch(data[,c("High","Low","Close")], nFastK = nFastK, nFastD = nFastD, nSlowD = nSlowD, maType = maType, bounded = bounded, smooth = smooth)[,'slowD']
    
    if(!is.null(from)){
      from <- paste("/",from,sep="")
      x.stoch[from,] <- NA
    }
    
    out <- iif(cross.dn(x.stoch, b), 1, iif(cross.up(x.stoch, 1 - b), -1, NA))
    #out <- iif(cross.up(x.stoch, b), 1, iif(cross.dn(x.stoch, 1 - b), -1, NA))
    
    out <- apply(out, 2, ifna.prev)
    out[is.na(out)] = 0
    out <- xts(out,order.by = index(data))
    out
  })
  
  out <- Reduce(cbind, datalist2)
  colnames(out) <- names(datalist)
  
  if(normalize){
    
    out <- normalizeSignal(out, normtype)
    
  }
  out
  
}


getOBV <- function(datalist, n = 20, type = 'binary', normalize = FALSE, normtype = 'norm', from = NULL, include.lag = TRUE){
  
  datalist2 <- lapply(datalist, function(data){
    
    obv <- OBV(data[,"Close"], data[,"Volume"])
    obv.ma <- EMA(obv, n = n)
    
    if(type == 'binary'){
      
      if(!is.null(from)){
        obv <- paste("/",from,sep="")
        obv.ma[from,] <- NA
      }
      
      out <- iif(cross.up(obv, obv.ma), 1, iif(cross.dn(obv, obv.ma), -1, NA))
      out <- apply(out, 2, ifna.prev)
      out[is.na(out)] = 0
      out <- xts(out,order.by = index(data))
    }
    
    if(type == 'sign'){
      
      if(!is.null(from)){
        obv <- paste("/",from,sep="")
        obv.ma[from,] <- NA
      }
      
      out <- sign(obv - obv.ma)

    }

    out
  })
  
  out <- Reduce(cbind, datalist2)
  colnames(out) <- names(datalist)
  
  if(normalize){
    
    out <- normalizeSignal(out, normtype)
    
  }
  
  if(include.lag) out <- mlag(out)
  
  out
}


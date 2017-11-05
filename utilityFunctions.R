
getWeeklyData <- function(xstdata, day = 3){
  
  days = as.POSIXlt(index(xstdata))$wday == day
  
  indx <- c(0, which(days))
  
  tmp <- period.apply(xstdata, INDEX=indx, FUN=last)
  
  return(tmp)
  
}

iif <- function
(
  cond,
  truepart,
  falsepart
)
{
  if(len(cond) == 1) { if(cond) truepart else falsepart }
  else {
    if(length(falsepart) == 1) {
      temp = falsepart
      falsepart = cond
      falsepart[] = temp
    }
    if(length(truepart) == 1)
      falsepart[cond] = truepart
    else {
      cond = ifna(cond,F)
      if(is.xts(truepart))
        falsepart[cond] = coredata(truepart)[cond]
      else
        falsepart[cond] = truepart[cond]
    }
    return(falsepart);
  }
}

len <- function(x){
  return(length(x))
}

ntop <- function(data, topn = 1, dirMaxMin = TRUE){
  temp = coredata(data)
  if(topn == ncol(data)) {
    index = is.na(temp)
    temp[index] = 0
    temp[!index] = 1
    out = data
    out[] = ifna(temp / rowSums(temp),0)
    return( out )
  }
  for( i in 1:nrow(data) ) {
    x = temp[i,]
    o = sort.list(x, na.last = TRUE, decreasing = dirMaxMin)
    index = which(!is.na(x))
    x[] = NA
    if(len(index)>0) {
      n = min(topn, len(index))
      x[o[1:n]] = 1/n
    }
    temp[i,] = x
  }
  temp[is.na(temp)] = 0
  out = data
  out[] = temp
  return( out )
}

bt.apply.matrix <- function(b, xfun=Cl, usecoredata = TRUE,...){
  out = b
  out[] = NA
  nsymbols = ncol(b)
  xfun = match.fun(xfun)
  
  if(usecoredata){
    b.data <- coredata(b)
  }else{
    b.data <- b
  }
  
  for( i in 1:nsymbols ) {
    msg = try( xfun( b.data[,i],... ) , silent=TRUE)
    if (class(msg)[1] == 'try-error')
      warning(i, msg, '\n')
    else
      out[,i] = msg
  }
  return(out)
}

mlag <- function(m, nlag = 1) {
  if( is.null(dim(m)) ) {
    n = len(m)
    if(nlag > 0) {
      m[(nlag+1):n] = m[1:(n-nlag)]
      m[1:nlag] = NA
    } else if(nlag < 0) {
      m[1:(n+nlag)] = m[(1-nlag):n]
      m[(n+nlag+1):n] = NA
    }
  } else {
    n = nrow(m)
    if(nlag > 0) {
      m[(nlag+1):n,] = m[1:(n-nlag),]
      m[1:nlag,] = NA
    } else if(nlag < 0) {
      m[1:(n+nlag),] = m[(1-nlag):n,]
      m[(n+nlag+1):n,] = NA
    }
  }
  return(m);
}


setupListMat <- function(charlist){
  lapply(charlist, function(x){coredata(x)})
  
}
setupListRange <- function(charlist, range, outtype = "mat"){
  
  lapply(charlist, function(x,range){
    if(outtype == "mat") return(coredata(x[range,]))
    if(outtype == "xts") return(x[range,])
    },range = range)
  
}


setupListColRange <- function(charlist, range){
  
  lapply(charlist, function(x,range){x[,range]},range = range)
  
}

volAdjustList <- function(charlist, vols, prices = NULL, nrisk = 1/100, type = "cvolvalue"){
  
  #if(type == "cvolvalue") vols[vols<0.1]<-0.1
  if(type == "riskparity") vols[vols<0.01]<-0.01
  out <- lapply(charlist, function(x){
    if(type == "cvolvalue") w <- (x*nrisk*prices)/vols
    if(type == "cvolweight") w <- (x*nrisk)/vols
    #if(type == "riskparity") w <- x*((1/vols)/rowSums(1/vols, na.rm=TRUE))
    if(type == "riskparity") w <- x*(1/vols)
    #w <- w/apply(w,1,norm_vec)
    #w <- apply(x,1,norm1_vec)*w/apply(w,1,norm1_vec)
    #w <- apply(x,1,norm2_vec)*w/apply(w,1,norm2_vec)
    #w <- w/apply(w,1,sd)
    w <- apply(x,1,norm1_vec)*t(apply(w,1,normalizeVec))
    
    w
  })
  out
}

repmat <- function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}

repmatxts <- function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  nx = dim(X)[2]
  Y <- matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
  #colnames(Y)<-paste('IR.',1:ncol(a),sep='')
  xts(Y, order.by = index(X))
}

getMacroFactors <- function(data,idx,symbol,k){
  out <- data[idx,symbol]
  out <- na.locf(out)
  out<-repmatxts(out,1,k)
}







risk.Weighted.strategy <- function(weight, vols, target = 10/100){
  vols[vols<0.01]<-0.01
  weight.risk <- target*weight / vols
  
  #weight.risk <- rowSums(abs(weight))*t(apply(weight.risk,1,normalizeVec))
  #weight.risk <- t(apply(weight.risk,1,normalizeVec))
  
  rs = rowSums(abs(weight.risk), na.rm=T)
  weight.risk = weight.risk / iif(rs > 3, rs/3, 1)
  weight.risk[is.na(weight.risk)] = 0
  
  weight.risk
}

risk.Weighted.strategy.bench <- function(weight, vols, target = 10/100){
  vols[vols<0.01]<-0.01
  weight.risk <- target*weight / vols
  weight.dollar <- ntop(vols, ncol(vols))
  weight.risk <- weight.risk*weight.dollar
  weight.risk
}

riskWeighted <- function(charlist, vols,target = 10/100){
  
  
  weight.dollar <- ntop(vols, ncol(vols))
  
  out <- lapply(charlist, function(x){
    #w <- weight.dollar*x
    w<-x
    w <- risk.Weighted.strategy(w, vols,target=target)
    
    #w/weight.dollar
  })
  out
}

Volatility.Weighted.strategy <- function(weight, vols){
  weight.risk <- weight / vols
  
  weight.risk <- t(apply(weight.risk,1,normalizeVec))
  weight.risk[is.na(weight.risk)] = 0
  
  return(weight.risk)
}

volatilityWeighted <- function(charlist, vols){
  vols[vols<0.01]<-0.01
  weight.dollar <- ntop(vols, ncol(vols))

  out <- lapply(charlist, function(x){
      w <- Volatility.Weighted.strategy(x, vols)
      w <- w/weight.dollar
      
      w[is.na(w)] = 0
    return(w)
  })
  names(out) <- paste(names(charlist),".vol",sep="")
  out
}

riskWeighted2 <- function(charlist, vols, CovList, targetvol = 10/100, from = 1){
  
  out <- lapply(charlist, function(x){
    #w <- weight.dollar*x
    w<-x
    w <- risk.Weighted.strategy2(w, vols,CovList, targetvol = targetvol, from = from)
    
    #w/weight.dollar
  })
  out
}



target.vol.strategy <- function(model, weight,
         target = 10/100,
         lookback.len = 21,
         max.portfolio.leverage = 100/100,
         annualperiods = 252
)
{
  ret = diff(log(model$equity))
  #hist.vol.model = sqrt(annualperiods) * runSD(ret, n = lookback.len)
  hist.vol.model = sqrt(annualperiods) * runsd(ret, k = lookback.len)
  hist.vol.model = as.vector(hist.vol.model)
  weight.target = weight * (target / hist.vol.model)
  
  
  
  rs = rowSums(abs(weight.target), na.rm=T)
  weight.target = weight.target / iif(rs > max.portfolio.leverage, rs/max.portfolio.leverage, 1)
  
  weight.target[is.na(weight.target)] = 0 
  return(weight.target)
}
targetVolatility <- function(charlist, prices, target = 5/100, lookback.len = 21, max.portfolio.leverage = 200/100, risk.weighted = FALSE, vols = NULL){
  
  
  data <- new.env()
  data$prices <- prices
  data$dates <- index(prices)
  
  naxts <- prices
  naxts[] <- NA
  
  data$weight <- naxts
  data$execution.price <- naxts
  
  weight.dollar <- ntop(prices, ncol(prices))
  
  out <- lapply(charlist, function(x){
    weight.timing <- weight.dollar * x
    
    if(risk.weighted) weight.timing <- risk.Weighted.strategy(weight.timing, vols)
    
    data$weight[] = NA
    data$weight = weight.timing
    models.timing = bt.run.share(data, clean.signal=F, silent = T)
    
    w <- target.vol.strategy(models.timing,weight.timing, target, lookback.len, max.portfolio.leverage)
    w <- w/weight.dollar
    w
  })
  out
}


targetVolatility2 <- function(charlist, prices, target = 5/100, lookback.len = 21, max.portfolio.leverage = 200/100, risk.weighted = TRUE, vols = NULL){
  
  
  data <- new.env()
  data$prices <- prices
  data$dates <- index(prices)
  
  naxts <- prices
  naxts[] <- NA
  
  data$weight <- naxts
  data$execution.price <- naxts
  
  weight.dollar <- ntop(prices, ncol(prices))
  
  vols[vols<0.01]<-0.01
  
  adj.vol <- 1/vols
  
  out <- lapply(charlist, function(x){
    
    weight.risk <- x*adj.vol/rowSums(adj.vol)
    
    data$weight[] = NA
    data$weight = weight.risk
    models.risk = bt.run.share(data, clean.signal=F, silent = T)
    
    w <- target.vol.strategy(models.risk, weight.risk, target, lookback.len, max.portfolio.leverage)
    w <- w/weight.dollar
    w
  })
  out
}


targetVolatility3 <- function(weightlist, prices, target = 5/100, lookback.len = 21, max.portfolio.leverage = 200/100,annualperiods = 252){
  
  
  data <- new.env()
  data$prices <- prices
  data$dates <- index(prices)
  
  naxts <- prices
  naxts[] <- NA
  
  data$weight <- naxts
  data$execution.price <- naxts
  
  weight.dollar <- ntop(prices, ncol(prices))
  
  out <- lapply(weightlist, function(x){
    weight.timing <- x*weight.dollar
    data$weight[] = NA
    data$weight = weight.timing
    models.timing = bt.run.share(data, clean.signal=F, silent = T)
    
    w <- target.vol.strategy(models.timing,weight.timing, target, lookback.len, max.portfolio.leverage,annualperiods=annualperiods)
    w <- w/weight.dollar
    w
  })
  out
}
eval_f <- function(w, covMat,signals,targetvol) {
  return( list( "objective" = -sum(log((w))), "gradient" = -1/w ) )
}


eval_g_ineq <- function(w, covMat,signals,targetvol) {
  
  comm <- as.numeric(sqrt(rbind(w)%*%covMat%*%cbind(w)))
  constr <- comm - targetvol
  grad <- (covMat%*%cbind(w))/comm
  
  return( list( "constraints"=constr, "jacobian"=grad ) )
}

eval_g_eq <- function(w, covMat, signals, targetvol) {
  
  constr <- sum(abs(w)) - 1
  grad <- sign(w)
  return( list( "constraints"=constr, "jacobian"=grad ) )
}
risk.parity.strategy0 <- function(signals, CovList, vols, targetvol = 10/100,from){
  vols[vols<0.01]<-0.01
  # w.opt <- vols
  # w.opt[] <- 0
  
  w.opt <- risk.Weighted.strategy(signals,vols)
  
  w.opt <- t(apply(w.opt,1,normalizeVec))
  
  for(i in from:nrow(vols)){
    
    idx <- which(signals[i,]!=0)
    
    if(length(idx) > 0){
      w0 <- 1/vols[i,idx]
      w0 <- w0/sum(w0)
      #  w0 <- w0*signals[i,idx]
      
      covMat <- CovList[[i]][idx,idx]
      
      rpw <- nloptr(x0=as.numeric(w0),
                    eval_f=eval_f,
                    eval_g_ineq = eval_g_ineq,
                    lb=c(rep(0,length(w0))),
                    opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1.0e-8), covMat = covMat, signals = signals[i,idx], targetvol = targetvol)
      
      w.opt[i,idx] <- rpw$solution
      
      w.opt[i,] <- signals[i,]*w.opt[i,]/sum(abs(w.opt[i,]))
      
      if(!is.null(targetvol)){
        rp.vol <- as.numeric(sqrt(rbind(as.numeric(w.opt[i,]))%*%CovList[[i]]%*%cbind(as.numeric(w.opt[i,]))))
        
        #rp.vol[rp.vol<0.01]<-0.01
        w.opt[i,] <- targetvol*w.opt[i,]/rp.vol
        
        
        rs = sum(abs(w.opt[i,]), na.rm=T)
        w.opt[i,] = w.opt[i,] / iif(rs > 3, rs/3, 1)
        
      }

    }
    
    
    
  }
  #w.opt <- rowSums(abs(signals))*w.opt/rowSums(abs(w.opt))
  w.opt[is.na(w.opt)] = 0
  
  w.opt <- w.opt*signals
  w.opt
}
risk.parity.strategy2 <- function(signals, CovList, vols, targetvol = 10/100,max.portfolio.leverage = 300/100, from){

  w.opt <- vols
  w.opt[] <- 0
  vols[vols<0.01]<-0.01
  #w.opt <- signals / vols
  #w.opt <- t(apply(w.opt,1,normalizeVec))
  #signals <- signals/vols
  
  signals <- t(apply(signals,1,normalizeVec))
  
  for(i in from:nrow(vols)){
    
    idx <- which(signals[i,]!=0)
    
    if(length(idx) > 2){
      w0 <- signals[i,idx]/vols[i,idx]
      
      #w0 <- signals[i,idx]
      w0 <- w0/sum(abs(w0))
      print(as.numeric(w0))
      #w0 <- w.opt[i,idx]
      covMat <- CovList[[i]][idx,idx]
      rpw <- nloptr(x0=as.numeric(w0),
                    eval_f=evalf,
                    eval_g_ineq = evalgineq,
                    opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1.0e-10), covMat = covMat, signals = signals[i,idx], targetvol = targetvol)
      
      
      w.opt[i,idx] <- rpw$solution
      
      w.opt[i,] <- w.opt[i,]/sum(abs(w.opt[i,]))
      
      if(!is.null(targetvol)){
        rp.vol <- as.numeric(sqrt(rbind(as.numeric(w.opt[i,]))%*%CovList[[i]]%*%cbind(as.numeric(w.opt[i,]))))
        

        w.opt[i,] <- targetvol*w.opt[i,]/rp.vol
        
        
        rs = sum(abs(w.opt[i,]), na.rm=T)
        w.opt[i,] = w.opt[i,] / iif(rs > max.portfolio.leverage, rs/max.portfolio.leverage, 1)
        
      }
    }
    
    
    
  }
  #w.opt <- w.opt/rowSums(abs(w.opt))
  w.opt[is.na(w.opt)] = 0
  w.opt
}

risk.Weighted.strategy2 <- function(weight, vols, CovList, targetvol = 10/100, max.portfolio.leverage = 300/100,from){
  vols[vols<0.01]<-0.01
  weight.risk <- weight / vols
  #weight.risk <- t(apply(weight.risk,1,normalizeVec))
  for(i in from:nrow(vols)){
    
    idx <- which(weight.risk[i,]!=0)
    
    if(length(idx) > 0){
      weight.risk[i,] <- weight.risk[i,]/sum(abs(weight.risk[i,]))
      
      rp.vol <- as.numeric(sqrt(rbind(as.numeric(weight.risk[i,]))%*%CovList[[i]]%*%cbind(as.numeric(weight.risk[i,]))))
      
      weight.risk[i,] <- targetvol*weight.risk[i,]/rp.vol
      
      
      rs = sum(abs(weight.risk[i,]), na.rm=T)
      weight.risk[i,] = weight.risk[i,] / iif(rs > max.portfolio.leverage, rs/max.portfolio.leverage, 1)
    }
    
  }
  
  
  weight.risk
}
target.Weighted.strategy2 <- function(weight, CovList, targetvol = 10/100, max.portfolio.leverage = 300/100,from){

  weight.risk <- weight

  for(i in from:nrow(weight)){
    
    idx <- which(weight.risk[i,]!=0)
    
    if(length(idx) > 0){
      weight.risk[i,] <- weight.risk[i,]/sum(abs(weight.risk[i,]))
      
      rp.vol <- as.numeric(sqrt(rbind(as.numeric(weight.risk[i,]))%*%CovList[[i]]%*%cbind(as.numeric(weight.risk[i,]))))
      
      weight.risk[i,] <- targetvol*weight.risk[i,]/rp.vol
      
      
      rs = sum(abs(weight.risk[i,]), na.rm=T)
      weight.risk[i,] = weight.risk[i,] / iif(rs > max.portfolio.leverage, rs/max.portfolio.leverage, 1)
    }
    
  }
  
  
  weight.risk
}

risk.parity.strategy <- function(signals, CovList, vols, targetvol,max.portfolio.leverage,from){
  vols[vols<0.01]<-0.01
  w.opt <- vols
  w.opt[] <- 0
  
  #w.opt <- risk.Weighted.strategy(signals,vols)
  
  #w.opt <- t(apply(w.opt,1,normalizeVec))
  
  for(i in from:nrow(vols)){
    
    idx <- which(signals[i,]!=0)
    
    if(length(idx) > 1){
      #w0 <- 1/vols[i,idx]
      #w0 <- w0/sum(abs(w0))
      #w0 <- signals[i,idx]*w0
      #w0 <- w.opt[i,idx]
      
      w0 <- signals[i,idx]/vols[i,idx]
      w0 <- w0/sum(abs(w0))
      covMat <- CovList[[i]][idx,idx]
      rpw <- nloptr(x0=as.numeric(w0),
                    eval_f=evalf,
                    eval_g_ineq = evalgineq,
                    eval_g_eq = evalgeq,
                    opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1.0e-8), covMat = covMat, signals = signals[i,idx]/vols[i,idx], targetvol = targetvol)
      
      w.opt[i,idx] <- rpw$solution
      
      if(!is.null(targetvol)){
        rp.vol <- as.numeric(sqrt(rbind(as.numeric(w.opt[i,]))%*%CovList[[i]]%*%cbind(as.numeric(w.opt[i,]))))
        

        w.opt[i,] <- targetvol*w.opt[i,]/rp.vol
        
        
        rs = sum(abs(w.opt[i,]), na.rm=T)
        w.opt[i,] = w.opt[i,] / iif(rs > max.portfolio.leverage, rs/max.portfolio.leverage, 1)
        
      }
      

    }
    
    
    
  }
  #w.opt <- w.opt/rowSums(abs(w.opt))
  w.opt[is.na(w.opt)] = 0
  w.opt
}

risk.parity.strategy.bench<- function(signals, CovList, vols, targetvol,from){
  vols[vols<0.01]<-0.01
  #w.opt <- vols
  #w.opt[] <- 0
  
  w.opt <- risk.Weighted.strategy(signals,vols)
  
  w.opt <- t(apply(w.opt,1,normalizeVec))
  
  for(i in from:nrow(vols)){
    
    idx <- which(signals[i,]!=0)
    
    if(length(idx) > 2){
      w0 <- 1/vols[i,idx]
      w0 <- w0/sum(abs(w0))
      w0 <- signals[i,idx]*w0
      #w0 <- w.opt[i,idx]
      covMat <- CovList[[i]][idx,idx]
      rpw <- nloptr(x0=as.numeric(w0),
                    eval_f=evalf,
                    eval_g_ineq = evalgineq,
                    eval_g_eq = evalgeq,
                    opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1.0e-8), covMat = covMat, signals = signals[i,idx], targetvol = targetvol)
      
      w.opt[i,idx] <- rpw$solution
    }
    
    
    
  }
  #w.opt <- w.opt/rowSums(abs(w.opt))
  w.opt[is.na(w.opt)] = 0
  w.opt
}


riskParity <- function(charlist, CovList, vols, targetvol = 10/100, from = 1, show = TRUE){
  vols[vols<0.01]<-0.01
  
  weight.dollar <- ntop(vols, ncol(vols))
  
  if(show){
    op <- pboptions(type="txt")
  }else {
    op <- pboptions(type="none")
  }
  
  out <- pblapply(charlist, function(x){
    
    w <- risk.parity.strategy(x, CovList, vols, targetvol = targetvol,from = from)
    w/weight.dollar
  })
  out
}
riskParity2 <- function(charlist, CovList, vols, targetvol = 10/100, from = 1, max.portfolio.leverage = 100/100, show = TRUE){
  vols[vols<0.01]<-0.01
  
  weight.dollar <- ntop(vols, ncol(vols))
  
  if(show){
    op <- pboptions(type="txt")
  }else {
    op <- pboptions(type="none")
  }
  
  out <- pblapply(charlist, function(x){
    
    w <- risk.parity.strategy2(x, CovList, vols, targetvol = targetvol,max.portfolio.leverage=max.portfolio.leverage,from = from)
    w/weight.dollar
  })
  out
}
riskParity0 <- function(charlist, CovList, vols, targetvol = 10/100, from = 1, show = TRUE){
  vols[vols<0.01]<-0.01
  
  weight.dollar <- ntop(vols, ncol(vols))
  
  if(show){
    op <- pboptions(type="txt")
  }else {
    op <- pboptions(type="none")
  }
  
  out <- pblapply(charlist, function(x){
    
    w <- risk.parity.strategy0(x, CovList, vols, targetvol = targetvol,from = from)
    w/weight.dollar
  })
  out
}

getCorrFactor <- function(signals, rets , lookback.len = 90, from){
  
  dates <- index(rets)
  if(is.character(from)) from <- which(dates ==  dates[dates>=from][1])

  corrFactor <- rets
  corrFactor[] <- 1
  

  
  for(i in from:nrow(rets)){
    
    


    idx <- which(signals[i,]!=0)

    if(length(idx) > 2){
      
      signedreturns <- t(t(rets[(i- lookback.len +1):i,]) * as.numeric(signals[i,]))
      
      corrMat <- cor(signedreturns[,idx], use="pairwise.complete.obs")
      #corrMat <-  cov2cor(tawny::cov.shrink(signedreturns[,idx]))
      #corrMat <-  corpcor::cor.shrink(signedreturns[,idx], verbose=FALSE)
      #corrMat <- cor(rets[(i- lookback.len +1):i,idx], use="pairwise.complete.obs")
      
      avcor <- mean(corrMat[lower.tri(corrMat)])
      
      n <- ncol(corrMat)
      
      avcori <- (colSums(corrMat) - 1)/(n - 1)
      
      #corrFactor[i,idx] <- n/sqrt(n + 2*(n - 1)*avcori +(n - 1)*(n - 2)*avcor)
      
      #corrFactor[i,idx] <- sqrt(n/(sum(signals[i,]*signals[i,])/n + (n - 1)*avcor))
      corrFactor[i,idx] <- sqrt(n/(1 + (n - 1)*avcor))
    }

  }
  
  corrFactor
}


risk.Weighted.Corr.strategy <- function(weight, vols, corrFactor,target = 10/100){
  vols[vols<0.01]<-0.01
  weight.risk <- target*weight*corrFactor / vols

  #weight.risk <- rowSums(abs(weight))*t(apply(weight.risk,1,normalizeVec))
  weight.risk
}


riskWeightedCorr <- function(charlist, rets , vols, target = 10/100, lookback.len = 90, from){
  
  
  weight.dollar <- ntop(rets, ncol(rets))
  
  out <- lapply(charlist, function(x){
    #w <- weight.dollar*x
    w <- x
    corrFactor <- getCorrFactor(x, rets, from = from)
    w <- risk.Weighted.Corr.strategy(w, vols, corrFactor,target=target)
    
    #w/weight.dollar
  })
  out
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

getRollingCov <- function(rets, lookback.len = 90, from = 1){
  
  windw <- embed(1:nrow(rets),lookback.len)
  windw <- windw[,lookback.len:1]
  
  l <- vector(mode="list", length = nrow(windw))
  
  for(i in from:nrow(rets)){
    if(i>=lookback.len){
      l[[i]] <- cov(rets[windw[i - lookback.len +1,],],use="pairwise.complete.obs")*252
      #l[[i]] <- tawny::cov.shrink(rets[windw[i - lookback.len +1,],])*252
      #l[[i]] <- corpcor::cov.shrink(rets[windw[i - lookback.len +1,],], verbose=FALSE)*252
    }
  }
  l
}
volatility <- function (OHLC, n = 10, calc = "close", N = 260, mean0 = FALSE, 
                        ...) 
{
  OHLC <- try.xts(OHLC, error = as.matrix)
  calc <- match.arg(calc, c("close", "absolute","garman.klass", "parkinson", 
                            "rogers.satchell", "gk.yz", "yang.zhang"))
  if (calc == "close") {
    if (NCOL(OHLC) == 1) {
      r <- ROC(OHLC[, 1], 1, ...)
    }
    else {
      r <- ROC(OHLC[, 4], 1, ...)
    }
    if (isTRUE(mean0)) {
      s <- sqrt(N) * sqrt(runSum(r^2, n - 1)/(n - 2))
    }
    else {
      s <- sqrt(N) * runSD(r, n - 1)
    }
  }
  if (calc == "absolute") {
    if (NCOL(OHLC) == 1) {
      r <- ROC(OHLC[, 1], 1, ...)
    }
    else {
      r <- ROC(OHLC[, 4], 1, ...)
    }
    
    s <- sqrt(N) * runMean(abs(r), n - 1)
    
  }
  if (calc == "garman.klass") {
    s <- sqrt(N/n * runSum(0.5 * log(OHLC[, 2]/OHLC[, 3])^2 - 
                             (2 * log(2) - 1) * log(OHLC[, 4]/OHLC[, 1])^2, n))
  }
  if (calc == "parkinson") {
    s <- sqrt(N/(4 * n * log(2)) * runSum(log(OHLC[, 2]/OHLC[, 
                                                             3])^2, n))
  }
  if (calc == "rogers.satchell") {
    s <- sqrt(N/n * runSum(log(OHLC[, 2]/OHLC[, 4]) * log(OHLC[, 
                                                               2]/OHLC[, 1]) + log(OHLC[, 3]/OHLC[, 4]) * log(OHLC[, 
                                                                                                                   3]/OHLC[, 1]), n))
  }
  if (calc == "gk.yz") {
    if (is.xts(OHLC)) {
      Cl1 <- lag.xts(OHLC[, 4])
    }
    else {
      Cl1 <- c(NA, OHLC[-NROW(OHLC), 4])
    }
    s <- sqrt(N/n * runSum(log(OHLC[, 1]/Cl1)^2 + 0.5 * log(OHLC[, 
                                                                 2]/OHLC[, 3])^2 - (2 * log(2) - 1) * log(OHLC[, 4]/OHLC[, 
                                                                                                                         1])^2, n))
  }
  if (calc == "yang.zhang") {
    if (is.xts(OHLC)) {
      Cl1 <- lag.xts(OHLC[, 4])
    }
    else {
      Cl1 <- c(NA, OHLC[-NROW(OHLC), 4])
    }
    dots <- list(...)
    if (is.null(dots$alpha)) {
      alpha <- 1.34
    }
    if (is.null(dots$k)) {
      k <- (alpha - 1)/(alpha + (n + 1)/(n - 1))
    }
    s2o <- N * runVar(log(OHLC[, 1]/Cl1), n = n)
    s2c <- N * runVar(log(OHLC[, 4]/OHLC[, 1]), n = n)
    s2rs <- volatility(OHLC = OHLC, n = n, calc = "rogers.satchell", 
                       N = N, ...)
    s <- sqrt(s2o + k * s2c + (1 - k) * (s2rs^2))
  }
  reclass(s, OHLC)
}
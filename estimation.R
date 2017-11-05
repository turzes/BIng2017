
# compute weight in PPP
computeWeight <- function(theta, charlist, BenchW){  
  #w <- lapply(1:length(theta),function(i){theta[i]*charlist[[i]]})
  #w <- Reduce("+", w)
  w <- matrix(0,nrow(BenchW),ncol(BenchW))
  for(i in 1:length(charlist)){
     w <- w + theta[i]*charlist[[i]]        
  }
  #Nt <- rowSums(BenchW != 0, na.rm=TRUE)
  Nt <- ncol(BenchW)
  w <- (1/Nt)*w
  w<-BenchW + w
  w
}

computeCombineSignals <- function(theta, charlist){  

  w <- lapply(1:length(theta),function(i){theta[i]*charlist[[i]]})
  w <- Reduce("+", w)
  w
}

getReturns <-function(theta, charlist, BenchW, RetData){
  #Nt <- rowSums(BenchW != 0)
  Nt <- ncol(BenchW)
  sumthetaXj <- lapply(1:length(theta),function(i,theta,x){theta[i]*x[[i]]},theta=theta,x=charlist)
  r <- rowSums((BenchW + Reduce("+", sumthetaXj)*(1/Nt))*RetData, na.rm=TRUE)
  r
}


getReturns2 <-function(theta, charlist, benchRet,invnTRet){
  
  sumthetaXj <- lapply(1:length(theta),function(i){theta[i]*charlist[[i]]*invnTRet})
  r <- rowSums(benchRet + Reduce("+", sumthetaXj), na.rm=TRUE)
  r
}



MeanVarianceRversion <- function(theta, charCovMat, covbenchVec, charMeanVec, gamma, invthetaW, lambda2){
  
  theta <- theta/invthetaW
  
  gamma*0.5*rbind(theta)%*%charCovMat%*%cbind(theta) + gamma*rbind(theta)%*%cbind(covbenchVec) - rbind(theta)%*%cbind(charMeanVec) + lambda2*sum(theta*theta)
}


MeanVarianceGradientRversion <- function(theta, charCovMat, covbenchVec, charMeanVec, gamma, invthetaW, lambda2){
  
  theta <- theta/invthetaW
  
  grvec <- gamma*charCovMat%*%cbind(theta) + gamma*covbenchVec - charMeanVec + 2*lambda2*theta;
  
  grvec/invthetaW
}

MeanVarianceCostRversion <- function(theta, charCovMat, covbenchVec, charMeanVec, CostData, charChangelist, costCharChangelist, BenchWChange, gamma, invthetaW, lambda2){
  
  theta <- theta/invthetaW
  
  wtheta <- lapply(1:length(theta),function(i){theta[i]*charChangelist[[i]]})
  wtheta <- Reduce("+", wtheta)
  wtheta <- wtheta + BenchWChange
  
  gamma*0.5*rbind(theta)%*%charCovMat%*%cbind(theta) + gamma*rbind(theta)%*%covbenchVec - rbind(theta)%*%charMeanVec + mean(rowSums(CostData*(sqrt(1e-10 + wtheta*wtheta) - 1e-05))) + lambda2*sum(theta*theta)
}


MeanVarianceGradCostRversion <- function(theta, charCovMat, covbenchVec, charMeanVec, CostData, charChangelist, costCharChangelist, BenchWChange, gamma, invthetaW, lambda2){
  
  theta <- theta/invthetaW
  
  wtheta <- lapply(1:length(theta),function(i){theta[i]*charChangelist[[i]]})
  wtheta <- Reduce("+", wtheta)
  wtheta <- wtheta + BenchWChange
  
  grvec = gamma*charCovMat%*%cbind(theta) + gamma*covbenchVec - charMeanVec + 2*lambda2*theta
  
  wtheta = wtheta/sqrt(1e-10 + wtheta*wtheta)
  
  for(i in 1:length(theta)){
    
    grvec[i] = (grvec[i] + mean(rowSums(wtheta*costCharChangelist[[i]])))/invthetaW[i];
  }
  
  grvec
}

  
softthresholdingRversion <- function(lambda1, lambda2, gamma, rhoj, deltaj, sigma2j, thetaW){
  
  a = -rhoj/(gamma*sigma2j + 2*lambda2)
  
  b = (deltaj + lambda1*thetaW)/(gamma*sigma2j + 2*lambda2)
  
  if(a > b){
    
    out <- a - b
  }
  else if(a < -b){
    
    out <- a + b
  }
  else{
    
    out <- 0
  }
  out
}



ccdPPPMeanVarianceRversion <- function(theta, charCovMat, covBenchVec, charMeanVec, gamma, lambda1, lambda2, invthetaW, tol, maxiter){
  
  pre_theta = theta
  
  sigma2j = diag(charCovMat)
  diag(charCovMat) <- 0
  
  rho = gamma*charCovMat%*%cbind(theta) + gamma*covBenchVec - charMeanVec
  
  for(j in 1:maxiter){
    
    for(i in 1:length(theta)){
      
      
      theta[i] = softthresholdingRversion(lambda1, lambda2,gamma, rho[i], 0, sigma2j[i], invthetaW[i]);
      
      rho = rho + (theta[i] - pre_theta[i])*gamma*charCovMat[,i];
    }
    
    
    if (norm(cbind(pre_theta - theta),"F") < tol) { break }
    pre_theta <- theta
    
  }
  
  list(theta)
}



ccdPPPMeanVarianceCostApproxRversion <- function(theta, charCovMat, covBenchVec, charMeanVec, charChangeCostMeanVec, gamma, lambda1, lambda2, invthetaW, tol, maxiter){
  
  pre_theta = theta
  
  sigma2j = diag(charCovMat)
  diag(charCovMat) <- 0
  
  rho = gamma*charCovMat%*%cbind(theta) + gamma*covBenchVec - charMeanVec
  
  for(j in 1:maxiter){
    
    for(i in 1:length(theta)){
      
      
      theta[i] = softthresholdingRversion(lambda1, lambda2,gamma, rho[i], charChangeCostMeanVec[i], sigma2j[i], invthetaW[i]);
      
      rho = rho + (theta[i] - pre_theta[i])*gamma*charCovMat[,i];
    }
    
    
    if (norm(cbind(pre_theta - theta),"F") < tol) { break }
    pre_theta <- theta
    
  }
  
  list(theta)
}



CRRA <-function(theta, charlist, BenchW, RetData, gamma) {
  U<-0
  #Nt <- rowSums(BenchW != 0) 
  Nt <- ncol(BenchW)
  sumthetaXj <- lapply(1:length(theta),function(i,theta,x){theta[i]*x[[i]]},theta=theta,x=charlist)
  r <- rowSums((BenchW + Reduce("+", sumthetaXj)*(1/Nt))*RetData, na.rm=TRUE)
  
  
  u<-(((1+r)^(1-gamma))/(1-gamma))
  
  U<--mean(u)
  
  U
}

computeMeanUtility <- function(retTheta, gamma){
  out <-  ((1+retTheta)^(1-gamma))/(1-gamma)
  -mean(out)
}

computeSharpeRatio <- function(retTheta, annualfactor = 252){
  
  -annualfactor*mean(retTheta, na.rm=TRUE)/(sqrt(annualfactor)*sd(retTheta, na.rm=TRUE))
}

computeSZ3 <- function(retTheta, gamma){
  
  sr <- mean(retTheta, na.rm=TRUE)/sd(retTheta, na.rm=TRUE)
  skew <- skewness(retTheta)
  
  sz <- (((1 + retTheta)^( - gamma))/(gamma*((1 + retTheta)^( -1 - gamma))))*((sr^2)/2 + gamma*(1 + gamma)*(sr^3)*skew/6)
  -mean(sz)
}

computeSZ4 <- function(retTheta, gamma){
  
  sr <- mean(retTheta, na.rm=TRUE)/sd(retTheta, na.rm=TRUE)
  skew <- skewness(retTheta, na.rm=TRUE)
  kurt <- kurtosis(retTheta, na.rm=TRUE)
  
  sz <- (((1 + retTheta)^( - gamma))/(gamma*((1 + retTheta)^( -1 - gamma))))*((sr^2)/2 + gamma*(1 + gamma)*(sr^3)*skew/6 - gamma*(1 + gamma)*(2 + gamma)*(sr^4)*(kurt - 3)/24)
  -mean(sz)
}

computeCEG <- function(retTheta, gamma, annualfactor = 252){
  
  u = (1 + retTheta)^(1 - gamma)
  out = mean(u, na.rm=TRUE)^(1/(1 - gamma)) - 1
  -annualfactor*out
  
}

getlambda1 <- function(charRetMat, benchRet, CostData, charChangelist, BenchWChange, costCharChangelist, startdate, enddate, gamma, alpha, minRatio, nLambda1, optcost){
  
  np <- ncol(charRetMat)
  thetaW <- numeric(np) + 1
  
  charCovMat <- cov(charRetMat[startdate:enddate,],use="pairwise.complete.obs")
  covbenchVec <- sapply(1:np, function(i){
    return(cov(benchRet[startdate:enddate],charRetMat[startdate:enddate,i]))
  })
  
  covbenchVec <- matrix(covbenchVec, ncol=1)
  
  charMeanVec <- colMeans(charRetMat[startdate:enddate,])
  charMeanVec <- matrix(charMeanVec, ncol=1)
  charChangelist.train <- setupListRange(charChangelist,startdate:enddate)
  costCharChangelist.train<- setupListRange(costCharChangelist,startdate:enddate)
  
  if(optcost == "EXACT_COST"){
    
    
    
     if(alpha != 0){
      ans.opt <- lbfgs::lbfgs(MeanVarianceCostRversion,MeanVarianceGradCostRversion, numeric(np),charCovMat=charCovMat, covbenchVec=covbenchVec, charMeanVec= charMeanVec,CostData=CostData[startdate:enddate,], charChangelist=charChangelist.train, costCharChangelist=costCharChangelist.train, BenchWChange=BenchWChange[startdate:enddate,], gamma=gamma, invthetaW = numeric(np) + 1, lambda2=0,invisible = 1,epsilon = 1e-05,linesearch_algorithm="LBFGS_LINESEARCH_BACKTRACKING")
      
      thetaW <- ans.opt$par
     }

      dg <-   MeanVarianceGradCostRversion(numeric(np), charCovMat, covbenchVec, charMeanVec, CostData[startdate:enddate,], charChangelist.train, costCharChangelist.train, BenchWChange[startdate:enddate,], gamma, numeric(np) + 1, 0)
  }else{
    if(alpha != 0){
      ans.opt <- lbfgs::lbfgs(MeanVarianceRversion,MeanVarianceGradientRversion, numeric(np),charCovMat=charCovMat, covbenchVec=covbenchVec, charMeanVec= charMeanVec, gamma=gamma, invthetaW = numeric(np) + 1, lambda2=0,invisible = 1,epsilon = 1e-05,linesearch_algorithm="LBFGS_LINESEARCH_BACKTRACKING")
      
      thetaW <- ans.opt$par
    }
    
    dg <-   MeanVarianceGradientRversion(numeric(np), charCovMat, covbenchVec, charMeanVec, gamma, numeric(np) + 1, 0)
    
  }

   maxLambda <- max(abs(dg)*(abs(thetaW) + 1/length(startdate:enddate))^alpha)

  maxLambda <- maxLambda
  if(is.null(minRatio)) minLambda <- 1e-18
  else minLambda <- maxLambda*minRatio
  lambda1 <- 10^(seq(log10(maxLambda), log10(minLambda),-(log10(maxLambda)-log10(minLambda))/(nLambda1-1)))
  lambda1
  
}


setupTrainData <- function(startdate,enddate,trainDataRatio){
  
  nobs <- enddate - startdate + 1
  trainstartdate <- startdate
  trainenddate <- trainstartdate + floor(nobs*trainDataRatio[1])
  teststartdate <- trainenddate + 1
  testenddate <- enddate
  c(trainstartdate, trainenddate, teststartdate, testenddate)
}



getPrepEst <- function(charlist, BenchW, RetData){
  
  np <- length(charlist)
  #nT <- rowSums(coredata(BenchW) != 0) 
  #
  nT <- ncol(BenchW)
  #nT = repmat(matrix(nT,ncol = 1),1,ncol(BenchW))
  
  benchRet <- coredata(BenchW)*coredata(mlag(RetData,-1))
  invnTRet <- (1/nT)*coredata(mlag(RetData,-1))
  
  benchRet <- rowSums(benchRet,na.rm=TRUE)
  
  charRetMat <- matrix(0,nrow = nrow(BenchW),ncol = np)
  
  
  for(i in 1:np){
    
    charRetMat[,i] <- rowSums(charlist[[i]]*invnTRet,na.rm=TRUE)
  }
  
  list(charRetMat,benchRet)
}

sortino.ratio = function(R,MAR = 0) {
  r = R[R < MAR]
  out <- (mean(R, na.rm=TRUE) - MAR) / sqrt(sum((MAR - r)^2, na.rm=TRUE) / length(r) ) 
  -out
}

omega.ratio = function(R,MAR = 0) {
  #sum(pmax(R - MAR, 0)) / sum(pmax(MAR - R, 0))
  r = R - MAR
  out <- - sum(r[r > 0], na.rm=TRUE) / sum(r[r < 0], na.rm=TRUE)
  -out
}



compute.mean <- function(x){
  
  temp = compute.annual.factor(x)
  x = as.vector(coredata(x))
  return(temp*mean(x))
  
}

compute.exposure <- function(weight,nsize = NULL)
{
  if(is.null(nsize)) nsize <- nrow(weight)
  
  sum( apply(weight, 1, function(x) sum(x != 0) ) != 0 ) / nsize
}

compute.risk<-function(x)
{
  temp = compute.annual.factor(x)
  x = as.vector(coredata(x))
  return( sqrt(temp)*sd(x) )
}
compute.sharpe <- function(x)
{
  temp = compute.annual.factor(x)
  
  x = as.vector(coredata(x))
  return(sqrt(temp) * mean(x)/sd(x) )
}

compute.drawdown <- function(x)
{
  return(x / cummax(c(1,x))[-1] - 1)
}

compute.max.drawdown <- function(x)
{
  as.double( min(compute.drawdown(x)) )
}

compute.avg.drawdown <- function(x)
{
  drawdown = c( 0, compute.drawdown(coredata(x)), 0 )
  dstart = which( drawdown == 0 & mlag(drawdown, -1) != 0 )
  dend = which(drawdown == 0 & mlag(drawdown, 1) != 0 )
  drawdowns = apply( cbind(dstart, dend), 1, function(x) min(drawdown[ x[1]:x[2] ], na.rm=T) )
  mean(drawdowns)
}

bt.summary.1 <- function(bt){
  
  
  out = list()
  out$Period = join( format( range(index.xts(bt$equity)), '%b%Y'), ' - ')
  #print('---')
  out$AveReturns = compute.mean(bt$ret)
  
  out$Skewness = skewness(bt$ret)
  out$Volatility = compute.risk(bt$ret)
  out$Sharpe = compute.sharpe(bt$ret)
  out$Ceg = bt$ceg
  out$MaxDD = compute.max.drawdown(bt$equity)
  out$AveDD = compute.avg.drawdown(bt$equity)
  
  
  month.ends = endpoints(bt$equity, 'months')
  mret = ROC(bt$equity[month.ends,], type = 'discrete')
  
  out$Win.Percent.Month = sum(mret > 0, na.rm = T) / len(mret)
  out$Best.Month = max(mret, na.rm = T)
  out$Worst.Month = min(mret, na.rm = T)
  out$AveNetLeverage = mean(rowSums(bt$weight))
  out$AveGrossLeverage = mean(rowSums(abs(bt$weight)))
  out$MaxNetLeverage = max(rowSums(bt$weight))
  out$MaxGrossLeverage = max(rowSums(abs(bt$weight)))
  #out$AveAnnualTurnover = compute.turnover(bt,b = NULL)
  
  #format(x,digits = 3)
  
  #out <- lapply(out, function(x) if(is.double(x)) round(x,3) else x)
  out <- lapply(out, function(x) if(is.double(x)) format(round(x, 3), nsmall = 2) else x)
  out
}





plot.ts.bt <- function(ts.bt,tilttxt = "Coefficient estimates"){
  exampleData <- data.frame(index=1:nrow(ts.bt$coefficients), ts.bt$coefficients)
  sectorMelt <- melt(exampleData, id = "index")
  colnames(sectorMelt) <- c("Time","Signal","Value")
  print(ggplot(data = sectorMelt, aes(x = Time, y = Value,group = Signal, colour = Signal))+geom_line(size=1) + theme_grey(base_size = 12) + ggtitle(tilttxt))
  
}

plotbtPPP <- function(bt,tilttxt = 'Peformance comparison'){
  
  layout(1:2)
  
  yrange <- range(bt$benchmark$equity,bt$ppp$equity)
  
  plota(bt$benchmark$equity,type = 'l',col = 'blue',  ylim = yrange, main = tilttxt)
  
  plota.lines(bt$ppp$equity, col='red')
  plota.legend(c('PPP','Benchmark'), 'red,blue')
  
  temp1 <- list2matrix(bt$ppp.summary )
  temp2 <- list2matrix(bt$benchmark.summary)
  temp <- cbind(temp1,temp2)
  temp = rbind(c('PPP','Benchmark'), temp)   # add header row
  plot.table(temp)
  
}

plotbt1 <- function(bt,tilttxt = 'Peformance comparison'){
  
  yrange<- range(bt[[1]]$benchmark$equity)
  legendnames <- 'Benchmark'
  
  temptable <- list2matrix(bt[[1]]$benchmark.summary)
  
  for(i in 1:length(bt)){
    
    yrange <- range(yrange,bt[[i]]$ppp$equity)
    legendnames <- c(legendnames,names(bt)[i])
    temptable <- cbind(temptable,list2matrix(bt[[i]]$ppp.summary))
  }
  temptable = rbind(legendnames, temptable)
  
  layout(1:2)
  
  plota(bt[[1]]$benchmark$equity,type = 'l',col = 1,ylim = yrange, main = tilttxt)
  
  for(i in 1:length(bt)){
    
    plota.lines(bt[[i]]$ppp$equity,col = i + 1)
  }
  plota.legend(legendnames, 1:(length(bt) + 1))
  
  plot.table(temptable)
  
}
#######################
bootPPPAENET <- function(charlist, rets, BenchW, charRetMat, benchRet, CostData, charChangelist, BenchWChange, costCharChangelist, charChangelist2, trainIndex,testIndex = NULL, gamma, lambda1, lambda2, alpha, optcost = "NO_COST", n.sim = 500, tol=1e-9, maxiter=2000){
  
  original.estimates <- simpleGS2(theta = numeric(length(charlist)), charlist, rets, BenchW, charRetMat, benchRet, CostData, charChangelist, BenchWChange, costCharChangelist, charChangelist2, trainIndex,NULL, gamma=gamma, lambda1=lambda1, lambda2=lambda2, alpha = alpha, tol=tol, maxiter=maxiter, optcost = optcost,compute.rets = FALSE)
  
  original.estimates <- as.numeric(original.estimates[[1]])
  set.seed(999)
  
  
  store.matrix <- matrix(NA, nrow=n.sim,ncol=length(charlist))
  
  for(i in 1:n.sim){
    
    trainIndex.sim <- sample(trainIndex,length(trainIndex),replace=T)
    ans<-simpleGS2(theta = numeric(length(charlist)), charlist, rets, BenchW, charRetMat, benchRet, CostData, charChangelist, BenchWChange, costCharChangelist, charChangelist2, trainIndex.sim,NULL, gamma=gamma, lambda1=lambda1, lambda2=lambda2, alpha = alpha, tol=tol, maxiter=maxiter, optcost = optcost,compute.rets = FALSE)
    store.matrix[i,] <- as.numeric(ans[[1]])
    
  }
  
  boot.means <- colMeans(store.matrix, na.rm=T)
  
  boot.medians <- apply(store.matrix,2,median, na.rm=T)
  
  boot.sds <- apply(store.matrix,2,sd, na.rm=T)
  
  boot.bias <- colMeans(store.matrix, na.rm=T) - original.estimates
  
  conf.mat <- matrix(apply(store.matrix, 2 ,quantile, c(0.025, 0.975), na.rm=T),
                     ncol=2, byrow=TRUE)
  colnames(conf.mat) <- c("95%-CI Lower", "95%-CI Upper")
  
  summary.frame <- data.frame(org=original.estimates,mean=boot.means, median=boot.medians,
                              sd=boot.sds, bias=boot.bias, "CI_lower"=conf.mat[,1], "CI_upper"=conf.mat[,2])
  return(summary.frame)
}

pppest <- function(thetatmp0, thetatmp, charCovMat, covbenchVec, charMeanVec, charChangeCostMeanVec, CostData, charChangelist, costCharChangelist,BenchWChange, trainIndex, gamma, lambda1, lambda2, alpha, invthetaW, optcost, tol, maxiter){
  
  np <- length(thetatmp0)
  
  if(optcost == "NO_COST"){
    
    if(alpha != 0){
      
      #ans.opt<-lbfgs::lbfgs(MeanVariance,MeanVarianceGrad, thetatmp0,charCovMat=charCovMat, covbenchVec=covbenchVec, charMeanVec= charMeanVec, gamma=gamma, invthetaW = numeric(np) + 1, lambda2=lambda2,invisible = 1,epsilon = 1e-09,linesearch_algorithm="LBFGS_LINESEARCH_BACKTRACKING",orthantwise_c = lambda1[i],orthantwise_start = 0,orthantwise_end=np)
      #invthetaW <- 1/((abs(as.numeric(ans.opt$par)) + 1/length(trainIndex))^alpha)
      #thetatmp0 <- as.numeric(ans.opt$par)
      
      ans.opt <- ccdPPPMeanVariance(thetatmp0, charCovMat, covbenchVec, charMeanVec, gamma, lambda1, lambda2, numeric(np) + 1, tol, maxiter)
      
      invthetaW <- 1/((abs(ans.opt[[1]]) + 1/length(trainIndex))^alpha)
      thetatmp0 <- as.numeric(ans.opt[[1]])
      
    }
    
    
    #fit <-lbfgs::lbfgs(MeanVariance,MeanVarianceGrad, thetatmp,charCovMat=charCovMat, covbenchVec=covbenchVec, charMeanVec= charMeanVec, gamma=gamma, invthetaW = invthetaW, lambda2=lambda2,invisible = 1,epsilon = 1e-09,linesearch_algorithm="LBFGS_LINESEARCH_BACKTRACKING",orthantwise_c = lambda1[i],orthantwise_start = 0,orthantwise_end=np)
    #thetatmp <- fit$par/invthetaW
    
    fit <- ccdPPPMeanVariance(thetatmp, charCovMat, covbenchVec, charMeanVec, gamma, lambda1, lambda2, invthetaW, tol, maxiter)
    
    thetatmp <- as.numeric(fit[[1]])
    
  }
  
  if(optcost == "EXACT_COST"){
    if(alpha != 0){
      
      ans.opt<-lbfgs::lbfgs(MeanVarianceCost, MeanVarianceGradCost, thetatmp0,charCovMat=charCovMat, covbenchVec=covbenchVec, charMeanVec= charMeanVec,CostData=CostData, charChangelist=charChangelist, costCharChangelist=costCharChangelist, BenchWChange=BenchWChange, gamma=gamma, invthetaW = numeric(np) + 1, lambda2=lambda2,invisible = 1,epsilon = 1e-05,linesearch_algorithm="LBFGS_LINESEARCH_BACKTRACKING",orthantwise_c = lambda1,orthantwise_start = 0,orthantwise_end=np)
      
      invthetaW <- 1/((abs(as.numeric(ans.opt$par)) + 1/length(trainIndex))^alpha)
      thetatmp0 <- as.numeric(ans.opt$par)
    }
    
    fit <- lbfgs::lbfgs(MeanVarianceCost,MeanVarianceGradCost, thetatmp,charCovMat=charCovMat, covbenchVec=covbenchVec, charMeanVec= charMeanVec,CostData=CostData, charChangelist=charChangelist, costCharChangelist=costCharChangelist,BenchWChange=BenchWChange, gamma=gamma, invthetaW = invthetaW, lambda2=lambda2,invisible = 1,epsilon = 1e-05,linesearch_algorithm="LBFGS_LINESEARCH_BACKTRACKING",orthantwise_c = lambda1,orthantwise_start = 0,orthantwise_end=np)
    thetatmp <- fit$par/invthetaW
    
    #thetatmp <- (1 + lambda2/length(trainIndex))*thetatmp
  }
  
  if(optcost == "EXTRA_COST"){
    if(alpha != 0){
      ans.opt <- ccdPPPMeanVarianceCostApprox(thetatmp0, charCovMat, covbenchVec, charMeanVec, charChangeCostMeanVec, gamma, lambda1, lambda2, numeric(np) + 1, tol, maxiter)
      
      invthetaW <- 1/((abs(ans.opt[[1]]) + 1/length(trainIndex))^alpha)
      thetatmp0 <- as.numeric(ans.opt[[1]])
    }
    
    
    fit <- ccdPPPMeanVarianceCostApprox(thetatmp, charCovMat, covbenchVec, charMeanVec, charChangeCostMeanVec, gamma, lambda1, lambda2, invthetaW, tol, maxiter)
    
    thetatmp <- as.numeric(fit[[1]])
    
  }
  
  list(thetatmp0, thetatmp)
}

pppestRversion <- function(thetatmp0, thetatmp, charCovMat, covbenchVec, charMeanVec, charChangeCostMeanVec, CostData, charChangelist, costCharChangelist,BenchWChange, trainIndex, gamma, lambda1, lambda2, alpha, invthetaW, optcost, tol, maxiter){
  
  np <- length(thetatmp0)
  
  if(optcost == "NO_COST"){
    
    if(alpha != 0){
      
      #ans.opt<-lbfgs::lbfgs(MeanVariance,MeanVarianceGrad, thetatmp0,charCovMat=charCovMat, covbenchVec=covbenchVec, charMeanVec= charMeanVec, gamma=gamma, invthetaW = numeric(np) + 1, lambda2=lambda2,invisible = 1,epsilon = 1e-09,linesearch_algorithm="LBFGS_LINESEARCH_BACKTRACKING",orthantwise_c = lambda1[i],orthantwise_start = 0,orthantwise_end=np)
      #invthetaW <- 1/((abs(as.numeric(ans.opt$par)) + 1/length(trainIndex))^alpha)
      #thetatmp0 <- as.numeric(ans.opt$par)
      
      ans.opt <- ccdPPPMeanVarianceRversion(thetatmp0, charCovMat, covbenchVec, charMeanVec, gamma, lambda1, lambda2, numeric(np) + 1, tol, maxiter)
      
      invthetaW <- 1/((abs(ans.opt[[1]]) + 1/length(trainIndex))^alpha)
      thetatmp0 <- as.numeric(ans.opt[[1]])
      
    }
    
    
    #fit <-lbfgs::lbfgs(MeanVariance,MeanVarianceGrad, thetatmp,charCovMat=charCovMat, covbenchVec=covbenchVec, charMeanVec= charMeanVec, gamma=gamma, invthetaW = invthetaW, lambda2=lambda2,invisible = 1,epsilon = 1e-09,linesearch_algorithm="LBFGS_LINESEARCH_BACKTRACKING",orthantwise_c = lambda1[i],orthantwise_start = 0,orthantwise_end=np)
    #thetatmp <- fit$par/invthetaW
    
    fit <- ccdPPPMeanVarianceRversion(thetatmp, charCovMat, covbenchVec, charMeanVec, gamma, lambda1, lambda2, invthetaW, tol, maxiter)
    
    thetatmp <- as.numeric(fit[[1]])
    
  }
  
  if(optcost == "EXACT_COST"){
    if(alpha != 0){
      
      ans.opt<-lbfgs::lbfgs(MeanVarianceCostRversion, MeanVarianceGradCostRversion, thetatmp0,charCovMat=charCovMat, covbenchVec=covbenchVec, charMeanVec= charMeanVec,CostData=CostData, charChangelist=charChangelist, costCharChangelist=costCharChangelist, BenchWChange=BenchWChange, gamma=gamma, invthetaW = numeric(np) + 1, lambda2=lambda2,invisible = 1,epsilon = 1e-05,linesearch_algorithm="LBFGS_LINESEARCH_BACKTRACKING",orthantwise_c = lambda1,orthantwise_start = 0,orthantwise_end=np)
      
      invthetaW <- 1/((abs(as.numeric(ans.opt$par)) + 1/length(trainIndex))^alpha)
      thetatmp0 <- as.numeric(ans.opt$par)
    }
    
    fit <- lbfgs::lbfgs(MeanVarianceCostRversion,MeanVarianceGradCostRversion, thetatmp,charCovMat=charCovMat, covbenchVec=covbenchVec, charMeanVec= charMeanVec,CostData=CostData, charChangelist=charChangelist, costCharChangelist=costCharChangelist,BenchWChange=BenchWChange, gamma=gamma, invthetaW = invthetaW, lambda2=lambda2,invisible = 1,epsilon = 1e-05,linesearch_algorithm="LBFGS_LINESEARCH_BACKTRACKING",orthantwise_c = lambda1,orthantwise_start = 0,orthantwise_end=np)
    thetatmp <- fit$par/invthetaW
    
    #thetatmp <- (1 + lambda2/length(trainIndex))*thetatmp
  }
  
  if(optcost == "EXTRA_COST"){
    if(alpha != 0){
      ans.opt <- ccdPPPMeanVarianceCostApproxRversion(thetatmp0, charCovMat, covbenchVec, charMeanVec, charChangeCostMeanVec, gamma, lambda1, lambda2, numeric(np) + 1, tol, maxiter)
      
      invthetaW <- 1/((abs(ans.opt[[1]]) + 1/length(trainIndex))^alpha)
      thetatmp0 <- as.numeric(ans.opt[[1]])
    }
    
    
    fit <- ccdPPPMeanVarianceCostApproxRversion(thetatmp, charCovMat, covbenchVec, charMeanVec, charChangeCostMeanVec, gamma, lambda1, lambda2, invthetaW, tol, maxiter)
    
    thetatmp <- as.numeric(fit[[1]])
    
  }
  
  list(thetatmp0, thetatmp)
}


simpleGS2 <- function(theta, charlist, rets, BenchW, charRetMat, benchRet, CostData, charChangelist, BenchWChange, costCharChangelist, charChangelist2, trainIndex,testIndex, gamma, lambda1, lambda2, alpha, tol, maxiter, optcost = "NO_COST", compute.rets = TRUE, usecfun = FALSE){

  np <- ncol(charRetMat)
  outMat <- matrix(0,nrow=length(lambda1), ncol = length(theta))
  
  thetatmp <- theta
  thetatmp0 <- theta
  
  if(compute.rets){
    
    outRet <- matrix(0,nrow=length(lambda1), ncol = length(testIndex))
    outWeight <- vector("list", length(lambda1))
  }
  
  
  invthetaW <- numeric(np) + 1

  
  charChangelistTrain <- setupListRange(charChangelist,trainIndex)
  charChangelistTest <- setupListRange(charChangelist,testIndex)
  
  #charRetSquareColMean <- colMeans(charRetSquareMat[trainIndex,])
  #charRetcharChangeCostColMean <- colMeans(charRetcharChangeCostMat[trainIndex,])
  #charChangeSquareColMean <- colMeans(charChangeSquareMat[trainIndex,])
  
  
  charCovMat <- cov(charRetMat[trainIndex,],use="pairwise.complete.obs")
  covbenchVec <- sapply(1:np, function(i){
    return(cov(benchRet[trainIndex],charRetMat[trainIndex,i]))
  })
  
  covbenchVec <- matrix(covbenchVec, ncol=1)
  
  charMeanVec <- colMeans(charRetMat[trainIndex,])
  charMeanVec <- matrix(charMeanVec, ncol=1)
  charChangelist.train <- setupListRange(charChangelist,trainIndex)
  costCharChangelist.train<- setupListRange(costCharChangelist,trainIndex)
  
  charChangeCostMeanVec <- colMeans(charChangelist2[trainIndex,])
  
  for(i in 1:length(lambda1)){
    
    if(usecfun){
      est <- pppest(thetatmp0, thetatmp, charCovMat, covbenchVec, charMeanVec, charChangeCostMeanVec, CostData[trainIndex,], charChangelist.train, costCharChangelist.train, BenchWChange[trainIndex,], trainIndex, gamma, lambda1[i], lambda2, alpha, invthetaW, optcost, 1e-8, 2000)
      
    }else{
      est <- pppestRversion(thetatmp0, thetatmp, charCovMat, covbenchVec, charMeanVec, charChangeCostMeanVec, CostData[trainIndex,], charChangelist.train, costCharChangelist.train, BenchWChange[trainIndex,], trainIndex, gamma, lambda1[i], lambda2, alpha, invthetaW, optcost, 1e-8, 2000)
   
    }
    
    
    thetatmp0 <- est[[1]]
    thetatmp <- est[[2]]
    
    outMat[i,] <- thetatmp
    
    
    if(compute.rets){
      if(is.null(CostData)){
        totalcost <- 0
      }else{
        wthetalist <- lapply(1:np,function(i){thetatmp[i]*charChangelistTest[[i]]})
        wtheta <- BenchWChange[testIndex,] + Reduce("+", wthetalist)
        totalcost <- rowSums(CostData[testIndex,]*abs(wtheta), na.rm=TRUE)
      }
      
      
      
      charlistTest <- setupListRange(charlist,testIndex)
      
      weight <- computeWeight(thetatmp, charlistTest, BenchW[testIndex,])
      
      #weight <- iif(rowSums(abs(weight)) > max.portfolio.leverage, max.portfolio.leverage,rowSums(abs(weight)))*weight/ rowSums(abs(weight))
    
      weight <- weight/ rowSums(abs(weight))
      
      outRet[i,] <- rowSums(weight*rets[(testIndex+1),], na.rm=TRUE)- totalcost
      
      # print(coredata(rets[(testIndex+1),1]))
      #  print(outRet[i,])

      
      outWeight[[i]] <- weight
    }

    
  }
  
  
  
  if(compute.rets) out <- list(outMat,outRet,matrix(benchRet[testIndex],nrow=1, ncol = length(testIndex)),outWeight)
  else out <- list(outMat)
  out
}


tsonefold2 <- function(folds, charlist, rets, BenchW, charRetMat, benchRet, CostData, charChangelist, BenchWChange, costCharChangelist, charChangelist2,gamma, lambda1, lambda2, alpha, tol, maxiter, optcost,annualfactor,ncores, show = FALSE, cl = NULL,compute.rets = TRUE, usecfun = FALSE){
  

  numtest <- length(folds$train)
  
  
  np <- ncol(charRetMat)
  if(ncores > 1){
    ans <- parLapply(cl, 1:numtest, function(i){
      trainIndex <- folds$train[[i]]
      testIndex <- folds$test[[i]]
      
      out <- simpleGS2(numeric(np),  charlist, rets, BenchW, charRetMat, benchRet, CostData, charChangelist, BenchWChange, costCharChangelist, charChangelist2, trainIndex, testIndex, gamma, lambda1, lambda2, alpha, tol, maxiter, optcost=optcost,compute.rets = compute.rets, usecfun=usecfun)
      
    })
    
  }else{
    
    if(show){
      loopf <- pblapply
    }else{
      loopf <- lapply
    }
    ans <- do.call(loopf,list(1:numtest,function(i){
      
      trainIndex <- folds$train[[i]]
      testIndex <- folds$test[[i]]
      
      
      
      out <- simpleGS2(numeric(np),  charlist, rets, BenchW, charRetMat, benchRet, CostData, charChangelist, BenchWChange, costCharChangelist, charChangelist2, trainIndex, testIndex, gamma, lambda1, lambda2, alpha, tol, maxiter, optcost=optcost,compute.rets = compute.rets, usecfun=usecfun)
      
      out
    }))
  }
  
  
  outThetaMean <- matrix(0, length(lambda1),np)
  outThetaSd <- matrix(0, length(lambda1),np)
  

  for(i in 1:np){
    
    outTheta <- sapply(ans,function(x) x[[1]][,i])
    if(length(outTheta)>1){
      outThetaMean[,i] <- matrix(rowMeans(outTheta),ncol=1)
      outThetaSd[,i] <- matrix(apply(outTheta,1,sd),ncol=1)
    }else{
      outThetaMean[,i] <- matrix(outTheta,ncol=1)
    }

  }
  
  if(compute.rets){
    
    retThetaMat <- lapply(ans,function(x) x[[2]])
    retThetaMat <- Reduce(cbind,retThetaMat)
    
    
    #mar <- max(0,mean(benchRetMat))
    mar <- 0
    
    outAveUtility <- matrix(apply(retThetaMat,1,computeMeanUtility,gamma),ncol = 1)
    outSR <- matrix(apply(retThetaMat,1,computeSharpeRatio,annualfactor),ncol = 1)
    outCEG <- matrix(apply(retThetaMat,1,computeCEG, gamma,annualfactor),ncol = 1)
    outSOR <- matrix(apply(retThetaMat,1,sortino.ratio,mar),ncol = 1)
    outOMG <- matrix(apply(retThetaMat,1,omega.ratio,mar),ncol = 1)  
    outSZ3 <- matrix(apply(retThetaMat,1,computeSZ3, gamma),ncol = 1)
    outSZ4 <- matrix(apply(retThetaMat,1,computeSZ4, gamma),ncol = 1)
    
    outMat <- cbind(outThetaMean,outThetaSd,outAveUtility,outSR,outCEG,outSOR,outOMG,outSZ3,outSZ4)
    
  }else{
    
    outMat <- cbind(outThetaMean,outThetaSd)
  }
 
  outMat
}





pppData <- setClass("pppData",representation(BenchW="xts",
                                             signallist = "list",
                                             prices = "xts",
                                             CostData = 'xts'
))

AENet <- setClass("AENet",representation(data="pppData",
                                         coefs = "matrix",
                                         gamma = "numeric",
                                         alpha = "numeric",
                                         lambda1 = "numeric",
                                         lambda2 = "numeric",
                                         optcost = 'character',
                                         opt.lambda1 = 'numeric',
                                         opt.lambda2 = 'numeric',
                                         startdate = "numeric",
                                         enddate = "numeric",
                                         from = "character",
                                         to = "character",
                                         name = "character"))



.AENetCV <- setClass("AENetCV",representation(nfolds = 'numeric',
                                              type = 'character',
                                              criterion.opt = "numeric",
                                              cv.results = "matrix",
                                              theta.opt = "numeric",
                                              theta.opt.sd="numeric",
                                              trainstartdate="numeric",
                                              trainenddate="numeric",
                                              teststartdate="numeric",
                                              testenddate="numeric",
                                              train.from="character",
                                              train.to="character",
                                              test.from="character",
                                              test.to="character"),
                     contains = 'AENet')


AENetCV <- function(superclass = NULL,...){
  
  if(is.null(superclass)){
    return(.AENetCV(...))
  }else if(is(superclass,"AENet")){
    
    return(.AENetCV(data = superclass@data, coefs = superclass@coefs, gamma = superclass@gamma, alpha = superclass@alpha, lambda1 = superclass@lambda1, lambda2 = superclass@lambda2 , opt.lambda1=superclass@opt.lambda1, opt.lambda2=superclass@opt.lambda2, optcost = superclass@optcost,startdate = superclass@startdate, enddate = superclass@enddate, from = superclass@from,to = superclass@to,...))
  }
}



.AENetBT <- setClass("AENetBT",representation(initialWindow = 'numeric',
                                              fixedWindow = 'logical',
                                              horizon = 'numeric',
                                              skip = 'numeric',
                                              oos = 'logical',
                                              criterion = "matrix",
                                              theta.mean = "numeric",
                                              theta.sd="numeric",
                                              trainstartdate="numeric",
                                              trainenddate="numeric",
                                              teststartdate="numeric",
                                              testenddate="numeric",
                                              ppp = "list",
                                              ppp.summary = "list",
                                              benchmark = "list",
                                              benchmark.summary = "list",
                                              train.from="character",
                                              train.to="character",
                                              test.from="character",
                                              test.to="character"),
contains = 'AENet')

AENetBT <- function(superclass = NULL,...){
  
  if(is.null(superclass)){
    return(.AENetBT(...))
  }else if(is(superclass,"AENet")){
    
    return(.AENetBT(data = superclass@data, coefs = superclass@coefs, gamma = superclass@gamma, alpha = superclass@alpha, lambda1 = superclass@lambda1, lambda2 = superclass@lambda2 , opt.lambda1=superclass@opt.lambda1, opt.lambda2=superclass@opt.lambda2, optcost = superclass@optcost, startdate = superclass@startdate, enddate = superclass@enddate, from = superclass@from,to = superclass@to,...))
  }
}

setGeneric(name = 'pppaenet', def = function(charlist,...){standardGeneric('pppaenet')})

setMethod('pppaenet',c(charlist = 'list'),
          function(charlist, BenchW, prices, CostData, gamma = 50, alpha = 0, startdate = 1, enddate = nrow(prices) - 1,trainDataRatio = c(0.8,0.2), lambda1 = NULL, minRatio = 1e-10, nLambda1 = 100, lambda2 = c(10^(-(0:7)),0), optcost = "NO_COST",tol = 1e-9, maxiter = 2000,ncores = 1, cfunfile = NULL,show = FALSE,...) {
            
            
            crossvalidate(charlist, BenchW, prices, CostData, gamma = gamma, alpha = alpha, startdate = startdate, enddate = enddate, nfolds = 1,lambda1 = lambda1, minRatio = minRatio, nLambda1 = nLambda1,lambda2 = lambda2,optcost = optcost,tol = tol, maxiter = maxiter,ncores = ncores, cfunfile = cfunfile,show = show,compute.rets = F,...)
          })

setMethod('pppaenet',c(charlist = 'pppData'),
          function(charlist,...) {
            
            pppaenet(charlist@signallist, charlist@BenchW, charlist@prices, charlist@CostData,...)
          })

setGeneric(name = 'crossvalidate', def = function(charlist,...){standardGeneric('crossvalidate')})

setMethod('crossvalidate',c(charlist = 'list'),
function(charlist, BenchW, prices, CostData, gamma = 50, alpha = 0, startdate = 1, enddate = nrow(rets) - 1,trainDataRatio = c(0.8,0.2), trainstartdate=NULL, trainenddate=NULL, teststartdate=NULL, testenddate=NULL, nfolds = 10, seed = 999, lambda1 = NULL, minRatio = 1e-10, nLambda1 = 100, lambda2 = c(10^(-(0:7)),0), optcost = "NO_COST", type = "OMG.out",tol = 1e-9, maxiter = 2000,ncores = 1, usecfun = FALSE, cfunfile = NULL, show = TRUE, compute.rets = TRUE,inference.type = "kfolds",name = "NULL",...) {
  
  if(usecfun){
    library(Rcpp)
    library(RcppArmadillo)
    
    sourceCpp(cfunfile)
  }
  
  
  if(is.numeric(CostData)&length(CostData)==1){
    CostData.tmp <- CostData
    CostData <- prices
    CostData[] <- NA
    CostData[] <- CostData.tmp
  }
  

  rets <- prices/mlag(prices) - 1
  
  data <- pppData(BenchW = BenchW, signallist = charlist, prices = prices, CostData = CostData)
  
  dates <- index(BenchW)
  if(is.character(startdate)) startdate <- which(dates ==  dates[dates>=startdate][1])
  if(is.character(enddate)) enddate <- which(dates ==  dates[dates<=enddate][length(dates[dates<=enddate])]) 
  
  np <- length(charlist)
  if(is.null(names(charlist))) charnames <- paste("Theta",1:np,sep="")
  else charnames <- names(charlist)
  
  
  
  if(!is.null(trainstartdate)&!is.null(trainenddate)&!is.null(teststartdate)&!is.null(testenddate)){
    
    if(is.character(trainstartdate)) trainstartdate <- which(dates ==  dates[dates>=trainstartdate][1])
    if(is.character(teststartdate)) teststartdate <- which(dates ==  dates[dates>=teststartdate][1])
    
    
    if(is.character(trainenddate)) trainenddate <- which(dates ==  dates[dates<=trainenddate][length(dates[dates<=trainenddate])]) 
    if(is.character(testenddate)) testenddate <- which(dates ==  dates[dates<=testenddate][length(dates[dates<=testenddate])]) 
    
    
  }else{
    
    dataIndex <- setupTrainData(startdate,enddate,trainDataRatio)
    trainstartdate <- dataIndex[1]
    trainenddate <- dataIndex[2]
    teststartdate <- dataIndex[3]
    testenddate <- dataIndex[4]
    
  }
  
  testid <- c(trainenddate, teststartdate, testenddate) - c(trainstartdate, trainenddate, teststartdate)
  if(any(testid <0 )) stop('train dates are not correctly set up')
  
  annualfactor <- compute.annual.factor(BenchW)
  
  allindex <- trainstartdate:trainenddate
  set.seed(seed)
  flds <- createFolds(allindex, k = nfolds, list = TRUE, returnTrain = FALSE)
  folds <- list()
  
  if(nfolds == 1){
    folds$train <- lapply(flds,function(x){allindex[x]})
    folds$test <- lapply(flds,function(x){allindex[-x]})
  }else{
    
    folds$test <- lapply(flds,function(x){allindex[x]})
    folds$train <- lapply(flds,function(x){allindex[-x]})
  }

  
  
  charlistTrain <- setupListRange(charlist,trainstartdate:trainenddate)
  prepEst <- getPrepEst(charlist, BenchW, rets)
  
  charRetMat <-  prepEst[[1]]
  benchRet <- prepEst[[2]]
  
  charChangelist2 <- NULL
  BenchWChange2 <- NULL
  
  charRetSquareMat <- charRetMat*charRetMat
  charRetcharChangeCostMat <- NULL
  charChangeSquareColMat <- NULL 
  
  if(!is.null(CostData)){
    #Nt <- rowSums(BenchW != 0) 
    Nt <- ncol(BenchW)
    BenchWChange <- BenchW - mlag(BenchW)
    #charChangelist <- lapply(1:length(charlist),function(i){coredata(charlist[[i]])/Nt - coredata(mlag(charlist[[i]]))/mlag(Nt)})
    charChangelist <- lapply(1:length(charlist),function(i){coredata(charlist[[i]])/Nt - coredata(mlag(charlist[[i]]))/(Nt)})
    CostData <- coredata(CostData)
    BenchWChange <- coredata(BenchWChange)
    
    charChangelist2 <- sapply(charChangelist, function(charChange){
      
      rowSums(CostData*abs(charChange))
    })
    
    BenchWChange2 <- rowSums(CostData*abs(BenchWChange))
    

    charRetcharChangeCostMat <- charRetMat*charChangelist2
    charChangeSquareMat <- charChangelist2*charChangelist2
    
    costCharChangelist <- lapply(charChangelist, function(charChange){
      
      CostData*charChange
    })
    
  }
  
  rets <- coredata(rets)
  BenchW <- coredata(BenchW)
  
  if(is.null(lambda1)){
    
    lambda1 <-getlambda1(charRetMat, benchRet, CostData, charChangelist, BenchWChange, costCharChangelist, trainstartdate, trainenddate, gamma, alpha, minRatio, nLambda1, optcost)
    
  }
  
  
  numtest <- length(folds$train)

  
  if(ncores > 1){
    
    cat("\n- Start multi cores estimation\n")
    cl <- makeMPIcluster(ncores)
    cfunfile <- cfunfile
    #print(cfunfile)
    clusterCall(cl,function(){
      library(xts)
      library(Rcpp)
      library(RcppArmadillo)
      library(lbfgs)
      sourceCpp(file = cfunfile)
    })
    
    clusterExport(cl, list("folds", "charlist", "rets", "BenchW","charRetMat","benchRet","CostData", "charChangelist", "BenchWChange","costCharChangelist", "charChangelist2","gamma","lambda1","lambda2","alpha","tol","maxiter","optcost","annualfactor","ncores","show","compute.rets","simpleGS2","setupListRange","computeWeight"),envir = environment())
    
    if(length(lambda2) == 1){
      cat("\nTimeseies cross validation on lambda1 grid, lambda2 is fixed (",lambda2,")\n")
      cat("NO. of loops/rolling windows to be estimated:",numtest,"\n")
      #outMat <- tsonefold(times, charRetMat, benchRet, gamma, lambda1, lambda2, invthetaW, alpha, tol, maxiter, type, annualfactor,ncores, cl=cl)
      out <- tsonefold2(folds,  charlist, rets, BenchW, charRetMat, benchRet, CostData, charChangelist, BenchWChange, costCharChangelist, charChangelist2, gamma, lambda1, lambda2, alpha, tol, maxiter, optcost, annualfactor,ncores,show = FALSE, cl=cl, compute.rets = compute.rets, usecfun=usecfun)
      
    }
    else{
      cat("\nTimeseies cross validation on lambda1 grid for each lambda2\n")
      cat("NO. of loops/rolling windows to be estimated:",numtest,"\n")
      if(show){
        op <- pboptions(type="txt")
      }else {
        op <- pboptions(type="none")
      }
      ans <- do.call(pblapply,list(1:length(lambda2),function(i){
      out <- tsonefold2(folds,  charlist, rets, BenchW, charRetMat, benchRet, CostData, charChangelist, BenchWChange, costCharChangelist, charChangelist2, gamma, lambda1, lambda2[i], alpha, tol, maxiter, optcost, annualfactor,ncores,show = FALSE, cl=cl, compute.rets = compute.rets, usecfun=usecfun)
        
        out
      }))
      outMat <- Reduce(rbind,(ans))
      #pboptions(op)
    }
    stopCluster(cl)
    cat("\n- End multi cores estimation\n")
    
    
    
  }else{
    
    if(length(lambda2) == 1){
      # cat("\nTimeseies simple cross validation on lambda1 grid, lambda2 is fixed (",lambda2,")\n")
      #cat("NO. of loops/rolling windows to be estimated:",numtest,"\n")
      outMat <- tsonefold2(folds,  charlist, rets, BenchW, charRetMat, benchRet, CostData, charChangelist, BenchWChange, costCharChangelist, charChangelist2, gamma, lambda1, lambda2, alpha, tol, maxiter, optcost, annualfactor,ncores,show = show, compute.rets = compute.rets, usecfun=usecfun)
      
    }else{
       cat("\nDouble cross validation on lambda1 grid for each lambda2\n")
      cat("NO. of loops/rolling windows to be estimated:",numtest,"\n")
      if(show){
        op <- pboptions(type="txt")
      }else {
        op <- pboptions(type="none")
      }
      

      ans <- do.call(pblapply,list(1:length(lambda2),function(i){
        out <- tsonefold2(folds,  charlist, rets, BenchW, charRetMat, benchRet, CostData, charChangelist, BenchWChange, costCharChangelist, charChangelist2, gamma, lambda1, lambda2[i], alpha, tol, maxiter, optcost, annualfactor,ncores,show = FALSE, compute.rets = compute.rets, usecfun=usecfun)
        
        out
      }))
      #pboptions(op)
      
      outMat <- Reduce(rbind,(ans))
    }
    
  }
  
  
  
  coefficients <- outMat[,1:np]
  
  if(length(lambda1) == 1) coefficients <- matrix(coefficients,nrow = 1)
  
  if(compute.rets){
    
    colnames(outMat) <-c(charnames, paste(charnames,".sd",sep=""), "Ave.Utility.out","SR.out","CEG.out","SOR.out","OMG.out","SZ3.out","SZ4.out")

    
    coefficients <- outMat[,1:np]
    minIndex <- which(outMat[,type] == min(outMat[,type],na.rm=TRUE))
    
    if(length(minIndex) > 1) minIndex <- minIndex[1]
    
    lambda1vec <- rep(lambda1,length(lambda2))
    lambda2vec <- rep(lambda2, each = length(lambda1))
    lambdaMat<-cbind(lambda1vec,lambda2vec)
    colnames(lambdaMat) <- c("lambda1","lambda2")
    
    opt.lambda1 <- as.numeric(lambdaMat[minIndex,1])
    opt.lambda2 <- as.numeric(lambdaMat[minIndex,2])
    
    
    criterion <- outMat[,c("Ave.Utility.out","SR.out","CEG.out","SOR.out","OMG.out")]
    criterion[,c("SR.out","CEG.out","SOR.out","OMG.out")] <- -criterion[,c("SR.out","CEG.out","SOR.out","OMG.out")]
    
    criterion.opt <- criterion[minIndex,]
    
    
    
    cv.results <- cbind(lambdaMat,coefficients,criterion)
    
    
    if(inference.type == "boot"){
        cat("Bootstrapping... \n")
        boot.results <- bootPPPAENET(charlist, rets, BenchW, charRetMat, benchRet, CostData, charChangelist, BenchWChange, costCharChangelist, charChangelist2, trainIndex=trainstartdate:trainenddate, gamma=gamma, alpha=alpha, lambda1=opt.lambda1, lambda2=opt.lambda2, optcost = optcost)
      
      theta.opt <- boot.results[,"mean"]
      theta.opt.sd <- boot.results[,"sd"]
    }else if(inference.type == "kfolds"){
      theta.opt <- outMat[minIndex,1:np]
      theta.opt.sd <- outMat[minIndex,(np + 1):(2*np)]
    }
    
    
    out <- AENetCV(data = data, coefs = coefficients, gamma = gamma, alpha = alpha, lambda1 = lambda1, lambda2 = lambda2 , optcost=optcost, opt.lambda1=opt.lambda1, opt.lambda2=opt.lambda2, startdate = startdate, enddate = enddate, from = as.character(index(rets)[startdate]),to = as.character(index(rets)[enddate]),name=name,
                   nfolds=nfolds, type=type, criterion.opt = criterion.opt,cv.results=cv.results,theta.opt=theta.opt,theta.opt.sd=theta.opt.sd,
                   trainstartdate=trainstartdate,trainenddate=trainenddate,teststartdate=teststartdate,testenddate=testenddate,
                   train.from=as.character(index(rets)[trainstartdate]),train.to=as.character(index(rets)[trainenddate]),
                   test.from=as.character(index(rets)[teststartdate]),test.to=as.character(index(rets)[testenddate]))
    
  }else{
    
    out <- AENet(data = data, coefs = coefficients, gamma = gamma, alpha = alpha, lambda1 = lambda1, lambda2 = lambda2 , optcost=optcost, opt.lambda1=0, opt.lambda2=0, startdate = startdate, enddate = enddate, from = as.character(index(prices)[startdate]),to = as.character(index(prices)[enddate]),name=name)
    
  }
  out
})


#setMethod('crossvalidate',c(charlist = 'pppData'),
#function(charlist, gamma = 50, alpha = 0, startdate = 1, enddate = nrow(rets) - 1,trainDataRatio = c(0.8,0.2), nfolds = 10, seed = 999, lambda1 = NULL, minRatio = 1e-10, nLambda1 = 100, lambda2 = c(10^(-(0:7)),0), optcost = TRUE, type = "CEG.out", max.portfolio.leverage = 400/100, tol = 1e-9, maxiter = 2000,ncores = 1, cfunfile = NULL, show = TRUE, compute.rets = TRUE,...) {

#  crossvalidate(charlist@signallist, charlist@BenchW, charlist@rets, CostData=charlist@CostData, gamma = gamma, alpha =alpha, startdate = startdate, enddate = enddate,trainDataRatio = trainDataRatio, nfolds = nfolds, seed = seed, lambda1 = lambda1, minRatio = minRatio, nLambda1 = nLambda1, lambda2 = lambda2, optcost = optcost, type = type, max.portfolio.leverage = max.portfolio.leverage, tol = tol, maxiter = maxiter,ncores = ncores, cfunfile = cfunfile, show = show, compute.rets = compute.rets,...)
  
#})

preparedata <- function(prices){
  data <- new.env()
  
  data$prices <- prices
  data$dates <- index(prices)
  
  naxts <- prices
  naxts[] <- NA
  
  data$weight <- naxts
  data$execution.price <- naxts
  
  data
}

setMethod('crossvalidate',c(charlist = 'pppData'),
function(charlist,...) {

  crossvalidate(charlist@signallist, charlist@BenchW, charlist@prices, charlist@CostData,...)

})

setGeneric(name = 'backtest', def = function(charlist, ...){standardGeneric('backtest')})

setMethod(f = 'backtest',
          signature = c(charlist = 'list'),
          definition = function(charlist, BenchW, prices, CostData, gamma = 50, alpha = 0, startdate = 1, enddate = nrow(rets) - 1,trainDataRatio =  c(0.8,0.2), trainstartdate=NULL, trainenddate=NULL, teststartdate=NULL, testenddate=NULL, lambda1 = 0.0001, lambda2 = 0.0001, optcost = TRUE, initialWindow = NULL, horizon = 1, skip = 0,fixedWindow = FALSE, oos = TRUE, max.portfolio.leverage = 400/100, tol = 1e-9, maxiter = 2000,usecfun = FALSE, cfunfile = NULL, show = FALSE, name = "NULL",run.bt = TRUE,...) {
            
            if(usecfun){
              library(Rcpp)
              library(RcppArmadillo)
              
              sourceCpp(cfunfile)
            }
            
            if(is.numeric(CostData)&length(CostData)==1){
              CostData.tmp <- CostData
              CostData <- prices
              CostData[] <- NA
              CostData[] <- CostData.tmp
            }
            
            rets <- prices/mlag(prices) - 1
            data <- pppData(BenchW = BenchW, signallist = charlist, prices = prices, CostData = CostData)
            
            dates <- index(BenchW)
            if(is.character(startdate)) startdate <- which(dates ==  dates[dates>=startdate][1])
            if(is.character(enddate)) enddate <- which(dates ==  dates[dates<=enddate][length(dates[dates<=enddate])]) 
            
            np <- length(charlist)
            if(is.null(names(charlist))) charnames <- paste("Theta",1:np,sep="")
            else charnames <- names(charlist)
            
            if(!is.null(trainstartdate)&!is.null(trainenddate)&!is.null(teststartdate)&!is.null(testenddate)){
              
              if(is.character(trainstartdate)) trainstartdate <- which(dates ==  dates[dates>=trainstartdate][1])
              if(is.character(teststartdate)) teststartdate <- which(dates ==  dates[dates>=teststartdate][1])
              
              
              if(is.character(trainenddate)) trainenddate <- which(dates ==  dates[dates<=trainenddate][length(dates[dates<=trainenddate])]) 
              if(is.character(testenddate)) testenddate <- which(dates ==  dates[dates<=testenddate][length(dates[dates<=testenddate])]) 
              
              
            }else{
              
              dataIndex <- setupTrainData(startdate,enddate,trainDataRatio)
              trainstartdate <- dataIndex[1]
              trainenddate <- dataIndex[2]
              teststartdate <- dataIndex[3]
              testenddate <- dataIndex[4]
              
            }
            
            testid <- c(trainenddate, teststartdate, testenddate) - c(trainstartdate, trainenddate, teststartdate)
            if(any(testid <0 )) stop('train dates are not correctly set up')
            
            
            
            
            
            annualfactor <- compute.annual.factor(BenchW)
            
            if(is.null(initialWindow)) initialWindow <- floor(sum(trainDataRatio)*(trainenddate - trainstartdate + 1)*trainDataRatio[1])
            
            
            trainstartdate1 <- trainstartdate
            trainenddate1 <- trainstartdate + initialWindow -1
            
            if(oos){
              
              times <- createTimeSlices(trainstartdate:testenddate, initialWindow = trainenddate - trainstartdate + 1, horizon = horizon,fixedWindow = fixedWindow, skip = skip)
              times$train <- lapply(times$train,function(x) x +trainstartdate-1)
              times$test <- lapply(times$test,function(x) x +trainstartdate-1)
              
              bt.type <- "OOS"
              
            }else{
              times <- createTimeSlices(trainstartdate:trainenddate, initialWindow = initialWindow, horizon = horizon,fixedWindow = fixedWindow, skip = skip)
              times$train <- lapply(times$train,function(x) x +trainstartdate-1)
              times$test <- lapply(times$test,function(x) x +trainstartdate-1)
              bt.type <- "IS"
            }
            
            
            #charlistTrain <- setupListRange(charlist,trainstartdate1:trainenddate1)
            
            prepEst <- getPrepEst(charlist, BenchW, rets)
            
            charRetMat <-  prepEst[[1]]
            benchRet <- prepEst[[2]]
            
            
            #charChangelist2 <- NULL
            #BenchWChange2 <- NULL
            
            #charRetSquareMat <- charRetMat*charRetMat
            #charRetcharChangeCostMat <- NULL
            #charChangeSquareColMat <- NULL 
            
            
              #Nt <- rowSums(BenchW != 0) 
              Nt <- ncol(BenchW)
              
              BenchWChange <- BenchW - mlag(BenchW)
              #charChangelist <- lapply(1:length(charlist),function(i){coredata(charlist[[i]])/Nt - coredata(mlag(charlist[[i]]))/mlag(Nt)})
              charChangelist <- lapply(1:length(charlist),function(i){coredata(charlist[[i]])/Nt - coredata(mlag(charlist[[i]]))/(Nt)})
              
              CostData <- coredata(CostData)
              BenchWChange <- coredata(BenchWChange)
              
              charChangelist2 <- sapply(charChangelist, function(charChange){
                
                rowSums(CostData*abs(charChange))
              })
              
              #BenchWChange2 <- rowSums(CostData*abs(BenchWChange))
              
              #charRetcharChangeCostMat <- charRetMat*charChangelist2
              #charChangeSquareMat <- charChangelist2*charChangelist2
              
              #charRetcharChangeCostMat <- charRetMat*charChangelist2
              #charChangeSquareMat <- charChangelist2*charChangelist2
              
              costCharChangelist <- lapply(charChangelist, function(charChange){
                
                CostData*charChange
              })
            
            
            
            numtest <- length(times$train)
            
            if(show){
              loopf <- pblapply
              op <- pboptions(type="txt")
            }else{
              loopf <- lapply
            }
            
            ans <- do.call(loopf,list(1:numtest,function(i){
              #trainstartdate1 <- times$train[[i]][1]
              #trainenddate1 <- times$train[[i]][length(times$train[[i]])]
              
              #teststartdate1 <- times$test[[i]][1]
              #testenddate1 <- times$test[[i]][length(times$test[[i]])]
              
              trainIndex <- times$train[[i]]
              testIndex <- times$test[[i]]
              
              invthetaW <- numeric(np) + 1
              
              thetatmp <- numeric(np)
              thetatmp0 <- numeric(np)
              
              
              
              charChangelistTrain <- setupListRange(charChangelist,trainIndex)
              charChangelistTest <- setupListRange(charChangelist,testIndex)
              
              
              #charRetSquareColMean <- colMeans(charRetSquareMat[trainstartdate1:trainenddate1,])
              #charRetcharChangeCostColMean <- colMeans(charRetcharChangeCostMat[trainstartdate1:trainenddate1,])
              #charChangeSquareColMean <- colMeans(charChangeSquareMat[trainstartdate1:trainenddate1,])
              
              
              charCovMat <- cov(charRetMat[trainIndex,],use="pairwise.complete.obs")
              covbenchVec <- sapply(1:np, function(i){
                return(cov(benchRet[trainIndex],charRetMat[trainIndex,i]))
              })
              
              covbenchVec <- matrix(covbenchVec, ncol=1)
              
              charMeanVec <- colMeans(charRetMat[trainIndex,])
              charMeanVec <- matrix(charMeanVec, ncol=1)
              
              charChangelist.train <- setupListRange(charChangelist,trainIndex)
              costCharChangelist.train<- setupListRange(costCharChangelist,trainIndex)
              
              charChangeCostMeanVec <- colMeans(charChangelist2[trainIndex,])
              
              if(usecfun){
                est <- pppest(thetatmp0, thetatmp, charCovMat, covbenchVec, charMeanVec, charChangeCostMeanVec, CostData[trainIndex,], charChangelist.train, costCharChangelist.train, BenchWChange[trainIndex,], trainIndex, gamma, lambda1, lambda2, alpha, invthetaW, optcost, 1e-8, 2000)
                
              }else{
                est <- pppestRversion(thetatmp0, thetatmp, charCovMat, covbenchVec, charMeanVec, charChangeCostMeanVec, CostData[trainIndex,], charChangelist.train, costCharChangelist.train, BenchWChange[trainIndex,], trainIndex, gamma, lambda1, lambda2, alpha, invthetaW, optcost, 1e-8, 2000)
                
              }
              
              thetatmp0 <- est[[1]]
              thetatmp <- est[[2]]
              
              charlistTest <- setupListRange(charlist,testIndex )
              
              weight <- computeWeight(thetatmp, charlistTest, BenchW[testIndex,])
              
              #weight <- iif(rowSums(abs(weight)) > max.portfolio.leverage, max.portfolio.leverage,rowSums(abs(weight)))*weight/ rowSums(abs(weight))
              
              weight <- weight/ rowSums(abs(weight))
              
              comsignals <- computeCombineSignals(thetatmp, charlistTest)
              
              out <- list(thetatmp,weight,comsignals)
              out
            }))
            
            
            outMat <- sapply(ans,function(x) (x[[1]]))
            outMat <- t(outMat)
            
            wout <- lapply(ans,function(x) (x[[2]]))
            wout <- do.call(rbind, wout)
            
            sout <- lapply(ans,function(x) (x[[3]]))
            sout <- do.call(rbind, sout)
            
            colnames(outMat) <-c(charnames)
            coefficients <- outMat[,1:np]
            
            
            trainIndex <- unique(as.numeric(unlist(times$train)))
            testIndex <- unique(as.numeric(unlist(times$test)))
           
            trainstartdate <- trainIndex[1]
            trainenddate <- trainIndex[length(trainIndex)]
            
              teststartdate <- testIndex[1]
              testenddate <- testIndex[length(testIndex)]
              
            timeindex <- index(prices)[testIndex]
            
            bench.weight <- BenchW[testIndex,]
            

            theta.mean <- colMeans(coefficients)
            theta.sd <- apply(coefficients, 2, sd)
            
            weight <- xts(wout, order.by = timeindex)
            comsignals <- xts(sout, order.by = timeindex)
            

              commission <- list(cps = 0.00, fixed =0, percentage =  mean(CostData))
              
              data.bt <- preparedata(prices[timeindex])
              data.bt$weight[] = NA
              data.bt$weight = weight
              ppp = bt.run.share(data.bt,commission = commission,silent=T)
              
              ppp$orgweight <- weight
              ppp$comsignals <- comsignals
              
              ppp.summary <- bt.detail.summary(ppp)
              
                
                data.bt$weight[] = NA
                data.bt$weight = (bench.weight)
                benchmark = bt.run.share(data.bt,commission = commission,silent=T)
                
                benchmark$orgweight <- bench.weight
                benchmark.summary <- bt.detail.summary(benchmark)
  
           

            out <- AENetBT(data = data, coefs = coefficients, gamma = gamma, alpha = alpha, lambda1 = lambda1, lambda2 = lambda2 , optcost=optcost, opt.lambda1=lambda1, opt.lambda2=lambda2, startdate = startdate, enddate = enddate, from = as.character(index(rets)[startdate]),to = as.character(index(rets)[enddate]),name = name,
                           initialWindow = initialWindow,fixedWindow = fixedWindow,horizon = horizon,skip = skip,oos = oos,theta.mean = theta.mean,theta.sd=theta.sd,trainstartdate=trainstartdate,trainenddate=trainenddate,teststartdate=teststartdate,testenddate=testenddate,ppp = ppp,ppp.summary = ppp.summary,benchmark = benchmark,benchmark.summary = benchmark.summary,train.from=as.character(index(rets)[trainstartdate]),train.to=as.character(index(rets)[trainenddate]),
                           test.from=as.character(index(rets)[teststartdate]),test.to=as.character(index(rets)[testenddate]))
            
            
            return(out)
          })






setMethod(f = 'backtest',
          signature = c(charlist = 'pppData'),
          definition = function(charlist, ...) {
            
            backtest(charlist@signallist, charlist@BenchW, charlist@prices, charlist@CostData, ...)
            
          })

setMethod(f = 'backtest',
          signature = c(charlist = 'AENet'),
          definition = function(charlist,...) {
            
            backtest(charlist@data@signallist,charlist@data@BenchW, charlist@data@prices, charlist@data@CostData, gamma = charlist@gamma, alpha =charlist@alpha, startdate = charlist@startdate, enddate = charlist@enddate,lambda1 = charlist@opt.lambda1, lambda2 = charlist@opt.lambda2, optcost = charlist@optcost, name=charlist@name,...)
            
          })


setMethod(f = 'backtest',
          signature = c(charlist = 'AENetCV'),
          definition = function(charlist,...) {
            
            backtest(charlist@data@signallist,charlist@data@BenchW, charlist@data@prices, charlist@data@CostData, gamma = charlist@gamma, alpha =charlist@alpha, startdate = charlist@startdate, enddate = charlist@enddate,lambda1 = charlist@opt.lambda1, lambda2 = charlist@opt.lambda2, optcost = charlist@optcost, name=charlist@name,trainstartdate=charlist@trainstartdate, trainenddate=charlist@trainenddate, teststartdate=charlist@teststartdate, testenddate=charlist@testenddate,...)
            
          })


setMethod('plot',c(x = 'AENetBT'),
          function(x,y=NULL, title = 'Peformance comparison',...) {
            models<-list()
            models$benchmark <- x@benchmark
            models$ppp <- x@ppp
            strategy.performance.snapshoot(models, T, title = title)  
          })

setMethod('plot',c(x = 'AENetCV'),
          function(x,y=NULL, xlabtxt =  expression(paste("log(",lambda,"1)",sep="")), ylabtxt = x@type, tilttxt = "Estimated model criterion",...) {
            lambda1vec <- rev(log(x@lambda1))
            criterion <- rev(as.numeric(x@cv.results[,x@type]))
            
            if(length(x@lambda2) == 1){
              exampleData <- data.frame(lambda1=lambda1vec, criterion = criterion)
              colnames(exampleData) <- c("Lambda1","Criterion")
              
              print(ggplot(data = exampleData, aes(x = Lambda1, y = Criterion)) + geom_line(colour="red", size=1) + ylab(ylabtxt) + xlab(xlabtxt) + theme_grey(base_size = 12) + ggtitle(tilttxt))
              
            }else{
              
              lambda2vec <- rev(log(x@lambda2))
              
              lambda1vec <- rep(lambda1vec,length(x@lambda2))
              lambda2vec <- rep(x@lambda2, each = length(x@lambda1))
              
              
              
              exampleData <- data.frame(lambda1=lambda1vec, lambda2=lambda2vec, criterion = criterion)
              colnames(exampleData) <- c("Lambda1","Lambda2","Criterion")
              
              ylabtxt <- expression(paste("log(",lambda,"2",sep=""))
              
              print(ggplot(data = exampleData, aes(x = Lambda1, y = Lambda2, z = Criterion))+stat_contour(aes(x = Lambda1, y = Lambda2, z = Criterion,colour=..level..) )+ ylab(ylabtxt) + xlab(xlabtxt) + theme_grey(base_size = 12) + ggtitle(tilttxt))
              
            }
          })



setMethod('print',c(x = 'AENet'),
          function(x,...) {
            
            cat("=========================================\n")
            cat("PPP AENET: path \n")
            cat("=========================================\n")
            
            cat("startdate: ", as.character(x@from), "\n")
            cat("enddate: ", as.character(x@to), "\n")
            cat("gamma: ", x@gamma, "\n")
            cat("alpha: ", x@alpha, "\n")
           
            cat("optcost: ", x@optcost, "\n")
            cat("lambda1 grids: ", paste("[",format(min(x@lambda1),digits = 3),",",format(max(x@lambda1),digits = 3),"]\n",sep=""))
            cat("lambda2 grids: ", paste("[",format(min(x@lambda2),digits = 3),",",format(max(x@lambda2),digits = 3),"]\n",sep=""))
            #cat("Coefficients: ", x@coefs, "\n")
            cat("=========================================\n")
            
          })


setMethod('print',c(x = 'AENetCV'),
          function(x,...) {
            
            cat("=========================================\n")
            cat(paste("PPP AENET: Cross Validation",x@name,"\n"))
            cat("=========================================\n")
            
            cat("startdate: ", as.character(x@train.from), "\n")
            cat("enddate: ", as.character(x@train.to), "\n")
            cat("gamma: ", x@gamma, "\n")
            cat("alpha: ", x@alpha, "\n")
            cat("type: ", x@type, "\n")
            cat("optcost: ", x@optcost, "\n")
            cat("lambda1 grids: ", paste("[",format(min(x@lambda1),digits = 3),",",format(max(x@lambda1),digits = 3),"]\n",sep=""))
            cat("lambda2 grids: ", paste("[",format(min(x@lambda2),digits = 3),",",format(max(x@lambda2),digits = 3),"]\n",sep=""))
            cat("=========================================\n")
            cat("resluts:\n")
            cat("=========================================\n")
            
            cat("Optimal lambda1: ",x@opt.lambda1,"\n")
            cat("Optimal lambda2: ",x@opt.lambda2,"\n")
            cat("=========================================\n")
            
          })

setMethod('print',c(x = 'AENetBT'),
          function(x,...) {
            
            cat("=========================================\n")
            cat(paste("PPP AENET:Backtest",x@name,"\n"))
            cat("=========================================\n")
            cat("startdate: ", as.character(x@test.from), "\n")
            cat("enddate: ", as.character(x@test.to), "\n")
            cat("gamma: ", x@gamma, "\n")
            cat("alpha: ", x@alpha, "\n")
            cat("optcost: ", x@optcost, "\n")
            cat("lambda1: ",x@opt.lambda1,"\n")
            cat("lambda2: ",x@opt.lambda2,"\n")
            cat("=========================================\n")
            
          })


bt.detail.summary <- function
(
  bt,
  trade.summary = NULL
)
{
  out.all = list()
  out = list()
  out$Period = join( format( range(index.xts(bt$equity)), '%b%Y'), ' - ')
  out$Cagr = compute.cagr(bt$equity)
  out$Sharpe = compute.sharpe(bt$ret) / 100
  out$DVR = compute.DVR(bt) / 100
  out$Volatility = compute.risk(bt$ret)
  out$MaxDD = compute.max.drawdown(bt$equity)
  out$AvgDD = compute.avg.drawdown(bt$equity)
  if( !is.null(trade.summary) ) {
    out$Profit.Factor = trade.summary$stats['profitfactor', 'All']
  }
  out$VaR = compute.var(bt$ret)
  out$CVaR = compute.cvar(bt$ret)
  out$Exposure = compute.exposure(bt$weight)
  out$AveNetLeverage = mean(rowSums(bt$weight)[-1])
  out$AveGrossLeverage = mean(rowSums(abs(bt$weight))[-1])
  out$MaxNetLeverage = max(rowSums(bt$weight)[-1])
  out$MaxGrossLeverage = max(rowSums(abs(bt$weight))[-1])
  out.all$System = lapply(out, function(x) if(is.double(x)) round(100*x,2) else x)
  if( !is.null(bt$trade.summary) ) trade.summary = bt$trade.summary
  out = list()
  if( !is.null(trade.summary) ) {
    out$Win.Percent = trade.summary$stats['win.prob', 'All']
    out$Avg.Trade = trade.summary$stats['avg.pnl', 'All']
    out$Avg.Win = trade.summary$stats['win.avg.pnl', 'All']
    out$Avg.Loss = trade.summary$stats['loss.avg.pnl', 'All']
    out = lapply(out, function(x) if(is.double(x)) round(100*x,1) else x)
    out$Best.Trade = max(as.double(trade.summary$trades[, 'return']))
    out$Worst.Trade = min(as.double(trade.summary$trades[, 'return']))
    out$WinLoss.Ratio = round( -trade.summary$stats['win.avg.pnl', 'All']/trade.summary$stats['loss.avg.pnl', 'All'] , 2)
    out$Avg.Len = round(trade.summary$stats['len', 'All'],2)
    out$Num.Trades = trade.summary$stats['ntrades', 'All']
  }
  out.all$Trade = out
  out = list()
  out$Win.Percent.Day = sum(bt$ret > 0, na.rm = T) / len(bt$ret)
  out$Best.Day = bt$best
  out$Worst.Day = bt$worst
  month.ends = endpoints(bt$equity, 'months')
  mret = ROC(bt$equity[month.ends,], type = 'discrete')
  out$Win.Percent.Month = sum(mret > 0, na.rm = T) / len(mret)
  out$Best.Month = max(mret, na.rm = T)
  out$Worst.Month = min(mret, na.rm = T)
  year.ends = endpoints(bt$equity, 'years')
  mret = ROC(bt$equity[year.ends,], type = 'discrete')
  out$Win.Percent.Year = sum(mret > 0, na.rm = T) / len(mret)
  out$Best.Year = max(mret, na.rm = T)
  out$Worst.Year = min(mret, na.rm = T)
  out.all$Period = lapply(out, function(x) if(is.double(x)) round(100*x,1) else x)
  return(out.all)
}

bt.detail.summary2 <- function
(
  bt,
  trade.summary = NULL
)
{
  out.all = list()
  out = list()
  out$Period = join( format( range(index.xts(bt$equity)), '%b%Y'), ' - ')
  out$AveReturns = compute.mean(bt$ret)
  
  out$Skewness = skewness(bt$ret)
  out$Volatility = compute.risk(bt$ret)
  out$Sharpe = compute.sharpe(bt$ret) / 100
  out$Cagr = compute.cagr(bt$equity)
  out$MaxDD = compute.max.drawdown(bt$equity)
  out$AveDD = compute.avg.drawdown(bt$equity)
  
  month.ends = endpoints(bt$equity, 'months')
  mret = ROC(bt$equity[month.ends,], type = 'discrete')
  
  out$Win.Percent.Month = sum(mret > 0, na.rm = T) / len(mret)
  out$Best.Month = max(mret, na.rm = T)
  out$Worst.Month = min(mret, na.rm = T)
  out$AveNetLeverage = mean(rowSums(bt$weight)[-1])
  out$AveGrossLeverage = mean(rowSums(abs(bt$weight))[-1])
  out$MaxNetLeverage = max(rowSums(bt$weight)[-1])
  out$MaxGrossLeverage = max(rowSums(abs(bt$weight))[-1])
  #out$AveAnnualTurnover = compute.turnover(bt,b = NULL)
  
  
  #out$VaR = compute.var(bt$ret)
  #out$CVaR = compute.cvar(bt$ret)
  #out$Exposure = compute.exposure(bt$weight)
  out.all$System = lapply(out, function(x) if(is.double(x)) round(100*x,2) else x)

  return(out.all)
}
plotbtlist <- function(x, title = 'Peformance comparison', showBenchmark = FALSE){
  
  models<-list()
  
  for(i in 1:length(x)){
    models[[i]]<- x[[i]]@ppp
  }
  
  names(models)<- names(x)
  
  if(showBenchmark){
    
    benchmark<- list(benchmark=x[[1]]@benchmark)
    models <- c(benchmark,models)
  }

  strategy.performance.snapshoot(models,one.page = T,title = title) 
  #models
}

setMethod('coef',c(object = 'AENetCV'),
          function(object,...) {
            coefnames <- names(object@data@signallist)
            idx <- which(nchar(coefnames)==0)
            if(length(idx) > 0){
              addnames <- paste("theta.",idx,sep="")
              coefnames[idx] <- addnames
            }
            
            out<- object@theta.opt
            names(out) <- coefnames
            return(out)
          })

setMethod('summary',c(object = 'AENetCV'),
          function(object,...) {
            print(object)
            cat("Estimated coefficients:\n")
            print(coef(object))
            cat("Standard errors:\n")
            print(coef.se(object))
            cat("=========================================\n")
          })


setGeneric(name = 'coef.se', def = function(object, ...){standardGeneric('coef.se')})

setMethod('coef.se',c(object = 'AENetCV'),
          function(object,...) {
            coefnames <- names(object@data@signallist)
            idx <- which(nchar(coefnames)==0)
            if(length(idx) > 0){
              addnames <- paste("theta.",idx,sep="")
              coefnames[idx] <- addnames
            }
            coefnames <- paste(coefnames,".se",sep="")
            out<- object@theta.opt.sd
            names(out) <- coefnames
            return(out)
          })


setGeneric(name = 'tuning.para', def = function(object, ...){standardGeneric('tuning.para')})

setMethod('tuning.para',c(object = 'AENet'),
          function(object,...) {
            
            out<- c(object@opt.lambda1,object@opt.lambda2)
            names(out) <- c("lamdba1","lambda2")
            return(out)
          })

setMethod('coef',c(object = 'AENetBT'),
          function(object,...) {
            coefnames <- names(object@data@signallist)
            idx <- which(nchar(coefnames)==0)
            if(length(idx) > 0){
              addnames <- paste("theta.",idx,sep="")
              coefnames[idx] <- addnames
            }
            
            out<- object@theta.mean
            names(out) <- coefnames
            return(out)
          })

setMethod('coef.se',c(object = 'AENetBT'),
          function(object,...) {
            coefnames <- names(object@data@signallist)
            idx <- which(nchar(coefnames)==0)
            if(length(idx) > 0){
              addnames <- paste("theta.",idx,sep="")
              coefnames[idx] <- addnames
            }
            coefnames <- paste(coefnames,".se",sep="")
            out<- object@theta.sd
            names(out) <- coefnames
            return(out)
          })


setGeneric(name = 'cv.results', def = function(object, ...){standardGeneric('cv.results')})
setMethod('cv.results',c(object = 'AENetCV'),
          function(object,...) {
            
            return(object@cv.results)
          })
setMethod('summary',c(object = 'AENetBT'),
          function(object,...) {
            print(object)
            cat("Estimated coefficients:\n")
            print(coef(object))
            cat("Standard errors:\n")
            print(coef.se(object))
            cat("=========================================\n")
            temp1 <- list2matrix(object@ppp.summary )
            temp2 <- list2matrix(object@benchmark.summary)
            temp <- cbind(temp2,temp1)
            temp = rbind(c("",'Benchmark',"",'PPP'), temp)   # add header row
            cat("Performance comparison:","\n")
            print.table(temp)
            cat("=========================================\n")
          })


setGeneric(name = 'plot.coef', def = function(object, ...){standardGeneric('plot.coef')})
setMethod('plot.coef',c(object = 'AENetBT'),
          function(object,tilttxt = "Coefficient estimates",coef.index = NULL,...) {
            
            if(is.null(coef.index)){
              exampleData <- data.frame(index=1:nrow(object@coefs), object@coefs)
            }else{
              exampleData <- data.frame(index=1:nrow(object@coefs[,coef.index]), object@coefs[,coef.index])
            }
            
            
            
            sectorMelt <- melt(exampleData, id = "index")
            colnames(sectorMelt) <- c("Time","Signal","Value")
            print(ggplot(data = sectorMelt, aes(x = Time, y = Value,group = Signal, colour = Signal))+geom_line(size=1) + theme_grey(base_size = 12) + ggtitle(tilttxt))
            
          })




setMethod('plot.coef',c(object = 'AENet'),
          function(object,xlabtxt =  expression(paste("log(",lambda,")",sep="")), ylabtxt = "Estimated Coefficients", tilttxt = expression(paste("Estimated ",theta," of PPP AENET",sep="")),...) {
            
            thetaMat <- apply(object@coefs, 2, rev)
            colnames(thetaMat) <- names(object@data@signallist)
            #rownames(thetaMat) <- NULL
            lambda <- rev(log(object@lambda1))
            exampleData <- data.frame(lambda1=lambda, thetaMat)
            sectorMelt <- melt(exampleData, id = "lambda1")
            colnames(sectorMelt) <- c("Lambda1","Signals","Value")
            print(ggplot(data = sectorMelt, aes(x = Lambda1, y = Value,group = Signals, colour = Signals))+geom_line(size=1)+ylab(ylabtxt)+xlab(xlabtxt) + theme_grey(base_size = 12) + ggtitle(tilttxt))
            
          })


bt.summary.3 <- function(bt){
  
  
  out = list()
  out$Period = join( format( range(index.xts(bt$equity)), '%b%Y'), ' - ')
  #print('---')
  out$AveReturns = compute.mean(bt$ret)
  
  out$Skewness = skewness(bt$ret)
  out$Volatility = compute.risk(bt$ret)
  out$Sharpe = compute.sharpe(bt$ret)/100
  out$MaxDD = compute.max.drawdown(bt$equity)
  out$AveDD = compute.avg.drawdown(bt$equity)
  
  
  month.ends = endpoints(bt$equity, 'months')
  mret = ROC(bt$equity[month.ends,], type = 'discrete')
  
  out$Win.Percent.Month = sum(mret > 0, na.rm = T) / len(mret)
  out$Best.Month = max(mret, na.rm = T)
  out$Worst.Month = min(mret, na.rm = T)
  
  #out <- lapply(out, function(x) if(is.double(x)) round(x,3) else x)
  #out <- lapply(out, function(x) if(is.double(x)) format(round(x, 3), nsmall = 2) else x)
  out <- lapply(out, function(x) if(is.double(x)) round(100*x,2) else x)
  out
}


get.bt.benchmarkIndex <- function(benchmark.ret, testIndex){
  
  benchmark.ret <- benchmark.ret[testIndex,]
  
  benchmark.ret[1] <- 0
  benchmark.equity <- cumprod(1 + benchmark.ret)
  
  bt.benchmark <- list()
  bt.benchmark$ret <- benchmark.ret
  bt.benchmark$equity <- benchmark.equity
  
  return(bt.benchmark)
}

plotbt2 <- function(bt,tilttxt = 'Peformance comparison'){
  
  yrange<- range(bt[[1]]$equity)
  legendnames <- names(bt)[1]
  
  temptable <- list2matrix(bt.summary.3(bt[[1]]))
  
  for(i in 2:length(bt)){
    
    yrange <- range(yrange,bt[[i]]$equity)
    legendnames <- c(legendnames,names(bt)[i])
    temptable <- cbind(temptable,list2matrix(bt.summary.3(bt[[i]])))
  }
  temptable = rbind(legendnames, temptable)
  
  layout(1:2)
  
  plota(bt[[1]]$equity,type = 'l',col = 1,ylim = yrange, main = tilttxt)
  
  for(i in 2:length(bt)){
    
    plota.lines(bt[[i]]$equity,col = i)
  }
  plota.legend(legendnames, 1:(length(bt) + 1))
  
  plot.table(temptable)
  
}

plotbt.custom.report.part1.2 <- function
(
  ...,
  dates = NULL,
  main = '',
  trade.summary = FALSE,
  x.highlight = NULL
)
{
  layout(1:3)
  models = variable.number.arguments( ... )
  model = models[[1]]
  plotbt(models, dates = dates, main = main, plotX = F, log = 'y', LeftMargin = 3, x.highlight = x.highlight)
  mtext('Cumulative Performance', side = 2, line = 1)
  plotbt(models, plottype = '12M', dates = dates, plotX = F, LeftMargin = 3, x.highlight = x.highlight)
  mtext('12 Month Rolling', side = 2, line = 1)
  plotbt(models, dates = dates, xfun = function(x) { 100 * compute.drawdown(x$equity) }, LeftMargin = 3, x.highlight = x.highlight)
  mtext('Drawdown', side = 2, line = 1)
}
####################################################################################
## Set of R functions for Sharpe screening
####################################################################################

.sharpeScreening = function(X, control = list()) {
  
  # process control
  ctr = processControl(control)
  
  # size of inputs and outputs
  T = nrow(X) 
  N = ncol(X)
  pval = dsharpe = tstat = matrix(data = NA, N, N)
  
  # determine which pairs can be compared (in a matrix way)
  Y = 1 * (!is.nan(X) & !is.na(X))
  YY = crossprod(Y) #YY = t(Y) %*% Y # row i indicates how many observations in common with column k
  YY[YY < ctr$minObs] = 0
  YY[YY > 0] = 1 
  liststocks = c(1:nrow(YY))[rowSums(YY) > ctr$minObsPi]  
  
  # determine bootstrap indices (do it before to speed up computations)
  bsids = bootIndices(T, ctr$nBoot, ctr$bBoot)
  
  if (length(liststocks) > 1){
    cl = makeCluster(c(rep("localhost", ctr$nCore)), type = "SOCK")
    
    liststocks = liststocks[1:(length(liststocks)-1)]
    
    z <- clusterApply(cl     = cl, 
                      x      = as.list(liststocks), 
                      fun    = sharpeScreeningi, 
                      rdata  = X,
                      T      = T, 
                      N      = N, 
                      nBoot  = ctr$nBoot, 
                      bsids  = bsids, 
                      minObs = ctr$minObs,
                      type   = ctr$type,
                      hac    = ctr$hac,
                      b      = ctr$bBoot,
                      ttype  = ctr$ttype,
                      pBoot  = ctr$pBoot)
    stopCluster(cl)
    
    for (i in 1:length(liststocks)) {
      out = z[[i]]  
      id  = liststocks[i]
      pval[id,id:N]    = pval[id:N,id] = out[[2]][id:N]
      dsharpe[id,id:N] = out[[1]][id:N]
      dsharpe[id:N,id] = -out[[1]][id:N]
      tstat[id,id:N]   = out[[3]][id:N]
      tstat[id:N,id]   = -out[[3]][id:N]
    }
  }
  
  # pi
  pi = computePi(pval = pval, dalpha = dsharpe, tstat = tstat, lambda = ctr$lambda, nBoot = ctr$nBoot)
  
  # info on the funds  
  info = infoFund(X)
  
  # form output
  out = list(n       = info$nObs, 
             npeer   = colSums(!is.na(pval)),
             sharpe  = info$sharpe, 
             dsharpe = dsharpe, 
             pval    = pval,
             tstat   = tstat,
             lambda  = pi$lambda,
             pizero  = pi$pizero,
             pipos   = pi$pipos,
             pineg   = pi$pineg)
  
  return(out)
}
sharpeScreening = compiler::cmpfun(.sharpeScreening)

## Sharpe ratio screening for fund i again its peers
.sharpeScreeningi = function(i, rdata, T, N, nBoot, bsids, minObs, type, hac, b, ttype, pBoot) {
  
  nPeer = N - i
  X = matrix(rdata[,i], nrow = T, ncol = nPeer)
  Y = matrix(rdata[, (i+1):N], nrow = T, ncol = nPeer)
  
  dXY = X - Y
  idx = (!is.nan(dXY)&!is.na(dXY))
  X[!idx] = NA
  Y[!idx] = NA    
  nObs = colSums(idx)
  
  pvali = dsharpei = tstati = rep(NA, N)
  
  k = 0
  for (j in (i + 1) : N) {
    k = k + 1
    if (nObs[k] < minObs) {
      next
    }
    rets = cbind(X[idx[,k],1], Y[idx[,k],k])
    
    if (type == 1) {
      tmp = sharpeTestAsymptotic(rets, hac, ttype)
    }
    else {
      tmp = sharpeTestBootstrap(rets, bsids, b, ttype, pBoot)
    }
    
    dsharpei[j] = tmp$dsharpe
    pvali[j]    = tmp$pval
    tstati[j]   = tmp$tstat 
  }
  
  out = list(dsharpei = dsharpei, pvali = pvali, tstati = tstati)
  return(out)   
}
sharpeScreeningi = compiler::cmpfun(.sharpeScreeningi)
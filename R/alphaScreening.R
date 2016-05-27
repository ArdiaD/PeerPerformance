####################################################################################
## Set of R function for alpha screening
####################################################################################

.alphaScreening = function(X, factors = NULL, control = list()) {  
  
  # process control
  ctr = processControl(control)
  
  T = nrow(X) 
  N = ncol(X)
  pval = dalpha = tstat = matrix(data = NA, N, N)
  
  # determine which pairs can be compared (in a matrix way)
  Y = 1 * (!is.nan(X) & !is.na(X))
  YY = crossprod(Y) #YY = t(Y) %*% Y # row i indicates how many observations in common with column k
  YY[YY < ctr$minObs] = 0
  YY[YY > 0] = 1 
  liststocks = c(1:nrow(YY))[rowSums(YY) > ctr$minObsPi]  
  
  if (length(liststocks) > 1){
    cl = makeCluster(c(rep("localhost", ctr$nCore)), type = "SOCK")
    
    if (ctr$hac){
      clusterEvalQ(cl, require("sandwich"))
      clusterEvalQ(cl, require("lmtest"))
    }
   
    liststocks = liststocks[1:(length(liststocks)-1)]
    
    z <- clusterApply(cl      = cl, 
                      x       = as.list(liststocks), 
                      fun     = alphaScreeningi, 
                      rdata   = X, 
                      factors = factors, 
                      T       = T, 
                      N       = N,
                      hac     = ctr$hac)
    
    stopCluster(cl)
    
    for (i in 1:length(liststocks)){
      out = z[[i]]  
      id  = liststocks[i]
      pval[id,id:N]   = pval[id:N,id] = out[[2]][id:N]
      dalpha[id,id:N] = out[[1]][id:N]
      dalpha[id:N,id] = -out[[1]][id:N]
      tstat[id,id:N]  =  out[[3]][id:N]
      tstat[id:N,id]  = -out[[3]][id:N]
    }
  }
  
  # pi
  pi = computePi(pval = pval, dalpha = dalpha, tstat = tstat, lambda = ctr$lambda, nBoot = ctr$nBoot)
  
  # info on the funds  
  info = infoFund(X, factors = factors)
  
  # form output
  out = list(n       = info$nObs, 
             npeer   = colSums(!is.na(pval)),
             alpha   = info$alpha, 
             dalpha  = dalpha, 
             pval    = pval,
             tstat   = tstat,
             lambda  = pi$lambda,
             pizero  = pi$pizero,
             pipos   = pi$pipos,
             pineg   = pi$pineg)
  
  return(out)
}
alphaScreening = compiler::cmpfun(.alphaScreening)

.alphaScreeningi = function(i, rdata, factors, T, N, hac){
  pvali = dalphai = tstati = rep(NA, N)
  
  nPeer = N - i
  X = matrix(rdata[,i], nrow = T, ncol = nPeer)
  Y = matrix(rdata[,(i+1):N], nrow = T, ncol = nPeer)
  dXY = X - Y 
  
  if (!hac){
    if (is.null(factors)){
      fit = lm(dXY ~ 1, na.action = na.omit)
    }
    else{
      fit = lm(dXY ~ 1 + factors, na.action = na.omit) 
    }
    sumfit = summary(fit)
    if (nPeer == 1){
      pvali[N]   = sumfit$coef[1,4]
      dalphai[N] = sumfit$coef[1,1]
      tstati[N]  = sumfit$coef[1,3]
    } 
    else{
      k = 1
      for (j in (i + 1) : N){
        pvali[j]   = sumfit[[k]]$coef[1,4]
        dalphai[j] = sumfit[[k]]$coef[1,1]
        tstati[j]  = sumfit[[k]]$coef[1,3]
        k = k + 1
      }
    }
  }
  else{
    if (nPeer == 1){
      if (is.null(factors)){
        fit = lm(dXY ~ 1, na.action = na.omit) 
      }
      else{
        fit = lm(dXY ~ 1 + factors, na.action = na.omit) 
      }
      sumfit = lmtest::coeftest(fit, vcov. = sandwich::vcovHAC(fit))
      pvali[N]   = sumfit[1,4]
      dalphai[N] = sumfit[1,1]
      tstati[N]  = sumfit[1,3]
    }
    else{
      k = 1
      for (j in (i + 1) : N){
        if (is.null(factors)){
          fit = lm(dXY[,k] ~ 1, na.action = na.omit) 
        }
        else{
          fit = lm(dXY[,k] ~ 1 + factors, na.action = na.omit) 
        }
        sumfit = coeftest(fit, vcov. = vcovHAC(fit))
        pvali[j]   = sumfit[1,4]
        dalphai[j] = sumfit[1,1]
        tstati[j]  = sumfit[1,3]
        k = k + 1
      }
    }
  }
  
  out = list(dalphai = dalphai, pvali = pvali, tstati = tstati)
  return(out)
} 
alphaScreeningi = compiler::cmpfun(.alphaScreeningi)
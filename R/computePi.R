## Set of R functions used to compute pi and lambda

# #' @name .computePi
# #' @title Compute pi0, pi+ and pi-
# #' @importFrom stats qnorm
# #' @import compiler
.computePi <- function(pval, dalpha, tstat, lambda = 0.5, nBoot = 500, 
                       bpos = 0.4, bneg = 0.6, adjust = TRUE) {
  if (!is.matrix(pval)) {
    pval <- matrix(pval, nrow = 1)
  }
  if (!is.matrix(dalpha)) {
    dalpha <- matrix(dalpha, nrow = 1)
  }
  if (!is.matrix(tstat)) {
    tstat <- matrix(tstat, nrow = 1)
  }
  
  m <- nrow(pval)
  # n1 = ncol(pval) # n + 1 funds
  if (length(lambda) == 1) {
    lambda <- rep(lambda, m)
  }
  
  pizero <- pipos <- pineg <- lambda_ <- rep(NA, m)
  for (i in 1:m) {
    pvali <- pval[i, ]
    dalphai <- dalpha[i, ]
    tstati <- tstat[i, ]
    if (all(is.na(pvali))) {
      next
    }
    if (is.null(lambda)) {
      lambdai <- computeOptLambda(pval = pvali, nBoot = nBoot, adjust = adjust)
    } else {
      lambdai <- lambda[i]
    }
    
    pizeroi <- computePizero(pvali, lambda = lambdai, adjust = adjust)
    idxOK <- !is.na(pvali) & !is.na(dalphai)
    n <- sum(idxOK)  # number of peers
    if (n <= 1) {
      next
    }
    
    ni0 <- pizeroi * n
    hn <- round(0.5 * n)
    piposi <- pinegi <- 0
    
    # idxOKnPizeroi = idxOK & pvali >= lambdai we need to fix this
    qpos <- stats::qnorm(p = bpos)
    qneg <- stats::qnorm(p = bneg)
    
    if (sum(dalphai[idxOK] >= 0) >= hn) {
      piposi <- (1/n) * min((n - ni0), max(sum(tstati[idxOK] >= qpos) - 
                                             ni0 * (1 - bpos), 0))
      pinegi <- min(max(1 - pizeroi - piposi, 0), 1)  # numerical stability
    } else {
      pinegi <- (1/n) * min((n - ni0), max(sum(tstati[idxOK] <= qneg) - 
                                             ni0 * bneg, 0))
      piposi <- min(max(1 - pizeroi - pinegi, 0), 1)  # numerical stability
    }
    
    pizero[i] <- pizeroi
    pipos[i] <- piposi
    pineg[i] <- pinegi
    lambda_[i] <- lambdai
  }
  
  out <- list(pizero = pizero, pipos = pipos, pineg = pineg, lambda = lambda_)
  return(out)
}
computePi <- compiler::cmpfun(.computePi)

# #' @name .computePizero
# #' @import compiler
.computePizero <- function(pval, lambda = 0.5, adjust = TRUE) {
  if (!is.matrix(pval)) {
    pval <- matrix(pval, nrow = 1)
  }
  
  n <- ncol(pval)
  pizero <- mean(pval >= lambda, na.rm = TRUE)
  pizero <- pizero * (1/(1 - lambda))
  pizero[pizero > 1] <- 1
  # adjust pi using truncated binomial
  if (adjust) {
    pizero <- adjustPi(pizero, n = n, lambda = lambda)
  }
  
  return(pizero)
}
computePizero <- compiler::cmpfun(.computePizero)

# #' @name .adjustPi
# #' @title Adjust estimated pi0 using quadratif fit
# #' @importFrom stats dnorm pnorm uniroot
# #' @import compiler
.adjustPi <- function(pi.hat, n = 100, lambda = 0.5) {
  
  asym.hatpi0 <- function(pi0) {
    npi0 <- pi0 * n
    nlambda <- npi0 * (1 - lambda)
    if (nlambda >= n) {
      hatpi0 <- pi0
    } else {
      s2 <- nlambda * (n - nlambda)/(n^3 * (1 - lambda)^2)
      s <- sqrt(s2)
      zcrit <- (1 - pi0)/s
      hatpi0 <- pi0 + s * (-stats::dnorm(zcrit) + (1 - stats::pnorm(zcrit)) * 
                             zcrit)
    }
    return(hatpi0)
  }
  
  asym.inverse <- function(f, lower = -100, upper = 100) {
    FUN <- function(y) stats::uniroot((function(x) f(x) - y), lower = lower, 
                                      upper = upper)$root
    return(FUN)
  }
  
  asym.invpi0 <- asym.inverse(asym.hatpi0, lower = 1e-05, upper = 1.5)
  
  m <- length(pi.hat)
  out <- vector("double", m)
  for (i in 1:m) {
    tmp <- pi.hat[i]
    test.it <- try({
      tmp <- asym.invpi0(pi.hat[i])
    }, silent = TRUE)
    if (class(test.it) == "try-error") {
      tmp <- pi.hat[i]
    }
    tmp[tmp > 1] <- 1
    tmp[tmp < 0] <- 0
    out[i] <- tmp
  }
  return(out)
}
adjustPi <- compiler::cmpfun(.adjustPi)

# #' @name .computeOptLambda
# #' @title Compute optimal lamba values
# #' @importFrom stats runif
# #' @import compiler
.computeOptLambda <- function(pval, nBoot = 500, adjust = TRUE) {
  if (!is.matrix(pval)) {
    pval <- matrix(pval, nrow = 1)
  }
  
  n <- nrow(pval)  # number of funds
  
  vlambda <- seq(0.3, 0.7, 0.02)
  nvlambda <- length(vlambda)
  
  mpizero <- matrix(data = NA, nrow = n, ncol = nvlambda)
  for (i in 1:nvlambda) {
    mpizero[, i] <- computePizero(pval, lambda = vlambda[i], adjust = adjust)
  }
  # pi0hat
  vminpizero <- apply(mpizero, 1, "min")
  
  idx <- !is.nan(pval) & !is.na(pval)
  nObs <- rowSums(idx)
  
  bsunif <- matrix(stats::runif(max(nObs) * nBoot), max(nObs), nBoot)
  optlambda <- rep(0.5, n)
  
  for (i in 1:n) {
    nObsi <- nObs[i]
    if (nObsi > 0) {
      bsidx <- ceiling(bsunif[1:nObsi, 1:nBoot] * nObsi)
      pvalb <- stats::na.omit(pval[i, ])
      pvalb <- matrix(pvalb[bsidx], nrow = nBoot)
      
      mpizerob <- matrix(data = NA, nrow = nBoot, ncol = nvlambda)
      for (j in 1:nvlambda) {
        mpizerob[, j] <- computePizero(pvalb, lambda = vlambda[j], 
                                       adjust = adjust)
      }
      
      vMSE <- colSums((mpizerob - vminpizero[i])^2)
      
      optlambda[i] <- vlambda[which.min(vMSE)]
    }
  }
  
  return(optlambda)
}
computeOptLambda <- compiler::cmpfun(.computeOptLambda)
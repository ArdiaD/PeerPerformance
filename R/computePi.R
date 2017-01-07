## Set of R functions used to compute pi and lambda

#' @name .computePi
#' @title Compute pi0, pi+ and pi-
#' @importFrom stats qnorm
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
    
    # if (sum(dalphai[idxOK] >= 0) >= hn){ piposi = (1 / n) * min((n -
    # ni0), max(sum(dalphai[idxOK] >= qpos) - ni0 * (1 - bpos), 0)) pinegi
    # = min(max(1 - pizeroi - piposi, 0), 1) # numerical stability } else{
    # pinegi = (1 / n) * min((n - ni0), max(sum(dalphai[idxOK] <= qneg) -
    # ni0 * bneg, 0)) piposi = min(max(1 - pizeroi - pinegi, 0), 1) #
    # numerical stability }
    
    pizero[i] <- pizeroi
    pipos[i] <- piposi
    pineg[i] <- pinegi
    lambda_[i] <- lambdai
  }
  
  out <- list(pizero = pizero, pipos = pipos, pineg = pineg, lambda = lambda_)
  return(out)
}
computePi <- compiler::cmpfun(.computePi)

#@name .computePizero
#@description Given a value of lambda, this function returns the pizero of the funds
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

#' @name .adjustPi
#' @title Adjust estimated pi0 using quadratif fit
#' @importFrom stats dnorm pnorm uniroot
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

#' @name .computeOptLambda
#' @title Compute optimal lamba values
#' @importFrom stats runif
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
      pvalb <- na.omit(pval[i, ])
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

# .coefAdjustPi = function(n, lambda){ ids = matrix(data =
# c(10,0.3,10,0.4,10,0.5,10,0.6,10,0.7,50,0.3,50,0.4,50,0.5,50,0.6,50,0.7,100,0.3,100,
# 0.4,100,0.5,100,0.6,100,0.7,500,0.3,500,0.4,500,0.5,500,0.6,500,0.7,1000,0.3,
# 1000,0.4,1000,0.5,1000,0.6,1000,0.7), nrow = 25, ncol = 2, byrow =
# TRUE) coef = matrix(data =
# c(12.4264397055725,-28.8462924603265,18.0061745439696,2.22784959599922,-4.33794978760025,
# 3.33846783477142,-3.5062288124517,9.94411690220581,-5.50860579817955,-9.59063438888672,
# 25.9944774753484,-16.0243234516758,-9.50704513388566,26.8509583669393,-17.2326763276237,
# 2.69973475372684,-5.18728038928751,3.56020419106916,2.90584925091467,-5.75369084260271,
# 3.94462528611671,3.10498649240963,-6.32430459839148,4.34712083046976,3.10341843532753,
# -6.459248111501,4.5207017397825,3.3201246167253,-7.13256919772734,5.03730846454779,
# 1.92615117037718,-3.36988186562355,2.48393491184823,2.23024137854889,-4.11257781958698,
# 2.93730040964901,2.53451778338116,-4.8778475515017,3.41778177692686,2.8712983587036,
# -5.72958870731188,3.95833690647981,2.65060599304652,-5.37075396882709,3.85060868153408,
# 0.782769992193963,-0.743978244550511,0.971372543293702,1.06310054658503,-1.37802062179133,
# 1.32960965871391,1.34679398014843,-2.02732476869674,1.7009379767174,1.79031480628225,
# -3.03536712430862,2.27395681782583,2.14237527425761,-3.87093079288463,2.76951593594632,
# 0.430884461964093,0.0411279494749582,0.533315651199973,0.659038682749762,-0.467637536214552,
# 0.816502791444037,0.902316149838157,-1.01323439030704,1.12219969436304,1.18316263219461,
# -1.64749879388121,1.48017617185269,1.52276073392502,-2.4266177990621,1.92694731845667),
# nrow = 25, ncol = 3, byrow = TRUE) absn = abs(n - ids[,1]) abslambda
# = abs(lambda - ids[,2]) idx = (absn == min(absn)) & abslambda ==
# min(abslambda) coef = coef[which(idx),] return(coef) } coefAdjustPi =
# cmpfun(.coefAdjustPi) .createArrayCoefAdjustPi = function(npi0 = 100,
# M = 1000){ n = c(10, 50, 100, 500, 1000) lambda = seq(0.3, 0.7, 0.1)
# n1 = length(n) n2 = length(lambda) ids = matrix(data = NA, nrow = n1
# * n2, ncol = 2) coef = matrix(data = NA, nrow = n1 * n2, ncol = 3) k
# = 1 for (i in 1 : n1){ for (j in 1 : n2){ tmp = computeCoefAdjustPi(n
# = n[i], lambda = lambda[j], npi0 = npi0, M = M) ids[k,] = c(n[i],
# lambda[j]) coef[k,] = tmp$coef k = k + 1 } } out = list(ids = ids,
# coef = coef) return(out) } createArrayCoefAdjustPi =
# cmpfun(.createArrayCoefAdjustPi) .computeCoefAdjustPi = function(n =
# 100, lambda = 0.5, npi0 = 100, M = 1000){ # true pi0 mesh lb = 0.85
# ub = 0.999 pi0 = seq(from = lb, to = ub, length.out = npi0) npi0 =
# length(pi0) pi0.unbiased = pi0.biased = vector('double', npi0) for (i
# in 1 : npi0){ n0 = floor(pi0[i] * n) # Monte Carlo estimation of pi0
# pval = matrix(data = 0, nrow = M, ncol = n) pval[1:M,1:n0] = runif(n
# = M * n0, min = 0, max = 1) tmp = rowSums(pval >= lambda) / (1 -
# lambda) pi0.unbiased[i] = mean(tmp) / n tmp[tmp > n] = n
# pi0.biased[i] = mean(tmp) / n } fit = lm(pi0 ~ 1 + pi0.biased +
# I(pi0.biased^2)) coef = as.numeric(fit$coef) out = list(coef = coef,
# pi0 = pi0, pi0.unbiased = pi0.unbiased, pi0.biased = pi0.biased, n =
# n, lambda = lambda, npi0 = npi0, M = M) return(out) }
# computeCoefAdjustPi = cmpfun(.computeCoefAdjustPi)
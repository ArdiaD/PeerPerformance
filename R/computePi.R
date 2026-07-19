## Set of R functions used to compute pi and lambda

# #' @name .computePi
# #' @title Compute pi0, pi+ and pi-
# #' @importFrom stats qnorm
# #' @import compiler
.computePi <- function(pval, dalpha, tstat, lambda = 0.5, nBoot = 499,
                       bpos = 0.4, bneg = 0.6, adjust = TRUE, fast = FALSE) {

  if (!is.matrix(pval) & !is.array(pval)) {
    pval <- matrix(pval, nrow = 1)
  }
  if (!is.matrix(dalpha) & !is.array(dalpha)) {
    dalpha <- matrix(dalpha, nrow = 1)
  }
  if (!is.matrix(tstat) & !is.array(tstat)) {
    tstat <- matrix(tstat, nrow = 1)
  }

  if(length(dim(pval))==2){
    m <- nrow(pval)
    factors_dim <- 1
    dim(pval) <- c(1, dim(pval))
    dim(dalpha) <- c(1, dim(dalpha))
    dim(tstat) <- c(1, dim(tstat))
  }else{
    m <- dim(pval)[2]
    factors_dim <- dim(pval)[1]
    }

  # n1 = ncol(pval) # n + 1 funds
  if (!is.null(lambda) && !(length(lambda) %in% c(1L, m))) {
    stop("'lambda' must be NULL, length 1, or length equal to the number of funds")
  }
  if (length(lambda) == 1) {
    lambda <- rep(lambda, m)
  }

  pizero <- pipos <- pineg <- lambda_ <- matrix(rep(NA, m*factors_dim), nrow=factors_dim)

  for(factor_dim in 1:factors_dim){
    for (i in 1:m) {
      pvali <- pval[factor_dim, i, ]
      dalphai <- dalpha[factor_dim, i, ]
      tstati <- tstat[factor_dim, i, ]
      if (all(is.na(pvali))) {
        next
      }
      if (is.null(lambda)) {
        lambdai <- computeOptLambda(pval = pvali, nBoot = nBoot,
                                    adjust = adjust, fast = fast)
      } else {
        lambdai <- lambda[i]
      }

      pizeroi <- computePizero(pvali, lambda = lambdai, adjust = adjust,
                               fast = fast)
      idxOK <- !is.na(pvali) & !is.na(dalphai) & !is.na(tstati)
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

      pizero[factor_dim, i] <- pizeroi
      pipos[factor_dim, i] <- piposi
      pineg[factor_dim, i] <- pinegi
      lambda_[factor_dim, i] <- lambdai
    }
  }

  if(factors_dim == 1){
  	pizero <- pizero[1, ]
  	pipos <- pipos[1, ]
  	pineg <- pineg[1, ]
  	lambda_ <- lambda_[1, ]
  }
  out <- list(pizero = pizero, pipos = pipos, pineg = pineg, lambda = lambda_)
  return(out)
}
computePi <- compiler::cmpfun(.computePi)

# #' @name .computePizero
# #' @import compiler
.computePizero <- function(pval, lambda = 0.5, adjust = TRUE, fast = FALSE) {
  if (!is.matrix(pval)) {
    pval <- matrix(pval, nrow = 1)
  }

  # Note: 'n' (number of trials in the truncated-binomial adjustment of adjustPi)
  # is the number of columns of 'pval'. This matches the original Ardia & Boudt
  # (2018) implementation; when some peers are NA it slightly over-counts the
  # effective number of peers. Kept as-is for consistency with published results.
  n <- ncol(pval)
  pizero <- apply(pval>=lambda, 1, mean, na.rm=TRUE)
  # pizero <- mean(pval >= lambda, na.rm = TRUE)
  pizero <- pizero * (1/(1 - lambda))
  pizero[pizero > 1] <- 1
  # adjust pi using truncated binomial
  if (adjust) {
    pizero <- adjustPi(pizero, n = n, lambda = lambda, fast = fast)
  }

  return(pizero)
}
computePizero <- compiler::cmpfun(.computePizero)

# #' @name .adjustPi
# #' @title Adjust estimated pi0 using quadratic fit
# #' @importFrom stats dnorm pnorm uniroot
# #' @import compiler
.adjustPi <- function(pi.hat, n = 100, lambda = 0.5, fast = FALSE) {

  # Fast path (control$fastAdjust): 'n' and 'lambda' are constant within a
  # call, so the same monotone map is inverted at every element. Instead of one
  # uniroot() per element, invert the whole vector at once by bisection. Forty
  # halvings of [1e-5, 1.5] give an absolute accuracy of about 1e-12, i.e. far
  # tighter than uniroot's default tolerance (.Machine$double.eps^0.25, about
  # 1.2e-4), so the fast path is if anything more accurate than the default one.
  if (isTRUE(fast)) {
    fwd <- function(pi0) {
      nlambda <- pi0 * n * (1 - lambda)
      out <- pi0
      i <- nlambda < n                       # else branch of asym.hatpi0
      if (any(i)) {
        s <- sqrt(nlambda[i] * (n - nlambda[i])/(n^3 * (1 - lambda)^2))
        zcrit <- (1 - pi0[i])/s
        out[i] <- pi0[i] + s * (-stats::dnorm(zcrit) +
                                  (1 - stats::pnorm(zcrit)) * zcrit)
      }
      out
    }
    m  <- length(pi.hat)
    lo <- rep(1e-05, m)
    hi <- rep(1.5, m)
    # outside the bracket uniroot() would fail; keep the input, as the
    # try-error fallback of the default path does
    inside <- pi.hat >= fwd(lo) & pi.hat <= fwd(hi)
    for (it in seq_len(40L)) {
      mid  <- (lo + hi)/2
      left <- fwd(mid) < pi.hat
      lo[left]  <- mid[left]
      hi[!left] <- mid[!left]
    }
    out <- (lo + hi)/2
    out[!inside] <- pi.hat[!inside]
    out[out > 1] <- 1
    out[out < 0] <- 0
    return(out)
  }

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
    if (inherits(test.it, "try-error")) {
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
.computeOptLambda <- function(pval, nBoot = 499, adjust = TRUE, fast = FALSE) {
  if (!is.matrix(pval)) {
    pval <- matrix(pval, nrow = 1)
  }

  n <- nrow(pval)  # number of funds

  # Grid for the data-driven lambda. Ardia & Boudt (2018) describe a finer grid
  # {0.30, 0.32, ..., 0.70}; the coarser 0.1 step is kept here for consistency
  # with the published package results (and is much faster).
  vlambda <- seq(0.3, 0.7, 0.1)
  nvlambda <- length(vlambda)

  mpizero <- matrix(data = NA, nrow = n, ncol = nvlambda)
  for (i in 1:nvlambda) {
    mpizero[, i] <- computePizero(pval, lambda = vlambda[i], adjust = adjust,
                                  fast = fast)
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
                                       adjust = adjust, fast = fast)
      }
      vMSE <- colSums((mpizerob - vminpizero[i])^2)

      optlambda[i] <- vlambda[which.min(vMSE)]
    }
  }

  return(optlambda)
}
computeOptLambda <- compiler::cmpfun(.computeOptLambda)

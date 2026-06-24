## Cross-group screening: each fund in group X against every fund in group Y.
## The single-fund-against-a-group case is simply nX = 1 (X a vector/one column).
## ----------------------------------------------------------------------------
## alpha
## ----------------------------------------------------------------------------

# #' @importFrom stats lm na.omit
# #' @importFrom lmtest coeftest
# #' @importFrom sandwich vcovHAC
.alphaScreeningXYi <- function(i, X, Y, factors, T, nY, hac, minObs,
                               screen_beta = FALSE) {
  if (screen_beta && !is.null(factors)) {
    ncoef <- 1L + ncol(factors)
    row_return <- 1:ncoef
  } else {
    ncoef <- 1L
    row_return <- 1L
  }
  pvali <- dalphai <- tstati <- matrix(NA, nrow = ncoef, ncol = nY)
  xi <- X[, i]
  for (j in 1:nY) {
    dxy <- xi - Y[, j]
    ok <- !is.na(dxy) & !is.nan(dxy)
    if (sum(ok) < minObs) {
      next
    }
    if (all(abs(dxy[ok]) < sqrt(.Machine$double.eps))) {
      next  # identical series (e.g. a fund compared with itself)
    }
    if (is.null(factors)) {
      fit <- stats::lm(dxy ~ 1, na.action = stats::na.omit)
    } else {
      fit <- stats::lm(dxy ~ 1 + factors, na.action = stats::na.omit)
    }
    if (!hac) {
      sm <- summary(fit)$coef
    } else {
      sm <- lmtest::coeftest(fit, vcov. = sandwich::vcovHAC(fit))
    }
    pvali[row_return, j]   <- sm[row_return, 4]
    dalphai[row_return, j] <- sm[row_return, 1]
    tstati[row_return, j]  <- sm[row_return, 3]
  }
  list(dalphai = dalphai, pvali = pvali, tstati = tstati)
}
alphaScreeningXYi <- compiler::cmpfun(.alphaScreeningXYi)

.alphaScreeningXY <- function(X, Y, factors = NULL, control = list(),
                              screen_beta = FALSE) {
  ctr <- processControl(control)
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  T <- nrow(X)
  if (nrow(Y) != T) {
    stop("'X' and 'Y' must have the same number of rows (time periods)")
  }
  nX <- ncol(X)
  nY <- ncol(Y)

  if (screen_beta && !is.null(factors)) {
    ncoef <- 1L + ncol(factors)
  } else {
    ncoef <- 1L
  }
  row_return <- 1:ncoef
  pval <- dalpha <- tstat <- array(NA, dim = c(ncoef, nX, nY))

  cl <- parallel::makeCluster(ctr$nCore)
  z <- parallel::clusterApplyLB(cl = cl, x = as.list(1:nX),
                                fun = alphaScreeningXYi, X = X, Y = Y,
                                factors = factors, T = T, nY = nY, hac = ctr$hac,
                                minObs = ctr$minObs, screen_beta = screen_beta)
  parallel::stopCluster(cl)

  for (i in 1:nX) {
    out <- z[[i]]
    dalpha[row_return, i, ] <- out$dalphai[row_return, ]
    pval[row_return, i, ]   <- out$pvali[row_return, ]
    tstat[row_return, i, ]  <- out$tstati[row_return, ]
  }

  pi <- computePi(pval = pval, dalpha = dalpha, tstat = tstat,
                  lambda = ctr$lambda, nBoot = ctr$nBoot,
                  bpos = ctr$gammaPos, bneg = ctr$gammaNeg)

  info <- infoFund(X, factors = factors, screen_beta = screen_beta)

  if (!screen_beta) {
    pval   <- matrix(pval[1, , ],   nrow = nX, ncol = nY)
    dalpha <- matrix(dalpha[1, , ], nrow = nX, ncol = nY)
    tstat  <- matrix(tstat[1, , ],  nrow = nX, ncol = nY)
    npeer  <- rowSums(!is.na(pval))
  } else {
    npeer <- apply(!is.na(pval), c(1, 2), sum)
    cn <- .coefNames(factors)
    rownames(pi$pizero) <- rownames(pi$pipos) <- rownames(pi$pineg) <- cn
    rownames(pi$lambda) <- cn
    rownames(info$alpha) <- cn
    rownames(npeer) <- cn
  }

  out <- list(n = info$nObs, npeer = npeer, ny = nY, alpha = info$alpha,
              dalpha = dalpha, pval = pval, tstat = tstat, lambda = pi$lambda,
              pizero = pi$pizero, pipos = pi$pipos, pineg = pi$pineg,
              cross = TRUE)
  class(out) <- "SCREENING"
  return(out)
}
alphaScreeningXY <- compiler::cmpfun(.alphaScreeningXY)

## ----------------------------------------------------------------------------
## Sharpe
## ----------------------------------------------------------------------------

.sharpeScreeningXYi <- function(i, X, Y, T, nY, nBoot, bsids, minObs, type,
                                hac, b, ttype, pBoot) {
  pvali <- dsharpei <- tstati <- rep(NA, nY)
  xi <- X[, i]
  for (j in 1:nY) {
    dxy <- xi - Y[, j]
    ok <- !is.na(dxy) & !is.nan(dxy)
    if (sum(ok) < minObs) {
      next
    }
    if (all(abs(dxy[ok]) < sqrt(.Machine$double.eps))) {
      next
    }
    rets <- cbind(xi[ok], Y[ok, j])
    if (type == 1) {
      tmp <- sharpeTestAsymptotic(rets, hac, ttype)
    } else {
      tmp <- sharpeTestBootstrap(rets, bsids, b, ttype, pBoot)
    }
    dsharpei[j] <- tmp$dsharpe
    pvali[j]    <- tmp$pval
    tstati[j]   <- tmp$tstat
  }
  list(dsharpei = dsharpei, pvali = pvali, tstati = tstati)
}
sharpeScreeningXYi <- compiler::cmpfun(.sharpeScreeningXYi)

.sharpeScreeningXY <- function(X, Y, control = list()) {
  ctr <- processControl(control)
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  T <- nrow(X)
  if (nrow(Y) != T) {
    stop("'X' and 'Y' must have the same number of rows (time periods)")
  }
  nX <- ncol(X)
  nY <- ncol(Y)

  bsids <- bootIndices(T, ctr$nBoot, ctr$bBoot)
  pval <- dsharpe <- tstat <- matrix(NA, nX, nY)

  cl <- parallel::makeCluster(ctr$nCore)
  z <- parallel::clusterApplyLB(cl = cl, x = as.list(1:nX),
                                fun = sharpeScreeningXYi, X = X, Y = Y, T = T,
                                nY = nY, nBoot = ctr$nBoot, bsids = bsids,
                                minObs = ctr$minObs, type = ctr$type,
                                hac = ctr$hac, b = ctr$bBoot, ttype = ctr$ttype,
                                pBoot = ctr$pBoot)
  parallel::stopCluster(cl)

  for (i in 1:nX) {
    out <- z[[i]]
    dsharpe[i, ] <- out$dsharpei
    pval[i, ]    <- out$pvali
    tstat[i, ]   <- out$tstati
  }

  pi <- computePi(pval = pval, dalpha = dsharpe, tstat = tstat,
                  lambda = ctr$lambda, nBoot = ctr$nBoot,
                  bpos = ctr$gammaPos, bneg = ctr$gammaNeg)
  info <- infoFund(X)

  out <- list(n = info$nObs, npeer = rowSums(!is.na(pval)), ny = nY,
              sharpe = info$sharpe, dsharpe = dsharpe, pval = pval,
              tstat = tstat, lambda = pi$lambda, pizero = pi$pizero,
              pipos = pi$pipos, pineg = pi$pineg, cross = TRUE)
  class(out) <- "SCREENING"
  return(out)
}
sharpeScreeningXY <- compiler::cmpfun(.sharpeScreeningXY)

## ----------------------------------------------------------------------------
## modified Sharpe
## ----------------------------------------------------------------------------

.msharpeScreeningXYi <- function(i, X, Y, level, T, nY, nBoot, bsids, minObs,
                                 na.neg, type, hac, b, ttype, pBoot) {
  pvali <- dmsharpei <- tstati <- rep(NA, nY)
  xi <- X[, i]
  for (j in 1:nY) {
    dxy <- xi - Y[, j]
    ok <- !is.na(dxy) & !is.nan(dxy)
    if (sum(ok) < minObs) {
      next
    }
    if (all(abs(dxy[ok]) < sqrt(.Machine$double.eps))) {
      next
    }
    rets <- cbind(xi[ok], Y[ok, j])
    if (type == 1) {
      tmp <- msharpeTestAsymptotic(rets, level, na.neg, hac, ttype)
    } else {
      tmp <- msharpeTestBootstrap(rets, level, na.neg, bsids, b, ttype, pBoot)
    }
    dmsharpei[j] <- tmp$dmsharpe
    pvali[j]     <- tmp$pval
    tstati[j]    <- tmp$tstat
  }
  list(dmsharpei = dmsharpei, pvali = pvali, tstati = tstati)
}
msharpeScreeningXYi <- compiler::cmpfun(.msharpeScreeningXYi)

.msharpeScreeningXY <- function(X, Y, level = 0.9, na.neg = TRUE, control = list()) {
  ctr <- processControl(control)
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  T <- nrow(X)
  if (nrow(Y) != T) {
    stop("'X' and 'Y' must have the same number of rows (time periods)")
  }
  nX <- ncol(X)
  nY <- ncol(Y)

  bsids <- bootIndices(T, ctr$nBoot, ctr$bBoot)
  pval <- dmsharpe <- tstat <- matrix(NA, nX, nY)

  cl <- parallel::makeCluster(ctr$nCore)
  z <- parallel::clusterApplyLB(cl = cl, x = as.list(1:nX),
                                fun = msharpeScreeningXYi, X = X, Y = Y,
                                level = level, T = T, nY = nY, nBoot = ctr$nBoot,
                                bsids = bsids, minObs = ctr$minObs,
                                na.neg = na.neg, type = ctr$type, hac = ctr$hac,
                                b = ctr$bBoot, ttype = ctr$ttype, pBoot = ctr$pBoot)
  parallel::stopCluster(cl)

  for (i in 1:nX) {
    out <- z[[i]]
    dmsharpe[i, ] <- out$dmsharpei
    pval[i, ]     <- out$pvali
    tstat[i, ]    <- out$tstati
  }

  pi <- computePi(pval = pval, dalpha = dmsharpe, tstat = tstat,
                  lambda = ctr$lambda, nBoot = ctr$nBoot,
                  bpos = ctr$gammaPos, bneg = ctr$gammaNeg)
  info <- infoFund(X, level = level, na.neg = na.neg)

  out <- list(n = info$nObs, npeer = rowSums(!is.na(pval)), ny = nY,
              msharpe = info$msharpe, dmsharpe = dmsharpe, pval = pval,
              tstat = tstat, lambda = pi$lambda, pizero = pi$pizero,
              pipos = pi$pipos, pineg = pi$pineg, cross = TRUE)
  class(out) <- "SCREENING"
  return(out)
}
msharpeScreeningXY <- compiler::cmpfun(.msharpeScreeningXY)

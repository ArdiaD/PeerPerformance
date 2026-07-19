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
    # usable rows: account for missing factors as lm(na.omit) would
    if (is.null(factors)) {
      ok <- !is.na(dxy) & !is.nan(dxy)
    } else {
      ok <- stats::complete.cases(dxy, factors)
    }
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
    sfit <- summary(fit)
    if (!is.finite(sfit$sigma) || sfit$sigma < sqrt(.Machine$double.eps)) {
      next  # (near) deterministic differential: skip to avoid spurious significance
    }
    if (!hac) {
      sm <- sfit$coef
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
  if (nX < 1L || nY < 1L) {
    stop("cross-group screening needs at least one fund in 'X' and one peer in 'Y'")
  }

  if (screen_beta && !is.null(factors)) {
    ncoef <- 1L + ncol(factors)
  } else {
    ncoef <- 1L
  }
  row_return <- 1:ncoef
  pval <- dalpha <- tstat <- array(NA, dim = c(ncoef, nX, nY))

  if (ctr$nCore == 1) {
    # serial path: no PSOCK cluster (avoids the per-call cluster overhead);
    # the anonymous wrapper avoids the clash with lapply's own 'X' argument
    z <- lapply(1:nX, function(ii)
      alphaScreeningXYi(ii, X = X, Y = Y, factors = factors, T = T, nY = nY,
                        hac = ctr$hac, minObs = ctr$minObs,
                        screen_beta = screen_beta))
  } else {
    cl <- parallel::makeCluster(ctr$nCore)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    z <- parallel::clusterApplyLB(cl = cl, x = as.list(1:nX),
                                  fun = alphaScreeningXYi, X = X, Y = Y,
                                  factors = factors, T = T, nY = nY, hac = ctr$hac,
                                  minObs = ctr$minObs, screen_beta = screen_beta)
  }

  for (i in 1:nX) {
    out <- z[[i]]
    dalpha[row_return, i, ] <- out$dalphai[row_return, ]
    pval[row_return, i, ]   <- out$pvali[row_return, ]
    tstat[row_return, i, ]  <- out$tstati[row_return, ]
  }

  pi <- computePi(pval = pval, dalpha = dalpha, tstat = tstat,
                  lambda = ctr$lambda, nBoot = ctr$nBoot,
                  bpos = ctr$gammaPos, bneg = ctr$gammaNeg,
                  fast = ctr$fastAdjust)

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
      # indices pre-generated in the master for this pair's length
      bs <- bsids[[as.character(sum(ok))]]
      if (is.null(bs)) {
        next  # block length exceeds this pair's sample size
      }
      tmp <- sharpeTestBootstrap(rets, bs, b, ttype, pBoot)
    }
    if (!is.finite(tmp$tstat)) {
      next  # degenerate pair (e.g. zero-variance series)
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
  if (ctr$bBoot == 0) {
    stop("'bBoot = 0' (data-driven block length) is not supported in screening; ",
         "set an explicit block length 'bBoot >= 1'.")
  }
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  T <- nrow(X)
  if (nrow(Y) != T) {
    stop("'X' and 'Y' must have the same number of rows (time periods)")
  }
  nX <- ncol(X)
  nY <- ncol(Y)
  if (nX < 1L || nY < 1L) {
    stop("cross-group screening needs at least one fund in 'X' and one peer in 'Y'")
  }

  # pre-generate the bootstrap indices in the master, one matrix per distinct
  # cross-pair complete-case length (workers select by length)
  # only needed for the bootstrap test, and only for pairs passing 'minObs',
  # so that the asymptotic path never depends on 'bBoot'
  pairLens <- crossprod(1 * (!is.na(X) & !is.nan(X)), 1 * (!is.na(Y) & !is.nan(Y)))
  bsids <- NULL
  if (ctr$type == 2) {
    bsids <- bootIndicesByLen(pairLens[pairLens >= ctr$minObs],
                              ctr$nBoot, ctr$bBoot)
  }
  pval <- dsharpe <- tstat <- matrix(NA, nX, nY)

  if (ctr$nCore == 1) {
    # serial path: no PSOCK cluster (avoids the per-call cluster overhead);
    # the anonymous wrapper avoids the clash with lapply's own 'X' argument
    z <- lapply(1:nX, function(ii)
      sharpeScreeningXYi(ii, X = X, Y = Y, T = T, nY = nY, nBoot = ctr$nBoot,
                         bsids = bsids, minObs = ctr$minObs, type = ctr$type,
                         hac = ctr$hac, b = ctr$bBoot, ttype = ctr$ttype,
                         pBoot = ctr$pBoot))
  } else {
    cl <- parallel::makeCluster(ctr$nCore)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    z <- parallel::clusterApplyLB(cl = cl, x = as.list(1:nX),
                                  fun = sharpeScreeningXYi, X = X, Y = Y, T = T,
                                  nY = nY, nBoot = ctr$nBoot, bsids = bsids,
                                  minObs = ctr$minObs, type = ctr$type,
                                  hac = ctr$hac, b = ctr$bBoot, ttype = ctr$ttype,
                                  pBoot = ctr$pBoot)
  }

  for (i in 1:nX) {
    out <- z[[i]]
    dsharpe[i, ] <- out$dsharpei
    pval[i, ]    <- out$pvali
    tstat[i, ]   <- out$tstati
  }

  pi <- computePi(pval = pval, dalpha = dsharpe, tstat = tstat,
                  lambda = ctr$lambda, nBoot = ctr$nBoot,
                  bpos = ctr$gammaPos, bneg = ctr$gammaNeg,
                  fast = ctr$fastAdjust)
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
      # indices pre-generated in the master for this pair's length
      bs <- bsids[[as.character(sum(ok))]]
      if (is.null(bs)) {
        next  # block length exceeds this pair's sample size
      }
      tmp <- msharpeTestBootstrap(rets, level, na.neg, bs, b, ttype, pBoot)
    }
    if (!is.finite(tmp$tstat)) {
      next  # degenerate pair or NA modified VaR (na.neg)
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
  if (ctr$bBoot == 0) {
    stop("'bBoot = 0' (data-driven block length) is not supported in screening; ",
         "set an explicit block length 'bBoot >= 1'.")
  }
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  T <- nrow(X)
  if (nrow(Y) != T) {
    stop("'X' and 'Y' must have the same number of rows (time periods)")
  }
  nX <- ncol(X)
  nY <- ncol(Y)
  if (nX < 1L || nY < 1L) {
    stop("cross-group screening needs at least one fund in 'X' and one peer in 'Y'")
  }

  # pre-generate the bootstrap indices in the master, one matrix per distinct
  # cross-pair complete-case length (workers select by length)
  # only needed for the bootstrap test, and only for pairs passing 'minObs',
  # so that the asymptotic path never depends on 'bBoot'
  pairLens <- crossprod(1 * (!is.na(X) & !is.nan(X)), 1 * (!is.na(Y) & !is.nan(Y)))
  bsids <- NULL
  if (ctr$type == 2) {
    bsids <- bootIndicesByLen(pairLens[pairLens >= ctr$minObs],
                              ctr$nBoot, ctr$bBoot)
  }
  pval <- dmsharpe <- tstat <- matrix(NA, nX, nY)

  if (ctr$nCore == 1) {
    # serial path: no PSOCK cluster (avoids the per-call cluster overhead);
    # the anonymous wrapper avoids the clash with lapply's own 'X' argument
    z <- lapply(1:nX, function(ii)
      msharpeScreeningXYi(ii, X = X, Y = Y, level = level, T = T, nY = nY,
                          nBoot = ctr$nBoot, bsids = bsids,
                          minObs = ctr$minObs, na.neg = na.neg,
                          type = ctr$type, hac = ctr$hac, b = ctr$bBoot,
                          ttype = ctr$ttype, pBoot = ctr$pBoot))
  } else {
    cl <- parallel::makeCluster(ctr$nCore)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    z <- parallel::clusterApplyLB(cl = cl, x = as.list(1:nX),
                                  fun = msharpeScreeningXYi, X = X, Y = Y,
                                  level = level, T = T, nY = nY, nBoot = ctr$nBoot,
                                  bsids = bsids, minObs = ctr$minObs,
                                  na.neg = na.neg, type = ctr$type, hac = ctr$hac,
                                  b = ctr$bBoot, ttype = ctr$ttype, pBoot = ctr$pBoot)
  }

  for (i in 1:nX) {
    out <- z[[i]]
    dmsharpe[i, ] <- out$dmsharpei
    pval[i, ]     <- out$pvali
    tstat[i, ]    <- out$tstati
  }

  pi <- computePi(pval = pval, dalpha = dmsharpe, tstat = tstat,
                  lambda = ctr$lambda, nBoot = ctr$nBoot,
                  bpos = ctr$gammaPos, bneg = ctr$gammaNeg,
                  fast = ctr$fastAdjust)
  info <- infoFund(X, level = level, na.neg = na.neg)

  out <- list(n = info$nObs, npeer = rowSums(!is.na(pval)), ny = nY,
              msharpe = info$msharpe, dmsharpe = dmsharpe, pval = pval,
              tstat = tstat, lambda = pi$lambda, pizero = pi$pizero,
              pipos = pi$pipos, pineg = pi$pineg, cross = TRUE)
  class(out) <- "SCREENING"
  return(out)
}
msharpeScreeningXY <- compiler::cmpfun(.msharpeScreeningXY)

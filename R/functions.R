## Set of various R functions

#@name processControl
#@title Control parameters processsing
processControl <- function(control) {
  if (!is.list(control) || length(control) == 0) {
    control <- list(type = 1, ttype = 2, hac = FALSE, nBoot = 499, 
                    bBoot = 1, pBoot = 1, nCore = 1, minObs = 10, minObsPi = 1, 
                    lambda = NULL)
  }
  nam <- names(control)
  if (!("type" %in% nam) || is.null(control$type)) {
    control$type <- 1
  }
  if (!("ttype" %in% nam) || is.null(control$ttype)) {
    control$ttype <- 2
  }
  if (!("hac" %in% nam) || is.null(control$hac)) {
    control$hac <- FALSE
  }
  if (!("nBoot" %in% nam) || is.null(control$nBoot)) {
    control$nBoot <- 499
  }
  if (!("bBoot" %in% nam) || is.null(control$bBoot)) {
    control$bBoot <- 1
  }
  if (!("pBoot" %in% nam) || is.null(control$pBoot)) {
    control$pBoot <- 1
  }
  if (!("nCore" %in% nam) || is.null(control$nCore)) {
    control$nCore <- 1
  }
  if (!("minObs" %in% nam) || is.null(control$minObs)) {
    control$minObs <- 10
  }
  if (!("minObsPi" %in% nam) || is.null(control$minObsPi)) {
    control$minObsPi <- 1
  }
  if (!("lambda" %in% nam) || is.null(control$lambda)) {
    control$lambda <- NULL
  }
  return(control)
}

#@name alphaFactor
#@title Compute alpha factor
alphaFactor <- function(X, factors = NULL) {
  fit <- stats::lm(X ~ 1 + factors)
  alpha <- as.vector(fit$coef[1, ])
  return(alpha)
}

#' @name sharpe
#' @title Compute Sharpe ratio
#' @description Function which computes the Sharpe ratio.
#' @details The Sharpe ratio (Sharpe 1992) is one industry standard for measuring the
#' absolute risk adjusted performance of hedge funds.
#' @param X vector (of lenght \eqn{T}) or matrix (of size \eqn{T \times
#' N}{TxN}) of returns for \eqn{N} funds. \code{NA} values are allowed.
#' @param na.rm a logical value indicating whether \code{NA} values should be
#' stripped before the computation. Default \code{na.rm = TRUE}
#' @return A scalar or a vector (of size \eqn{N}) with the Sharpe ratios.
#' @note Please cite the package in publications. Use
#' \code{citation('PeerPerformance')}.
#' @author David Ardia and Kris Boudt.
#' @seealso \code{\link{sharpeTesting}}, \code{\link{sharpeScreening}} and
#' \code{\link{msharpe}}.
#' @references 
#' Ardia, D., Boudt, K. (2015).  Testing equality of modified
#' Sharpe ratios \emph{Finance Research Letters} \bold{13}, pp.97--104.
#' 
#' Ardia, D., Boudt, K. (2016).  \emph{The Peer Performance Ratios of Hedge
#' Funds}.  \url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2000901}
#' 
#' Sharpe, W. F. (1994).  The Sharpe ratio.  \emph{Journal of Portfolio
#' Management} Fall, pp.49--58.
#' @keywords htest
#' @examples
#' ## Load the data
#' data('hfdata')
#' 
#' ## Compute the Sharpe ratio
#' out = sharpe(hfdata)
#' print(out)
#' 
#' out = sharpe(hfdata, na.rm = FALSE)
#' print(out)
#' @export
sharpe <- function(X, na.rm = TRUE) {
  X <- as.matrix(X)
  N <- ncol(X)
  tmp <- .sharpe(X, na.rm)
  out <- tmp$mu.hat/tmp$sig.hat
  return(out)
}

#@name .sharpe
#@title Compute Sharpe ratio
.sharpe <- function(X, na.rm) {
  nObs <- colSums(!is.nan(X), na.rm = na.rm)
  mu.hat <- colMeans(X, na.rm = na.rm)
  X_ <- sweep(x = X, MARGIN = 2, STATS = mu.hat, FUN = "-")
  sig.hat <- sqrt(colSums(X_^2, na.rm = na.rm)/(nObs - 1))
  out <- list(mu.hat = mu.hat, sig.hat = sig.hat)
  return(out)
}

#' @name msharpe
#' @title Compute modified Sharpe ratio
#' @description Function which computes the modified Sharpe ratio
#' @details The modified Sharpe ratio (Favre and Galeano 2002) is one industry
#' standard for measuring the absolute risk adjusted performance of hedge
#' funds.
#' 
#' @param X vector (of lenght \eqn{T}) or matrix (of size \eqn{T \times
#' N}{TxN}) of returns.  \code{NA} values are allowed.
#' @param level modified Value-at-Risk level. Default: \code{level = 0.90}.
#' @param na.rm a logical value indicating whether \code{NA} values should be
#' stripped before the computation. Default \code{na.rm = TRUE}.
#' @param na.neg a logical value indicating whether \code{NA} values should be
#' returned if a negative modified Value-at-Risk is obtained.  Default
#' \code{na.neg = TRUE}.
#' @return scalar or a vector (of size \eqn{N}) with the modified Sharpe
#' ratios.
#' @note Please cite the package in publications. Use
#' \code{citation('PeerPerformance')}.
#' @author David Ardia and Kris Boudt.
#' @seealso \code{\link{msharpeTesting}}, \code{\link{msharpeScreening}} and
#' \code{\link{sharpe}}.
#' @references 
#' Ardia, D., Boudt, K. (2015).  Testing equality of modified
#' Sharpe ratios \emph{Finance Research Letters} \bold{13}, pp.97--104.
#' 
#' Ardia, D., Boudt, K. (2016).  \emph{The Peer Performance Ratios of Hedge
#' Funds}.  \url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2000901}
#' 
#' Favre, L., Galeano, J.A. (2002).  Mean-modified Value-at-Risk with Hedge
#' Funds.  \emph{The Journal of Alternative Investments} \bold{5}, pp.21--25.
#' 
#' Gregoriou, G. N., Gueyie, J.-P. (2003).  Risk-adjusted performance of funds
#' of hedge funds using a modified Sharpe ratio.  \emph{The Journal of Wealth
#' Management} \bold{Winter}, pp.77--83.
#' @keywords htest
#' @examples
#' ## Load the data (randomized data of monthly hedge fund returns)
#' data('hfdata')
#' 
#' out = msharpe(hfdata)
#' print(out)
#' 
#' out = msharpe(hfdata, na.rm = FALSE)
#' print(out)
#' @export
#' @importFrom stats qnorm
msharpe <- function(X, level = 0.9, na.rm = TRUE, na.neg = TRUE) {
  X <- as.matrix(X)
  N <- ncol(X)
  tmp <- .msharpe(X, level, na.rm, na.neg)
  out <- tmp$m1/tmp$mVaR
  return(out)
}

.msharpe <- function(X, level, na.rm, na.neg) {
  m1 <- colMeans(X, na.rm = na.rm)
  X_ <- sweep(x = X, MARGIN = 2, STATS = m1, FUN = "-")
  m2 <- colMeans(X_^2, na.rm = na.rm)
  m3 <- colMeans(X_^3, na.rm = na.rm)
  m4 <- colMeans(X_^4, na.rm = na.rm)
  za <- stats::qnorm(1 - level)
  skew <- m3/m2^(3/2)
  kurt <- (m4/m2^2) - 3
  mVaR <- -m1 + sqrt(m2) * (-za - (1/6) * (za^2 - 1) * skew - (1/24) * 
                              (za^3 - 3 * za) * kurt + (1/36) * (2 * za^3 - 5 * za) * skew^2)
  if (na.neg) {
    mVaR[mVaR < 0] <- NA
  }
  out <- list(m1 = m1, mVaR = mVaR)
  return(out)
}

# #' @name .infoFund
# #' @import compiler
.infoFund <- function(X, factors = NULL, level = NULL, na.rm = TRUE, na.neg = TRUE) {
  X <- as.matrix(X)
  N <- ncol(X)
  nObs <- colSums(!is.nan(X))
  muX <- colMeans(X, na.rm = na.rm)
  rX <- sweep(x = X, MARGIN = 2, STATS = muX, FUN = "-")
  sigX <- sqrt(colSums(rX^2, na.rm = na.rm)/(nObs - 1))
  sharpe_ <- muX/sigX
  
  if (is.null(factors)) {
    fit <- stats::lm(X ~ 1)
  } else {
    fit <- stats::lm(X ~ 1 + factors)
  }
  alpha_ <- as.vector(fit$coef[1, ])
  
  msharpe_ <- NULL
  if (!is.null(level)) {
    msharpe_ <- msharpe(X, level = level, na.rm = na.rm, na.neg = na.neg)
  }
  
  out <- list(nObs = nObs, mu = muX, sig = sigX, sharpe = sharpe_, alpha = alpha_, 
              msharpe = msharpe_)
  return(out)
}
infoFund <- cmpfun(.infoFund)

# #' @name .bootIndices
# #' @import compiler
.bootIndices <- function(T, nBoot, bBoot) {
  idsBoot <- matrix(data = NA, nrow = T, ncol = nBoot)
  if (bBoot == 1) {
    idsBoot <- matrix(sample.int(T, size = T * nBoot, replace = TRUE), 
                      nrow = T, ncol = nBoot)
  } else {
    for (i in 1:nBoot) {
      l <- floor(T/bBoot)
      ids <- c(1:T, 1:bBoot)
      seqb <- vector("integer", T)
      start.points <- sample.int(T, size = l, replace = TRUE)
      for (j in (1:l)) {
        start <- start.points[j]
        seqb[((j - 1) * bBoot + 1):(j * bBoot)] <- ids[start:(start + 
                                                                bBoot - 1)]
      }
      idsBoot[, i] <- seqb
    }
  }
  return(idsBoot)
}
bootIndices <- compiler::cmpfun(.bootIndices)
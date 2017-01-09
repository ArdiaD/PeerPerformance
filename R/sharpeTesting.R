## Set of R functions for Sharpe ratio testing

#@name .sharpeTesting
#@title See sharpeTesting
.sharpeTesting <- function(x, y, control = list()) {
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  # process control parameters
  ctr <- processControl(control)
  
  # check if enough data are available for testing
  dxy <- x - y
  idx <- (!is.nan(dxy) & !is.na(dxy))
  rets <- cbind(x[idx], y[idx])
  T <- sum(idx)
  if (T < ctr$minObs) {
    stop("intersection of 'x' and 'y' is shorter than 'minObs'")
  }
  
  # sharpe testing
  if (ctr$type == 1) {
    # ==> asymptotic approach
    tmp <- sharpeTestAsymptotic(rets, ctr$hac, ctr$ttype)
  } else {
    # ==> bootstrap approach (iid and circular block bootstrap)
    if (ctr$bBoot == 0) {
      ctr$bBoot <- sharpeBlockSize(x, y, ctr)
    }
    bsids <- bootIndices(T, ctr$nBoot, ctr$bBoot)
    tmp <- sharpeTestBootstrap(rets, bsids, ctr$bBoot, ctr$ttype, ctr$pBoot)
  }
  
  # info on the funds
  info <- infoFund(rets)
  
  ## form output
  out <- list(n = T, sharpe = info$sharpe, dsharpe = as.vector(tmp$dsharpe), 
              tstat = as.vector(tmp$tstat), pval = as.vector(tmp$pval))
  return(out)
}

#' @name sharpeTesting
#' @title Testing the difference of Sharpe ratios
#' @description Function which performs the testing of the difference of Sharpe ratios.
#' @details The Sharpe ratio (Sharpe 1992) is one industry standard for measuring the
#' absolute risk adjusted performance of hedge funds. This function performs
#' the testing of Sharpe ratio difference for two funds using the approach by
#' Ledoit and Wolf (2002).
#' 
#' For the testing, only the intersection of non-\code{NA} observations for the
#' two funds are used.
#' 
#' The argument \code{control} is a list that can supply any of the following
#' components:
#' \itemize{
#' \item \code{'type'} asymptotic approach (\code{type = 1}) or
#' studentized circular bootstrap approach (\code{type = 2}). Default:
#' \code{type = 1}.
#' \item \code{'ttype'} test based on ratio (\code{type = 1})
#' or product (\code{type = 2}). Default: \code{type = 2}.
#' \item \code{'hac'} heteroscedastic-autocorrelation consistent standard
#' errors. Default: \code{hac = FALSE}.
#' \item \code{'nBoot'} number of boostrap replications for computing the p-value. Default: \code{nBoot =
#' 499}. 
#' \item \code{'bBoot'} block length in the circular bootstrap. Default:
#' \code{bBoot = 1}, i.e. iid bootstrap. \code{bBoot = 0} uses optimal
#' block-length. 
#' \item \code{'pBoot'} symmetric p-value (\code{pBoot = 1}) or
#' asymmetric p-value (\code{pBoot = 2}). Default: \code{pBoot = 1}.
#' }
#' @param x vector (of lenght \eqn{T}) of returns for the first fund. \code{NA}
#' values are allowed.
#' @param y vector (of lenght \eqn{T}) returns for the second fund. \code{NA}
#' values are allowed.
#' @param control control parameters (see *Details*).
#' @return % A list with the following components:\cr
#' 
#' \code{n}: number of non-\code{NA} concordant observations.\cr
#' 
#' \code{sharpe}: vector (of length 2) of unconditional Sharpe ratios.\cr
#' 
#' \code{dsharpe}: Sharpe ratios difference.\cr
#' 
#' \code{tstat}: t-stat of Sharpe ratios differences.\cr
#' 
#' \code{pval}: pvalues of test of Sharpe ratios differences.
#' @note Further details on the methdology with an application to the hedge
#' fund industry is given in in Ardia and Boudt (2016). 
#' 
#' Some internal functions where adapted from Wolf's R code.
#' 
#' Please cite the package in publications. Use
#' \code{citation('PeerPerformance')}.
#' @author David Ardia and Kris Boudt.
#' @seealso \code{\link{sharpe}}, \code{\link{sharpeScreening}} and
#' \code{\link{msharpeTesting}}.
#' @references 
#' Ardia, D., Boudt, K. (2015).  Testing equality of modified
#' Sharpe ratios \emph{Finance Research Letters} \bold{13}, pp.97--104.
#' 
#' Ardia, D., Boudt, K. (2016).  \emph{The Peer Performance Ratios of Hedge
#' Funds}.  \url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2000901}
#' 
#' Ledoit, O., Wolf, M. (2008). Robust performance hypothesis testing with the
#' Sharpe ratio.  \emph{Journal of Empirical Finance} \bold{15}, pp.850--859.
#' 
#' Sharpe, W. F. (1994).  The Sharpe ratio.  \emph{Journal of Portfolio
#' Management} Fall, pp.49--58.
#' @keywords htest
#' @examples
#' ## Load the data (randomized data of monthly hedge fund returns)
#' data("hfdata")
#' x = hfdata[,1]
#' y = hfdata[,2]
#' 
#' ## Run Sharpe testing (asymptotic)
#' ctr = list(type = 1)
#' out = sharpeTesting(x, y, control = ctr)
#' print(out)
#'   
#' ## Run Sharpe testing (asymptotic hac)
#' ctr = list(type = 1, hac = TRUE)
#' out = sharpeTesting(x, y, control = ctr)
#' print(out)
#'   
#' ## Run Sharpe testing (iid bootstrap)
#' set.seed(1234)
#' ctr = list(type = 2, nBoot = 250)
#' out = sharpeTesting(x, y, control = ctr)
#' print(out)
#'   
#' ## Run Sharpe testing (circular bootstrap)
#' set.seed(1234)
#' ctr = list(type = 2, nBoot = 250, bBoot = 5)
#' out = sharpeTesting(x, y, control = ctr)
#' print(out)
#' @export
#' @importFrom compiler cmpfun
sharpeTesting <- compiler::cmpfun(.sharpeTesting)

#@name .sharpe.ratio.diff
#@title Difference of sharpe ratios
.sharpe.ratio.diff <- function(X, Y, ttype) {
  
  if (is.null(Y)) {
    Y <- X[, 2, drop = FALSE]
    X <- X[, 1, drop = FALSE]
  }
  n <- nrow(X)
  mu1.hat <- colMeans(X)
  mu2.hat <- colMeans(Y)
  X_ <- sweep(x = X, MARGIN = 2, STATS = mu1.hat, FUN = "-")
  Y_ <- sweep(x = Y, MARGIN = 2, STATS = mu2.hat, FUN = "-")
  sig1.hat <- sqrt(colSums(X_^2)/(n - 1))
  sig2.hat <- sqrt(colSums(Y_^2)/(n - 1))
  
  if (ttype == 1) {
    SR1.hat <- mu1.hat/sig1.hat
    SR2.hat <- mu2.hat/sig2.hat
  } else {
    SR1.hat <- mu1.hat * sig2.hat
    SR2.hat <- mu2.hat * sig1.hat
  }
  diff <- SR1.hat - SR2.hat
  return(diff)
}
sharpe.ratio.diff <- compiler::cmpfun(.sharpe.ratio.diff)

#' @name .sharpeTestAsymptotic
#' @title Asymptotic Sharpe test
#' @importFrom stats pnorm
.sharpeTestAsymptotic <- function(rets, hac, ttype) {
  
  dsharpe <- sharpe.ratio.diff(rets, Y = NULL, ttype)
  se <- se.sharpe.asymptotic(rets, hac, ttype)
  tstat <- dsharpe/se
  pval <- 2 * stats::pnorm(-abs(tstat))  # asymptotic normal p-value
  out <- list(dsharpe = dsharpe, tstat = tstat, se = se, pval = pval)
  return(out)
}
sharpeTestAsymptotic <- compiler::cmpfun(.sharpeTestAsymptotic)

#' @name .se.sharpe.asymptotic
#' @title Asymptotic standard error
#' @importFrom stats cov ar
.se.sharpe.asymptotic <- function(X, hac, ttype) {
  
  # estimation of (robust) Psi function; see Ledoit Wolf paper
  compute.Psi.hat <- function(V.hat, hac) {
    
    if (hac) {
      T <- length(V.hat[, 1])
      alpha.hat <- compute.alpha.hat(V.hat)
      S.star <- 2.6614 * (alpha.hat * T)^0.2
      Psi.hat <- compute.Gamma.hat(V.hat, 0)
      j <- 1
      while (j < S.star) {
        Gamma.hat <- compute.Gamma.hat(V.hat, j)
        Psi.hat <- Psi.hat + kernel.Parzen(j/S.star) * (Gamma.hat + 
                                                          t(Gamma.hat))
        j <- j + 1
      }
      Psi.hat <- (T/(T - 4)) * Psi.hat
    } else {
      Psi.hat <- stats::cov(V.hat)
    }
    return(Psi.hat)
    
  }
  
  # Parzen kernel
  kernel.Parzen <- function(x) {
    
    if (abs(x) <= 0.5) 
      result <- 1 - 6 * x^2 + 6 * abs(x)^3 else if (abs(x) <= 1) 
        result <- 2 * (1 - abs(x))^3 else result <- 0
        return(result)
        
  }
  
  compute.alpha.hat <- function(V.hat) {
    
    p <- ncol(V.hat)
    num <- den <- 0
    for (i in 1:p) {
      fit <- stats::ar(V.hat[, i], 0, 1, method = "ols")
      rho.hat <- as.numeric(fit[2])
      sig.hat <- sqrt(as.numeric(fit[3]))
      num <- num + 4 * rho.hat^2 * sig.hat^4/(1 - rho.hat)^8
      den <- den + sig.hat^4/(1 - rho.hat)^4
    }
    return(num/den)
    
  }
  
  compute.Gamma.hat <- function(V.hat, j) {
    
    T <- nrow(V.hat)
    p <- ncol(V.hat)
    Gamma.hat <- matrix(0, p, p)
    if (j >= T) 
      stop("j must be smaller than the row dimension!")
    for (i in ((j + 1):T)) {
      Gamma.hat <- Gamma.hat + tcrossprod(V.hat[i, ], V.hat[i - j, 
                                                            ])
    }
    Gamma.hat <- Gamma.hat/T
    return(Gamma.hat)
    
  }
  
  T <- nrow(X)
  
  if (ttype == 1) {
    mu.hat <- colMeans(X)
    gamma.hat <- colMeans(X^2)
    gradient <- vector("double", 4)
    gradient[1] <- gamma.hat[1]/(gamma.hat[1] - mu.hat[1]^2)^1.5
    gradient[2] <- -gamma.hat[2]/(gamma.hat[2] - mu.hat[2]^2)^1.5
    gradient[3] <- -0.5 * mu.hat[1]/(gamma.hat[1] - mu.hat[1]^2)^1.5
    gradient[4] <- 0.5 * mu.hat[2]/(gamma.hat[2] - mu.hat[2]^2)^1.5
    V.hat <- matrix(NA, T, 4)
    V.hat[, 1:2] <- sweep(x = X, MARGIN = 2, STATS = mu.hat, FUN = "-")
    V.hat[, 3:4] <- sweep(x = X^2, MARGIN = 2, STATS = gamma.hat, FUN = "-")
  } else {
    m1 <- colMeans(X)
    X_ <- sweep(x = X, MARGIN = 2, STATS = m1, FUN = "-")
    m2 <- colMeans(X_^2)
    g2 <- m2 + m1^2
    dm1i <- c(1, 0, 0, 0)
    dm1j <- c(0, 0, 1, 0)
    dsigi <- 1/(2 * sqrt(m2[1])) * c(-2 * m1[1], 1, 0, 0)
    dsigj <- 1/(2 * sqrt(m2[2])) * c(0, 0, -2 * m1[2], 1)
    tmp1 <- dm1i * sqrt(m2[2]) + dsigj * m1[1]
    tmp2 <- dm1j * sqrt(m2[1]) + dsigi * m1[2]
    gradient <- tmp1 - tmp2
    V.hat <- matrix(NA, T, 4)
    V.hat[, c(1, 3)] <- sweep(x = X, MARGIN = 2, STATS = m1, FUN = "-")
    V.hat[, c(2, 4)] <- sweep(x = X^2, MARGIN = 2, STATS = g2, FUN = "-")
  }
  
  Psi.hat <- compute.Psi.hat(V.hat, hac)
  se <- as.numeric(sqrt(crossprod(gradient, Psi.hat %*% gradient)/T))
  return(se)
  
}
se.sharpe.asymptotic <- compiler::cmpfun(.se.sharpe.asymptotic)

#@name .sharpeTestBootstrap
#@title Test Sharpe difference using circular studentized boostrap of Ledoit and Wolf
.sharpeTestBootstrap <- function(rets, bsids, b, ttype, pBoot, d = 0) {
  
  T <- nrow(rets)
  x <- rets[, 1, drop = FALSE]
  y <- rets[, 2, drop = FALSE]
  dsharpe <- as.numeric(sharpe.ratio.diff(x, y, ttype) - d)
  se <- se.sharpe.bootstrap(x, y, b, ttype)
  # se = se.sharpe.asymptotic(X = cbind(x, y), hac = TRUE, ttype = ttype)
  
  # bootstrap indices
  nBoot <- ncol(bsids)
  bsidx <- 1 + bsids%%T  # ensure that the bootstrap indices match the length of the time series
  bsX <- matrix(x[bsidx], T, nBoot)
  bsY <- matrix(y[bsidx], T, nBoot)
  
  bsdsharpe <- sharpe.ratio.diff(bsX, bsY, ttype)
  bsse <- se.sharpe.bootstrap(bsX, bsY, b, ttype)
  tstat <- dsharpe/se
  
  if (pBoot == 1) {
    # first type p-value calculation
    bststat <- abs(bsdsharpe - dsharpe)/bsse
    pval <- (sum(bststat > abs(tstat)) + 1)/(nBoot + 1)
    # pval = sum(bststat > abs(tstat)) / nBoot
  } else {
    # second type p-value calculation (as in Barras)
    bststat <- (bsdsharpe - dsharpe)/bsse
    pval <- 2 * min(sum(bststat > tstat) + 1, sum(bststat < tstat) + 
                      1)/(nBoot + 1)
    # pval = 2 * min(sum(bststat > tstat), sum(bststat < tstat)) / nBoot
  }
  
  out <- list(dsharpe = dsharpe, tstat = tstat, se = se, bststat = bststat, 
              pval = pval)
  return(out)
}
sharpeTestBootstrap <- compiler::cmpfun(.sharpeTestBootstrap)

#' @name .se.sharpe.bootstrap
#' @title Bootstrap standard error
#' @importFrom stats cov
.se.sharpe.bootstrap <- function(X, Y, b, ttype) {
  
  ## Compute Psi with two approaches: 1) iid bootstrap, 2) circular block
  ## bootstrap
  compute.Psi.hat <- function(V.hat, b) {
    
    T <- length(V.hat[, 1])
    if (b == 1) {
      # ==> standard estimation
      Psi.hat <- stats::cov(V.hat)
    } else {
      # ==> block estimation
      l <- floor(T/b)
      Psi.hat <- matrix(0, 4, 4)
      for (j in (1:l)) {
        zeta <- b^0.5 * colMeans(V.hat[((j - 1) * b + 1):(j * b), 
                                       , drop = FALSE])
        Psi.hat <- Psi.hat + tcrossprod(zeta)
      }
      Psi.hat <- Psi.hat/l
    }
    return(Psi.hat)
    
  }
  
  T <- nrow(X)
  N <- ncol(Y)
  if (ttype == 1) {
    mu1.hat <- colMeans(X)
    mu2.hat <- colMeans(Y)
    gamma1.hat <- colMeans(X^2)
    gamma2.hat <- colMeans(Y^2)
    gradient <- array(NA, c(4, 1, N))
    gradient[1, 1, ] <- gamma1.hat/(gamma1.hat - mu1.hat^2)^1.5
    gradient[2, 1, ] <- -gamma2.hat/(gamma2.hat - mu2.hat^2)^1.5
    gradient[3, 1, ] <- -0.5 * mu1.hat/(gamma1.hat - mu1.hat^2)^1.5
    gradient[4, 1, ] <- 0.5 * mu2.hat/(gamma2.hat - mu2.hat^2)^1.5
    V.hat <- array(NA, c(T, 4, N))
    V.hat[, 1, ] <- sweep(x = X, MARGIN = 2, STATS = mu1.hat, FUN = "-")
    V.hat[, 2, ] <- sweep(x = Y, MARGIN = 2, STATS = mu2.hat, FUN = "-")
    V.hat[, 3, ] <- sweep(x = X^2, MARGIN = 2, STATS = gamma1.hat, 
                          FUN = "-")
    V.hat[, 4, ] <- sweep(x = Y^2, MARGIN = 2, STATS = gamma2.hat, 
                          FUN = "-")
  } else {
    m1X <- colMeans(X)
    m1Y <- colMeans(Y)
    X_ <- sweep(x = X, MARGIN = 2, STATS = m1X, FUN = "-")
    Y_ <- sweep(x = Y, MARGIN = 2, STATS = m1Y, FUN = "-")
    m2X <- colMeans(X_^2)
    m2Y <- colMeans(Y_^2)
    g2X <- m2X + m1X^2
    g2Y <- m2Y + m1Y^2
    
    cst1X <- 1/(2 * sqrt(m2X))
    cst1Y <- 1/(2 * sqrt(m2Y))
    
    dm1X <- matrix(rep(c(1, 0, 0, 0), N), 4, N, FALSE)
    dm1Y <- matrix(rep(c(0, 0, 1, 0), N), 4, N, FALSE)
    dsigX <- rbind(-2 * cst1X * m1X, cst1X, 0, 0)
    dsigY <- rbind(0, 0, -2 * cst1Y * m1Y, cst1Y)
    
    # matrix form
    m1X_ <- matrix(m1X, nrow = 4, ncol = N, byrow = TRUE)
    m1Y_ <- matrix(m1Y, nrow = 4, ncol = N, byrow = TRUE)
    m2X_ <- matrix(m2X, nrow = 4, ncol = N, byrow = TRUE)
    m2Y_ <- matrix(m2Y, nrow = 4, ncol = N, byrow = TRUE)
    
    dm1X_ <- matrix(dm1X, nrow = 4, ncol = N, byrow = FALSE)
    dm1Y_ <- matrix(dm1Y, nrow = 4, ncol = N, byrow = FALSE)
    dsigX_ <- matrix(dsigX, nrow = 4, ncol = N, byrow = FALSE)
    dsigY_ <- matrix(dsigY, nrow = 4, ncol = N, byrow = FALSE)
    
    cst2X_ <- sqrt(m2X_)
    cst2Y_ <- sqrt(m2Y_)
    
    # gradient
    tmp1 <- dm1X_ * cst2Y_ + dsigY_ * m1X_
    tmp2 <- dm1Y_ * cst2X_ + dsigX_ * m1Y_
    
    # =======
    gradient <- array(NA, c(4, 1, N))
    gradient[1:4, 1, ] <- tmp1 - tmp2
    V.hat <- array(NA, c(T, 4, N))
    V.hat[, 1, ] <- sweep(x = X, MARGIN = 2, STATS = m1X, FUN = "-")
    V.hat[, 3, ] <- sweep(x = Y, MARGIN = 2, STATS = m1Y, FUN = "-")
    V.hat[, 2, ] <- sweep(x = X^2, MARGIN = 2, STATS = g2X, FUN = "-")
    V.hat[, 4, ] <- sweep(x = Y^2, MARGIN = 2, STATS = g2Y, FUN = "-")
  }
  Psi.hat <- array(apply(X = V.hat, MARGIN = 3, FUN = compute.Psi.hat, 
                         b = b), c(4, 4, N))
  se <- vector("double", N)
  for (i in 1:N) {
    se[i] <- sqrt(crossprod(gradient[, , i], Psi.hat[, , i] %*% gradient[, 
                                                                         , i])/T)
  }
  return(se)
  
}
se.sharpe.bootstrap <- compiler::cmpfun(.se.sharpe.bootstrap)

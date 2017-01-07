## Set of R functions for the optimal block length computation for the
## modified Sharpe ratio

#' @name .msharpeBlockSize
#' @title See msharpeBlockSize
#' @importFrom stats rgeom
.msharpeBlockSize <- function(x, y, level = 0.9, na.neg = TRUE, control = list(), 
                              b.vec = c(1, 3, 6, 10), alpha = 0.05, M = 199, K = 500, b.av = 5, T.start = 50) {
  
  sb.sequence <- function(T, b.av, length = T) {
    
    index.sequence <- c(1:T, 1:T)
    sequence <- rep.int(0, length + T)
    current <- 0
    while (current < length) {
      start <- sample(1:T, 1)
      b <- stats::rgeom(1, 1/b.av) + 1
      sequence[(current + 1):(current + b)] <- index.sequence[start:(start + 
                                                                       b - 1)]
      current <- current + b
    }
    out <- sequence[1:length]
    return(out)
    
  }
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  rets <- cbind(x, y)
  
  # process control parameters
  ctr <- processControl(control)
  
  b.len <- length(b.vec)
  emp.reject.probs <- rep.int(0, b.len)
  d <- msharpe.ratio.diff(x, y, level, na.neg, ctr$ttype)
  T <- length(x)
  Var.data <- matrix(data = 0, nrow = T.start + T, ncol = 2)
  Var.data[1, ] <- rets[1, ]
  fit1 <- lm(x[2:T] ~ x[1:(T - 1)] + y[1:(T - 1)])
  fit2 <- lm(y[2:T] ~ x[1:(T - 1)] + y[1:(T - 1)])
  coef1 <- as.numeric(fit1$coef)
  coef2 <- as.numeric(fit2$coef)
  resid.mat <- cbind(as.numeric(fit1$resid), as.numeric(fit2$resid))
  for (k in 1:K) {
    ids <- sb.sequence(T - 1, b.av, T.start + T - 1)
    resid.mat.star <- rbind(c(0, 0), resid.mat[ids, ])
    for (t in 2:(T.start + T)) {
      Var.data[t, 1] <- coef1[1] + coef1[2] * Var.data[t - 1, 1] + 
        coef1[3] * Var.data[t - 1, 2] + resid.mat.star[t, 1]
      Var.data[t, 2] <- coef2[1] + coef2[2] * Var.data[t - 1, 2] + 
        coef2[3] * Var.data[t - 1, 2] + resid.mat.star[t, 2]
    }
    Var.data.trunc <- Var.data[(T.start + 1):(T.start + T), ]
    for (j in 1:b.len) {
      bsids <- bootIndices(T, M, b.vec[j])
      tmp <- msharpeTestBootstrap(Var.data.trunc, level, na.neg, 
                                  bsids, b.vec[j], ctr$ttype, pBoot = 1, d)
      if (tmp$pval <= alpha) {
        emp.reject.probs[j] <- emp.reject.probs[j] + 1
      }
    }
  }
  emp.reject.probs <- emp.reject.probs/K
  b.order <- order(abs(emp.reject.probs - alpha))
  b.opt <- b.vec[b.order[1]]
  return(b.opt)
  
}

#' @name msharpeBlockSize
#' @title Optimal block length for bootstrap test of difference of modified Sharpe
#' ratios
#' @description Function which computes the optimal block length for testing of the
#' difference of modified Sharpe ratios
#' @details The argument \code{control} is a list that can supply any of the following
#' components:
#' \itemize{
#' \item \code{'type'} asymptotic approach (\code{type = 1}) or
#' studentized circular bootstrap approach (\code{type = 2}). Default:
#' \code{type = 1}.
#' \item \code{'ttype'} test based on ratio (\code{type = 1})
#' or product (\code{type = 2}). Default: \code{type = 2}.
#' \item \code{'hac'} heteroscedastic-autocorrelation consistent standard
#' errors. Default: \code{hac = FALSE}.
#' \item \code{'minObs'} minimum number
#' of concordant observations to compute the ratios. Default: \code{minObs =
#' 10}. 
#' \item \code{'nBoot'} number of boostrap replications for computing the
#' p-value. Default: \code{nBoot = 499}.
#' \item \code{'bBoot'} block length in
#' the circular bootstrap. Default: \code{bBoot = 1}, i.e. iid bootstrap.
#' \code{bBoot = 0} uses optimal block-length. 
#' \item \code{'pBoot'} symmetric p-value (\code{pBoot = 1}) or asymmetric p-value (\code{pBoot = 2}).
#' Default: \code{pBoot = 1}.
#' }
#' @param x vector (of lenght \eqn{T}) of returns for the first fund. \code{NA}
#' values are allowed.
#' @param y vector (of lenght \eqn{T}) returns for the second fund. \code{NA}
#' values are allowed.
#' @param level modified Value-at-Risk level. Default: \code{level = 0.90}.
#' @param na.neg a logical value indicating whether \code{NA} values should be
#' returned if a negative modified Value-at-Risk is obtained.  Default
#' \code{na.neg = TRUE}.
#' @param control control parameters (see *Details*).
#' @param b.vec vector of block to be tested.
#' @param alpha significance level.
#' @param M bootstrap replications.
#' @param K number of cross-validation.
#' @param b.av average block length in the stationary bootstrap.
#' @param T.start starting point for bootstrap
#' @return The optimal block length.
#' @note Please cite the package in publications. Use
#' \code{citation('PeerPerformance')}.
#' @author David Ardia and Kris Boudt.
#' @seealso \code{\link{msharpe}}, \code{\link{msharpeScreening}} and
#' \code{\link{sharpeTesting}}.
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
#' 
#' Ledoit, O., Wolf, M. (2008).  Robust performance hypothesis testing with the
#' Sharpe ratio.  \emph{Journal of Empirical Finance} \bold{15}, pp.850--859.
#' 
#' Sharpe, W. F. (1994).  The Sharpe ratio.  \emph{Journal of Portfolio
#' Management} Fall, pp.49--58.
#' @keywords htest
#' @examples
#' # !!! DA ADD EXAMPLES HERE !!!
#' @export
msharpeBlockSize <- compiler::cmpfun(.msharpeBlockSize)

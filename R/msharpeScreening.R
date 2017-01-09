## Set of R functions for modified Sharpe screening

#' @name .msharpeScreening
#' @title See msharpeScreening
#' @importFrom snow makeCluster clusterApply stopCluster
.msharpeScreening <- function(X, level = 0.9, na.neg = TRUE, control = list()) {
  
  # process control
  ctr <- processControl(control)
  
  # size of inputs and outputs
  T <- nrow(X)
  N <- ncol(X)
  pval <- dmsharpe <- tstat <- matrix(data = NA, N, N)
  
  # determine which pairs can be compared (in a matrix way)
  Y <- 1 * (!is.nan(X) & !is.na(X))
  YY <- crossprod(Y)  #YY = t(Y) %*% Y # row i indicates how many observations in common with column k
  YY[YY < ctr$minObs] <- 0
  YY[YY > 0] <- 1
  liststocks <- c(1:nrow(YY))[rowSums(YY) > ctr$minObsPi]
  
  # determine bootstrap indices (do it before to speed up computations)
  bsids <- bootIndices(T, ctr$nBoot, ctr$bBoot)
  
  if (length(liststocks) > 1) {
    cl <- snow::makeCluster(c(rep("localhost", ctr$nCore)), type = "SOCK")
    
    liststocks <- liststocks[1:(length(liststocks) - 1)]
    
    z <- snow::clusterApply(cl = cl, x = as.list(liststocks), fun = msharpeScreeningi, 
                            rdata = X, level = level, T = T, N = N, na.neg = na.neg, nBoot = ctr$nBoot, 
                            bsids = bsids, minObs = ctr$minObs, type = ctr$type, hac = ctr$hac, 
                            b = ctr$bBoot, ttype = ctr$ttype, pBoot = ctr$pBoot)
    snow::stopCluster(cl)
    
    for (i in 1:length(liststocks)) {
      out <- z[[i]]
      id <- liststocks[i]
      pval[id, id:N] <- pval[id:N, id] <- out[[2]][id:N]
      dmsharpe[id, id:N] <- out[[1]][id:N]
      dmsharpe[id:N, id] <- -out[[1]][id:N]
      tstat[id, id:N] <- out[[3]][id:N]
      tstat[id:N, id] <- -out[[3]][id:N]
    }
  }
  
  # pi
  pi <- computePi(pval = pval, dalpha = dmsharpe, tstat = tstat, lambda = ctr$lambda, 
                  nBoot = ctr$nBoot)
  
  # info on the funds
  info <- infoFund(X, level = level, na.neg = na.neg)
  
  # form output
  out <- list(n = info$nObs, npeer = colSums(!is.na(pval)), msharpe = info$msharpe, 
              dmsharpe = dmsharpe, pval = pval, tstat = tstat, lambda = pi$lambda, 
              pizero = pi$pizero, pipos = pi$pipos, pineg = pi$pineg)
  
  return(out)
}

#' @name msharpeScreening
#' @title Screening using the modified Sharpe outperformance ratio
#' @description Function which performs the screening of a universe of returns, and
#' compute the modified Sharpe outperformance ratio.
#' @details The modified Sharpe ratio (Favre and Galeano 2002, Gregoriou and Gueyie
#' 2003) is one industry standard for measuring the absolute risk adjusted
#' performance of hedge funds. We propose to complement the modified Sharpe
#' ratio with the fund's outperformance ratio, defined as the percentage number
#' of funds that have a significantly lower modified Sharpe ratio. In a
#' pairwise testing framework, a fund can have a significantly higher modified
#' Sharpe ratio because of luck. We correct for this by applying the false
#' discovery rate approach by (Storey 2002).
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
#' \item \code{'nBoot'} number of
#' boostrap replications for computing the p-value. Default: \code{nBoot =
#' 499}.
#' \item \code{'bBoot'} block length in the circular bootstrap. Default:
#' \code{bBoot = 1}, i.e. iid bootstrap. \code{bBoot = 0} uses optimal
#' block-length.
#' \item \code{'pBoot'} symmetric p-value (\code{pBoot = 1}) or
#' asymmetric p-value (\code{pBoot = 2}). Default: \code{pBoot = 1}.
#' \item \code{'nCore'} number of cores to be used. Default: \code{nCore = 1}.
#' \item \code{'minObs'} minimum number of concordant observations to compute
#' the ratios. Default: \code{minObs = 10}.
#' \item \code{'minObsPi'} minimum number of observations to compute pi0. Default: \code{minObsPi = 1}.
#' \item \code{'lambda'} threshold value to compute pi0. Default: \code{lambda
#' = NULL}, i.e. data driven choice.
#' }
#' @param X matrix \eqn{(T \times N)}{(TxN)} of \eqn{T} returns for the \eqn{N}
#' funds. \code{NA} values are allowed.
#' @param level modified Value-at-Risk level. Default: \code{level = 0.90}.
#' @param na.neg a logical value indicating whether \code{NA} values should be
#' returned if a negative modified Value-at-Risk is obtained.  Default
#' \code{na.neg = TRUE}.
#' @param control control parameters (see *Details*).
#' @return % A list with the following components:\cr
#' 
#' \code{n}: vector (of length \eqn{N}) of number of non-\code{NA}
#' observations.\cr
#' 
#' \code{npeer}: vector (of length \eqn{N}) of number of available peers.\cr
#' 
#' \code{msharpe}: vector (of length \eqn{N}) of unconditional modified Sharpe
#' ratios.\cr
#' 
#' \code{dmsharpe}: matrix (of size \eqn{N \times N}{NxN}) of modified Sharpe
#' ratios differences.\cr
#' 
#' \code{pval}: matrix (of size \eqn{N \times N}) of p-values of test for
#' modified Sharpe ratios differences.\cr
#' 
#' \code{lambda}: vector (of length \eqn{N}) of lambda values.\cr
#' 
#' \code{pizero}: vector (of length \eqn{N}) of probability of equal
#' performance.\cr
#' 
#' \code{pipos}: vector (of length \eqn{N}) of probability of outperformance
#' performance.\cr
#' 
#' \code{pineg}: vector (of length \eqn{N}) of probability of underperformance
#' performance.
#' @note % Further details on the methdology with an application to the hedge
#' fund industry is given in in Ardia and Boudt (2016). The file
#' \file{ThePeerPerformanceOfHedgeFunds.txt} in the \file{/doc} package's
#' folder allows the reprodution of the steps followed in the article. See also
#' the presentation by Kris Boudt at the R/Finance conference 2012 at
#' \url{http://www.rinfinance.com}.
#' 
#' Some internal functions where adapted from Wolf's R code.
#' 
#' Application of the false discovery rate approach applied to the mutual fund
#' industry has been presented in Barraz, Scaillet and Wermers (2010).
#' 
#' Please cite the package in publications. Use
#' \code{citation('PeerPerformance')}.
#' @author David Ardia and Kris Boudt.
#' @seealso \code{\link{msharpe}}, \code{\link{msharpeTesting}},
#' \code{\link{sharpeScreening}} and \code{\link{alphaScreening}}.
#' @references 
#' Ardia, D., Boudt, K. (2016).  \emph{The Peer Performance
#' Ratios of Hedge Funds}.
#' \url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2000901}
#' 
#' Barras, L., Scaillet, O., Wermers, R. (2010).  False discoveries in mutual
#' fund performance: Measuring luck in estimated alphas.  \emph{Journal of
#' Finance} \bold{5}, pp.179--216.
#' 
#' Favre, L., Galeano, J.A. (2002).  Mean-modified Value-at-Risk with Hedge
#' Funds \emph{The Journal of Alternative Investments} \bold{5}, pp.21--25.
#' 
#' Gregoriou, G. N., Gueyie, J.-P. (2003).  Risk-adjusted performance of funds
#' of hedge funds using a modified Sharpe ratio \emph{The Journal of Wealth
#' Management} \bold{Winter}, pp.77--83.
#' 
#' Ledoit, O., Wolf, M. (2008). Robust performance hypothesis testing with the
#' Sharpe ratio.  \emph{Journal of Empirical Finance} \bold{15}, pp.850--859.
#' 
#' Sharpe, W. F. (1994).  The Sharpe ratio.  \emph{Journal of Portfolio
#' Management} Fall, pp.49--58.
#' 
#' Storey, J. (2002).  A direct approach to false discovery rates.
#' \emph{Journal of the Royal Statistical Society B} \bold{64}, pp.479--498.
#' @keywords htest
#' @examples
#' ## Load the data (randomized data of monthly hedge fund returns)
#' data("hfdata")
#' rets = hfdata[,1:10]
#' 
#' ## Modified Sharpe screening 
#' msharpeScreening(rets, control = list(nCore = 1))
#' 
#' ## Modified Sharpe screening with bootstrap and HAC standard deviation
#' msharpeScreening(rets, control = list(nCore = 1, type = 2, hac = TRUE))
#' @export
#' @importFrom compiler cmpfun
msharpeScreening <- compiler::cmpfun(.msharpeScreening)

#@name .msharpeScreeningi
#@title Sharpe ratio screening for fund i again its peers
.msharpeScreeningi <- function(i, rdata, level, T, N, nBoot, bsids, minObs, 
                               na.neg, type, hac, b, ttype, pBoot) {
  
  nPeer <- N - i
  X <- matrix(rdata[, i], nrow = T, ncol = nPeer)
  Y <- matrix(rdata[, (i + 1):N], nrow = T, ncol = nPeer)
  
  dXY <- X - Y
  idx <- (!is.nan(dXY) & !is.na(dXY))
  X[!idx] <- NA
  Y[!idx] <- NA
  nObs <- colSums(idx)
  
  pvali <- dmsharpei <- tstati <- rep(NA, N)
  
  k <- 0
  for (j in (i + 1):N) {
    k <- k + 1
    if (nObs[k] < minObs) {
      next
    }
    rets <- cbind(X[idx[, k], 1], Y[idx[, k], k])
    
    if (type == 1) {
      tmp <- msharpeTestAsymptotic(rets, level, na.neg, hac, ttype)
    } else {
      tmp <- msharpeTestBootstrap(rets, level, na.neg, bsids, b, 
                                  ttype, pBoot)
    }
    
    dmsharpei[j] <- tmp$dmsharpe
    pvali[j] <- tmp$pval
    tstati[j] <- tmp$tstat
  }
  
  out <- list(dmsharpei = dmsharpei, pvali = pvali, tstati = tstati)
  return(out)
}
msharpeScreeningi <- compiler::cmpfun(.msharpeScreeningi)
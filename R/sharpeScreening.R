## Set of R functions for Sharpe screening

# #' @name .sharpeScreening
# #' @import compiler
.sharpeScreening <- function(X, control = list()) {
  
  # process control
  ctr <- processControl(control)
  
  # size of inputs and outputs
  T <- nrow(X)
  N <- ncol(X)
  pval <- dsharpe <- tstat <- matrix(data = NA, N, N)
  
  # determine which pairs can be compared (in a matrix way)
  Y <- 1 * (!is.nan(X) & !is.na(X))
  YY <- crossprod(Y)  #YY = t(Y) %*% Y # row i indicates how many observations in common with column k
  YY[YY < ctr$minObs] <- 0
  YY[YY > 0] <- 1
  liststocks <- c(1:nrow(YY))[rowSums(YY) > ctr$minObsPi]
  
  # determine bootstrap indices (do it before to speed up computations)
  bsids <- bootIndices(T, ctr$nBoot, ctr$bBoot)
  
  if (length(liststocks) > 1) {
    cl <- makeCluster(c(rep("localhost", ctr$nCore)), type = "SOCK")
    
    liststocks <- liststocks[1:(length(liststocks) - 1)]
    
    z <- clusterApply(cl = cl, x = as.list(liststocks), fun = sharpeScreeningi, 
                      rdata = X, T = T, N = N, nBoot = ctr$nBoot, bsids = bsids, 
                      minObs = ctr$minObs, type = ctr$type, hac = ctr$hac, b = ctr$bBoot, 
                      ttype = ctr$ttype, pBoot = ctr$pBoot)
    stopCluster(cl)
    
    for (i in 1:length(liststocks)) {
      out <- z[[i]]
      id <- liststocks[i]
      pval[id, id:N] <- pval[id:N, id] <- out[[2]][id:N]
      dsharpe[id, id:N] <- out[[1]][id:N]
      dsharpe[id:N, id] <- -out[[1]][id:N]
      tstat[id, id:N] <- out[[3]][id:N]
      tstat[id:N, id] <- -out[[3]][id:N]
    }
  }
  
  # pi
  pi <- computePi(pval = pval, dalpha = dsharpe, tstat = tstat, lambda = ctr$lambda, 
                  nBoot = ctr$nBoot)
  
  # info on the funds
  info <- infoFund(X)
  
  # form output
  out <- list(n = info$nObs, npeer = colSums(!is.na(pval)), sharpe = info$sharpe, 
              dsharpe = dsharpe, pval = pval, tstat = tstat, lambda = pi$lambda, 
              pizero = pi$pizero, pipos = pi$pipos, pineg = pi$pineg)
  
  return(out)
}

#' @name sharpeScreening
#' @title Screening using the Sharpe outperformance ratio
#' @description Function which performs the screening of a universe of returns, and
#' compute the Sharpe outperformance ratio.
#' @details The Sharpe ratio (Sharpe 1992) is one industry standard for measuring the
#' absolute risk adjusted performance of hedge funds. We propose to complement
#' the Sharpe ratio with the fund's outperformance ratio, defined as the
#' percentage number of funds that have a significantly lower Sharpe ratio. In
#' a pairwise testing framework, a fund can have a significantly higher Sharpe
#' ratio because of luck. We correct for this by applying the false discovery
#' rate approach by (Storey 2002).
#' 
#' For the testing, only the intersection of non-\code{NA} observations for the
#' two funds are used.
#' 
#' The methodology proceeds as follows: 
#' \itemize{
#' \item (1) compute all
#' pairwise tests of Sharpe differences using the bootstrap approach of Ledoit
#' and Wolf (2002). This means that for a universe of \eqn{N} funds, we perform
#' \eqn{N(N-1)/2}{N*(N-1)/2} tests. The algorithm has been parallelized and the
#' computational burden can be slip across several cores. The number of cores
#' can be defined in \code{control}, see below.
#' \item (2) for each fund, the
#' false discovery rate approach by Storey (2002) is used to determine the
#' proportions over, equal, and underperfoming funds, in terms of Sharpe ratio,
#' in the database.
#' }
#' The argument \code{control} is a list that can supply any of the following
#' components:
#' \itemize{
#' \item \code{'type'} Asymptotic approach (\code{type = 1}) or
#' studentized circular bootstrap approach (\code{type = 2}). Default:
#' \code{type = 1}.
#' \item \code{'ttype'} Test based on ratio (\code{type = 1})
#' or product (\code{type = 2}). Default: \code{type = 2}.
#' \item \code{'hac'} Heteroscedastic-autocorrelation consistent standard
#' errors. Default: \code{hac = FALSE}.
#' \item \code{'nBoot'} Number of boostrap replications for computing the p-value. Default: \code{nBoot =
#' 499}. 
#' \item \code{'bBoot'} Block length in the circular bootstrap. Default:
#' \code{bBoot = 1}, i.e. iid bootstrap. \code{bBoot = 0} uses optimal
#' block-length. 
#' \item \code{'pBoot'} Symmetric p-value (\code{pBoot = 1}) or
#' asymmetric p-value (\code{pBoot = 2}). Default: \code{pBoot = 1}.
#' \item \code{'nCore'} Number of cores to be used. Default: \code{nCore = 1}.
#' \item \code{'minObs'} Minimum number of concordant observations to compute
#' the ratios. Default: \code{minObs = 10}.
#' \item \code{'minObsPi'} Minimum
#' number of observations to compute pi0. Default: \code{minObsPi = 1}.
#' \item \code{'lambda'} Threshold value to compute pi0. Default: \code{lambda
#' = NULL}, i.e. data driven choice.
#' }
#' @param X Matrix \eqn{(T \times N)}{(TxN)} of \eqn{T} returns for the \eqn{N}
#' funds. \code{NA} values are allowed.
#' @param control Control parameters (see *Details*).
#' @return A list with the following components:\cr
#' 
#' \code{n}: Vector (of length \eqn{N}) of number of non-\code{NA}
#' observations.\cr
#' 
#' \code{npeer}: Vector (of length \eqn{N}) of number of available peers.\cr
#' 
#' \code{sharpe}: Vector (of length \eqn{N}) of unconditional Sharpe ratios.\cr
#' 
#' \code{dsharpe}: Matrix (of size \eqn{N \times N}{NxN}) of Sharpe ratios
#' differences.\cr
#' 
#' \code{pval}: Matrix (of size \eqn{N \times N}) of pvalues of test for Sharpe
#' ratios differences.\cr
#' 
#' \code{lambda}: vector (of length \eqn{N}) of lambda values.\cr
#' 
#' \code{pizero}: vector (of length \eqn{N}) of probability of equal
#' performance.\cr
#' 
#' \code{pipos}: vector (of length \eqn{N}) of probability of outperformance
#' performance.\cr
#' 
#' \code{pineg}: Vector (of length \eqn{N}) of probability of underperformance
#' performance.
#' @note Further details on the methdology with an application to the hedge
#' fund industry is given in in Ardia and Boudt (2016).
#' 
#' Some internal functions where adapted from Michael Wolf MATLAB code.
#' 
#' Application of the false discovery rate approach applied to the mutual fund
#' industry has been presented in Barraz, Scaillet and Wermers (2010).
#' @author David Ardia and Kris Boudt.
#' @seealso \code{\link{sharpe}}, \code{\link{sharpeTesting}},
#' \code{\link{msharpeScreening}} and \code{\link{alphaScreening}}.
#' @references 
#' Ardia, D., Boudt, K. (2015).  
#' Testing equality of modified Sharpe ratios.
#' \emph{Finance Research Letters} \bold{13}, pp.97--104. 
#' \doi{10.1016/j.frl.2015.02.008}
#' 
#' Ardia, D., Boudt, K. (2016).  
#' The Peer Ratios Performance of Hedge Funds. 
#' \emph{Working paper}.
#' \doi{10.2139/ssrn.2000901}
#' 
#' Barras, L., Scaillet, O., Wermers, R. (2010).  
#' False discoveries in mutual fund performance: Measuring luck in estimated alphas.  
#' \emph{Journal of Finance} \bold{65}(1), pp.179--216.
#' \doi{10.1111/j.1540-6261.2009.01527.x}
#' 
#' Sharpe, W.F. (1994).  
#' The Sharpe ratio.  
#' \emph{Journal of Portfolio Management} \bold{21}(1), pp.49--58.
#' \doi{10.3905/jpm.1994.409501}
#' 
#' Ledoit, O., Wolf, M. (2008). 
#' Robust performance hypothesis testing with the Sharpe ratio.  
#' \emph{Journal of Empirical Finance} \bold{15}(5), pp.850--859.
#' \doi{10.1016/j.jempfin.2008.03.002}
#' 
#' Storey, J. (2002).  
#' A direct approach to false discovery rates.
#' \emph{Journal of the Royal Statistical Society B} \bold{64}(3), pp.479--498.
#' \doi{10.1111/1467-9868.00346}
#' @keywords htest
#' @examples
#' ## Load the data (randomized data of monthly hedge fund returns)
#' data("hfdata")
#' rets = hfdata[,1:10]
#' 
#' ## Sharpe screening 
#' sharpeScreening(rets, control = list(nCore = 1))
#' 
#' ## Sharpe screening with bootstrap and HAC standard deviation
#' sharpeScreening(rets, control = list(nCore = 1, type = 2, hac = TRUE))
#' @export
#' @import compiler
sharpeScreening <- compiler::cmpfun(.sharpeScreening)

#@name .sharpeScreeningi
#@title Sharpe ratio screening for fund i again its peers
.sharpeScreeningi <- function(i, rdata, T, N, nBoot, bsids, minObs, type, 
                              hac, b, ttype, pBoot) {
  
  nPeer <- N - i
  X <- matrix(rdata[, i], nrow = T, ncol = nPeer)
  Y <- matrix(rdata[, (i + 1):N], nrow = T, ncol = nPeer)
  
  dXY <- X - Y
  idx <- (!is.nan(dXY) & !is.na(dXY))
  X[!idx] <- NA
  Y[!idx] <- NA
  nObs <- colSums(idx)
  
  pvali <- dsharpei <- tstati <- rep(NA, N)
  
  k <- 0
  for (j in (i + 1):N) {
    k <- k + 1
    if (nObs[k] < minObs) {
      next
    }
    rets <- cbind(X[idx[, k], 1], Y[idx[, k], k])
    
    if (type == 1) {
      tmp <- sharpeTestAsymptotic(rets, hac, ttype)
    } else {
      tmp <- sharpeTestBootstrap(rets, bsids, b, ttype, pBoot)
    }
    
    dsharpei[j] <- tmp$dsharpe
    pvali[j] <- tmp$pval
    tstati[j] <- tmp$tstat
  }
  
  out <- list(dsharpei = dsharpei, pvali = pvali, tstati = tstati)
  return(out)
}
sharpeScreeningi <- compiler::cmpfun(.sharpeScreeningi)
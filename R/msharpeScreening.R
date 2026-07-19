## Set of R functions for modified Sharpe screening

# #' @name .msharpeScreening
# #' @title See msharpeScreening
# #' @importFrom parallel makeCluster clusterApplyLB stopCluster
# #' @import compiler
.msharpeScreening <- function(X, level = 0.9, na.neg = TRUE, control = list(),
                              Y = NULL) {

  # cross-group screening: each fund in X against every fund in group Y
  if (!is.null(Y)) {
    return(.msharpeScreeningXY(X, Y, level = level, na.neg = na.neg,
                               control = control))
  }

  # process control
  ctr <- processControl(control)
  if (ctr$bBoot == 0) {
    stop("'bBoot = 0' (data-driven block length) is not supported in screening; ",
         "set an explicit block length 'bBoot >= 1'.")
  }

  # size of inputs and outputs
  X <- as.matrix(X)
  T <- nrow(X)
  N <- ncol(X)
  if (N < 2L) {
    stop("within-group screening needs at least two funds in 'X'; supply 'Y' to screen a single fund against a peer group")
  }
  pval <- dmsharpe <- tstat <- matrix(data = NA, N, N)

  # determine which pairs can be compared (in a matrix way)
  Y <- 1 * (!is.nan(X) & !is.na(X))
  YY <- crossprod(Y)  #YY = t(Y) %*% Y # row i indicates how many observations in common with column k
  pairLens <- YY[upper.tri(YY)]  # pairwise complete-case lengths (before thresholding)
  YY[YY < ctr$minObs] <- 0
  YY[YY > 0] <- 1
  liststocks <- c(1:nrow(YY))[rowSums(YY) > ctr$minObsPi]

  # pre-generate the bootstrap indices in the master, one matrix per distinct
  # pair length (workers select by length; see bootIndicesByLen)
  bsids <- bootIndicesByLen(pairLens, ctr$nBoot, ctr$bBoot)

  if (length(liststocks) > 1) {
    liststocks <- liststocks[1:(length(liststocks) - 1)]

    if (ctr$nCore == 1) {
      # serial path: no PSOCK cluster (avoids the per-call cluster overhead)
      z <- lapply(as.list(liststocks), msharpeScreeningi,
                  rdata = X, level = level, T = T, N = N, na.neg = na.neg,
                  nBoot = ctr$nBoot, bsids = bsids, minObs = ctr$minObs,
                  type = ctr$type, hac = ctr$hac, b = ctr$bBoot,
                  ttype = ctr$ttype, pBoot = ctr$pBoot)
    } else {
      cl <- parallel::makeCluster(ctr$nCore)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      z <- parallel::clusterApplyLB(cl = cl, x = as.list(liststocks), fun = msharpeScreeningi,
                                    rdata = X, level = level, T = T, N = N, na.neg = na.neg, nBoot = ctr$nBoot,
                                    bsids = bsids, minObs = ctr$minObs, type = ctr$type, hac = ctr$hac,
                                    b = ctr$bBoot, ttype = ctr$ttype, pBoot = ctr$pBoot)
    }

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
                  nBoot = ctr$nBoot, bpos = ctr$gammaPos, bneg = ctr$gammaNeg,
                  fast = ctr$fastAdjust)

  # info on the funds
  info <- infoFund(X, level = level, na.neg = na.neg)

  # form output
  out <- list(n = info$nObs, npeer = colSums(!is.na(pval)), msharpe = info$msharpe,
              dmsharpe = dmsharpe, pval = pval, tstat = tstat, lambda = pi$lambda,
              pizero = pi$pizero, pipos = pi$pipos, pineg = pi$pineg)
  class(out) <- "SCREENING"

  return(out)
}

#' @name msharpeScreening
#' @title Screening using the modified Sharpe outperformance ratio
#' @description Function which performs the screening of a universe of returns, and
#' computes the modified Sharpe outperformance ratio.
#' @details The modified Sharpe ratio (Favre and Galeano 2002, Gregoriou and Gueyie
#' 2003) is one industry standard for measuring the absolute risk adjusted
#' performance of hedge funds. We propose to complement the modified Sharpe
#' ratio with the fund's outperformance ratio, defined as the percentage number
#' of funds that have a significantly lower modified Sharpe ratio. In a
#' pairwise testing framework, a fund can have a significantly higher modified
#' Sharpe ratio because of luck. We correct for this by applying the false
#' discovery rate approach by Storey (2002).
#'
#' For the testing, only the intersection of non-\code{NA} observations for the
#' two funds are used.
#'
#' The argument \code{control} is a list that can supply any of the following
#' components:
#' \itemize{
#' \item \code{'type'} Asymptotic approach (\code{type = 1}) or
#' studentized circular bootstrap approach (\code{type = 2}). Default:
#' \code{type = 1}.
#' \item \code{'ttype'} Test based on ratio (\code{type = 1})
#' or product (\code{type = 2}). Default: \code{type = 2}.
#' \item \code{'hac'} heteroscedastic-autocorrelation consistent standard
#' errors. Default: \code{hac = FALSE}.
#' \item \code{'nBoot'} Number of
#' bootstrap replications for computing the p-value. Default: \code{nBoot =
#' 499}.
#' \item \code{'bBoot'} Block length in the circular bootstrap. Default:
#' \code{bBoot = 1}, i.e. iid bootstrap. (The data-driven choice
#' \code{bBoot = 0} is only available in \code{\link{msharpeTesting}}, not in
#' screening.)
#' \item \code{'pBoot'} Symmetric p-value (\code{pBoot = 1}) or
#' asymmetric p-value (\code{pBoot = 2}). Default: \code{pBoot = 1}.
#' \item \code{'nCore'} Number of cores to be used. Default: \code{nCore = 1}.
#' \item \code{'minObs'} Minimum number of concordant observations to compute
#' the ratios. Default: \code{minObs = 10}.
#' \item \code{'minObsPi'} Minimum number of observations to compute pi0. Default: \code{minObsPi = 1}.
#' \item \code{'lambda'} Threshold value to compute pi0. Default: \code{lambda
#' = NULL}, i.e. data driven choice.
#' \item \code{'gammaPos'} One-sided quantile level (of the standard Normal
#' distribution) used as the critical value for counting outperformed peers:
#' a peer counts as outperformed when the pairwise t-statistic exceeds
#' \code{qnorm(gammaPos)} (a \emph{negative} threshold for
#' \code{gammaPos < 0.5}), and the expected fraction \code{1 - gammaPos} of
#' false positives among the equal-performing peers is then subtracted.
#' Default: \code{gammaPos = 0.4} (the value recommended in Ardia and Boudt, 2018).
#' \item \code{'gammaNeg'} Mirror image of \code{gammaPos} for the peers that
#' outperform the focal fund: the count uses \code{tstat <= qnorm(gammaNeg)}
#' and subtracts the expected fraction \code{gammaNeg} of false positives.
#' Default: \code{gammaNeg = 0.6}.
#' \item \code{'fastAdjust'} Use a fast vectorised inversion in the
#' truncated-normal bias correction of \eqn{\pi^0}{pi0} instead of one
#' \code{uniroot} call per value. This is the dominant cost when \code{lambda}
#' is data driven and gives a large speed-up on big universes; the result agrees
#' with the default path to about 1e-12, i.e. well inside \code{uniroot}'s own
#' tolerance. Default: \code{fastAdjust = FALSE}, i.e. the original code path,
#' kept as default for exact reproducibility of published results.
#' }
#' @param X Matrix \eqn{(T \times N)}{(TxN)} of \eqn{T} returns for the \eqn{N}
#' funds. \code{NA} values are allowed.
#' @param level Modified Value-at-Risk level. Default: \code{level = 0.90}.
#' @param na.neg A logical value indicating whether \code{NA} values should be
#' returned if a negative modified Value-at-Risk is obtained.  Default
#' \code{na.neg = TRUE}.
#' @param control Control parameters (see *Details*).
#' @param Y Optional matrix \eqn{(T \times M)}{(TxM)} of returns for a second
#' (peer) group of \eqn{M} funds. When supplied, the ratios are computed for
#' each fund in \code{X} against the funds in \code{Y} (cross-group screening)
#' instead of against the other funds in \code{X}; a single focal fund versus a
#' peer group corresponds to \code{X} being a vector. Columns of \code{Y}
#' identical to the focal fund are automatically excluded. Default:
#' \code{Y = NULL}, i.e. within-group screening.
#' @return A list with the following components:\cr
#'
#' \code{n}: Vector (of length \eqn{N}) of number of non-\code{NA}
#' observations.\cr
#'
#' \code{npeer}: Vector (of length \eqn{N}) of number of available peers.\cr
#'
#' \code{msharpe}: Vector (of length \eqn{N}) of unconditional modified Sharpe
#' ratios.\cr
#'
#' \code{dmsharpe}: Matrix (of size \eqn{N \times N}{NxN}) of modified Sharpe
#' ratios differences.\cr
#'
#' \code{tstat}: Matrix (of size \eqn{N \times N}{NxN}) of t-statistics.\cr
#'
#' \code{pval}: Matrix (of size \eqn{N \times N}{NxN}) of p-values of test for
#' modified Sharpe ratios differences.\cr
#'
#' \code{lambda}: Vector (of length \eqn{N}) of lambda values.\cr
#'
#' \code{pizero}: Vector (of length \eqn{N}) of probability of equal
#' performance.\cr
#'
#' \code{pipos}: Vector (of length \eqn{N}) of probability of outperformance
#' performance.\cr
#'
#' \code{pineg}: Vector (of length \eqn{N}) of probability of underperformance
#' performance.
#' @note Further details on the methodology with an application to the hedge
#' fund industry is given in Ardia and Boudt (2018).
#'
#' Some internal functions where adapted from Michael Wolf MATLAB code.
#'
#' Application of the false discovery rate approach applied to the mutual fund
#' industry has been presented in Barras, Scaillet and Wermers (2010).
#' @author David Ardia and Kris Boudt.
#' @seealso \code{\link{msharpe}}, \code{\link{msharpeTesting}},
#' \code{\link{sharpeScreening}} and \code{\link{alphaScreening}}.
#' @references
#' Ardia, D., Boudt, K. (2015).
#' Testing equality of modified Sharpe ratios.
#' \emph{Finance Research Letters} \bold{13}, 97--104.
#'
#' Ardia, D., Boudt, K. (2018).
#' The peer performance ratios of hedge funds.
#' \emph{Journal of Banking and Finance} \bold{87}, 351--368.
#'
#' Barras, L., Scaillet, O., Wermers, R. (2010).
#' False discoveries in mutual fund performance: Measuring luck in estimated alphas.
#' \emph{Journal of Finance} \bold{65}(1), 179--216.
#'
#' Favre, L., Galeano, J.A. (2002).
#' Mean-modified Value-at-Risk Optimization with Hedge Funds.
#' \emph{Journal of Alternative Investments} \bold{5}(2), 21--25.
#'
#' Gregoriou, G. N., Gueyie, J.-P. (2003).
#' Risk-adjusted performance of funds of hedge funds using a modified Sharpe ratio.
#' \emph{Journal of Wealth Management} \bold{6}(3), 77--83.
#'
#' Ledoit, O., Wolf, M. (2008).
#' Robust performance hypothesis testing with the Sharpe ratio.
#' \emph{Journal of Empirical Finance} \bold{15}(5), 850--859.
#'
#' Storey, J. (2002).
#' A direct approach to false discovery rates.
#' \emph{Journal of the Royal Statistical Society B} \bold{64}(3), 479--498.
#' @keywords htest
#' @examples
#' ## Load the data (randomized data of monthly hedge fund returns)
#' data("hfdata")
#' rets = hfdata[,1:4]
#'
#' ## Modified Sharpe screening
#' msharpeScreening(rets, control = list(nCore = 1))
#'
#' ## Modified Sharpe screening with the studentized circular bootstrap
#' ## (the 'hac' flag is ignored when type = 2)
#' msharpeScreening(rets, control = list(nCore = 1, type = 2))
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
    rets <- cbind(X[idx[, k], k], Y[idx[, k], k])

    if (type == 1) {
      tmp <- msharpeTestAsymptotic(rets, level, na.neg, hac, ttype)
    } else {
      # indices pre-generated in the master for this pair's length
      bs <- bsids[[as.character(nObs[k])]]
      if (is.null(bs)) {
        next  # block length exceeds this pair's sample size
      }
      tmp <- msharpeTestBootstrap(rets, level, na.neg, bs, b,
                                  ttype, pBoot)
    }
    if (!is.finite(tmp$tstat)) {
      next  # degenerate pair or NA modified VaR (na.neg)
    }

    dmsharpei[j] <- tmp$dmsharpe
    pvali[j] <- tmp$pval
    tstati[j] <- tmp$tstat
  }

  out <- list(dmsharpei = dmsharpei, pvali = pvali, tstati = tstati)
  return(out)
}
msharpeScreeningi <- compiler::cmpfun(.msharpeScreeningi)

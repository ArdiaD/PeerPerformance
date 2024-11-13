## Set of R function for alpha screening

#@name .alphaScreening
#@description See alphaScreening
.alphaScreening <- function(X, factors = NULL, control = list(), screen_beta=FALSE) {

  # process control
  ctr <- processControl(control)

  T <- nrow(X)
  N <- ncol(X)

  if (screen_beta & !is.null(factors)) {
    row_return <- 1:(1 + ncol(factors))
    pval <- dalpha <- tstat <- array(rep(NA, N * N * (1 + ncol(factors))), dim = c((1 + ncol(factors)), N, N))
  } else {
    row_return <- 1
    pval <- dalpha <- tstat <- array(rep(NA, N*N), dim = c(1, N, N))
  }
  # pval <- dalpha <- tstat <- matrix(data = NA, N, N)

  # determine which pairs can be compared (in a matrix way)
  Y <- 1 * (!is.nan(X) & !is.na(X))
  YY <- crossprod(Y)  #YY = t(Y) %*% Y # row i indicates how many observations in common with column k
  YY[YY < ctr$minObs] <- 0
  YY[YY > 0] <- 1
  liststocks <- c(1:nrow(YY))[rowSums(YY) > ctr$minObsPi]

  if (length(liststocks) > 1) {
    cl <- parallel::makeCluster(ctr$nCore)

    liststocks <- liststocks[1:(length(liststocks) - 1)]

  	z <- parallel::clusterApplyLB(cl = cl, x = as.list(liststocks), fun = alphaScreeningi,
                                    rdata = X, factors = factors, T = T, N = N, hac = ctr$hac, screen_beta)

    parallel::stopCluster(cl)

    for (i in 1:length(liststocks)) {
      out <- z[[i]]
      id <- liststocks[i]
      pval[row_return, id, id:N] <- pval[row_return, id:N, id] <- out[[2]][row_return, id:N]
      dalpha[row_return, id, id:N] <- out[[1]][row_return, id:N]
      dalpha[row_return, id:N, id] <- -out[[1]][row_return, id:N]
      tstat[row_return, id, id:N] <- out[[3]][row_return, id:N]
      tstat[row_return, id:N, id] <- -out[[3]][row_return, id:N]
    }
  }

  # pi
  pi <- computePi(pval = pval, dalpha = dalpha, tstat = tstat, lambda = ctr$lambda,
                  nBoot = ctr$nBoot)


  # info on the funds
  info <- infoFund(X, factors = factors, screen_beta = screen_beta)

  if (screen_beta == FALSE) {
  	pval <- pval[1, , ]
  	dalpha <- dalpha[1, , ]
  	tstat <- tstat[1, , ]
  	npeer <- colSums(!is.na(pval))
  } else {
    npeer <- apply(!is.na(pval), c(1, 3), sum)
  }

  # form output
  out <- list(n = info$nObs, npeer = npeer, alpha = info$alpha,
              dalpha = dalpha, pval = pval, tstat = tstat, lambda = pi$lambda,
              pizero = pi$pizero, pipos = pi$pipos, pineg = pi$pineg)
  class(out) <- "SCREENING"

  return(out)
}

#' @name alphaScreening
#' @title Screening using the alpha outperformance ratio
#' @description Function which performs the screening of a universe of returns, and
#' computes the alpha outperformance ratio.
#' @details The alpha measure (Treynor and Black 1973, Carhart 1997, Fung and Hsieh
#' 2004) is one industry standard for measuring the absolute risk adjusted
#' performance of hedge funds. We propose to complement the alpha measure with
#' the fund's alpha outperformance ratio, defined as the percentage number of
#' funds that have a significantly lower alpha. In a pairwise testing
#' framework, a fund can have a significantly higher alpha because of luck. We
#' correct for this by applying the false discovery rate approach by Storey (2002).
#'
#' The methodology proceeds as follows:
#' \itemize{
#' \item (1) compute all pairwise tests of alpha differences. This means that for a universe of
#' \eqn{N} funds, we perform \eqn{N(N-1)/2}{N*(N-1)/2} tests. The algorithm has
#' been parallelized and the computational burden can be split across several
#' cores. The number of cores can be defined in \code{control}, see below.
#' \item (2) for each fund, the false discovery rate approach by Storey (2002)
#' is used to determine the proportions of over, equal, and underperforming
#' funds, in terms of alpha, in the database.}
#' The argument \code{control} is a list that can supply any of the following
#' components:
#' \itemize{
#' \item \code{'hac'} Heteroscedastic-autocorrelation consistent
#' standard errors. Default: \code{hac = FALSE}.
#' \item \code{'minObs'} Minimum number of concordant observations to compute the ratios. Default:
#' \code{minObs = 10}.
#' \item \code{'minObsPi'} Minimum number of observations
#' for computing the p-values). Default: \code{minObsPi = 1}.
#' \item \code{'nCore'} Number of cores used to perform the screening. Default:
#' \code{nCore = 1}.
#' \item \code{'lambda'} Threshold value to compute pi0.
#' Default: \code{lambda = NULL}, i.e. data driven choice.
#' }
#' @param X Matrix \eqn{(T \times N)}{(TxN)} of \eqn{T} returns for the \eqn{N}
#' funds. \code{NA} values are allowed.
#' @param factors Matrix \eqn{(T \times K)}{(TxK)} of \eqn{T} returns for the
#' \eqn{K} factors. \code{NA} values are allowed.
#' @param control Control parameters (see *Details*).
#' @param screen_beta Boolean to screen all factors' coefficients (beta).
#' Default: \code{screen_beta=FALSE} (i.e. only outputs the alpha).
#' If \code{screen_beta=TRUE}, each element of the returned list will have a new first dimension
#' representing each coefficient (the first one being alpha)
#' @return A list with the following components:\cr
#'
#' \code{n}: Vector (of length \eqn{N}) of number of non-\code{NA}
#' observations.\cr
#'
#' \code{npeer}: Vector (of length \eqn{N}) of number of available peers.\cr
#'
#' \code{alpha}: Vector (of length \eqn{N}) of unconditional alpha.\cr
#'
#' \code{dalpha}: Matrix (of size \eqn{N \times N}{NxN}) of alpha
#' differences.\cr
#'
#' \code{tstat}: Matrix (of size \eqn{N \times N}{NxN}) of t-statistics.\cr
#'
#' \code{pval}: Matrix (of size \eqn{N \times N}) of p-values of test for alpha
#' differences.\cr
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
#' Application of the false discovery rate approach applied to the mutual fund
#' industry has been presented in Barras, Scaillet and Wermers (2010).
#'
#' Currently, the HAC asymptotic and studentized circular block bootstrap
#' presented in Ledoit and Wolf (2008) are not supported by the
#' \code{alphaScreening} function.
#' @author David Ardia and Kris Boudt.
#' @seealso \code{\link{sharpeScreening}} and \code{\link{msharpeScreening}}.
#' @references
#' Ardia, D., Boudt, K. (2015).
#' Testing equality of modified Sharpe ratios.
#' \emph{Finance Research Letters} \bold{13}, pp.97--104.
#' \doi{10.1016/j.frl.2015.02.008}
#'
#' Ardia, D., Boudt, K. (2018).
#' The peer performance ratios of hedge funds.
#' \emph{Journal of Banking and Finance} \bold{87}, pp.351-.368.
#' \doi{10.1016/j.jbankfin.2017.10.014}
#'
#' Barras, L., Scaillet, O., Wermers, R. (2010).
#' False discoveries in mutual fund performance: Measuring luck in estimated alphas.
#' \emph{Journal of Finance} \bold{65}(1), pp.179--216.
#'
#' Carhart, M. (1997).
#' On persistence in mutual fund performance.
#' \emph{Journal of Finance} \bold{52}(1), pp.57--82.
#'
#' Fama, E., French, K. (2010).
#' Luck versus skill in the cross-section of mutual fund returns.
#' \emph{Journal of Finance} \bold{65}(5), pp.1915--1947.
#'
#' Fung, W., Hsieh, D. (2004).
#' Hedge fund benchmarks: A risk based approach.
#' \emph{Financial Analysts Journal} \bold{60}(5), pp.65--80.
#'
#' Storey, J. (2002).
#' A direct approach to false discovery rates.
#' \emph{Journal of the Royal Statistical Society B} \bold{64}(3), pp.479--498.
#'
#' Treynor, J. L., Black, F. (1973).
#' How to use security analysis to improve portfolio selection.
#' \emph{Journal of Business} \bold{46}(1), pp.66--86.
#' @keywords htest
#' @examples
#' ## Load the data (randomized data of monthly hedge fund returns)
#' data("hfdata")
#' rets = hfdata[,1:5]
#'
#' ## Run alpha screening
#' ctr = list(nCore = 1)
#' alphaScreening(rets, control = ctr)
#'
#' ## Run alpha screening with HAC standard deviation
#' ctr = list(nCore = 1, hac = TRUE)
#' alphaScreening(rets, control = ctr)
#' @export
#' @importFrom parallel makeCluster clusterApplyLB stopCluster
#' @importFrom compiler cmpfun
alphaScreening <- compiler::cmpfun(.alphaScreening)

# #' @name .alphaScreeningi
# #' @title Screening for fund i again its peers
# #' @importFrom stats lm na.omit
# #' @importFrom lmtest coeftest
# #' @importFrom sandwich vcovHAC
.alphaScreeningi <- function(i, rdata, factors, T, N, hac, screen_beta=FALSE) {


  if(screen_beta & !is.null(factors)){
    row_return <- 1:(1+ncol(factors))
    pvali <- dalphai <- tstati <- matrix(rep(NA, N*(1+ncol(factors))), ncol=N)
  }else{
    row_return <- 1
    pvali <- dalphai <- tstati <- matrix(rep(NA, N), ncol=N)
  }

  nPeer <- N - i
  X <- matrix(rdata[, i], nrow = T, ncol = nPeer)
  Y <- matrix(rdata[, (i + 1):N], nrow = T, ncol = nPeer)
  dXY <- X - Y

  # Additional filter: Nonoverlapping observations
  # Iterate over selIds
  D <- !is.na(dXY)
  # See that it has no shared obs with the second one.
  selId.in  <- which(colSums(D) != 0)
  selId.out <- selId.in + i

  if (nPeer == 1) {
    if (is.null(factors)) {
      fit <- stats::lm(dXY ~ 1, na.action = stats::na.omit)
    } else {
      fit <- stats::lm(dXY ~ 1 + factors, na.action = stats::na.omit)
    } # end of factors/no factors

    # HAC within loop.
    if (!hac) {
      sumfit <- summary(fit)
      pvali[row_return, N] <- sumfit$coef[row_return, 4]
      dalphai[row_return, N] <- sumfit$coef[row_return, 1]
      tstati[row_return, N] <- sumfit$coef[row_return, 3]
    } else {
      sumfit <- lmtest::coeftest(fit, vcov. = sandwich::vcovHAC(fit))
      pvali[row_return, N] <- sumfit[row_return, 4]
      dalphai[row_return, N] <- sumfit[row_return, 1]
      tstati[row_return, N] <- sumfit[row_return, 3]
    }
  } else {
    # end of nPeer == 1

    # k selects the columns in dXY
    #   k will match with redefined selId.in
    # j plugs them into the correct list, but it uses an efficient allocation "(i + 1):N"
    #   j matches with selId.out
    # We make a correction for k in (20190004)

    # k <- 1
    # for (j in (i + 1):N) {
    for (idx in 1:length(selId.in)) {

      # proper indices
      k <- selId.in[idx]
      j <- selId.out[idx]

      if (is.null(factors)) {
        fit <- stats::lm(dXY[, k] ~ 1, na.action = stats::na.omit)
      } else {
        fit <- stats::lm(dXY[, k] ~ 1 + factors, na.action = stats::na.omit)
      } # end of factors/no factors

      # HAC within loop.
      if (!hac) {
        sumfit <- summary(fit)
        pvali[row_return, j] <- sumfit$coef[row_return, 4]
        dalphai[row_return, j] <- sumfit$coef[row_return, 1]
        tstati[row_return, j] <- sumfit$coef[row_return, 3]
      } else{
        sumfit <- lmtest::coeftest(fit, vcov. = sandwich::vcovHAC(fit))
        pvali[row_return, j] <- sumfit[row_return, 4]
        dalphai[row_return, j] <- sumfit[row_return, 1]
        tstati[row_return, j] <- sumfit[row_return, 3]
      }

      # k <- k + 1
    }
  }

  out <- list(dalphai = dalphai, pvali = pvali, tstati = tstati)
  return(out)
}
alphaScreeningi <- compiler::cmpfun(.alphaScreeningi)

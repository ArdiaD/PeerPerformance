## Set of R function for alpha screening

#@name .alphaScreening
#@description See alphaScreening
.alphaScreening <- function(X, factors = NULL, control = list()) {
  
  # process control
  ctr <- processControl(control)
  
  T <- nrow(X)
  N <- ncol(X)
  pval <- dalpha <- tstat <- matrix(data = NA, N, N)
  
  # determine which pairs can be compared (in a matrix way)
  Y <- 1 * (!is.nan(X) & !is.na(X))
  YY <- crossprod(Y)  #YY = t(Y) %*% Y # row i indicates how many observations in common with column k
  YY[YY < ctr$minObs] <- 0
  YY[YY > 0] <- 1
  liststocks <- c(1:nrow(YY))[rowSums(YY) > ctr$minObsPi]
  
  if (length(liststocks) > 1) {
    #cl <- snow::makeCluster(c(rep("localhost", ctr$nCore)), type = "SOCK")
    cl <- parallel::makeCluster(ctr$nCore)
    
    # if (ctr$hac){ snow::clusterEvalQ(cl, require('sandwich'))
    # snow::clusterEvalQ(cl, require('lmtest')) }
    
    liststocks <- liststocks[1:(length(liststocks) - 1)]
    
    #z <- snow::clusterApply(cl = cl, x = as.list(liststocks), fun = alphaScreeningi, 
    #                        rdata = X, factors = factors, T = T, N = N, hac = ctr$hac)
    
    z <- parallel::clusterApplyLB(cl = cl, x = as.list(liststocks), fun = alphaScreeningi,
                                  rdata = X, factors = factors, T = T, N = N, hac = ctr$hac)
    
    #snow::stopCluster(cl)
    parallel::stopCluster(cl)
    
    for (i in 1:length(liststocks)) {
      out <- z[[i]]
      id <- liststocks[i]
      pval[id, id:N] <- pval[id:N, id] <- out[[2]][id:N]
      dalpha[id, id:N] <- out[[1]][id:N]
      dalpha[id:N, id] <- -out[[1]][id:N]
      tstat[id, id:N] <- out[[3]][id:N]
      tstat[id:N, id] <- -out[[3]][id:N]
    }
  }
  
  # pi
  pi <- computePi(pval = pval, dalpha = dalpha, tstat = tstat, lambda = ctr$lambda, 
                  nBoot = ctr$nBoot)
  
  # info on the funds
  info <- infoFund(X, factors = factors)
  
  # form output
  out <- list(n = info$nObs, npeer = colSums(!is.na(pval)), alpha = info$alpha, 
              dalpha = dalpha, pval = pval, tstat = tstat, lambda = pi$lambda, 
              pizero = pi$pizero, pipos = pi$pipos, pineg = pi$pineg)
  
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
#' been parallelized and the computational burden can be splitted across several
#' cores. The number of cores can be defined in \code{control}, see below.
#' 
#' \item (2) for each fund, the false discovery rate approach by Storey (2002)
#' is used to determine the proportions of over, equal, and underperfoming
#' funds, in terms of alpha, in the database.}
#' 
#' The argument \code{control} is a list that can supply any of the following
#' components:
#' \itemize{
#' \item \code{'hac'} Heteroscedastic-autocorrelation consistent
#' standard errors. Default: \code{hac = FALSE}.
#' \item \code{'minObs'} Minimum number of concordant observations to compute the ratios. Default:
#' \code{minObs = 10}.
#' \item \code{'minObsPi'} Minimum number of observations
#' for computing the p-values). Default: \code{minObsPi = 1}.
#' \item \code{'nCore'} Number of cores used to perform the screeing. Default:
#' \code{nCore = 1}.
#' \item \code{'lambda'} Threshold value to compute pi0.
#' Default: \code{lambda = NULL}, i.e. data driven choice.
#' }
#' @param X Matrix \eqn{(T \times N)}{(TxN)} of \eqn{T} returns for the \eqn{N}
#' funds. \code{NA} values are allowed.
#' @param factors Matrix \eqn{(T \times K)}{(TxK)} of \eqn{T} returns for the
#' \eqn{K} factors. \code{NA} values are allowed.
#' @param control Control parameters (see *Details*).
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
#' @note Further details on the methdology with an application to the hedge
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
#' The Peer Ratios Performance of Hedge Funds. 
#' \emph{Journal of Banking and Finance} \bold{87}, pp.351-.368.
#' \doi{10.1016/j.jbankfin.2017.10.014}
#' 
#' Barras, L., Scaillet, O., Wermers, R. (2010).  
#' False discoveries in mutual fund performance: Measuring luck in estimated alphas.  
#' \emph{Journal of Finance} \bold{65}(1), pp.179--216.
#' \doi{10.1111/j.1540-6261.2009.01527.x}
#' 
#' Carhart, M. (1997).  
#' On persistence in mutual fund performance.
#' \emph{Journal of Finance} \bold{52}(1), pp.57--82.
#' \doi{10.1111/j.1540-6261.1997.tb03808.x}
#' 
#' Fama, E., French, K. (2010).  
#' Luck versus skill in the cross-section of mutual fund returns.  
#' \emph{Journal of Finance} \bold{65}(5), pp.1915--1947.
#' \doi{10.1111/j.1540-6261.2010.01598.x}
#' 
#' Fung, W., Hsieh, D. (2004).  
#' Hedge fund benchmarks: A risk based approach.
#' \emph{Financial Analysts Journal} \bold{60}(5), pp.65--80.
#' \doi{10.2469/faj.v60.n5.2657}
#' 
#' Storey, J. (2002).  
#' A direct approach to false discovery rates.
#' \emph{Journal of the Royal Statistical Society B} \bold{64}(3), pp.479--498.
#' \doi{10.1111/1467-9868.00346}
#' 
#' Treynor, J. L., Black, F. (1973).  
#' How to use security analysis to improve portfolio selection.  
#' \emph{Journal of Business} \bold{46}(1), pp.66--86.
#' \doi{10.1002/9781119196679.ch60}
#' @keywords htest
#' @examples
#' ## Load the data (randomized data of monthly hedge fund returns)
#' data("hfdata")
#' rets = hfdata[,1:10]
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
.alphaScreeningi <- function(i, rdata, factors, T, N, hac) {
  pvali <- dalphai <- tstati <- rep(NA, N)
  
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
  
  # #####################################
  # # Correct structure of if/else.  
  # if (!hac) {
  #   
  #   #####################################
  #   # Correct this function. Will be slower (#20190002)
  #   
  #   # if (is.null(factors)) {
  #   #   fit <- stats::lm(dXY ~ 1, na.action = stats::na.omit) 
  #   # } else {
  #   #   fit <- stats::lm(dXY ~ 1 + factors, na.action = stats::na.omit)
  #   # }
  #   # sumfit <- summary(fit)
  #   # if (nPeer == 1) {
  #   #   pvali[N] <- sumfit$coef[1, 4]
  #   #   dalphai[N] <- sumfit$coef[1, 1]
  #   #   tstati[N] <- sumfit$coef[1, 3]
  #   # } else {
  #   #   k <- 1
  #   #   for (j in (i + 1):N) {
  #   #     pvali[j] <- sumfit[[k]]$coef[1, 4]
  #   #     dalphai[j] <- sumfit[[k]]$coef[1, 1]
  #   #     tstati[j] <- sumfit[[k]]$coef[1, 3]
  #   #     k <- k + 1
  #   #   }
  #   # }
  #   
  #   #####################################
  #   # New proposal that's in line with hac = TRUE
  #   # Doesn't use lists. 
  #   
  #   if (nPeer == 1) {
  #     if (is.null(factors)) {
  #       fit <- stats::lm(dXY ~ 1, na.action = stats::na.omit)
  #     } else {
  #       fit <- stats::lm(dXY ~ 1 + factors, na.action = stats::na.omit)
  #     }
  #     sumfit <- summary(fit)
  #     pvali[N] <- sumfit$coef[1, 4]
  #     dalphai[N] <- sumfit$coef[1, 1]
  #     tstati[N] <- sumfit$coef[1, 3]
  #   } else { # end of nPeer == 1
  #     k <- 1
  #     for (j in (i + 1):N) {
  #       if (is.null(factors)) {
  #         fit <- stats::lm(dXY[, k] ~ 1, na.action = stats::na.omit)
  #       } else {
  #         fit <- stats::lm(dXY[, k] ~ 1 + factors, na.action = stats::na.omit)
  #       }
  #       sumfit <- summary(fit)
  #       pvali[j] <- sumfit$coef[1, 4]
  #       dalphai[j] <- sumfit$coef[1, 1]
  #       tstati[j] <- sumfit$coef[1, 3]
  #       k <- k + 1
  #     }
  #   }
  #   
  #   #####################################
  #   #
  # } else { # end of HAC = FALSE
  #   if (nPeer == 1) {
  #     if (is.null(factors)) {
  #       fit <- stats::lm(dXY ~ 1, na.action = stats::na.omit)
  #     } else {
  #       fit <- stats::lm(dXY ~ 1 + factors, na.action = stats::na.omit)
  #     }
  #     sumfit <- lmtest::coeftest(fit, vcov. = sandwich::vcovHAC(fit))
  #     pvali[N] <- sumfit[1, 4]
  #     dalphai[N] <- sumfit[1, 1]
  #     tstati[N] <- sumfit[1, 3]
  #   } else {
  #     k <- 1
  #     for (j in (i + 1):N) {
  #       if (is.null(factors)) {
  #         fit <- stats::lm(dXY[, k] ~ 1, na.action = stats::na.omit)
  #       } else {
  #         fit <- stats::lm(dXY[, k] ~ 1 + factors, na.action = stats::na.omit)
  #       }
  #       sumfit <- lmtest::coeftest(fit, vcov. = sandwich::vcovHAC(fit))
  #       pvali[j] <- sumfit[1, 4]
  #       dalphai[j] <- sumfit[1, 1]
  #       tstati[j] <- sumfit[1, 3]
  #       k <- k + 1
  #     }
  #   }
  # }
  
  #####################################
  # Updated function. 
  
  # Now will not iterate from 1:N; but it will iterate over selIds. 
  # This means that potentially the invertion of the matrix doesn't work? 
  # It does; dXY is properly defined. 
  
  if (nPeer == 1) {
    if (is.null(factors)) {
      fit <- stats::lm(dXY ~ 1, na.action = stats::na.omit)
    } else {
      fit <- stats::lm(dXY ~ 1 + factors, na.action = stats::na.omit)
    } # end of factors/no factors
    
    # HAC within loop.
    if (!hac) {
      sumfit <- summary(fit)
      pvali[N] <- sumfit$coef[1, 4]
      dalphai[N] <- sumfit$coef[1, 1]
      tstati[N] <- sumfit$coef[1, 3]
    } else {
      sumfit <- lmtest::coeftest(fit, vcov. = sandwich::vcovHAC(fit))
      pvali[N] <- sumfit[1, 4]
      dalphai[N] <- sumfit[1, 1]
      tstati[N] <- sumfit[1, 3]
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
        pvali[j] <- sumfit$coef[1, 4]
        dalphai[j] <- sumfit$coef[1, 1]
        tstati[j] <- sumfit$coef[1, 3]
      } else{
        sumfit <- lmtest::coeftest(fit, vcov. = sandwich::vcovHAC(fit))
        pvali[j] <- sumfit[1, 4]
        dalphai[j] <- sumfit[1, 1]
        tstati[j] <- sumfit[1, 3]
      }
      
      # k <- k + 1
    }
  }
  
  out <- list(dalphai = dalphai, pvali = pvali, tstati = tstati)
  return(out)
}
alphaScreeningi <- compiler::cmpfun(.alphaScreeningi)
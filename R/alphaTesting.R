## Set of R functions for Sharpe ratio testing

#@name .alphaTesting
#@description See alphaTesting
# .alphaTesting <- function(x, y, factors=NULL, control=list()){
# 	X <- cbind(x, y)
# 	screen <- alphaScreening(X, factors, control)
# 	out <- list(n = screen$n[[1]], alpha = screen$alpha, dalpha = screen$dalpha[1,2],
# 				tstat = screen$tstat[1,2], pval = screen$pval[1,2])
# 	class(out) <- "TESTING"
# 	return(out)
# }

.alphaTesting <- function(x, y, factors=NULL, control=list(), screen_beta=FALSE){
	# process control
	ctr <- processControl(control)
	hac <- ctr$hac
	X <- as.matrix(x)
	Y <- as.matrix(y)
	dXY <- X - Y
	nObs <- colSums(is.finite(dXY))
	N <- 2


	if (screen_beta & !is.null(factors)) {
		n_coef <- ncol(factors) + 1
		row_return <- 1:n_coef
	} else {
		row_return <- 1
	}

	if (is.null(factors)) {
		fit <- stats::lm(dXY ~ 1, na.action = stats::na.omit)
		fitX <- stats::lm(X ~ 1, na.action = stats::na.omit)
		fitY <- stats::lm(Y ~ 1, na.action = stats::na.omit)
	} else {
		beta <- factors
		fit <- stats::lm(dXY ~ 1 + beta, na.action = stats::na.omit)
		fitX <- stats::lm(X ~ 1 + beta, na.action = stats::na.omit)
		fitY <- stats::lm(Y ~ 1 + beta, na.action = stats::na.omit)
	} # end of factors/no factors

	# HAC within loop.
	if (!hac) {
		sumfit <- summary(fit)
		sumfitX <- summary(fitX)
		sumfitY <- summary(fitY)
		alpha <- cbind(sumfitX$coef[row_return,1], sumfitY$coef[row_return,1])
		pval <- sumfit$coef[row_return, 4]
		dalpha <- sumfit$coef[row_return, 1]
		tstat <- sumfit$coef[row_return, 3]
	} else {
		sumfit <- lmtest::coeftest(fit, vcov. = sandwich::vcovHAC(fit))
		sumfitX <- lmtest::coeftest(fitX, vcov. = sandwich::vcovHAC(fitX))
		sumfitY <- lmtest::coeftest(fitY, vcov. = sandwich::vcovHAC(fitY))
		alpha <- c(sumfitX[row_return,1], sumfitY[row_return,1])
		pval <- sumfit[row_return, 4]
		dalpha <- sumfit[row_return, 1]
		tstat <- sumfit[row_return, 3]
	}

	out <- list(n = nObs,
				alpha = alpha,
				dalpha = dalpha,
				tstat = tstat,
				pval = pval,
				screen_beta = screen_beta)
	class(out) <- "TESTING"
	return(out)
}

#' @name alphaTesting
#' @title Testing the difference of alpha outperformance ratios
#' @description Function which performs the testing of the difference of alpha outperformance ratios.
#' @details The alpha measure (Treynor and Black 1973, Carhart 1997, Fung and Hsieh
#' 2004) is one industry standard for measuring the absolute risk adjusted
#' performance of hedge funds. This function performs
#' the testing of alpha outperformance ratio difference for two funds.
#'
#' For the testing, only the intersection of non-\code{NA} observations for the
#' two funds are used.
#'
##' The argument \code{control} is a list that can supply any of the following
#' components:
#' \itemize{
#' \item \code{'hac'} Heteroscedastic-autocorrelation consistent
#' standard errors. Default: \code{hac = FALSE}.
#' }
#' @param x Vector (of length \eqn{T}) of returns for the first fund. \code{NA}
#' values are allowed.
#' @param y Vector (of length \eqn{T}) returns for the second fund. \code{NA}
#' values are allowed.
#' @param factors Matrix \eqn{(T \times K)}{(TxK)} of \eqn{T} returns for the
#' \eqn{K} factors. \code{NA} values are allowed.
#' @param screen_beta Boolean to screen all factors' coefficients (beta).
#' Default: \code{screen_beta=FALSE} (i.e. only outputs the alpha).
#' If \code{screen_beta=TRUE}, each element of the returned list will have a new first dimension
#' representing each coefficient (the first one being alpha)
#'
#' @param control Control parameters (see *Details*).
#' @return A list with the following components:\cr
#'
#' \code{n}: Number of non-\code{NA} concordant observations.\cr
#'
#' \code{alpha}: Vector (of length 2) of unconditional alpha outperformance ratios.\cr
#'
#' \code{dalpha}: alpha outperformance ratios difference.\cr
#'
#' \code{tstat}: t-stat of alpha outperformance ratios differences.\cr
#'
#' \code{pval}: pvalues of test of alpha outperformance ratios differences.
#' @note Further details on the methodology with an application to the hedge
#' fund industry is given in in Ardia and Boudt (2018).
#'
#' Some internal functions where adapted from Michael Wolf MATLAB code.
#' @author David Ardia and Kris Boudt.
#' @seealso \code{\link{alphaScreening}}..
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
#' Sharpe, W.F. (1994).
#' The Sharpe ratio.
#' \emph{Journal of Portfolio Management} \bold{21}(1), 49--58.
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
#' x = hfdata[,1]
#' y = hfdata[,2]
#'
#' ## Run alpha testing
#' alphaTesting(x, y)

#' @export
#' @import compiler
alphaTesting <- compiler::cmpfun(.alphaTesting)

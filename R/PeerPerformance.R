#' @name PeerPerformance
#' @docType package
#' @title PeerPerformance
#' @description \code{PeerPerformance} is an \R package for the peer-performance evaluation with 
#' luck-correction, useful in the financial industry. In particular, it implements the peer performance ratios of Ardia and Boudt 
#' (2016) which measure the percentage of peers a focal hedge fund outperforms and underperforms, after 
#' correction for luck. In addition, it implements the testing framework for the Sharpe and modified Sharpe ratios, described in 
#' Ledoit and Wolf (2008) and Ardia and Boudt (2015).
#' 
#' @section Functions:
#' \itemize{
#' \item Sharpe ratio: \code{\link{sharpe}}, \code{\link{sharpeTesting}} and \code{\link{sharpeScreening}}\cr
#' \item Modified Share ratio: \code{\link{msharpe}}, \code{\link{msharpeTesting}} and \code{\link{msharpeScreening}}\cr
#' \item Screening function: \code{\link{alphaScreening}}, \code{\link{sharpeScreening}} and \code{\link{msharpeScreening}}.
#' }
#' 
#' @section Update:
#' The latest version of the package is available at \url{https://github.com/ArdiaD/PeerPerformance}
#' @author David Ardia and Kris Boudt.
#' @references 
#' Ardia, D., Boudt, K. (2015).  Testing equality of modified
#' Sharpe ratios \emph{Finance Research Letters} \bold{13}, pp.97--104.
#' 
#' Ardia, D., Boudt, K. (2016).  \emph{The Peer Ratios Performance of Hedge
#' Funds}. \url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2000901}
#' 
#' Barras, L., Scaillet, O., Wermers, R. (2010).  False discoveries in mutual
#' fund performance: Measuring luck in estimated alphas.  \emph{Journal of
#' Finance} \bold{5}, pp.179--216.
#' 
#' Ledoit, O., Wolf, M. (2008). Robust performance hypothesis testing with the
#' Sharpe ratio.  \emph{Journal of Empirical Finance} \bold{15}, pp.850--859.
#' 
#' Storey, J. (2002).  A direct approach to false discovery rates.
#' \emph{Journal of the Royal Statistical Society B} \bold{64}, pp.479--498.
#' @import lmtest
#' @import sandwich
#' @import compiler
NULL
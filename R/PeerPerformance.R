#' @name PeerPerformance
#' @docType package
#' @title : PeerPerformance: Luck-Corrected Peer Performance Analysis in R
#' @description \code{PeerPerformance} is an \R package for the peer-performance evaluation of financial investments with 
#' luck-correction, useful in the financial industry. In particular, it implements the peer performance ratios of Ardia and Boudt 
#' (2016) which measure the percentage of peers a focal (hedge) fund outperforms and underperforms, after 
#' correction for luck. It is useful for fund or portfolio managers to benchmark their investments or screen a universe of new funds. 
#' In addition, the package implements the testing framework for the Sharpe and modified Sharpe ratios, described in 
#' Ledoit and Wolf (2008) and Ardia and Boudt (2015).
#' @section Functions:
#' \itemize{
#' \item Sharpe ratio: \code{\link{sharpe}}, \code{\link{sharpeTesting}} and \code{\link{sharpeScreening}};
#' \item Modified Share ratio: \code{\link{msharpe}}, \code{\link{msharpeTesting}} and \code{\link{msharpeScreening}};
#' \item Screening function: \code{\link{alphaScreening}}, \code{\link{sharpeScreening}} and \code{\link{msharpeScreening}}.
#' }
#' @section Update:
#' The latest version of the package is available at \url{https://github.com/ArdiaD/PeerPerformance}
#' @author David Ardia and Kris Boudt.
#' @note Full description of the methodologies implemented in the various functions is available 
#' in Ledoit and Wolf (2008) and Ardia and Boudt (2015, 2016).
#' @references 
#' Ardia, D., Boudt, K. (2015).  
#' Testing equality of modified Sharpe ratios.
#' \emph{Finance Research Letters} \bold{13}, pp.97--104. 
#' \doi{10.1016/j.frl.2015.02.008}
#' 
#' Ardia, D., Boudt, K. (2016).  
#' The Peer Ratios Performance of Hedge Funds. 
#' \emph{Working paper}.
#' \url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2000901}
#' 
#' Barras, L., Scaillet, O., Wermers, R. (2010).  
#' False discoveries in mutual fund performance: Measuring luck in estimated alphas.  
#' \emph{Journal of Finance} \bold{65}(1), pp.179--216.
#' \doi{10.1111/j.1540-6261.2009.01527.x}
#' 
#' Favre, L., Galeano, J.A. (2002).  
#' Mean-modified Value-at-Risk Optimization with Hedge Funds.  
#' \emph{Journal of Alternative Investments} \bold{5}(2), pp.21--25.
#' \doi{10.3905/jai.2002.319052}
#' 
#' Gregoriou, G. N., Gueyie, J.-P. (2003).  
#' Risk-adjusted performance of funds of hedge funds using a modified Sharpe ratio.  
#' \emph{Journal of Wealth Management} \bold{6}(3), pp.77--83.
#' 
#' Ledoit, O., Wolf, M. (2008). 
#' Robust performance hypothesis testing with the Sharpe ratio.  
#' \emph{Journal of Empirical Finance} \bold{15}(5), pp.850--859.
#' \doi{10.1016/j.jempfin.2008.03.002}
#' 
#' Sharpe, W.F. (1994).  
#' The Sharpe ratio.  
#' \emph{Journal of Portfolio Management} \bold{21}(1), pp.49--58.
#' \doi{10.3905/jpm.1994.409501}
#' 
#' Storey, J. (2002).  
#' A direct approach to false discovery rates.
#' \emph{Journal of the Royal Statistical Society B} \bold{64}(3), pp.479--498.
#' \doi{10.1111/1467-9868.00346}
#' @import lmtest
#' @import sandwich
#' @import compiler
NULL
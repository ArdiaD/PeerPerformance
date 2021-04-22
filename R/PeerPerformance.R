#' @name PeerPerformance
#' @docType package
#' @title PeerPerformance: Luck-corrected peer performance analysis in R
#' @description \code{PeerPerformance} (Ardia and Boudt, 20xx) is an \R package for the peer-performance evaluation of financial investments with 
#' luck-correction, useful in the financial industry. In particular, it implements the peer performance ratios of Ardia and Boudt 
#' (2018) which measure the percentage of peers a focal (hedge) fund outperforms and underperforms, after 
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
#' @note By using \code{PeerPerformance} you agree to the following rules: (1) You must cite Ardia and Boudt (2018) in 
#' working papers and published papers that use \code{PeerPerformance} (use \code{citation("PeerPerformance")}), (2) you 
#' must place the URL \url{https://CRAN.R-project.org/package=PeerPerformance} in a footnote to help 
#' others find \code{PeerPerformance}, and (3) you assume all risk for the use of \code{PeerPerformance}.
#' @note Full description of the methodologies implemented in the various functions is available 
#' in Ledoit and Wolf (2008) and Ardia and Boudt (2015, 2018).
#' @references 
#' Ardia, D., Boudt, K. (2015).  
#' Testing equality of modified Sharpe ratios.
#' \emph{Finance Research Letters} \bold{13}, pp.97--104. 
#' \doi{10.1016/j.frl.2015.02.008}
#' 
#' Ardia, D., Boudt, K. (2018).  
#' The peer performance ratios of hedge Funds. 
#' \emph{Journal of Banking and Finance} \bold{87}, pp.351-.368.
#' \doi{10.1016/j.jbankfin.2017.10.014}
#' 
#' Ardia, D., Boudt, K. (20xx).    
#' \emph{PeerPerformance: Luck-corrected peer performance analysis in R}.    
#' R package.   
#' \url{https://github.com/ArdiaD/PeerPerformance}
#' 
#' Barras, L., Scaillet, O., Wermers, R. (2010).  
#' False discoveries in mutual fund performance: Measuring luck in estimated alphas.  
#' \emph{Journal of Finance} \bold{65}(1), pp.179--216.
#' 
#' Favre, L., Galeano, J.A. (2002).  
#' Mean-modified Value-at-Risk Optimization with Hedge Funds.  
#' \emph{Journal of Alternative Investments} \bold{5}(2), pp.21--25.
#' 
#' Gregoriou, G. N., Gueyie, J.-P. (2003).  
#' Risk-adjusted performance of funds of hedge funds using a modified Sharpe ratio.  
#' \emph{Journal of Wealth Management} \bold{6}(3), pp.77--83.
#' 
#' Ledoit, O., Wolf, M. (2008). 
#' Robust performance hypothesis testing with the Sharpe ratio.  
#' \emph{Journal of Empirical Finance} \bold{15}(5), pp.850--859.
#' 
#' Sharpe, W.F. (1994).  
#' The Sharpe ratio.  
#' \emph{Journal of Portfolio Management} \bold{21}(1), pp.49--58.
#' 
#' Storey, J. (2002).  
#' A direct approach to false discovery rates.
#' \emph{Journal of the Royal Statistical Society B} \bold{64}(3), pp.479--498.
#' @import lmtest
#' @import sandwich
#' @import compiler
NULL
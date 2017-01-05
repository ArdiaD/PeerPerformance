\name{alphaScreening}
\alias{alphaScreening}
\concept{Screeing}
\concept{Outperformance ratio}
\concept{False discovery rate}
\concept{Alpha}
\title{Screening using the Alpha outperformance ratio}
\description{%
  Function which performs the screening of a universe of returns, and compute the Alpha outperformance ratio}
\usage{%
  alphaScreening(X, factors = NULL, control = list())
}
\arguments{%
  \item{X}{matrix \eqn{(T \times N)}{(TxN)} of \eqn{T} returns for the \eqn{N} funds. \code{NA} values are allowed.}
  \item{factors}{matrix \eqn{(T \times K)}{(TxK)} of \eqn{T} returns for the \eqn{K} factors. \code{NA} values are allowed.}
  \item{control}{control parameters (see *Details*).}
}
\details{%
  The alpha measure (Treynor and Black 1973, Carhart 1997, Fung and Hsieh 2004) is one industry standard for measuring 
  the absolute risk adjusted performance of hedge funds. We propose to complement the alpha measure with the fund's alpha outperformance
  ratio, defined as the percentage number of funds that have a significantly lower alpha. In a pairwise testing framework, a fund can have a significantly
  higher alpha because of luck. We correct for this by applying the false discovery rate approach by (Storey 2002).

  The methodology proceeds as follows:
  \describe{%
    \item{(1)}{compute all pairwise tests of alpha differences. This means
       that for a universe of \eqn{N} funds, we perform \eqn{N(N-1)/2}{N*(N-1)/2} tests. The algorithm has been parallelized and the computational
       burden can be slip across several cores. The number of cores can be defined in \code{control}, see below.}
    \item{(2)}{for each fund, the false discovery rate approach by Storey (2002) is used to determine the proportions of over, equal, and
       underperfoming funds, in terms of alpha, in the database.}
  }
  
  The argument \code{control} is a list that can supply any of the following components:
  
  \describe{%
    \item{\code{hac}}{heteroscedastic-autocorrelation consistent standard errors. Default: \code{hac = FALSE}.}
    \item{\code{minObs}}{minimum number of concordant observations to compute the ratios. Default: \code{minObs = 10}.}
    \item{\code{minObsPi}}{minimum number of observations for computing the p-values). Default: \code{minObsPi = 1}.}
    \item{\code{nCore}}{number of cores used to perform the screeing. Default: \code{nCore = 1}.}
    \item{\code{lambda}}{threshold value to compute pi0. Default: \code{lambda = NULL}, i.e. data driven choice.}
  }
}
\value{%
  A list with the following components:\cr

  \code{n}: vector (of length \eqn{N}) of number of non-\code{NA} observations.\cr

  \code{npeer}: vector (of length \eqn{N}) of number of available peers.\cr

  \code{alpha}: vector (of length \eqn{N}) of unconditional alpha.\cr
  
  \code{dalpha}: matrix (of size \eqn{N \times N}{NxN}) of alpha differences.\cr
  
  \code{pval}: matrix (of size \eqn{N \times N}) of p-values of test for alpha differences.\cr

  \code{lambda}: vector (of length \eqn{N}) of lambda values.\cr

  \code{pizero}: vector (of length \eqn{N}) of probability of equal performance.\cr

  \code{pipos}: vector (of length \eqn{N}) of probability of outperformance performance.\cr

  \code{pineg}: vector (of length \eqn{N}) of probability of underperformance performance.
}
\note{%
  Further details on the methdology with an application to the hedge fund industry is given in 
  in Ardia and Boudt (2016). The file \file{ThePeerPerformanceOfHedgeFunds.txt} in the \file{/doc} package's folder 
  allows the reprodution of the steps followed in the article. See also the presentation by Kris Boudt
  at the R/Finance conference 2012 at \url{http://www.rinfinance.com}.

  Application of the false discovery rate approach applied to the mutual fund industry has been presented in 
  Barras, Scaillet and Wermers (2010). 

  Currently, the hac asymptotic and studentized circular block bootstrap presented in Ledoit and Wolf (2008) are not supported by
  the \code{alphaScreening} function.

  Please cite the package in publications. Use \code{citation("PeerPerformance")}.
}
\author{%
  David Ardia and Kris Boudt.
}
\references{%
   Ardia, D., Boudt, K. (2015).
  Testing equality of modified Sharpe ratios
  \emph{Finance Research Letters} \bold{13}, pp.97--104.

  Ardia, D., Boudt, K. (2016).
  \emph{The Peer Ratios Performance of Hedge Funds}.
  \url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2000901}

  Barras, L., Scaillet, O., Wermers, R. (2010). 
  False discoveries in mutual fund performance: Measuring luck in estimated alphas. 
  \emph{Journal of Finance} \bold{5}, pp.179--216.

  Carhart, M. (1997). 
  On persistence in mutual fund performance. 
  \emph{Journal of Finance} \bold{52}, pp.57--82.
  
  Fama, E., French, K. (2010). 
  Luck versus skill in the cross-section of mutual fund returns. 
  \emph{Journal of Finance} \bold{65}, pp.1915--1947.

  Fung, W., Hsieh, D. (2004). 
  Hedge fund benchmarks: A risk based approach. 
  \emph{Financial Analysts Journal} \bold{60}, pp.65--80.

  Storey, J. (2002). 
  A direct approach to false discovery rates. 
  \emph{Journal of the Royal Statistical Society B} \bold{64}, pp.479--498.
  
  Treynor, J. L., Black, F. (1973). 
  How to use security analysis to improve portfolio selection. 
  \emph{Journal of Business 46} \bold{1}, pp.66--86.
}
\seealso{%
  \code{\link{sharpeScreening}} and \code{\link{msharpeScreening}}.
}
\keyword{htest}

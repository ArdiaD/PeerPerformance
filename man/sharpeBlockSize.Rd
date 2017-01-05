\name{sharpeBlockSize}
\alias{sharpeBlockSize}
\concept{Testing}
\concept{Sharpe}
\title{Optimal block length for bootstrap test of difference of Sharpe ratios}
\description{%
  Function which computes the optimal block length for testing of the difference of Sharpe ratios}
\usage{%
  sharpeBlockSize(x, y, control = list(),
                  b.vec = c(1, 3, 6, 10), alpha = 0.05, M = 199, K = 500, b.av = 5, T.start = 50)
}
\arguments{%
  \item{x}{vector (of lenght \eqn{T}) of returns for the first fund. \code{NA} values are allowed.}
  \item{y}{vector (of lenght \eqn{T}) returns for the second fund. \code{NA} values are allowed.}
  \item{control}{control parameters (see *Details*).}
  \item{b.vec}{vector of block to be tested.}
  \item{alpha}{significance level.}
  \item{M}{bootstrap replications.}
  \item{K}{number of cross-validation.}
  \item{b.av}{average block length in the stationary bootstrap.}
  \item{T.start}{starting point for bootstrap}
}
\details{%
  The argument \code{control} is a list that can supply any of the following components:
  
  \describe{%
    \item{\code{type}}{asymptotic approach (\code{type = 1}) or studentized circular bootstrap approach (\code{type = 2}). Default: \code{type = 1}.}
    \item{\code{ttype}}{test based on ratio (\code{type = 1}) or product (\code{type = 2}). Default: \code{type = 2}.}
    \item{\code{hac}}{heteroscedastic-autocorrelation consistent standard errors. Default: \code{hac = FALSE}.}
    \item{\code{minObs}}{minimum number of concordant observations to compute the ratios. Default: \code{minObs = 10}.}
    \item{\code{nBoot}}{number of boostrap replications for computing the p-value. Default: \code{nBoot = 499}.}
    \item{\code{bBoot}}{block length in the circular bootstrap. Default: \code{bBoot = 1}, i.e. iid bootstrap. \code{bBoot = 0} uses optimal block-length.}
    \item{\code{pBoot}}{symmetric p-value (\code{pBoot = 1}) or asymmetric p-value (\code{pBoot = 2}). Default: \code{pBoot = 1}.}
  }
}
\value{%
  The optimal block length.
}
\note{%
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
  \emph{The Peer Performance Ratios of Hedge Funds}.
  \url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2000901}

  Ledoit, O., Wolf, M. (2008). Robust performance hypothesis testing with the Sharpe ratio.
  \emph{Journal of Empirical Finance} \bold{15}, pp.850--859.

  Sharpe, W. F. (1994). 
  The Sharpe ratio. 
  \emph{Journal of Portfolio Management} Fall, pp.49--58.
}
\seealso{%
  \code{\link{sharpe}}, \code{\link{sharpeScreening}} and \code{\link{msharpeTesting}}.
}
\keyword{htest}
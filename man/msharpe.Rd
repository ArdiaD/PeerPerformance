\name{msharpe}
\alias{msharpe}
\concept{Estimation}
\concept{Sharpe}
\concept{Modified Sharpe}
\title{Compute Modified Sharpe ratio}
\description{%
  Function which computes the modified Sharpe ratio}
\usage{%
  msharpe(X, level = 0.90, na.rm = TRUE, na.neg = TRUE)
}
\arguments{%
  \item{X}{vector (of lenght \eqn{T}) or matrix (of size \eqn{T \times N}{TxN}) of returns. 
  \code{NA} values are allowed.}
  \item{level}{modified Value-at-Risk level. Default: \code{level = 0.90}.}
  \item{na.rm}{a logical value indicating whether \code{NA} values should be stripped before the computation. Default \code{na.rm = TRUE}.}
  \item{na.neg}{a logical value indicating whether \code{NA} values should be returned if a negative modified Value-at-Risk is obtained.
  Default \code{na.neg = TRUE}.}
}
\details{%
  The modified Sharpe ratio (Favre and Galeano 2002) is one industry standard for measuring the absolute risk adjusted performance
  of hedge funds.
}
\value{%
  A scalar or a vector (of size \eqn{N}) with the modified Sharpe ratios.
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

  Ardia, D., Boudt, K. (2015).
  \emph{The Peer Performance of Hedge Funds}.
  \url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2000901}

  Favre, L., Galeano, J.A. (2002).
  Mean-modified Value-at-Risk with Hedge Funds.
  \emph{The Journal of Alternative Investments} \bold{5}, pp.21--25.

  Gregoriou, G. N., Gueyie, J.-P. (2003).
  Risk-adjusted performance of funds of hedge funds using a modified Sharpe ratio.
  \emph{The Journal of Wealth Management} \bold{Winter}, pp.77--83.
}
\seealso{%
  \code{\link{msharpeTesting}}, \code{\link{msharpeScreening}} and \code{\link{sharpe}}.
}
\examples{%
  ## Load the data (randomized data of monthly hedge fund returns)
  data("hfdata")

  out = msharpe(hfdata)
  print(out)

  out = msharpe(hfdata, na.rm = FALSE)
  print(out)
}
\keyword{htest}
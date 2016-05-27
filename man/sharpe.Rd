\name{sharpe}
\alias{sharpe}
\concept{Estimation}
\concept{Sharpe}
\title{Compute Sharpe ratio}
\description{%
  Function which computes the Sharpe ratio}
\usage{%
  sharpe(X, na.rm = TRUE)
}
\arguments{%
  \item{X}{vector (of lenght \eqn{T}) or matrix (of size \eqn{T \times N}{TxN}) of returns 
  for \eqn{N} funds. \code{NA} values are allowed.}
  \item{na.rm}{a logical value indicating whether \code{NA} values should be stripped before the computation. Default \code{na.rm = TRUE}}
}
\details{%
  The Sharpe ratio (Sharpe 1992) is one industry standard for measuring the absolute risk adjusted performance of hedge funds.
}
\value{%
  A scalar or a vector (of size \eqn{N}) with the Sharpe ratios.
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

  Sharpe, W. F. (1994). 
  The Sharpe ratio. 
  \emph{Journal of Portfolio Management} Fall, pp.49--58.
}
\seealso{%
  \code{\link{sharpeTesting}}, \code{\link{sharpeScreening}} and \code{\link{msharpe}}.
}
\examples{%
  ## Load the data (randomized data of monthly hedge fund returns)
  data("hfdata")

  out = sharpe(hfdata)
  print(out)

  out = sharpe(hfdata, na.rm = FALSE)
  print(out)
}
\keyword{htest}
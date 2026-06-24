## Targeted peer-performance screening for a selected subset of funds.
## This is a thin convenience wrapper over the cross-group screening
## (`*Screening(X[, funds], Y = X)`): each focal fund is screened against the
## whole universe `X` (its own column is excluded automatically). Originally
## contributed by Murilo Andre Peres Pereira.

.targetPeerPerformance <- function(X, funds,
                                   method = c("alpha", "sharpe", "msharpe"),
                                   factors = NULL, level = 0.9, na.neg = TRUE,
                                   control = list()) {
  method <- match.arg(method)
  X <- as.matrix(X)
  N <- ncol(X)
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("Fund ", seq_len(N))
  }
  fund_names <- colnames(X)

  # resolve 'funds' to column indices (preserving the user-supplied order)
  idx <- if (is.numeric(funds)) {
    if (anyNA(funds) || any(!is.finite(funds)) || any(funds != as.integer(funds))) {
      stop("numeric 'funds' must be whole, finite, positive indices")
    }
    as.integer(funds)
  } else {
    m <- match(funds, fund_names)
    if (anyNA(m)) stop("unknown fund(s): ", paste(funds[is.na(m)], collapse = ", "))
    m
  }
  if (any(idx < 1L | idx > N)) stop("some 'funds' indices are out of range")
  if (anyDuplicated(idx)) stop("duplicate entries in 'funds' are not allowed")
  if (method != "alpha" && !is.null(factors)) {
    stop("'factors' is only used when method = 'alpha'")
  }

  Xf <- X[, idx, drop = FALSE]

  # each focal fund vs the whole universe (self column excluded by the
  # cross-group machinery)
  out <- switch(method,
                alpha   = alphaScreening(Xf, factors = factors, control = control, Y = X),
                sharpe  = sharpeScreening(Xf, control = control, Y = X),
                msharpe = msharpeScreening(Xf, level = level, na.neg = na.neg,
                                           control = control, Y = X))

  # label the focal funds (rows / vectors) and the peers (columns)
  fn <- fund_names[idx]
  dname <- switch(method, alpha = "dalpha", sharpe = "dsharpe", msharpe = "dmsharpe")
  for (el in c("n", "npeer", "lambda", "pizero", "pipos", "pineg", method)) {
    v <- out[[el]]
    if (is.null(v)) next
    if (is.matrix(v)) {
      # screen_beta: coefficient-by-focal-fund matrices
      if (ncol(v) == length(fn)) colnames(out[[el]]) <- fn
    } else if (length(v) == length(fn)) {
      names(out[[el]]) <- fn
    }
  }
  for (el in c("pval", "tstat", dname)) {
    if (is.matrix(out[[el]])) {
      rownames(out[[el]]) <- fn
      colnames(out[[el]]) <- fund_names
    }
  }
  out
}

#' @name targetPeerPerformance
#' @title Targeted peer-performance screening for selected funds
#' @description Performs peer-performance screening for a user-selected subset of
#' funds (the \dQuote{focal} funds) against the whole universe, returning
#' detailed outputs for those focal funds only. It is a convenience wrapper that
#' calls the corresponding screening function with the cross-group argument
#' \code{Y} set to the full matrix \code{X}; equivalently,
#' \code{targetPeerPerformance(X, funds, "alpha")} is
#' \code{alphaScreening(X[, funds], Y = X)}.
#' @param X Matrix \eqn{(T \times N)}{(TxN)} of \eqn{T} returns for \eqn{N}
#' funds. \code{NA} values are allowed.
#' @param funds Integer indices or character column names identifying the focal
#' funds. The output preserves the order provided in \code{funds}.
#' @param method Screening method: \code{"alpha"} (default), \code{"sharpe"}, or
#' \code{"msharpe"}.
#' @param factors Optional matrix of factor returns (only when
#' \code{method = "alpha"}).
#' @param level Modified Value-at-Risk level (only when \code{method = "msharpe"}).
#' @param na.neg Logical; return \code{NA} for a negative modified
#' Value-at-Risk (only when \code{method = "msharpe"}).
#' @param control Control parameters (see \code{\link{alphaScreening}}).
#' @return A \code{SCREENING} object (see \code{\link{alphaScreening}}) whose
#' fund-level components (\code{pizero}, \code{pipos}, \code{pineg}, the
#' performance measure, \code{npeer}, \ldots) cover the focal funds only, and
#' whose \code{pval}/\code{tstat}/difference matrices have one row per focal fund
#' and one column per universe fund. Works with the \code{print}, \code{plot},
#' \code{summary}, and \code{as.data.frame} methods.
#' @author Murilo Andre Peres Pereira, David Ardia and Kris Boudt.
#' @seealso \code{\link{alphaScreening}}, \code{\link{sharpeScreening}},
#' \code{\link{msharpeScreening}}.
#' @keywords htest
#' @examples
#' data("hfdata")
#' rets <- hfdata[, 1:10]
#'
#' ## focal Sharpe screening (funds by index) against the whole universe
#' out <- targetPeerPerformance(rets, funds = c(2, 5, 7), method = "sharpe",
#'                              control = list(nCore = 1))
#' out$sharpe
#' @export
#' @importFrom compiler cmpfun
targetPeerPerformance <- compiler::cmpfun(.targetPeerPerformance)

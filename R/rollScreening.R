## Rolling-window peer performance screening (time series of ratios).

## internal: cross-sectional aggregation of a SCREENING result
.rollAgg <- function(sc) {
  pz <- sc$pizero; pp <- sc$pipos; pn <- sc$pineg
  if (is.matrix(pz)) {
    cn <- rownames(pz)
    if (is.null(cn)) cn <- paste0("coef", seq_len(nrow(pz)))
    data.frame(coefficient = cn,
               pizero = as.vector(rowMeans(pz, na.rm = TRUE)),
               pipos  = as.vector(rowMeans(pp, na.rm = TRUE)),
               pineg  = as.vector(rowMeans(pn, na.rm = TRUE)),
               stringsAsFactors = FALSE)
  } else {
    data.frame(coefficient = NA_character_,
               pizero = mean(pz, na.rm = TRUE),
               pipos  = mean(pp, na.rm = TRUE),
               pineg  = mean(pn, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }
}

#' @name rollScreening
#' @title Rolling-window peer performance screening
#' @description Runs a peer performance screening on a rolling time window and
#' returns the time series of cross-sectionally averaged peer performance
#' ratios. This is the design underlying the dynamic analyses of Ardia et al.
#' (2022, 2023): on each window the screening is computed and the ratios are
#' averaged across funds (per coefficient when \code{screen_beta = TRUE}),
#' yielding a series of (factor exposure) heterogeneity measures.
#' @param X Matrix \eqn{(T \times N)}{(TxN)} of \eqn{T} returns for the \eqn{N}
#' funds.
#' @param factors Matrix \eqn{(T \times K)}{(TxK)} of factor returns (for
#' \code{screen = "alpha"}). Default: \code{NULL}.
#' @param Y Optional matrix \eqn{(T \times M)}{(TxM)} of returns for a peer group
#' (cross-group screening, see \code{\link{alphaScreening}}). Default:
#' \code{NULL}.
#' @param screen Performance measure to screen on: \code{"alpha"} (default),
#' \code{"sharpe"}, or \code{"msharpe"}.
#' @param width Window length (number of observations). Default: \code{36}.
#' @param by Step (in observations) between successive windows. Default:
#' \code{1}.
#' @param level Modified Value-at-Risk level (for \code{screen = "msharpe"}).
#' Default: \code{0.9}.
#' @param na.neg A logical value passed to \code{\link{msharpeScreening}}
#' (for \code{screen = "msharpe"}) indicating whether a negative modified
#' Value-at-Risk yields \code{NA}. Default: \code{TRUE}.
#' @param control Control parameters passed to the screening function (see
#' \code{\link{alphaScreening}}); set \code{control = list(screen_beta = TRUE)}
#' to obtain per-factor heterogeneity series. With \code{nCore > 1} each
#' window spins up its own PSOCK cluster, which is wasteful for many small
#' windows; the default \code{nCore = 1} runs serially without any cluster.
#' @param dates Optional vector of length \eqn{T} (e.g. \code{Date}) used to
#' label the windows by their end date. Default: \code{NULL}.
#' @return A \code{data.frame} of class \code{rollScreening} with one row per
#' window (or per window/coefficient when \code{screen_beta = TRUE}), with
#' columns \code{window}, \code{index} (end-of-window row), optionally
#' \code{date}, \code{coefficient}, the averaged \code{pizero}/\code{pipos}/
#' \code{pineg}, and \code{heterogeneity} (\eqn{1-\bar\pi^0}).
#' @author David Ardia and Kris Boudt.
#' @references
#' Ardia, D., Bluteau, K., Tran, D. (2022).
#' How easy is it for investment managers to deploy their talent in green and brown stocks?
#' \emph{Finance Research Letters} \bold{48}, 102992.
#'
#' Ardia, D., Bluteau, K., Lortie-Cloutier, G., Tran, D. (2023).
#' Factor exposure heterogeneity in green and brown stocks.
#' \emph{Finance Research Letters} \bold{55}, Part A, 103900.
#' @seealso \code{\link{alphaScreening}}, \code{\link{exposureHeterogeneity}}.
#' @keywords htest
#' @examples
#' \donttest{
#' data("hfdata")
#' set.seed(1234)
#' roll <- rollScreening(hfdata, screen = "alpha", width = 36, by = 6,
#'                       control = list(nCore = 1))
#' plot(roll)
#' }
#' @export
rollScreening <- function(X, factors = NULL, Y = NULL,
                          screen = c("alpha", "sharpe", "msharpe"),
                          width = 36L, by = 1L, level = 0.9, na.neg = TRUE,
                          control = list(), dates = NULL) {
  screen <- match.arg(screen)
  X <- as.matrix(X)
  T <- nrow(X)
  if (width < 2L || width > T) {
    stop("'width' must be between 2 and nrow(X)")
  }
  if (length(by) != 1L || !is.finite(by) || by < 1L) {
    stop("'by' must be a single positive integer")
  }
  by <- as.integer(by)
  if (!is.null(dates) && length(dates) != T) {
    stop("'dates' must have length nrow(X)")
  }
  ctr <- processControl(control)
  # screen_beta only applies to the alpha screening with factors
  screen_beta <- isTRUE(ctr$screen_beta) && screen == "alpha" && !is.null(factors)
  Fmat <- if (!is.null(factors)) as.matrix(factors) else NULL
  Ymat <- if (!is.null(Y)) as.matrix(Y) else NULL

  starts <- seq.int(1L, T - width + 1L, by = by)
  frags <- vector("list", length(starts))
  for (k in seq_along(starts)) {
    idx <- starts[k]:(starts[k] + width - 1L)
    Xw <- X[idx, , drop = FALSE]
    Fw <- if (!is.null(Fmat)) Fmat[idx, , drop = FALSE] else NULL
    Yw <- if (!is.null(Ymat)) Ymat[idx, , drop = FALSE] else NULL
    sc <- switch(screen,
                 alpha   = alphaScreening(Xw, factors = Fw, control = control, Y = Yw),
                 sharpe  = sharpeScreening(Xw, control = control, Y = Yw),
                 msharpe = msharpeScreening(Xw, level = level, na.neg = na.neg,
                                            control = control, Y = Yw))
    frag <- .rollAgg(sc)
    frag <- cbind(window = k, index = idx[length(idx)], frag)
    frags[[k]] <- frag
  }
  out <- do.call(rbind, frags)
  if (!is.null(dates)) {
    out <- cbind(out[, 1:2, drop = FALSE], date = dates[out$index],
                 out[, -(1:2), drop = FALSE])
  }
  out$heterogeneity <- 1 - out$pizero
  rownames(out) <- NULL
  class(out) <- c("rollScreening", "data.frame")
  attr(out, "screen") <- screen
  attr(out, "screen_beta") <- screen_beta
  out
}

#' @name plot.rollScreening
#' @title Plot a rolling peer performance screening
#' @description Time-series plot of the rolling peer performance ratios produced
#' by \code{\link{rollScreening}}.
#' @param x A \code{rollScreening} object.
#' @param what Either \code{"heterogeneity"} (\eqn{1-\bar\pi^0}) or
#' \code{"ratios"} (the averaged \eqn{\bar\pi^+,\bar\pi^0,\bar\pi^-}). Default:
#' \code{"heterogeneity"} for a beta screening, \code{"ratios"} otherwise.
#' @param legend A logical value indicating whether to draw a legend. Default:
#' \code{TRUE}.
#' @param ... Further graphical arguments.
#' @return Invisibly returns \code{x}.
#' @author David Ardia and Kris Boudt.
#' @seealso \code{\link{rollScreening}}.
#' @keywords hplot
#' @examples
#' \donttest{
#'   data("hfdata")
#'   roll <- rollScreening(hfdata[, 1:15], width = 36, by = 12,
#'                         control = list(nCore = 1))
#'   plot(roll)
#' }
#' @export
#' @importFrom graphics matplot legend lines
#' @importFrom grDevices hcl.colors
plot.rollScreening <- function(x, what = NULL, legend = TRUE, ...) {
  beta <- isTRUE(attr(x, "screen_beta"))
  if (is.null(what)) {
    what <- if (beta) "heterogeneity" else "ratios"
  }
  hasdate <- !is.null(x$date)

  if (what == "ratios" && !beta) {
    ord <- order(x$index)
    xs <- (if (hasdate) x$date else x$index)[ord]
    graphics::matplot(xs, cbind(x$pipos[ord], x$pizero[ord], x$pineg[ord]),
                      type = "l", lty = 1, lwd = 2,
                      col = c("black", "grey60", "grey40"),
                      xlab = "", ylab = "average ratio", ylim = c(0, 1), ...)
    if (legend) {
      graphics::legend("topright", bty = "n", lwd = 2,
                       col = c("black", "grey60", "grey40"),
                       legend = c(expression(bar(pi)^"+"),
                                  expression(bar(pi)^0),
                                  expression(bar(pi)^"-")))
    }
  } else {
    coefs <- if (beta) unique(x$coefficient) else NA
    cols <- grDevices::hcl.colors(max(length(coefs), 1), "Dark3")
    ylab <- if (what == "heterogeneity") expression(1 - bar(pi)^0) else what
    first <- TRUE
    for (i in seq_along(coefs)) {
      sub <- if (beta) x[x$coefficient == coefs[i], , drop = FALSE] else x
      ord <- order(sub$index)
      xs <- (if (hasdate) sub$date else sub$index)[ord]
      ys <- sub[[what]][ord]
      if (first) {
        graphics::plot(xs, ys, type = "l", lwd = 2, col = cols[i],
                       ylim = c(0, 1), xlab = "", ylab = ylab, ...)
        first <- FALSE
      } else {
        graphics::lines(xs, ys, lwd = 2, col = cols[i])
      }
    }
    if (legend && beta) {
      graphics::legend("topright", bty = "n", lwd = 2, col = cols, legend = coefs)
    }
  }
  invisible(x)
}

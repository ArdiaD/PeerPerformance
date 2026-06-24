## S3 methods for the 'TESTING' and 'SCREENING' objects

## internal helper: identify the performance measure carried by a SCREENING
## object and return its label and the (possibly matrix) values
.screeningPerf <- function(x) {
  if (!is.null(x$alpha)) {
    return(list(value = x$alpha, label = "alpha"))
  } else if (!is.null(x$sharpe)) {
    return(list(value = x$sharpe, label = "Sharpe ratio"))
  } else if (!is.null(x$msharpe)) {
    return(list(value = x$msharpe, label = "modified Sharpe ratio"))
  }
  stop("unknown performance measure in the 'SCREENING' object")
}

## internal: single focal fund (one stacked bar of pi+/pi0/pi-)
.plotScreeningSingle <- function(pipos, pizero, pineg, ny, colorset, label, ...) {
  mat <- matrix(100 * c(pipos, pizero, pineg), ncol = 1)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::par(mar = c(4, 1, 4, 1))
  main <- if (is.null(ny)) {
    paste0("Peer performance (", label, ")")
  } else {
    paste0("Focal fund vs peer group of ", ny, " (", label, ")")
  }
  graphics::barplot(mat, horiz = TRUE, xlim = c(0, 100), col = colorset,
                    border = NA, names.arg = "", axes = FALSE, main = main, ...)
  graphics::axis(side = 1, at = seq(0, 100, by = 25),
                 labels = paste0(seq(0, 100, by = 25), "%"))
  graphics::box()
  graphics::legend("top", horiz = TRUE, bty = "n", inset = c(0, -0.16),
                   xpd = TRUE, fill = colorset,
                   legend = c(expression(hat(pi)^"+"), expression(hat(pi)^0),
                              expression(hat(pi)^"-")))
  invisible(t(mat / 100))
}

#' @name print.TESTING
#' @title Print method for the 'TESTING' object
#' @description Print method for the \code{TESTING} object returned by
#' \code{\link{sharpeTesting}}, \code{\link{msharpeTesting}} and
#' \code{\link{alphaTesting}}.
#' @param x A \code{TESTING} object.
#' @param ... Further arguments passed to or from other methods.
#' @return Invisibly returns \code{x}.
#' @author David Ardia and Kris Boudt.
#' @seealso \code{\link{sharpeTesting}}, \code{\link{msharpeTesting}} and
#' \code{\link{alphaTesting}}.
#' @keywords htest
#' @export
print.TESTING <- function(x, ...) {
  if (!is.null(x$sharpe)) {
    label <- "Sharpe ratio"; vals <- x$sharpe; dval <- x$dsharpe
  } else if (!is.null(x$msharpe)) {
    label <- "modified Sharpe ratio"; vals <- x$msharpe; dval <- x$dmsharpe
  } else {
    label <- "alpha"; vals <- x$alpha; dval <- x$dalpha
  }
  vals <- as.vector(vals); dval <- as.vector(dval)[1]
  cat("\nPeer performance test (", label, ")\n", sep = "")
  cat("  Number of (concordant) observations: ", x$n[1], "\n", sep = "")
  cat("  Measure (fund x, fund y) : ",
      formatC(vals[1], format = "f", digits = 4), ", ",
      formatC(vals[2], format = "f", digits = 4), "\n", sep = "")
  cat("  Difference (x - y)       : ",
      formatC(dval, format = "f", digits = 4), "\n", sep = "")
  cat("  t-statistic              : ",
      formatC(as.vector(x$tstat)[1], format = "f", digits = 4), "\n", sep = "")
  cat("  p-value                  : ",
      formatC(as.vector(x$pval)[1], format = "f", digits = 4), "\n", sep = "")
  invisible(x)
}

#' @name print.SCREENING
#' @title Print method for the 'SCREENING' object
#' @description Print method for the \code{SCREENING} object returned by
#' \code{\link{alphaScreening}}, \code{\link{sharpeScreening}} and
#' \code{\link{msharpeScreening}}.
#' @param x A \code{SCREENING} object.
#' @param ... Further arguments passed to or from other methods.
#' @return Invisibly returns \code{x}.
#' @author David Ardia and Kris Boudt.
#' @seealso \code{\link{alphaScreening}}, \code{\link{sharpeScreening}},
#' \code{\link{msharpeScreening}} and \code{\link{plot.SCREENING}}.
#' @keywords htest
#' @export
print.SCREENING <- function(x, ...) {
  perf <- .screeningPerf(x)
  pipos <- x$pipos; pizero <- x$pizero; pineg <- x$pineg
  if (is.matrix(pipos)) {
    pipos <- pipos[1, ]; pizero <- pizero[1, ]; pineg <- pineg[1, ]
  }
  N <- length(pipos)
  cat("\nPeer performance screening (", perf$label, ")\n", sep = "")
  cat("  Number of funds : ", N, "\n", sep = "")
  cat("  Average ratios  : ",
      "pi+ = ", formatC(100 * mean(pipos,  na.rm = TRUE), format = "f", digits = 1), "%",
      "  pi0 = ", formatC(100 * mean(pizero, na.rm = TRUE), format = "f", digits = 1), "%",
      "  pi- = ", formatC(100 * mean(pineg,  na.rm = TRUE), format = "f", digits = 1), "%",
      "\n", sep = "")
  cat("  Use plot() to display the peer performance screening plot.\n")
  invisible(x)
}

#' @name as.data.frame.SCREENING
#' @title Coerce a 'SCREENING' object to a data frame
#' @description Returns a tidy \code{data.frame} with one row per fund (or, when
#' the screening was run with \code{screen_beta = TRUE}, one row per
#' fund/coefficient), suitable for downstream filtering, ranking, and plotting.
#' @param x A \code{SCREENING} object.
#' @param row.names,optional Currently ignored.
#' @param ... Further arguments (ignored).
#' @return A \code{data.frame} with columns \code{fund}, the performance measure
#' (\code{alpha}, \code{sharpe}, or \code{msharpe}), \code{pipos}, \code{pizero},
#' \code{pineg}, \code{npeer}, and \code{lambda}; a \code{coefficient} column is
#' added when \code{screen_beta = TRUE}.
#' @author David Ardia and Kris Boudt.
#' @seealso \code{\link{alphaScreening}}.
#' @keywords manip
#' @export
as.data.frame.SCREENING <- function(x, row.names = NULL, optional = FALSE, ...) {
  short <- if (!is.null(x$alpha)) "alpha" else if (!is.null(x$sharpe)) "sharpe" else "msharpe"
  measure <- .screeningPerf(x)$value
  pipos <- x$pipos; pizero <- x$pizero; pineg <- x$pineg
  lambda <- x$lambda; npeer <- x$npeer

  if (!is.matrix(pizero)) {
    df <- data.frame(fund = seq_along(pizero))
    df[[short]] <- as.vector(measure)
    df$pipos <- as.vector(pipos); df$pizero <- as.vector(pizero)
    df$pineg <- as.vector(pineg)
    df$npeer <- as.vector(npeer)
    df$lambda <- if (is.null(lambda)) NA_real_ else as.vector(lambda)
  } else {
    ncoef <- nrow(pizero); N <- ncol(pizero)
    cn <- rownames(pizero); if (is.null(cn)) cn <- paste0("coef", seq_len(ncoef))
    df <- data.frame(coefficient = rep(cn, times = N),
                     fund = rep(seq_len(N), each = ncoef),
                     stringsAsFactors = FALSE)
    df[[short]] <- as.vector(measure)
    df$pipos <- as.vector(pipos); df$pizero <- as.vector(pizero)
    df$pineg <- as.vector(pineg)
    df$npeer <- as.vector(npeer)
    df$lambda <- if (is.null(lambda)) NA_real_ else as.vector(lambda)
  }
  rownames(df) <- NULL
  df
}

#' @name summary.SCREENING
#' @title Summary method for the 'SCREENING' object
#' @description Produces a ranked summary of a screening: distribution of the
#' performance measure, the peer performance ratios, optional win/loss counts
#' from the pairwise tests, and the top/bottom funds (overall and by
#' out-/underperformance ratio). Originally contributed by Murilo Andre Peres Pereira.
#' @param object A \code{SCREENING} object.
#' @param coef For a \code{screen_beta} screening, which coefficient row to use
#' (default 1, the alpha).
#' @param top Number of best/worst funds to display. Default: 5.
#' @param p_level Significance level for the win/loss counts (used only when the
#' object carries square \code{pval}/\code{tstat} matrices). Default: 0.05.
#' @param ... Further arguments (ignored).
#' @return An object of class \code{summary.SCREENING}.
#' @author Murilo Andre Peres Pereira, David Ardia and Kris Boudt.
#' @seealso \code{\link{alphaScreening}}.
#' @keywords htest
#' @export
#' @importFrom stats quantile median
#' @importFrom utils head tail
summary.SCREENING <- function(object, coef = 1L, top = 5L, p_level = 0.05, ...) {
  x <- object
  measure <- if (!is.null(x$alpha)) "alpha" else if (!is.null(x$sharpe)) "sharpe" else if (!is.null(x$msharpe)) "msharpe" else stop("no recognized measure in the 'SCREENING' object")

  est <- x[[measure]]
  if (is.matrix(est)) est <- est[coef, ]
  est <- as.numeric(est)

  funds <- names(x$pizero)
  if (is.null(funds)) funds <- names(x$n)
  if (is.null(funds) || length(funds) != length(est)) {
    funds <- paste0("Fund ", seq_along(est))
  }
  pick <- function(v) {
    if (is.null(v)) return(rep(NA_real_, length(est)))
    if (is.matrix(v)) v <- v[coef, ]
    as.numeric(v)
  }

  tab <- data.frame(fund = funds, n = pick(x$n), npeer = pick(x$npeer),
                    estimate = est, lambda = pick(x$lambda),
                    pizero = pick(x$pizero), pipos = pick(x$pipos),
                    pineg = pick(x$pineg), stringsAsFactors = FALSE)

  ## win/loss counts only when pval/tstat are square (within-group screening)
  if (is.matrix(x$pval) && is.matrix(x$tstat) &&
      nrow(x$pval) == length(funds) && ncol(x$pval) == length(funds)) {
    tab$wins   <- rowSums((x$pval < p_level) & (x$tstat > 0), na.rm = TRUE)
    tab$losses <- rowSums((x$pval < p_level) & (x$tstat < 0), na.rm = TRUE)
    tab$net    <- tab$wins - tab$losses
  }

  tab$rank <- rank(-tab$estimate, ties.method = "first")
  tab <- tab[order(tab$rank), , drop = FALSE]

  est_stats <- c(min  = min(tab$estimate, na.rm = TRUE),
                 q25  = unname(stats::quantile(tab$estimate, 0.25, na.rm = TRUE)),
                 med  = stats::median(tab$estimate, na.rm = TRUE),
                 mean = mean(tab$estimate, na.rm = TRUE),
                 q75  = unname(stats::quantile(tab$estimate, 0.75, na.rm = TRUE)),
                 max  = max(tab$estimate, na.rm = TRUE))

  top <- max(1L, min(as.integer(top), nrow(tab)))
  ord_plus  <- order(-tab$pipos, -tab$estimate)
  ord_minus <- order(-tab$pineg, -tab$estimate)

  res <- list(measure = measure, coef = coef, p_level = p_level,
              n_funds = nrow(tab), stats = est_stats, table = tab,
              top = utils::head(tab, top), bottom = utils::tail(tab, top),
              top_pi_plus  = utils::head(tab[ord_plus, , drop = FALSE], top),
              top_pi_minus = utils::head(tab[ord_minus, , drop = FALSE], top))
  class(res) <- "summary.SCREENING"
  res
}

#' @name print.summary.SCREENING
#' @title Print method for the 'summary.SCREENING' object
#' @description Print method for the object returned by
#' \code{\link{summary.SCREENING}}.
#' @param x A \code{summary.SCREENING} object.
#' @param ... Further arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @author Murilo Andre Peres Pereira, David Ardia and Kris Boudt.
#' @keywords htest
#' @export
print.summary.SCREENING <- function(x, ...) {
  fmt <- function(z) formatC(as.numeric(z), format = "f", digits = 4)
  round_df <- function(df) {
    cols <- intersect(c("estimate", "pizero", "pipos", "pineg", "lambda"), names(df))
    for (cc in cols) df[[cc]] <- fmt(df[[cc]])
    df
  }
  cat("\nPeer performance screening summary (", x$measure, ")\n", sep = "")
  cat("  Funds: ", x$n_funds, "\n\n", sep = "")
  s <- x$stats
  cat("  ", x$measure, " distribution: min ", fmt(s["min"]), " | 25% ", fmt(s["q25"]),
      " | med ", fmt(s["med"]), " | mean ", fmt(s["mean"]), " | 75% ", fmt(s["q75"]),
      " | max ", fmt(s["max"]), "\n\n", sep = "")
  cat("Top funds (by ", x$measure, "):\n", sep = "")
  print(round_df(x$top), row.names = FALSE, right = TRUE)
  cat("\nTop funds (by outperformance ratio pi+):\n")
  print(round_df(x$top_pi_plus), row.names = FALSE, right = TRUE)
  invisible(x)
}

#' @name exposureHeterogeneity
#' @title Factor exposure heterogeneity from a beta screening
#' @description Aggregates the per-fund equal-performance ratios of a
#' \code{screen_beta = TRUE} \code{\link{alphaScreening}} into the
#' factor-by-factor heterogeneity measure of Ardia et al. (2023). For each
#' coefficient \eqn{k} (the alpha and each factor beta) it reports the average
#' equal-exposure ratio \eqn{\pi^0_k = \frac1N\sum_i \pi^0_{i,k}} and the
#' heterogeneity \eqn{1-\pi^0_k}: the share of peers that are significantly
#' differentiated on coefficient \eqn{k}. A value close to one indicates large
#' heterogeneity (much room to differentiate); close to zero, homogeneity.
#' @param object A \code{SCREENING} object produced with
#' \code{screen_beta = TRUE}.
#' @return A \code{data.frame} of class \code{exposureHeterogeneity} with columns
#' \code{coefficient}, \code{equalExposure} (\eqn{\pi^0_k}), and
#' \code{heterogeneity} (\eqn{1-\pi^0_k}).
#' @author David Ardia and Kris Boudt.
#' @references
#' Ardia, D., Bluteau, K., Lortie-Cloutier, G., Tran, D. (2023).
#' Factor exposure heterogeneity in green and brown stocks.
#' \emph{Finance Research Letters} \bold{55}, Part A, 103900.
#' @seealso \code{\link{alphaScreening}}.
#' @keywords htest
#' @examples
#' \donttest{
#' data("hfdata")
#' set.seed(1234)
#' fac <- matrix(rnorm(nrow(hfdata) * 2), ncol = 2,
#'               dimnames = list(NULL, c("MKT", "SMB")))
#' sc <- alphaScreening(hfdata[, 1:20], factors = fac, screen_beta = TRUE,
#'                      control = list(nCore = 1))
#' exposureHeterogeneity(sc)
#' }
#' @export
exposureHeterogeneity <- function(object) {
  if (!inherits(object, "SCREENING")) {
    stop("'object' must be a 'SCREENING' object")
  }
  if (!is.matrix(object$pizero)) {
    stop("'exposureHeterogeneity' requires a screening run with screen_beta = TRUE")
  }
  pizero <- object$pizero
  cn <- rownames(pizero)
  if (is.null(cn)) cn <- paste0("coef", seq_len(nrow(pizero)))
  eq <- rowMeans(pizero, na.rm = TRUE)
  out <- data.frame(coefficient = cn,
                    equalExposure = as.vector(eq),
                    heterogeneity = as.vector(1 - eq),
                    stringsAsFactors = FALSE)
  rownames(out) <- NULL
  class(out) <- c("exposureHeterogeneity", "data.frame")
  out
}

#' @name plot.exposureHeterogeneity
#' @title Plot factor exposure heterogeneity
#' @description Bar plot of the heterogeneity measure \eqn{1-\pi^0_k} returned by
#' \code{\link{exposureHeterogeneity}}.
#' @param x An \code{exposureHeterogeneity} object.
#' @param ... Further arguments passed to \code{\link[graphics]{barplot}}.
#' @return Invisibly returns \code{x}.
#' @author David Ardia and Kris Boudt.
#' @seealso \code{\link{exposureHeterogeneity}}.
#' @keywords hplot
#' @export
#' @importFrom graphics barplot
plot.exposureHeterogeneity <- function(x, ...) {
  h <- x$heterogeneity
  names(h) <- x$coefficient
  graphics::barplot(h, ylim = c(0, 1), ylab = expression(1 - hat(pi)^0),
                    main = "Factor exposure heterogeneity", ...)
  invisible(x)
}

#' @name plot.SCREENING
#' @title Peer performance screening plot
#' @description Plot method for the \code{SCREENING} object returned by
#' \code{\link{alphaScreening}}, \code{\link{sharpeScreening}} and
#' \code{\link{msharpeScreening}}. It reproduces the peer performance
#' screening plot of Ardia and Boudt (2018): funds are sorted by their
#' performance measure and, for each fund, the estimated outperformance
#' (\eqn{\hat\pi^+}), equal-performance (\eqn{\hat\pi^0}) and underperformance
#' (\eqn{\hat\pi^-}) ratios are displayed as a horizontal stacked bar. The
#' dashed diagonal lines depict the naive percentile-rank benchmark
#' (\eqn{\hat\pi^0 = 0}); the gap between the black \eqn{\hat\pi^+} area and the
#' diagonal visualizes the luck correction.
#' @param x A \code{SCREENING} object.
#' @param nblock Optional number of equally-sized blocks into which the funds
#' (sorted by performance) are aggregated by averaging, as in Ardia and Boudt
#' (2018). Default: \code{nblock = NULL}, i.e. one bar per fund.
#' @param reference A logical value indicating whether the dashed
#' percentile-rank reference lines should be drawn. Default:
#' \code{reference = NULL}, i.e. drawn for within-group screening and omitted
#' for cross-group (\code{Y}-based) screening, where the within-group
#' percentile-rank benchmark does not apply. A single focal fund (cross-group
#' with one fund in \code{X}) is shown as a single stacked bar.
#' @param band Half-width (in percentage points) of the reference band drawn
#' around the central diagonal. Default: \code{band = 27.5}.
#' @param colorset Vector of three colors for the \eqn{\hat\pi^+},
#' \eqn{\hat\pi^0} and \eqn{\hat\pi^-} areas. Default:
#' \code{c(gray(0), gray(0.8), gray(0.5))}.
#' @param ... Further graphical arguments passed to \code{\link[graphics]{barplot}}.
#' @return Invisibly returns the (sorted, possibly aggregated) matrix of
#' ratios that is plotted.
#' @details If the screening was run with \code{screen_beta = TRUE}, the alpha
#' (first) coefficient is used.
#' @author David Ardia and Kris Boudt.
#' @references
#' Ardia, D., Boudt, K. (2018).
#' The peer performance ratios of hedge funds.
#' \emph{Journal of Banking and Finance} \bold{87}, pp.351--368.
#' \doi{10.1016/j.jbankfin.2017.10.014}
#' @seealso \code{\link{alphaScreening}}, \code{\link{sharpeScreening}} and
#' \code{\link{msharpeScreening}}.
#' @keywords hplot
#' @examples
#' \donttest{
#' data("hfdata")
#' set.seed(1234)
#' sc <- alphaScreening(hfdata[, 1:30], control = list(nCore = 1))
#' plot(sc)
#' }
#' @export
#' @importFrom graphics barplot axis box grid par plot.default
#' @importFrom grDevices gray
plot.SCREENING <- function(x, nblock = NULL, reference = NULL, band = 27.5,
                           colorset = c(grDevices::gray(0), grDevices::gray(0.8),
                                        grDevices::gray(0.5)), ...) {
  cross <- isTRUE(x$cross)
  # the percentile-rank diagonal is a within-group benchmark: off for cross-group
  if (is.null(reference)) {
    reference <- !cross
  }
  perf  <- .screeningPerf(x)
  value <- perf$value
  pipos <- x$pipos; pizero <- x$pizero; pineg <- x$pineg

  # screen_beta = TRUE => matrices; use the alpha (first) coefficient
  if (is.matrix(pipos)) {
    pipos <- pipos[1, ]; pizero <- pizero[1, ]; pineg <- pineg[1, ]
  }
  if (is.matrix(value)) {
    value <- value[1, ]
  }
  value <- as.vector(value)

  # keep funds with complete information
  ok <- !is.na(pipos) & !is.na(pizero) & !is.na(pineg) & !is.na(value)
  value <- value[ok]; pipos <- pipos[ok]; pizero <- pizero[ok]; pineg <- pineg[ok]
  if (length(value) == 0) {
    stop("no fund with non-missing ratios to plot")
  }

  # a single focal fund (typically X-vs-Y with one fund): one stacked bar
  if (length(value) == 1) {
    return(.plotScreeningSingle(pipos, pizero, pineg, x$ny, colorset,
                                perf$label, ...))
  }

  # sort by performance (decreasing): best fund first (drawn at the bottom)
  ord <- order(value, decreasing = TRUE)
  value <- value[ord]; pipos <- pipos[ord]; pizero <- pizero[ord]; pineg <- pineg[ord]
  n <- length(value)

  # optional aggregation into 'nblock' equally-sized blocks
  if (!is.null(nblock) && nblock >= 1 && nblock < n) {
    grp <- ceiling(seq_len(n) / (n / nblock))
    grp[grp > nblock] <- nblock
    value <- as.vector(tapply(value, grp, mean))
    m <- cbind(as.vector(tapply(pipos,  grp, mean)),
               as.vector(tapply(pizero, grp, mean)),
               as.vector(tapply(pineg,  grp, mean)))
    m <- m / rowSums(m)  # renormalize so the triple sums to one
    pipos <- m[, 1]; pizero <- m[, 2]; pineg <- m[, 3]
    n <- length(value)
  }

  mat <- rbind(pipos, pizero, pineg)

  # two-panel layout: performance | peer performance ratios
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::par(mfrow = c(1, 2))

  # left panel: performance measure, best fund at the bottom
  graphics::plot.default(value, seq_len(n), type = "b", pch = 20, las = 1,
                         axes = FALSE, xlab = "", ylab = "fund (sorted by performance)",
                         main = perf$label)
  graphics::box(); graphics::grid()
  graphics::axis(side = 1); graphics::axis(side = 2, las = 1)

  # right panel: stacked outperformance / equal / underperformance ratios
  mainstr <- expression(hat(pi)^"+" * " / " * hat(pi)^0 * " / " * hat(pi)^"-")
  graphics::barplot(100 * mat, horiz = TRUE, space = 0, names.arg = NULL,
                    col = colorset, border = NA, xlim = c(0, 100), axes = FALSE,
                    main = mainstr, ...)
  graphics::axis(side = 1, at = seq(0, 100, by = 25),
                 labels = paste0(seq(0, 100, by = 25), "%"))
  graphics::box()

  # dashed percentile-rank reference lines (naive benchmark, pi0 = 0)
  if (reference) {
    yy <- seq(0, 100, length.out = 200)
    xx <- seq(104, -4, length.out = 200)
    for (shift in c(0, -band, band)) {
      graphics::par(new = TRUE)
      graphics::plot.default(xx + shift, yy, type = "l", lty = "dashed", lwd = 1,
                             xlim = c(0, 100), ylim = c(0, 100), axes = FALSE,
                             xlab = "", ylab = "")
    }
  }

  invisible(t(mat))
}

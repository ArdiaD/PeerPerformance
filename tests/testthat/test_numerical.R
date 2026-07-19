context("Numerical regression")

data("hfdata")

test_that("alphaScreening peer ratios are stable and sum to one", {
  rets <- hfdata[, 15:20]
  res  <- alphaScreening(rets, control = list(nCore = 1))
  ## pinned values (verified against the article example)
  expect_equal(unname(round(res$pipos[1], 3)), 0.6)
  expect_equal(unname(round(res$pineg[1], 3)), 0.4)
  expect_equal(unname(res$npeer), rep(5, 6))
  ## the triple sums to one for every fund
  s <- res$pizero + res$pipos + res$pineg
  expect_true(all(abs(s - 1) < 1e-8))
})

test_that("gammaPos / gammaNeg are honoured and only affect the +/- split", {
  rets <- hfdata[, 1:30]
  ## fix lambda so that pizero is deterministic (the data-driven lambda uses an
  ## unseeded bootstrap); pizero must then be invariant to the gamma thresholds
  a <- alphaScreening(rets, control = list(nCore = 1, lambda = 0.5,
                                           gammaPos = 0.2, gammaNeg = 0.8))
  b <- alphaScreening(rets, control = list(nCore = 1, lambda = 0.5,
                                           gammaPos = 0.6, gammaNeg = 0.4))
  ## pizero is invariant to the gamma thresholds
  expect_equal(a$pizero, b$pizero, tolerance = 1e-10)
  ## the out/under split does respond
  expect_false(isTRUE(all.equal(a$pipos, b$pipos)))
  ## ratios remain valid
  for (res in list(a, b)) {
    expect_true(all(res$pipos  >= 0 & res$pipos  <= 1, na.rm = TRUE))
    expect_true(all(res$pineg  >= 0 & res$pineg  <= 1, na.rm = TRUE))
    expect_true(all(res$pizero >= 0 & res$pizero <= 1, na.rm = TRUE))
  }
})

test_that("alphaTesting returns a coherent htest-like list", {
  x <- hfdata[, 1]; y <- hfdata[, 2]
  out <- alphaTesting(x, y)
  expect_true(out$pval >= 0 && out$pval <= 1)
  expect_equal(length(out$alpha), 2)
  ## dalpha equals the difference of the two fund alphas
  expect_equal(unname(out$dalpha), unname(out$alpha[1] - out$alpha[2]), tolerance = 1e-8)
})

test_that("alphaScreening with screen_beta returns one row per coefficient", {
  set.seed(1)
  rets <- hfdata[, 1:5]
  fmat <- matrix(rnorm(nrow(rets)), ncol = 1)        # single factor
  res  <- alphaScreening(rets, factors = fmat, control = list(nCore = 1),
                         screen_beta = TRUE)
  ## first dimension = intercept (alpha) + 1 factor = 2
  expect_equal(nrow(res$pizero), 2L)
  expect_equal(ncol(res$pizero), 5L)
})

test_that("NA values are handled without error", {
  rets <- hfdata[, 1:6]
  rets[1:5, 1] <- NA
  res <- alphaScreening(rets, control = list(nCore = 1))
  expect_equal(length(res$pizero), 6L)
  expect_true(all(is.finite(res$npeer)))
})

test_that("audit-v2: output naming and control validation", {
  X <- hfdata[, 1:10]; colnames(X) <- paste0("HF", 1:10)
  ## as.data.frame keeps fund names; targetPeerPerformance labels its outputs
  tp <- targetPeerPerformance(X, funds = c("HF3", "HF7"), method = "sharpe",
                              control = list(nCore = 1))
  expect_equal(as.data.frame(tp)$fund, c("HF3", "HF7"))
  ## rollScreening: screen_beta requires factors (alpha) -> attr FALSE without
  suppressWarnings(
    rb <- rollScreening(hfdata[, 1:12], screen = "alpha", width = 40, by = 20,
                        control = list(nCore = 1, screen_beta = TRUE)))
  expect_false(isTRUE(attr(rb, "screen_beta")))
  ## input / control validation
  expect_error(targetPeerPerformance(X, funds = 1.9, control = list(nCore = 1)))
  expect_error(alphaScreening(X[, 1:4], control = list(nCore = 1, hac = c(TRUE, FALSE))))
  expect_error(alphaScreening(X[, 1:4], control = list(nCore = 1, lambda = c(0.4, 0.5))))
})

test_that("sharpe/msharpe screening on an unbalanced panel has no NaN (PR #14)", {
  set.seed(1)
  Tn <- 80
  X <- matrix(rnorm(Tn * 4, 0.01, 0.05), Tn, 4)
  X[1:25, 2] <- NA   # the FIRST peer of fund 1 is missing early
  ss <- sharpeScreening(X, control = list(nCore = 1))
  ms <- msharpeScreening(X, control = list(nCore = 1))
  ## with the old 'X[idx[, k], 1]' indexing the focal returns were NA-contaminated
  ## for some pairs, producing NaN p-values
  expect_false(any(is.nan(ss$pval)))
  expect_false(any(is.nan(ms$pval)))
  expect_true(all(ss$pizero + ss$pipos + ss$pineg - 1 < 1e-8, na.rm = TRUE))
})

test_that("audit fixes: degenerate inputs are handled", {
  ## A1: screen_beta = TRUE without factors warns and falls back (no crash)
  expect_warning(r <- alphaScreening(hfdata[, 1:5], control = list(nCore = 1, screen_beta = TRUE)))
  expect_false(is.matrix(r$pizero))

  ## A2: bBoot = 0 in screening errors cleanly
  expect_error(sharpeScreening(hfdata[, 1:4], control = list(nCore = 1, type = 2, bBoot = 0)))
  expect_error(msharpeScreening(hfdata[, 1:4], control = list(nCore = 1, type = 2, bBoot = 0)))

  ## A3: a peer differing from the focal fund by a constant is excluded
  Y <- hfdata[, 11:14]; Y[, 2] <- hfdata[, 1] + 0.5
  s <- alphaScreening(hfdata[, 1], Y = Y, control = list(nCore = 1))
  expect_equal(s$npeer, 3L)
  expect_false(any(s$pval == 0, na.rm = TRUE))

  ## A4: rollScreening only stamps screen_beta for the alpha screen
  rb <- rollScreening(hfdata[, 1:20], screen = "sharpe", width = 40, by = 20,
                      control = list(nCore = 1, screen_beta = TRUE))
  expect_false(isTRUE(attr(rb, "screen_beta")))

  ## B5: invalid control values are rejected
  expect_error(alphaScreening(hfdata[, 1:4], control = list(nCore = 1, type = 3)))
  expect_error(alphaScreening(hfdata[, 1:4], control = list(nCore = 1, gammaPos = 1.5)))
})

test_that("cross-group screening (X vs Y) works and excludes self", {
  ## single focal fund against a peer group
  s <- alphaScreening(hfdata[, 1], Y = hfdata[, 11:30], control = list(nCore = 1))
  expect_equal(length(s$pizero), 1L)
  expect_equal(s$npeer, 20L)
  expect_true(abs(s$pizero + s$pipos + s$pineg - 1) < 1e-8)

  ## group X against group Y: ratios over the nY peers
  g <- alphaScreening(hfdata[, 1:5], Y = hfdata[, 11:30], control = list(nCore = 1))
  expect_equal(dim(g$pval), c(5L, 20L))
  expect_equal(g$ny, 20L)
  expect_true(all(abs(g$pizero + g$pipos + g$pineg - 1) < 1e-8))

  ## self-exclusion when X is a subset of Y
  s2 <- alphaScreening(hfdata[, 1], Y = hfdata[, 1:30], control = list(nCore = 1))
  expect_equal(s2$npeer, 29L)

  ## Sharpe and modified Sharpe cross-group
  sh <- sharpeScreening(hfdata[, 1:3], Y = hfdata[, 11:30], control = list(nCore = 1))
  expect_equal(dim(sh$pval), c(3L, 20L))
  ms <- msharpeScreening(hfdata[, 1:3], Y = hfdata[, 11:30], level = 0.95,
                         control = list(nCore = 1))
  expect_equal(dim(ms$pval), c(3L, 20L))
})

test_that("as.data.frame and exposureHeterogeneity work", {
  set.seed(1)
  w <- alphaScreening(hfdata[, 1:8], control = list(nCore = 1))
  df <- as.data.frame(w)
  expect_equal(nrow(df), 8L)
  expect_true(all(c("fund", "alpha", "pipos", "pizero", "pineg", "npeer") %in% names(df)))

  set.seed(1)
  fac <- matrix(rnorm(nrow(hfdata) * 2), ncol = 2,
                dimnames = list(NULL, c("MKT", "SMB")))
  scb <- alphaScreening(hfdata[, 1:15], factors = fac, screen_beta = TRUE,
                        control = list(nCore = 1))
  expect_equal(rownames(scb$pizero), c("alpha", "MKT", "SMB"))
  dfb <- as.data.frame(scb)
  expect_equal(nrow(dfb), 15L * 3L)
  expect_true("coefficient" %in% names(dfb))

  eh <- exposureHeterogeneity(scb)
  expect_s3_class(eh, "exposureHeterogeneity")
  expect_equal(nrow(eh), 3L)
  expect_equal(eh$heterogeneity, 1 - eh$equalExposure, tolerance = 1e-12)
  expect_true(all(eh$heterogeneity >= 0 & eh$heterogeneity <= 1))
})

test_that("screen_beta can be set via control and the argument overrides it", {
  set.seed(1)
  fac <- matrix(rnorm(nrow(hfdata) * 2), ncol = 2,
                dimnames = list(NULL, c("MKT", "SMB")))
  a <- alphaScreening(hfdata[, 1:8], factors = fac,
                      control = list(nCore = 1, screen_beta = TRUE))
  expect_true(is.matrix(a$pizero))
  expect_equal(rownames(a$pizero), c("alpha", "MKT", "SMB"))
  ## explicit argument wins over the control value
  b <- alphaScreening(hfdata[, 1:8], factors = fac, screen_beta = FALSE,
                      control = list(nCore = 1, screen_beta = TRUE))
  expect_false(is.matrix(b$pizero))
})

test_that("rollScreening returns a tidy time series", {
  set.seed(1234)
  roll <- rollScreening(hfdata[, 1:30], screen = "alpha", width = 36, by = 8,
                        control = list(nCore = 1))
  expect_s3_class(roll, "rollScreening")
  expect_true(all(c("window", "index", "pizero", "pipos", "pineg",
                    "heterogeneity") %in% names(roll)))
  expect_equal(roll$heterogeneity, 1 - roll$pizero, tolerance = 1e-12)
  expect_equal(nrow(roll), length(seq.int(1, nrow(hfdata) - 36 + 1, by = 8)))

  ## screen_beta -> one row per window/coefficient
  set.seed(1)
  fac <- matrix(rnorm(nrow(hfdata) * 2), ncol = 2,
                dimnames = list(NULL, c("MKT", "SMB")))
  rb <- rollScreening(hfdata[, 1:30], factors = fac, width = 36, by = 12,
                      control = list(nCore = 1, screen_beta = TRUE))
  ## set comparison (locale-independent: avoids C vs UTF-8 sort order)
  expect_setequal(unique(rb$coefficient), c("alpha", "MKT", "SMB"))
  pf <- tempfile(fileext = ".pdf"); pdf(pf); plot(rb); dev.off(); unlink(pf)
})

test_that("plot works on a single focal fund (cross-group)", {
  s <- alphaScreening(hfdata[, 1], Y = hfdata[, 11:30], control = list(nCore = 1))
  pf <- tempfile(fileext = ".pdf"); pdf(pf)
  out <- plot(s)
  dev.off(); unlink(pf)
  expect_equal(length(out), 3L)
})

test_that("targetPeerPerformance equals screening with Y = X", {
  rets <- hfdata[, 1:10]
  focals <- c(2, 5, 7)
  tp <- targetPeerPerformance(rets, funds = focals, method = "alpha",
                              control = list(nCore = 1, lambda = 0.5))
  yx <- alphaScreening(rets[, focals], Y = rets,
                       control = list(nCore = 1, lambda = 0.5))
  expect_equal(unname(tp$pizero), unname(yx$pizero), tolerance = 1e-10)
  expect_equal(unname(tp$pipos),  unname(yx$pipos),  tolerance = 1e-10)
  expect_equal(unname(tp$npeer),  rep(9L, 3))
  expect_s3_class(tp, "SCREENING")
  ## selection by name works and rows are labelled
  colnames(rets) <- paste0("F", 1:10)
  tp2 <- targetPeerPerformance(rets, funds = c("F2", "F5"), method = "sharpe",
                               control = list(nCore = 1))
  expect_equal(rownames(tp2$pval), c("F2", "F5"))
})

test_that("summary.SCREENING produces a ranked table", {
  set.seed(1)
  sc <- alphaScreening(hfdata[, 1:12], control = list(nCore = 1))
  s <- summary(sc)
  expect_s3_class(s, "summary.SCREENING")
  expect_equal(nrow(s$table), 12L)
  expect_true(all(c("estimate", "pipos", "pineg", "wins", "losses") %in% names(s$table)))
  expect_identical(print(s), s)
})

test_that("print and plot methods dispatch and return invisibly", {
  set.seed(1234)
  sc <- alphaScreening(hfdata[, 1:12], control = list(nCore = 1))
  expect_s3_class(sc, "SCREENING")
  expect_identical(print(sc), sc)

  tt <- msharpeTesting(hfdata[, 1], hfdata[, 2], level = 0.95)
  expect_s3_class(tt, "TESTING")
  expect_identical(print(tt), tt)

  ## plot to a throwaway device; returns the plotted ratio matrix invisibly
  pf <- tempfile(fileext = ".pdf"); pdf(pf)
  out <- plot(sc, nblock = 6)
  dev.off(); unlink(pf)
  expect_equal(ncol(out), 3L)
  expect_true(all(abs(rowSums(out) - 1) < 1e-8))
})

test_that("seeded bootstrap tests are reproducible", {
  x <- hfdata[, 1]; y <- hfdata[, 2]
  set.seed(321); p1 <- sharpeTesting(x, y, control = list(type = 2, nBoot = 200))$pval
  set.seed(321); p2 <- sharpeTesting(x, y, control = list(type = 2, nBoot = 200))$pval
  expect_identical(p1, p2)

  set.seed(321); q1 <- msharpeTesting(x, y, control = list(type = 2, nBoot = 200))$pval
  set.seed(321); q2 <- msharpeTesting(x, y, control = list(type = 2, nBoot = 200))$pval
  expect_identical(q1, q2)
})

test_that("confint.SCREENING brackets the point estimate and is valid", {
  rets <- hfdata[, 1:12]
  sc   <- alphaScreening(rets, control = list(nCore = 1))
  set.seed(42)
  ci <- confint(sc, parm = "pipos", nBoot = 200)
  ## shape and naming
  expect_equal(nrow(ci), ncol(rets))
  expect_equal(ncol(ci), 2L)
  expect_equal(rownames(ci), colnames(rets))
  est <- attr(ci, "estimate")
  ok  <- !is.na(ci[, 1]) & !is.na(ci[, 2])
  ## bounds are ordered, inside [0, 1], and contain the point estimate
  expect_true(all(ci[ok, 1] <= ci[ok, 2]))
  expect_true(all(ci[ok, ] >= 0 & ci[ok, ] <= 1))
  expect_true(all(ci[ok, 1] <= est[ok] + 1e-8 & est[ok] <= ci[ok, 2] + 1e-8))
  ## all three ratios are supported: valid, ordered bounds that bracket the estimate
  for (p in c("pizero", "pineg")) {
    set.seed(42)
    cp  <- confint(sc, parm = p, nBoot = 200)
    ep  <- attr(cp, "estimate")
    okp <- !is.na(cp[, 1]) & !is.na(cp[, 2])
    expect_equal(dim(cp), c(ncol(rets), 2L))
    expect_true(all(cp[okp, 1] <= cp[okp, 2]))
    expect_true(all(cp[okp, ] >= 0 & cp[okp, ] <= 1))
    expect_true(all(cp[okp, 1] <= ep[okp] + 1e-8 & ep[okp] <= cp[okp, 2] + 1e-8))
  }
  ## screen_beta screenings are rejected
  scb <- alphaScreening(rets, factors = hfdata[, 50, drop = FALSE],
                        control = list(nCore = 1, screen_beta = TRUE))
  expect_error(confint(scb), "screen_beta")
})

test_that("alphaTesting screen_beta returns a (K+1) x 2 alpha matrix (HAC and not)", {
  x <- hfdata[, 1]; y <- hfdata[, 2]
  fac <- hfdata[, 50:51]                       # K = 2 factors
  for (use_hac in c(FALSE, TRUE)) {
    res <- alphaTesting(x, y, factors = fac,
                        control = list(hac = use_hac), screen_beta = TRUE)
    expect_equal(dim(res$alpha), c(ncol(fac) + 1L, 2L))   # rows = coef, cols = x/y
    expect_length(res$dalpha, ncol(fac) + 1L)
    ## print shows the two funds' *alphas* (row 1), not a beta
    expect_output(print(res), "Peer performance test")
  }
})

test_that("control validation rejects non-whole and out-of-range values", {
  rets <- hfdata[, 1:4]
  expect_error(alphaScreening(rets, control = list(nCore = 1, nBoot = 2.5)), "whole number")
  expect_error(alphaScreening(rets, control = list(nCore = 1, minObs = -1)))
  ## the block length is only relevant to the bootstrap test (type = 2)
  expect_error(sharpeScreening(rets, control = list(nCore = 1, type = 2,
                                                    bBoot = nrow(rets) + 1)),
               "cannot exceed")
  ## degenerate within-group input
  expect_error(alphaScreening(hfdata[, 1], control = list(nCore = 1)), "at least two funds")
})

test_that("bootstrap screening is valid and seeded-reproducible on unbalanced panels", {
  ## unbalanced panel: distinct complete-case lengths across pairs
  X <- hfdata[, 1:6]
  X[1:10, 2]  <- NA
  X[1:20, 3]  <- NA
  X[41:60, 4] <- NA
  ctr <- list(nCore = 1, type = 2, bBoot = 3, nBoot = 99)

  set.seed(99); s1 <- sharpeScreening(X, control = ctr)
  set.seed(99); s2 <- sharpeScreening(X, control = ctr)
  expect_identical(s1$pval, s2$pval)                    # all RNG in the master
  ok <- !is.na(s1$pval)
  expect_true(any(ok))
  expect_true(all(s1$pval[ok] >= 0 & s1$pval[ok] <= 1))
  s <- s1$pizero + s1$pipos + s1$pineg
  expect_true(all(abs(s[!is.na(s)] - 1) < 1e-8))

  ## cross-group path with the same panel
  set.seed(99); x1 <- sharpeScreening(X[, 1:2], Y = X[, 3:6], control = ctr)
  set.seed(99); x2 <- sharpeScreening(X[, 1:2], Y = X[, 3:6], control = ctr)
  expect_identical(x1$pval, x2$pval)
  okx <- !is.na(x1$pval)
  expect_true(any(okx))
  expect_true(all(x1$pval[okx] >= 0 & x1$pval[okx] <= 1))

  ## a block length larger than a short pair's sample leaves that pair NA
  ## (fund 7 overlaps fund 8 on only 5 observations < bBoot = 8)
  Z <- hfdata[, 7:9]
  Z[1:55, 2] <- NA
  set.seed(7)
  expect_warning(
    sz <- sharpeScreening(Z, control = list(nCore = 1, type = 2, bBoot = 8,
                                            nBoot = 99, minObs = 3)),
    "left untested")
  expect_true(is.na(sz$pval[1, 2]))                     # untestable pair
  expect_false(is.na(sz$pval[1, 3]))                    # full-length pair still tested
})

test_that("serial (nCore = 1) and cluster (nCore = 2) paths give identical results", {
  skip_on_cran()   # keep CRAN runs light; 2 cores are exercised locally/CI
  rets <- hfdata[, 1:8]
  ## fix lambda so the comparison does not depend on the data-driven bootstrap
  ctr1 <- list(nCore = 1, lambda = 0.5)
  ctr2 <- list(nCore = 2, lambda = 0.5)
  a1 <- alphaScreening(rets, control = ctr1)
  a2 <- alphaScreening(rets, control = ctr2)
  expect_equal(a1$pval,  a2$pval,  tolerance = 1e-12)
  expect_equal(a1$pipos, a2$pipos, tolerance = 1e-12)
  ## bootstrapped Sharpe screening: indices are drawn in the master, so the
  ## result must not depend on the number of workers
  set.seed(3); b1 <- sharpeScreening(rets, control = c(ctr1, type = 2, bBoot = 2, nBoot = 99))
  set.seed(3); b2 <- sharpeScreening(rets, control = c(ctr2, type = 2, bBoot = 2, nBoot = 99))
  expect_equal(b1$pval, b2$pval, tolerance = 1e-12)
})

test_that("control$fastAdjust matches the default path and is faster", {
  ## adjustPi: the fast vectorised inversion agrees with the uniroot path
  ## across the (n, lambda) grid, well inside uniroot's own default tolerance
  set.seed(1)
  ph <- runif(200, 0.3, 1.0)
  for (nn in c(5, 9, 29, 99)) {
    for (lam in c(0.3, 0.4, 0.5, 0.6, 0.7)) {
      a <- PeerPerformance:::adjustPi(ph, n = nn, lambda = lam, fast = FALSE)
      b <- PeerPerformance:::adjustPi(ph, n = nn, lambda = lam, fast = TRUE)
      expect_equal(a, b, tolerance = 1e-3)       # uniroot tol is ~1.2e-4
      expect_true(all(b >= 0 & b <= 1))
    }
  }
  ## the flag must not change semantics, only speed: NA/NaN inputs and targets
  ## outside the uniroot bracket must behave exactly as the default path
  edge <- c(0.5, NA, NaN, 0.9, -0.5, 0, 1.4, 3.0)
  e0 <- PeerPerformance:::adjustPi(edge, n = 99, lambda = 0.5, fast = FALSE)
  e1 <- PeerPerformance:::adjustPi(edge, n = 99, lambda = 0.5, fast = TRUE)
  expect_equal(is.na(e0), is.na(e1))                    # same NA pattern
  expect_equal(e0[!is.na(e0)], e1[!is.na(e1)], tolerance = 1e-3)
  ## NaN can reach adjustPi through an all-NA p-value row; must not error
  pv <- matrix(c(NA, NA, NA, 0.2, 0.6, 0.9), nrow = 2, byrow = TRUE)
  expect_silent(PeerPerformance:::computePizero(pv, lambda = 0.5,
                                                adjust = TRUE, fast = TRUE))

  ## end to end: a screening with a fixed lambda is unchanged by the flag
  rets <- hfdata[, 1:12]
  s0 <- alphaScreening(rets, control = list(nCore = 1, lambda = 0.5))
  s1 <- alphaScreening(rets, control = list(nCore = 1, lambda = 0.5,
                                            fastAdjust = TRUE))
  expect_equal(s0$pizero, s1$pizero, tolerance = 1e-3)
  expect_equal(s0$pipos,  s1$pipos,  tolerance = 1e-3)
  ## the flag is validated like the other logical controls
  expect_error(alphaScreening(rets, control = list(nCore = 1,
                                                   fastAdjust = c(TRUE, FALSE))))

  ## round-trip: the fast path really solves the inversion (this is the strong
  ## check; agreement with the default is capped by uniroot's own tolerance)
  fwd <- function(pi0, n, lambda) {
    nlambda <- pi0 * n * (1 - lambda)
    out <- pi0
    i <- nlambda < n
    s <- sqrt(nlambda[i] * (n - nlambda[i])/(n^3 * (1 - lambda)^2))
    z <- (1 - pi0[i])/s
    out[i] <- pi0[i] + s * (-dnorm(z) + (1 - pnorm(z)) * z)
    out
  }
  for (nn in c(9, 99)) {
    for (lam in c(0.3, 0.5, 0.7)) {
      r <- PeerPerformance:::adjustPi(ph, n = nn, lambda = lam, fast = TRUE)
      ## only entries that were actually inverted: the target must lie inside
      ## the bracket [fwd(1e-5), fwd(1.5)] (otherwise uniroot fails and both
      ## paths fall back to the input), and the result must be unclamped
      brack <- ph >= fwd(1e-05, nn, lam) & ph <= fwd(1.5, nn, lam)
      inner <- brack & r > 0 & r < 1
      expect_true(any(inner))
      expect_equal(fwd(r[inner], nn, lam), ph[inner], tolerance = 1e-9)
    }
  }
})

test_that("asymptotic screening does not depend on the bootstrap block length", {
  ## bBoot is irrelevant when type = 1; an unbalanced panel with short pairs
  ## must not warn or error just because bBoot exceeds some pair length
  X <- hfdata[, 1:6]
  X[1:52, 2] <- NA                            # one pair with only 8 obs
  expect_silent(s1 <- sharpeScreening(X, control = list(nCore = 1, type = 1,
                                                        bBoot = 20, minObs = 5)))
  expect_silent(m1 <- msharpeScreening(X, control = list(nCore = 1, type = 1,
                                                         bBoot = 20, minObs = 5)))
  expect_true(any(!is.na(s1$pval)))
  expect_true(any(!is.na(m1$pval)))
})

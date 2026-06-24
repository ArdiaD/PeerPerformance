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
  expect_equal(sort(unique(rb$coefficient)), c("alpha", "MKT", "SMB"))
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

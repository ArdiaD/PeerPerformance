## ---------------------------------------------------------------------------
## Extended Monte-Carlo validation of the PeerPerformance estimators
##
## The companion script validation.R checks the estimators under i.i.d. Gaussian
## returns with independent funds. That is the setting in which the luck
## correction is easiest; it is not the setting the package exists for. This
## script covers the two questions that regime cannot answer:
##
##   (D) What does the estimator return under an exact null at the design of a
##       realistic application (N = 100 funds, T = 60 months)? Without this the
##       empirical mean pi0 of about 0.78 on hfdata has no benchmark: a reader
##       cannot tell whether 22% non-ties is more than the estimator reports
##       when nothing is there.
##
##   (E) How does it behave when the assumptions are violated in the ways fund
##       returns actually violate them -- fat tails, serial correlation, and
##       cross-sectional correlation and volatility heterogeneity (hfdata has
##       a mean pairwise correlation of about 0.25 and a 16-fold spread in
##       standard deviations)?
##
##   (F) Does the modified Sharpe test keep its size under non-normality and
##       autocorrelation, and does the studentized circular bootstrap
##       with prespecified blocks (type = 2; bBoot = 1 for i.i.d. returns and
##       bBoot = 5 for AR(1)) repair it where the asymptotic test fails? This
##       experiment validates the fixed-block bootstrap.
##
## This script is slower than validation.R (tens of minutes). Run with:
##   source(system.file("scripts", "validation-extended.R",
##                      package = "PeerPerformance"))
## ---------------------------------------------------------------------------

library("PeerPerformance")
set.seed(1234)

## fastAdjust only changes the pi0 bias correction by a few 1e-5 (see NEWS) and
## makes the data-driven lambda affordable at these sample sizes.
ctr <- list(nCore = 1, fastAdjust = TRUE)

RS <- 50L    # replications for the screening scenarios
RT <- 400L   # replications for the test-size scenarios

## helper: mean of the triple over funds, for one screening
ratios <- function(X) {
  sc <- alphaScreening(X, control = ctr)
  c(pi0 = mean(sc$pizero, na.rm = TRUE),
    pip = mean(sc$pipos,  na.rm = TRUE),
    pim = mean(sc$pineg,  na.rm = TRUE))
}

## ===========================================================================
## (D) Null calibration across the design grid
## ===========================================================================
cat("\n(D) Exact null: what the estimator returns when nothing is there\n")
cat("    (true pi0 = 1; the non-tie mass pi+ + pi- is the false-discovery floor)\n\n")
cat(sprintf("    %5s %5s   %14s %14s\n", "N", "T", "pi0 (s.e.)", "pi+ + pi- (s.e.)"))
gridD <- expand.grid(N = c(30L, 100L), T = c(60L, 120L))
outD <- vector("list", nrow(gridD))
for (g in seq_len(nrow(gridD))) {
  N <- gridD$N[g]; TT <- gridD$T[g]
  m <- matrix(NA_real_, RS, 2)
  for (r in seq_len(RS)) {
    X <- matrix(stats::rnorm(TT * N, 0.005, 0.04), TT, N)
    z <- ratios(X)
    m[r, ] <- c(z["pi0"], z["pip"] + z["pim"])
  }
  outD[[g]] <- c(N = N, T = TT, colMeans(m),
                 apply(m, 2, stats::sd) / sqrt(RS))
  cat(sprintf("    %5d %5d   %.3f (%.3f)   %.3f (%.3f)\n", N, TT,
              mean(m[, 1]), stats::sd(m[, 1])/sqrt(RS),
              mean(m[, 2]), stats::sd(m[, 2])/sqrt(RS)))
}
cat("\n    For comparison, the empirical mean pi0 on hfdata (N = 100, T = 60)\n")
data("hfdata")
set.seed(1234)
emp <- mean(alphaScreening(hfdata, control = ctr)$pizero, na.rm = TRUE)
cat(sprintf("    is %.3f, i.e. a non-tie mass of %.3f.\n", emp, 1 - emp))

## ===========================================================================
## (E) Departures from the i.i.d. Gaussian null, at the empirical design
## ===========================================================================
cat("\n(E) Exact null at N = 100, T = 60 under violated assumptions\n\n")
N <- 100L; TT <- 60L
data("hfdata")
hfSd <- apply(hfdata, 2, stats::sd, na.rm = TRUE)

genGauss <- function() matrix(stats::rnorm(TT * N, 0.005, 0.04), TT, N)
genT5    <- function() {
  ## i.i.d. Student-t(5) returns, standardized to the Gaussian reference
  ## variance while preserving the common expected return of 0.005.
  matrix(0.005 + 0.04 * stats::rt(TT * N, df = 5)/sqrt(5/3), TT, N)
}
genAR1   <- function(rho = 0.3) {
  e <- matrix(stats::rnorm(TT * N, 0, 0.04 * sqrt(1 - rho^2)), TT, N)
  X <- matrix(NA_real_, TT, N)
  X[1, ] <- stats::rnorm(N, 0, 0.04)
  for (t in 2:TT) X[t, ] <- rho * X[t - 1, ] + e[t, ]
  X + 0.005
}
genFac   <- function(rc = 0.25) {
  ## Heterogeneous standardized loadings leave (beta_i - beta_j) f_t in
  ## pairwise differences, while preserving a marginal volatility of 0.04.
  beta <- stats::runif(N, 0.5, 1.5)
  f <- stats::rnorm(TT)
  X <- matrix(stats::rnorm(TT * N), TT, N)
  X <- sweep(X, 2, sqrt(1 - rc * beta^2), "*") +
       tcrossprod(sqrt(rc) * f, beta)
  0.005 + 0.04 * X
}
genVol   <- function() {
  ## Preserve the empirical cross-sectional distribution of volatility while
  ## imposing the exact null of a common expected return.
  0.005 + sweep(matrix(stats::rnorm(TT * N), TT, N), 2, hfSd, "*")
}
## AR(1) idiosyncratic returns plus an AR(1) common factor with heterogeneous
## standardized loadings and the empirical marginal volatilities. This matches
## the serial dependence, average cross-sectional dependence, and volatility
## distribution of hfdata under the exact null of a common expected return.
genBoth  <- function(rho = 0.20, rc = 0.25) {
  ar1 <- function(n, s) {
    e <- stats::rnorm(n, 0, s * sqrt(1 - rho^2))
    x <- numeric(n); x[1] <- stats::rnorm(1, 0, s)
    for (t in 2:n) x[t] <- rho * x[t - 1] + e[t]
    x
  }
  beta <- stats::runif(N, 0.5, 1.5)
  f <- ar1(TT, 1)
  X <- matrix(0, TT, N)
  for (j in seq_len(N)) {
    X[, j] <- hfSd[j] * (sqrt(rc) * beta[j] * f +
                          sqrt(1 - rc * beta[j]^2) * ar1(TT, 1))
  }
  0.005 + X
}
gens <- list("i.i.d. Gaussian (reference)" = genGauss,
             "i.i.d. standardized t(5)"    = genT5,
             "AR(1), rho = 0.2"            = function() genAR1(0.20),
             "AR(1), rho = 0.3"            = genAR1,
             "heterogeneous loadings"       = genFac,
             "heterogeneous volatility"     = genVol,
             "calibrated null"              = genBoth)
cat(sprintf("    %-28s %14s %14s\n", "return process", "pi0 (s.e.)", "pi+ + pi- (s.e.)"))
for (nm in names(gens)) {
  m <- matrix(NA_real_, RS, 2)
  for (r in seq_len(RS)) {
    z <- ratios(gens[[nm]]())
    m[r, ] <- c(z["pi0"], z["pip"] + z["pim"])
  }
  cat(sprintf("    %-28s %.3f (%.3f)   %.3f (%.3f)\n", nm,
              mean(m[, 1]), stats::sd(m[, 1])/sqrt(RS),
              mean(m[, 2]), stats::sd(m[, 2])/sqrt(RS)))
}

## ===========================================================================
## (E1) Does the package's own HAC option repair the inflated floor?
##
## The article recommends hac = TRUE for autocorrelated series, so that
## recommendation should rest on evidence. Both columns are computed on the
## SAME replications, so they differ only in the hac setting.
## ===========================================================================
cat("\n(E1) False-discovery floor with and without HAC standard errors\n\n")
hacGens <- list("i.i.d. Gaussian"           = genGauss,
                "AR(1), rho = 0.2"          = function() genAR1(0.20),
                "AR(1), rho = 0.3"          = genAR1,
                "calibrated null"            = genBoth)
RS_HAC <- 40L   # the HAC panel is twice the cost per replication
cat(sprintf("    %-26s %18s %18s\n", "null process",
            "floor hac = FALSE", "floor hac = TRUE"))
for (nm in names(hacGens)) {
  m <- matrix(NA_real_, RS_HAC, 2)
  for (r in seq_len(RS_HAC)) {
    X <- hacGens[[nm]]()
    m[r, 1] <- 1 - mean(alphaScreening(X, control = ctr)$pizero, na.rm = TRUE)
    m[r, 2] <- 1 - mean(alphaScreening(X, control = c(ctr, list(hac = TRUE)))$pizero,
                        na.rm = TRUE)
  }
  cat(sprintf("    %-26s   %.3f (%.3f)     %.3f (%.3f)\n", nm,
              mean(m[, 1]), stats::sd(m[, 1])/sqrt(RS_HAC),
              mean(m[, 2]), stats::sd(m[, 2])/sqrt(RS_HAC)))
}
## the same contrast on the real data
data("hfdata")
for (h in c(FALSE, TRUE)) {
  set.seed(1234)
  sc <- alphaScreening(hfdata, control = c(ctr, list(hac = h)))
  cat(sprintf("    hfdata, hac = %-5s: non-tie mass %.3f\n", h,
              1 - mean(sc$pizero, na.rm = TRUE)))
}

## ===========================================================================
## (E2) The maximum pi+ over the cross-section, under the calibrated null
##
## Picking the best fund out of N is an extreme-value operation: max_i pi+_i
## must be judged against the distribution of that maximum under the null,
## not against zero. Without this, the single most striking number in any
## applied screening has no benchmark at all.
## ===========================================================================
cat("\n(E2) Largest pi+ across the cross-section, calibrated null\n",
    "     (AR(1) 0.2 + factor 0.25 + heterogeneous loadings and volatility)\n\n")
mx <- numeric(RS)
for (r in seq_len(RS)) {
  sc <- alphaScreening(genBoth(), control = ctr)
  mx[r] <- max(sc$pipos, na.rm = TRUE)
}
data("hfdata")
set.seed(1234)
empmax <- max(alphaScreening(hfdata, control = ctr)$pipos, na.rm = TRUE)
cat(sprintf("    null max pi+ : mean %.3f (s.d. %.3f), 95th pct %.3f\n",
            mean(mx), stats::sd(mx), stats::quantile(mx, 0.95)))
cat(sprintf("    hfdata max pi+ : %.3f -> %.0f%% of null replications exceed it\n",
            empmax, 100 * mean(mx >= empmax)))

## Dependence actually present in hfdata, for reference
ac1 <- apply(hfdata, 2, function(z) {
  z <- z[is.finite(z)]; stats::cor(z[-1], z[-length(z)])
})
cc <- stats::cor(hfdata, use = "pairwise.complete.obs")
cat(sprintf("    hfdata dependence: mean lag-1 AC %.3f (%.0f%% above 0.2),",
            mean(ac1, na.rm = TRUE), 100 * mean(ac1 > 0.2, na.rm = TRUE)))
cat(sprintf(" mean pairwise cor %.3f\n", mean(cc[upper.tri(cc)], na.rm = TRUE)))
cat(sprintf("    hfdata volatility: min %.3f, max %.3f (%.1f-fold spread)\n",
            min(hfSd), max(hfSd), max(hfSd) / min(hfSd)))

## ===========================================================================
## (F) Size of the modified Sharpe test: asymptotic vs studentized bootstrap
## ===========================================================================
cat("\n(F) Empirical size of the modified Sharpe equality test (nominal 5%)\n")
cat("    asymptotic (type = 1, the default) vs studentized circular bootstrap\n")
cat("    (type = 2; block 1 for i.i.d. returns, block 5 for AR(1))\n\n")
set.seed(1234)
binse <- function(p, R) sqrt(p * (1 - p)/R)
sizeSummary <- function(p) {
  ok <- is.finite(p)
  n <- sum(ok)
  rate <- if (n > 0L) mean(p[ok] < 0.05) else NA_real_
  c(rate = rate,
    se = if (n > 0L) binse(rate, n) else NA_real_)
}
TT <- 120L
pairGauss <- function() cbind(stats::rnorm(TT, 0.006, 0.04),
                              stats::rnorm(TT, 0.006, 0.04))
pairT5    <- function() cbind(0.006 + 0.04 * stats::rt(TT, 5)/sqrt(5/3),
                              0.006 + 0.04 * stats::rt(TT, 5)/sqrt(5/3))
pairAR1   <- function(rho = 0.3) {
  sim <- function() as.numeric(stats::arima.sim(list(ar = rho), TT,
                                                sd = 0.04 * sqrt(1 - rho^2))) + 0.006
  cbind(sim(), sim())
}
# The t(5) design is a heavy-tail stress test outside the finite-eighth-moment
# condition of the moment-based asymptotics. The t(10) design satisfies it.
pairT10   <- function() cbind(0.006 + 0.04 * stats::rt(TT, 10)/sqrt(10/8),
                              0.006 + 0.04 * stats::rt(TT, 10)/sqrt(10/8))
pairs <- list("i.i.d. Gaussian" = pairGauss,
              "i.i.d. standardized t(5)" = pairT5,
              "AR(1), rho = 0.3" = pairAR1,
              "i.i.d. standardized t(10)" = pairT10)
bootBlock <- c("i.i.d. Gaussian" = 1L,
               "i.i.d. standardized t(5)" = 1L,
               "AR(1), rho = 0.3" = 5L,
               "i.i.d. standardized t(10)" = 1L)
cat(sprintf("    %-28s %5s %18s %18s\n", "return process", "block",
            "asymptotic (s.e.)", "bootstrap (s.e.)"))
for (nm in names(pairs)) {
  p1 <- p2 <- rep(NA_real_, RT)
  for (r in seq_len(RT)) {
    xy <- pairs[[nm]]()
    p1[r] <- msharpeTesting(xy[, 1], xy[, 2], level = 0.90)$pval
    p2[r] <- msharpeTesting(xy[, 1], xy[, 2], level = 0.90,
                            control = list(type = 2, bBoot = bootBlock[[nm]],
                                           nBoot = 199))$pval
  }
  s1 <- sizeSummary(p1)
  s2 <- sizeSummary(p2)
  cat(sprintf("    %-28s %5d %.3f (%.3f)      %.3f (%.3f)\n", nm,
              bootBlock[[nm]], s1["rate"], s1["se"],
              s2["rate"], s2["se"]))
}

cat("\nsettings: replications", RS, "(screening) and", RT,
    "(test size); seed 1234; fastAdjust = TRUE\n")

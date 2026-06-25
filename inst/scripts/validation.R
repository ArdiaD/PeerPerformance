## ---------------------------------------------------------------------------
## Monte-Carlo validation of the PeerPerformance estimators
##
## Reproducible checks that the implementation behaves as the methodology of
## Ardia and Boudt (2018, JBF) and Ardia and Boudt (2015, FRL) predicts:
##
##   (A) Under a null of equal performance, a naive percentile rank spuriously
##       labels ~half of each fund's peers as out-performed (pure luck). The
##       luck correction collapses this: the estimated out-performance ratio
##       pi+ is small and most of the mass is (correctly) assigned to the
##       equal-performance ratio pi0. This is the correction the package
##       delivers.
##
##   (B) Under a planted skill structure, the estimator is deliberately
##       conservative: only statistically detectable differences are flagged,
##       so undetected true differences are pooled into pi0 and the estimated
##       pi+/pi- are lower bounds on their population values. It errs toward
##       "equal performance" rather than toward false discoveries.
##
##   (C) The modified Sharpe ratio equality test has approximately correct
##       size under the null, and power against a genuine difference.
##
## Run with:
##   source(system.file("scripts", "validation.R", package = "PeerPerformance"))
## Increase the replication counts below for publication-quality figures; the
## defaults are chosen to run in a couple of minutes on a laptop.
## ---------------------------------------------------------------------------

library("PeerPerformance")
set.seed(1234)

ctr <- list(nCore = 1)          # single core for reproducibility

## ===========================================================================
## (A) Equal performance: pi0 should be ~1, naive rank ~0.5 outperformance
## ===========================================================================
N  <- 30L                       # funds
TT <- 120L                      # months
RA <- 100L                      # replications

pi0A <- posA <- rankA <- numeric(RA)
for (r in seq_len(RA)) {
  X  <- matrix(stats::rnorm(TT * N, mean = 0.005, sd = 0.04), TT, N)
  sc <- alphaScreening(X, control = ctr)
  pi0A[r]  <- mean(sc$pizero, na.rm = TRUE)
  posA[r]  <- mean(sc$pipos,  na.rm = TRUE)
  ## naive percentile rank: fraction of peers each fund beats on mean return,
  ## averaged over funds (ignores statistical significance)
  mu <- colMeans(X)
  rankA[r] <- mean(sapply(seq_len(N), function(i) mean(mu[i] > mu[-i])))
}

cat("\n(A) Equal-performance null (true pi0 = 1, true out-performance = 0)\n")
cat(sprintf("    naive rank out-performance   : %.3f  <- spurious 'skill' (pure luck)\n",
            mean(rankA)))
cat(sprintf("    luck-corrected pi+           : %.3f  <- collapses toward 0\n",
            mean(posA)))
cat(sprintf("    luck-corrected pi0           : %.3f  (sd %.3f) <- mass restored to 'equal'\n",
            mean(pi0A), stats::sd(pi0A)))

## ===========================================================================
## (B) Planted skill: one third skilled, one third neutral, one third poor
## ===========================================================================
RB <- 100L
mu_true <- c(rep(0.012, N %/% 3), rep(0.005, N %/% 3),
             rep(-0.002, N - 2L * (N %/% 3)))
## true population ratios (per fund), from the planted means (ties -> pi0)
truePos <- mean(sapply(seq_len(N), function(i) mean(mu_true[i] > mu_true[-i])))
trueNeg <- mean(sapply(seq_len(N), function(i) mean(mu_true[i] < mu_true[-i])))

posB <- negB <- numeric(RB)
for (r in seq_len(RB)) {
  X  <- sapply(mu_true, function(m) stats::rnorm(TT, mean = m, sd = 0.04))
  sc <- alphaScreening(X, control = ctr)
  posB[r] <- mean(sc$pipos, na.rm = TRUE)
  negB[r] <- mean(sc$pineg, na.rm = TRUE)
}

cat("\n(B) Planted skill structure (universe-average ratios; estimates are\n",
    "    conservative lower bounds -- undetected differences fall into pi0)\n", sep = "")
cat(sprintf("    out-performance pi+ : estimated %.3f  vs true %.3f\n",
            mean(posB), truePos))
cat(sprintf("    under-performance pi-: estimated %.3f  vs true %.3f\n",
            mean(negB), trueNeg))

## ===========================================================================
## (C) Modified Sharpe equality test: size (null) and power (alternative)
## ===========================================================================
RC    <- 500L
alpha <- 0.05

## size: two funds drawn from the SAME distribution -> H0 true
rej_size <- 0L
for (r in seq_len(RC)) {
  a <- stats::rnorm(TT, 0.006, 0.04)
  b <- stats::rnorm(TT, 0.006, 0.04)
  rej_size <- rej_size + (msharpeTesting(a, b, level = 0.90)$pval < alpha)
}

## power: clearly different risk/return -> H0 false
rej_pow <- 0L
for (r in seq_len(RC)) {
  a <- stats::rnorm(TT, 0.010, 0.030)
  b <- stats::rnorm(TT, 0.002, 0.055)
  rej_pow <- rej_pow + (msharpeTesting(a, b, level = 0.90)$pval < alpha)
}

cat("\n(C) Modified Sharpe equality test (asymptotic, nominal 5%)\n")
cat(sprintf("    empirical size (H0 true) : %.3f  (target %.3f)\n",
            rej_size / RC, alpha))
cat(sprintf("    empirical power (H0 false): %.3f\n", rej_pow / RC))
cat("\nNote: the asymptotic test can be mildly oversized in short samples;\n",
    "use control = list(type = 2) for the studentized bootstrap.\n", sep = "")

## ===========================================================================
## Compact summary with Monte-Carlo standard errors (for reporting)
## ===========================================================================
binse <- function(p, R) sqrt(p * (1 - p) / R)   # binomial s.e. for a rate
cat("\n--- summary (estimate (MC s.e.)) ----------------------------------\n")
cat(sprintf("(A) naive-rank pi+        : %.3f (%.3f)\n", mean(rankA), stats::sd(rankA)/sqrt(RA)))
cat(sprintf("(A) luck-corrected pi+    : %.3f (%.3f)   [true 0]\n",     mean(posA), stats::sd(posA)/sqrt(RA)))
cat(sprintf("(A) luck-corrected pi0    : %.3f (%.3f)   [true 1]\n",     mean(pi0A), stats::sd(pi0A)/sqrt(RA)))
cat(sprintf("(B) pi+ (planted skill)   : %.3f (%.3f)   [true %.3f]\n",  mean(posB), stats::sd(posB)/sqrt(RB), truePos))
cat(sprintf("(B) pi- (planted skill)   : %.3f (%.3f)   [true %.3f]\n",  mean(negB), stats::sd(negB)/sqrt(RB), trueNeg))
cat(sprintf("(C) mSR test size         : %.3f (%.3f)   [target 0.05]\n", rej_size/RC, binse(rej_size/RC, RC)))
cat(sprintf("(C) mSR test power        : %.3f (%.3f)\n",                 rej_pow/RC,  binse(rej_pow/RC, RC)))
cat(sprintf("settings: N=%d, T=%d, reps A/B=%d, reps C=%d, seed=1234\n", N, TT, RA, RC))

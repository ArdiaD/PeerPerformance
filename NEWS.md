# Version 2.4.0 (DA)
- New optional `control$fastAdjust` (default `FALSE`): the truncated-normal
  bias correction of `pi0` inverts its monotone map for the whole vector at
  once (vectorised bisection) instead of one `uniroot()` call per value. This
  is the dominant cost of a screening with a data-driven `lambda` -- a
  100-fund `alphaScreening()` drops from about 19s to about 4.5s (4.2x). The
  bisection locates the root to about 1e-12; because `uniroot` stops at its own
  tolerance (about 1.2e-4), the two paths typically differ by a few 1e-5, the
  fast path being the more accurate. The original code path remains the default
  so that published results reproduce exactly. Bootstrap indices are now built
  only for the bootstrap test (`type = 2`), so the asymptotic path no longer
  depends on `bBoot`
- Bootstrap screening on unbalanced panels fixed: the bootstrap indices are now
  pre-generated in the master for every distinct pairwise complete-case length
  (previously a single full-length index matrix was remapped by modulo inside
  the workers, which biased the iid resampling and broke the circular-block
  structure for `bBoot > 1`); pairs shorter than the block length are left
  untested; balanced panels are unaffected. All randomness stays in the master,
  so seeded results do not depend on `nCore`
- `nCore = 1` (the default) now runs serially without creating a PSOCK
  cluster, removing the per-call cluster overhead (noticeable in
  `rollScreening`'s window loop); results are identical to the cluster path
- The symmetric bootstrap p-value of the Sharpe test now uses `>=`, consistent
  with the modified Sharpe test (ties have probability zero for continuous
  returns)
- Documentation: HAC availability in `alphaScreening` clarified; the
  `gammaPos`/`gammaNeg` counting rules are stated explicitly; misleading
  "bootstrap and HAC" example headers fixed (`hac` is ignored when
  `type = 2`); `summary` documents that win/loss counts are within-group only
- Robustness (pre-submission audit): PSOCK clusters are now closed with
  `on.exit()` so workers are not leaked on error; `processControl` requires the
  count-like controls (`nBoot`, `bBoot`, `nCore`, `minObs`, `minObsPi`) to be
  whole numbers and rejects a bootstrap block length exceeding the sample size;
  within-group screening and cross-group screening now stop with a clear
  message on degenerate inputs (a single fund, an empty peer group)
- `alphaScreening`/`alphaTesting` now enforce `minObs` on the factor
  complete-case sample (factor `NA`s were previously ignored), and
  `alphaTesting(screen_beta = TRUE, hac = TRUE)` now returns the `alpha`
  component as a coefficient-by-fund matrix, consistent with the non-HAC path
  (the `print` method reports the alpha row)
- `confint` is now tested for all three ratios (`pipos`/`pizero`/`pineg`)
- Added a `confint` method for `SCREENING` objects: nonparametric peer
  (pairwise) bootstrap confidence intervals for the peer performance ratios
  (`pipos`/`pizero`/`pineg`)
- Added a reproducible Monte-Carlo validation script
  (`system.file("scripts", "validation.R", package = "PeerPerformance")`):
  checks the near-unbiasedness of `pizero` under the equal-performance null and
  the size/power of the modified Sharpe equality test
- Added a package vignette ("Luck-Corrected Peer Performance Analysis with
  PeerPerformance"), a `pkgdown` configuration, and a package `CITATION` entry
- Robustness: `processControl` now validates that scalar control values are
  single finite numbers / logicals; `computePi` checks the `lambda` length;
  `targetPeerPerformance` rejects non-whole `funds`; `rollScreening` validates
  `by` and only flags `screen_beta` for the alpha screen with factors; output
  fund names are preserved in `as.data.frame` and `targetPeerPerformance`
- Fixed a bug in `sharpeScreening`/`msharpeScreening` on unbalanced panels: the
  focal fund's returns were indexed with the first peer's missing-value mask
  (`X[idx[, k], 1]`) instead of the current pair's (`X[idx[, k], k]`), which
  could inject `NA`s and yield `NaN` p-values for some pairs (reported by
  GitHub user NenoJo)
- Added `targetPeerPerformance()` (contributed by Murilo Andre Peres Pereira): screens a
  selected subset of funds against the whole universe; a convenience wrapper
  over the cross-group screening (`*Screening(X[, funds], Y = X)`)
- Added a `summary` method for `SCREENING` objects (contributed by Murilo Andre Peres Pereira):
  ranked table, distribution of the measure, win/loss counts and top funds
- Added `rollScreening()`: rolling-window screening returning the time series of
  cross-sectionally averaged ratios (per factor when `screen_beta = TRUE`), with
  a `plot` method -- the dynamic design of Ardia et al. (2022, 2023)
- `plot` on a cross-group screening omits the (within-group) percentile-rank
  diagonal; a single focal fund is shown as one stacked bar
- `screen_beta` can now also be set through `control` (e.g.
  `control = list(screen_beta = TRUE)`); the function argument still works and
  takes precedence
- Cross-group screening: `alphaScreening`, `sharpeScreening` and
  `msharpeScreening` gain a `Y` argument to screen each fund in `X` against a
  second peer group `Y` (a single focal fund versus a group is `X` a vector);
  columns of `Y` identical to the focal fund are excluded automatically
- Added `as.data.frame` method for `SCREENING` objects (tidy, one row per fund,
  or per fund/coefficient with `screen_beta = TRUE`)
- `screen_beta = TRUE` output now labels the coefficient rows (alpha + factor
  names), and `exposureHeterogeneity()` aggregates them into the factor
  exposure heterogeneity measure of Ardia et al. (2023), with a `plot` method
- Added `print` methods for the `TESTING` and `SCREENING` objects, and a
  `plot` method for the `SCREENING` object that reproduces the peer performance
  screening plot of Ardia and Boudt (2018)
- `gammaPos` and `gammaNeg` (default 0.4 and 0.6) are now exposed in the `control`
  list of the screening functions, controlling the one-sided thresholds used for
  the out- and underperformance counts
- Fixed the default number of bootstrap replications (`nBoot = 499`) when an empty
  `control` list is supplied
- Fixed the VAR(1) data-generating process in the optimal block-length routines
  (`sharpeBlockSize`, `msharpeBlockSize`): the lagged cross term in the second
  equation now uses the correct series
- `sharpe()` now counts observations with `is.finite()`, consistent with the
  other moment computations, when `NA`/`NaN` are present
- Moved run-time dependencies from `Depends` to `Imports`
- Documentation fixes (Sharpe 1994 reference, `alphaTesting` return values)

# Version 2.3.2 (DA)
- Several fixes in documentation and good practices
  
# Version 2.3.1 (DA)
- Doc fixed

# Version 2.3.0 (DA,SL)
- Doc fixed
- alphaScreening now also outputs the betas
- Bug in lambda resampling fixed

# Version 2.2.5 (DA)
- Doc fixed

# Version 2.2.3 (DA)
- url fixed

# Version 2.2.1 (DA)
- Switch to parallel package
- References updated

# Version 2.1.4 (DA)
- Small fix in counting NA

# Version 2.1.3 (DA)
- Documentation fixes

# Version 2.1.2 (DA)
- Documentation fixes

# Version 2.1.1 (DA)
- Documentation fixes
- First CRAN release

# Version 2.1.0 (DA)
- New PeerPerformance documentation
- Compiler imported directly within function

# Version 2.0.11 (DA)
- Examples added
- Block length not exported anymore

# Version 2.0.10 (DA)
- Roxygen documentation
- testthat added
- Format of code

# Version 2.0.9 (DA)
- Update CITATION and DESCRIPTION

# Version 2.0.8 (DA)
- Small improvements with compiler
- Citations updated

# Version 2.0.6 and 2.0.7 (DA)
- Various improvements
- Small fix in documentation

# Version 2.0.5 (DA)
- Fix in documentation for modified Sharpe testing
- Small fix in pvalue computation by bootstrap (symmetric)

# Version 2.0.4 (DA)
- Bug fix for alpha screening when NA are in the dataset
- Risk-free rate removed
- tstat used for attribution

# Version 2.0.3 (DA)
- Major functions contain risk-free rates (zero by default)
- Documentation updated
- Adjustment factor robustified

# Version 2.0.2 (DA)
- alphaScreening fixed
- alphaScreening now encompasses hac estimation with sandwich and lmtest
- citation file updated

# Version 2.0.1 (DA)
- New package's name
- New package's number 
- New package's structure

# Version 1-00.15 (DA)
- pi+ fixed
- Default settings for lambda = NULL

# Version 1-00.14 (DA)
- Fix of errors in examples

# Version 1-00.13 (DA)
- Control parameters for lambda data driven (NULL)
- Documentation updated
- Function for optimal lambda corrected and enhanced
- Function pizero and pi corrected

# Version 1-00.01 (DA)
- First release
- Package includes (parallel) alpha and sharpe screening algorithms

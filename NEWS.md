# Version 2.4.0 (DA)
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

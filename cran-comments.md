## Submission summary

This is a bug-fix release. Nothing here adds features. It carries two sets of
corrections: the defect that originally motivated the version, and a second set
found in a code audit while the release was still unsubmitted.

Two of the fixes are user-visible, so they are called out individually below.

## The defect that motivated the release

`bootIndices()` generated the circular block bootstrap by drawing
`floor(T / bBoot)` blocks. When the block length did not divide the number of
concordant observations, the last `T %% bBoot` indices of every bootstrap
sample kept their initial value of zero, and because the callers remap indices
with `1 + ids %% T`, each of those zeros selected the first observation. At
`T = 50` and `bBoot = 6` this over-sampled observation 1 by a factor of about
three and biased the bootstrap null distribution, with no error, warning or
`NA` to signal it. The routine now draws `ceiling(T / bBoot)` blocks and keeps
the first `T` indices.

The default is `bBoot = 1` (i.i.d. resampling), which was never affected, and
results are unchanged whenever `bBoot` divides `T`.

A fund with no usable observation also no longer aborts a screening:
`infoFund()` called `lm()` on every column, so an all-`NA` fund stopped the run
with "0 (non-NA) cases" even though the pairwise stage handles such a column
correctly. Its summary statistics are now returned as `NA`. Relatedly,
non-finite returns are now normalised to `NA` on entry, so an `Inf` no longer
counts as an observation and then fails the fit.

## The two user-visible changes from the audit

**Studentized bootstrap centring.** `sharpeTesting()` and `msharpeTesting()`
replaced the observed difference by the difference minus the hypothesised value
and then centred the bootstrap draws on that quantity, so the reference
distribution was shifted rather than centred on the estimate. The test
statistic itself was correct, so the default `d = 0` was never affected. The
automatic block-length selectors, however, pass the observed difference as the
hypothesised value, so their calibration computed every rejection rate from a
mis-centred distribution. Users who called `bBoot = 0` will see different block
lengths, and therefore different p-values; `bBoot = 1`, the default, is
unchanged.

**Returned field semantics.** With the default `ttype = 2` the screening
functions stored the cross-multiplied statistic in `dsharpe` and `dmsharpe`,
fields named as differences, while the pairwise testing functions returned the
raw ratio difference — so the same field name meant two different things. The
screenings now store raw ratio differences as well. The sign is unchanged, so
the ratios and win/loss counts these feed are unaffected; only the reported
differential changes scale.

## The narrower audit fixes

The true-null bias correction is now scaled by the number of valid comparisons
rather than by the number of columns, so appending peer columns that are
entirely missing no longer changes the estimate. `minObsPi` is enforced on each
focal fund's final count of valid peers in the cross-group path as well as the
within-group one. And `alphaTesting()` fits both reported alphas on the same
concordant sample as the difference it tests, so the two reported alphas now
differ by the reported `dalpha`.

`NEWS.md` lists all of the above. All fixes have regression tests.

## Test environments

* Local: macOS 27.0.1 (aarch64-apple-darwin23), R 4.6.1

## R CMD check results

0 errors | 0 warnings | 0 notes

Checked with `--as-cran --run-donttest` from outside a synced folder:
`Status: OK`.

## Reverse dependencies

None on CRAN.

## Submission summary

This is a bug-fix release, submitted shortly after 2.4.0. It carries two sets of
corrections: the one that motivated the version, and a second found in a code
audit while the release was still unsubmitted. Nothing here adds features.

The second set contains one correctness fix that changes results, and one that
changes the meaning of a returned field, so both are called out below.

**Studentized bootstrap centring.** `sharpeTesting()` and `msharpeTesting()`
replaced the observed difference by the difference minus the hypothesised value
and then centred the bootstrap draws on that quantity, so the reference
distribution was shifted rather than centred on the estimate. The test
statistic was correct, so the default `d = 0` was never affected. The automatic
block-length selectors, however, pass the observed difference as the
hypothesised value, so their calibration computed every rejection rate from a
mis-centred distribution and `bBoot = 0` returned an essentially arbitrary
block length. Users who called `bBoot = 0` will see different block lengths and
therefore different p-values; `bBoot = 1`, the default, is unchanged.

**Returned field semantics.** With the default `ttype = 2` the screening
functions stored the cross-multiplied statistic in `dsharpe`/`dmsharpe`, fields
named as differences, while the pairwise tests returned the raw ratio
difference. The screenings now store raw ratio differences too. The sign is
unchanged, so the ratios and counts these feed are unaffected; only the
reported differential changes scale.

The remaining audit fixes are narrower: the true-null bias correction is now
scaled by the number of valid comparisons rather than the column count, so
appending entirely missing peers no longer moves the estimate; `minObsPi` is
enforced in the cross-group path; and `alphaTesting()` fits both reported alphas
on the same sample as the difference it tests. NEWS.md lists all of them.

## The defect that motivated the release Apologies for the
quick turnaround: 2.4.0 introduced no new defect, but a pre-existing one was
found immediately after its acceptance and it silently biases a statistical
result, so it seemed better to correct it than to leave it in place.

`bootIndices()` generated the circular block bootstrap by drawing
`floor(T / bBoot)` blocks. When the block length did not divide the number of
concordant observations, the last `T %% bBoot` indices of every bootstrap
sample kept their initial value of zero, and because the callers remap indices
with `1 + ids %% T` each of those zeros selected the first observation. At
`T = 50` and `bBoot = 6` this over-sampled observation 1 by a factor of about
three and biased the bootstrap null distribution, with no error, warning or
`NA` to signal it. The routine now draws `ceiling(T / bBoot)` blocks and keeps
the first `T` indices.

The default is `bBoot = 1` (i.i.d. resampling), which was never affected, and
results are unchanged whenever `bBoot` divides `T`.

A second fix is included: a fund with no usable observation no longer aborts a
screening (`infoFund()` called `lm()` on every column, so an all-`NA` fund
stopped the run with "0 (non-NA) cases"); its summary statistics are now
returned as `NA`.

Both fixes have regression tests.

## Test environments

* Local: macOS 26.5 (aarch64-apple-darwin20), R 4.5.2

## R CMD check results

0 errors | 0 warnings | 2 notes

```
* checking CRAN incoming feasibility ... NOTE
Days since last update: 1
```
Explained above.

```
* checking for future file timestamps ... NOTE
unable to verify current time
```
Local clock/network issue, not a package issue.

## Reverse dependencies

None on CRAN.

## Test environments

macOS 15 (arm64), R 4.6.1: `R CMD check --as-cran --run-donttest`
returns **Status: OK** -- 0 errors, 0 warnings, 0 notes.

# PeerPerformance: Luck-Corrected Peer Performance Analysis in R

`PeerPerformance` is an R package for the peer-performance evaluation of financial investments with
luck-correction. In particular, it implements the peer performance ratios 
of [Ardia and Boudt (2018)](https://doi.org/10.1016/j.jbankfin.2017.10.014) which measure the percentage of peers a focal fund outperforms and underperforms, after
correction for luck. It is useful for fund or portfolio managers to 
benchmark their investments or screen a universe of new funds. 
In addition, it implements the testing framework for the Sharpe and modified Sharpe ratios, described 
in [Ledoit and Wolf (2008)](https://doi.org/10.1016/j.jempfin.2008.03.002) 
and [Ardia and Boudt (2015)](https://doi.org/10.1016/j.frl.2015.02.008). See also Ardia et al. (2022,2023) for applications in sustainable finance.

## Installation

The stable version is on [CRAN](https://CRAN.R-project.org/package=PeerPerformance):

```r
install.packages("PeerPerformance")
```

The development version can be installed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("ArdiaD/PeerPerformance")
```

## Features

- **Peer performance screening** of a fund universe with luck correction:
  `alphaScreening()` (risk-adjusted alphas, optionally with factor exposures
  via `screen_beta = TRUE`), `sharpeScreening()` and `msharpeScreening()`
  (Sharpe / modified Sharpe). Each returns the out-/equal-/under-performance
  ratios (pi+, pi0, pi-).
- **Cross-group screening**: the `Y` argument screens each fund (or a single
  focal fund) against a *separate* peer group; `targetPeerPerformance()` is a
  convenience wrapper for screening a chosen subset against the whole universe.
- **Pairwise testing**: `alphaTesting()`, `sharpeTesting()`, `msharpeTesting()`.
- **Methods** for screening results: `print()`, `summary()`, `plot()` (the
  Ardia and Boudt 2018 screening plot), `confint()` (bootstrap confidence
  intervals for the ratios), and `as.data.frame()` (tidy output).
- **Dynamic and factor analyses**: `rollScreening()` (rolling-window ratios)
  and `exposureHeterogeneity()` (factor exposure heterogeneity of Ardia et al.
  2023).
- A vignette (`vignette("PeerPerformance")`) and a reproducible Monte-Carlo
  validation script
  (`system.file("scripts", "validation.R", package = "PeerPerformance")`).

## Quick start

```r
library("PeerPerformance")
data("hfdata")

## screen a universe of funds, luck-corrected
sc <- alphaScreening(hfdata[, 1:30], control = list(nCore = 1))
summary(sc)                    # ranked table with win/loss counts
plot(sc)                       # peer performance screening plot
confint(sc, parm = "pipos")    # bootstrap CIs for the outperformance ratios
```

## Please cite the package in publications!

By using `PeerPerformance` you agree to the following rules: 

1) You must cite [Ardia and Boudt (2018)](https://doi.org/10.1016/j.jbankfin.2017.10.014) in working papers and published papers that use `PeerPerformance`.
2) You must place the following URL in a footnote to help others find `PeerPerformance`: [https://CRAN.R-project.org/package=PeerPerformance](https://CRAN.R-project.org/package=PeerPerformance) 
3) You assume all risk for the use of `PeerPerformance`.

Ardia, D., Boudt, K. (2018).      
The peer performance ratios of hedge funds.      
_Journal of Banking and Finance_, 87, 351-368.    
[https://doi.org/10.1016/j.jbankfin.2017.10.014](https://doi.org/10.1016/j.jbankfin.2017.10.014)  
[https://doi.org/10.2139/ssrn.2000901](https://doi.org/10.2139/ssrn.2000901)

## Other references

Ardia, D., Boudt, K. (2015).
Testing equality of modified Sharpe ratios.  
_Finance Research Letters_, 13, 97-104.   
[https://doi.org/10.1016/j.frl.2015.02.008](https://doi.org/10.1016/j.frl.2015.02.008)   
[https://doi.org/10.2139/ssrn.2516591](https://doi.org/10.2139/ssrn.2516591)

Ardia, D., Bluteau, K., Tran, D. (2022).
How easy is it for investment managers to deploy their talent in green and brown stocks?
_Finance Research Letters_, 48, 102992.
[https://doi.org/10.1016/j.frl.2022.102992](https://doi.org/10.1016/j.frl.2022.102992)   
[https://doi.org/10.2139/ssrn.4009286](https://doi.org/10.2139/ssrn.4009286)

Ardia, D., Bluteau, K., Lortie-Cloutier, G., Tran, D. (2023).
Factor exposure heterogeneity in green and brown stocks.
_Finance Research Letters_, 55, Part A, pp.103900.
[https://doi.org/10.1016/j.frl.2023.103900](https://doi.org/10.1016/j.frl.2023.103900)   
[https://doi.org/10.2139/ssrn.4362696](https://doi.org/10.2139/ssrn.4362696)

Ledoit, O., Wolf, M. (2008).   
Robust performance hypothesis testing with the Sharpe ratio.    
_Journal of Empirical Finance_, 15(5), 850-859.

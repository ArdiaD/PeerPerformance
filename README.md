---
title: 'PeerPerformance: Peer performance analysis in R'
date: "13 January 2017"
tags:
- performance
- peer
- alpha
- Sharpe
- modified Sharpe
- alpha
- screening
- false discovery rate
authors:
- affiliation: Institute of Financial Analysis - University of Neuchâtel
  name: David Ardia
  orcid: 0000-0003-2823-782X
- affiliation: Solvay Business School - Vrije Universiteit Brussel
  name: Kris Boudt
  orcid: null
---

[![CRAN](http://www.r-pkg.org/badges/version/PeerPerformance)](https://cran.r-project.org/package=PeerPerformance) 
[![Downloads](http://cranlogs.r-pkg.org/badges/PeerPerformance?color=brightgreen)](http://www.r-pkg.org/pkg/PeerPerformance)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/PeerPerformance?color=brightgreen)](http://www.r-pkg.org/pkg/PeerPerformance)

# Introduction
`PeerPerformance` is an R package for the peer-performance evaluation of financial investments with
luck-correction. In particular, it implements the peer performance ratios 
of [Ardia and Boudt (2016)](http://dx.doi.org/10.2139/ssrn.2000901) which measure the percentage of peers a focal fund outperforms and underperforms, after
correction for luck. It is useful for fund or portfolio managers to 
benchmark their investments or screen a universe of new funds. 
In addition, it implements the testing framework for the Sharpe and modified Sharpe ratios, described 
in [Ledoit and Wolf (2008)](http://dx.doi.org/10.1016/j.jempfin.2008.03.002) 
and [Ardia and Boudt (2015)](http://dx.doi.org/10.1016/j.frl.2015.02.008).

# Installation
To install the package, run the following commands in R:

R> install.packages("devtools")

R> devtools::install_github("ArdiaD/PeerPerformance", dependencies = TRUE)

Then check the help of the various files, and run the examples:

R> library("PeerPerformance")

R> ?PeerPerformance
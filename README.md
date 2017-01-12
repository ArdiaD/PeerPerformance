---
title: 'PeerPerformance: Peer performance analysis in R'
date: "09 January 2017"
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

# Introduction
`PeerPerformance` is an R package for the peer-performance evaluation with
luck-correction. In particular, it implements the peer performance ratios of Ardia and Boudt
(2016) which measure the percentage of peers a focal fund outperforms and underperforms, after
correction for luck. In addition, it implements the testing framework for the Sharpe and modified Sharpe ratios.

# Installation
To install the package, run the following commands in R:

R> install.packages("devtools")

R> devtools::install_github("ArdiaD/PeerPerformance", dependencies = TRUE)

Then check the help of the various files, and run the examples:

R> library("PeerPerformance")

R> ?PeerPerformance
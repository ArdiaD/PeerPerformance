---
title: 'PeerPerformance: Luck--corrected peer performance analysis in R'
date: "2 February 2017"
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

# PeerPerformance
`PeerPerformance` (Ardia and Boudt, 2017) is an R package for the peer-performance evaluation of financial investments with
luck-correction. In particular, it implements the peer performance ratios 
of [Ardia and Boudt (2016)](http://dx.doi.org/10.2139/ssrn.2000901) which measure the percentage of peers a focal fund outperforms and underperforms, after
correction for luck. It is useful for fund or portfolio managers to 
benchmark their investments or screen a universe of new funds. 
In addition, it implements the testing framework for the Sharpe and modified Sharpe ratios, described 
in [Ledoit and Wolf (2008)](http://dx.doi.org/10.1016/j.jempfin.2008.03.002) 
and [Ardia and Boudt (2015)](http://dx.doi.org/10.1016/j.frl.2015.02.008).

# Installation
To install the latest stable version of `PeerPerformance`, run the following commands in R:

    R> install.packages("PeerPerformance")

To install the development version of `PeerPerformance`, run the following commands in R:

    R> install.packages("devtools")

    R> devtools::install_github("ArdiaD/PeerPerformance")

Then check the help of the various files and run the examples:

    R> library("PeerPerformance")

    R> ?PeerPerformance
    
    
Please cite `PeerPerformance` in publications:

    R> citation("PeerPerformance")


# References

Ardia, D., Boudt, K. (2015)  
Testing equality of modified Sharpe ratios.  
_Finance Research Letters_ **13**, pp.97--104.   
http://dx.doi.org/10.1016/j.frl.2015.02.008

Ardia, D., Boudt, K. (2016).    
_The peer performance ratios of hedge funds_.    
Working paper.  
http://dx.doi.org/10.2139/ssrn.2000901

Ardia, D., Boudt, K. (2017).    
_PeerPerformance: Peer performance analysis in R_.    
R package.   
https://github.com/ArdiaD/PeerPerformance  

Ledoit, O., Wolf, M. (2008).   
Robust performance hypothesis testing with the Sharpe ratio.    
_Journal of Empirical Finance_ **15**(5), pp.850--859.  
http://dx.doi.org/10.1016/j.jempfin.2008.03.002  

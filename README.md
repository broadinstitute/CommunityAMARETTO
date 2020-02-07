[//]: # (TODO: Bioconductor support?)
[//]: # (TODO: Some examples)

<p align="center">
  <a href="https://github.com/broadinstitute/CommunityAMARETTO/">
    <img height="400" src="https://github.com/broadinstitute/CommunityAMARETTO/blob/develop/inst/extdata/CommunityAMARETTO_logo.png">
  </a>
  <h1 align="center"></h1>
</p>


[![Version](https://img.shields.io/badge/version-0.99.1-lightgrey.svg)]()

The goal of the CommunityAMARETTO algorithm (Champion et al., EBioMedicine 2018) is to identify cell circuits and their drivers that are shared and distinct across biological systems. Specifically, Community-AMARETTO takes as input multiple regulatory networks inferred using the AMARETTO algorithm that are based on multi-omics and imaging data fusion. Next, Community-AMARETTO learns communities or subnetworks, in particular, regulatory modules comprising of cell circuits and their drivers, that are shared and distinct across multiple regulatory networks derived from multiple cohorts, diseases, or biological systems more generally, using the Girvan-Newman "edge betweenness community detection" algorithm (Girvan and Newman, Physical Review E. 2004).

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Running CommunityAMARETTO](#running communityamaretto)
- [References](#references)

## Introduction

Uncovering gene regulatory mechanisms underlying diseases and cancer
systems have long time been sought by many researchers. 
This led to the development of many novel computational algorithms for
regulatory network inference using multiomics such as genetics,
epigenetics and transcriptomics. These algorithms are designed to infer
system-specific gene regulatory network, leading to the understanding
of underlying regulatory mechanisms of the system and facilitating
biomarker discovery.

It is also believed that the integration of regulatory networks across
multiple systems provides key information about cross-systems shared
and distinct mechanisms. To this end, we developed *Community-AMARETTO*
algorithm to integrate multiple regulatory networks inferred by the
AMARETTO algorithm [1]. Specifically, Community-AMARETTO algorithm consists of
1) constructing a master network composed of multiple regulatory networks
followed by 2) detecting groups (communities) of circuits that are shared
across systems as well as highliting circuits that are system-specific
and distinct. 

## Installation

Install from the GitHub repository using devtools:

``` r
library(devtools)
install_github("broadinstitute/CommunityAMARETTO")
```

## Running CommunityAMARETTO

* The vignettes contains an example R script for a typical AMARETTO analysis. Please try!

* Detailed information on `CommunityAMARETTO` package functions can be obtained in the help files. For example, to view the help file for the function `CommunityAMARETTO` in a R session, use `?CommunityAMARETTO`.

## References


1.	Gevaert, O., Villalobos, V., Sikic, B. I. & Plevritis, S. K. Identification of ovarian cancer driver genes by using module network integration of multi-omics data. Interface Focus 3, 20130013–20130013 (2013).
2.	Gevaert, O. MethylMix: an R package for identifying DNA methylation-driven genes. Bioinformatics 31, 1839–1841 (2015).
3. AMARETTO package in Bioconductor.

## Useful Links

GitHub:<br/>
AMARETTO: https://github.com/gevaertlab/AMARETTO<br/>
Community-AMARETTO: https://github.com/broadinstitute/CommunityAMARETTO

GenePattern: <under development><br/>
AMARETTO: https://beta.genepattern.org/gp/pages/index.jsf?lsid=urn:lsid:broad.mit.edu:cancer.software.genepattern.module.analysis:00378:0.52<br/>
Community-AMARETTO: https://beta.genepattern.org/gp/pages/index.jsf?lsid=urn:lsid:broad.mit.edu:cancer.software.genepattern.module.analysis:00380:999999999
  
Docker: <under development><br/>
AMARETTO: https://hub.docker.com/r/genepattern/docker-amaretto<br/>
Community-AMARETTO: https://hub.docker.com/r/genepattern/community-amaretto


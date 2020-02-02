[//]: # (TODO: Bioconductor support?)
[//]: # (TODO: Some examples)

<p align="center">
  <a href="https://github.com/broadinstitute/CommunityAMARETTO/">
    <img height="400" src="https://github.com/broadinstitute/CommunityAMARETTO/blob/develop/inst/extdata/CommunityAMARETTO_logo.png">
  </a>
  <h1 align="center"></h1>
</p>




# CommunityAMARETTO

The goal of the CommunityAMARETTO algorithm (Champion et al., EBioMedicine 2018) is to identify cell circuits and their drivers that are shared and distinct across biological systems. Specifically, Community-AMARETTO takes as input multiple regulatory networks inferred using the AMARETTO algorithm that are based on multi-omics and imaging data fusion. Next, Community-AMARETTO learns communities or subnetworks, in particular, regulatory modules comprising of cell circuits and their drivers, that are shared and distinct across multiple regulatory networks derived from multiple cohorts, diseases, or biological systems more generally, using the Girvan-Newman "edge betweenness community detection" algorithm (Girvan and Newman, Physical Review E. 2004).


## Installation

Install from the GitHub repository using devtools:

``` r
library(devtools)
install_github("broadinstitute/CommunityAMARETTO")
```

## Example

* The vignette provides a small case study, demonstrating how to implement the workflow of CommunityAMARETTO package in R. Please try!

``` r
## basic example code
```


## References

1.	Champion, M. et al. Module Analysis Captures Pancancer Genetically and Epigenetically Deregulated Cancer Driver Genes for Smoking and Antiviral Response. EBioMedicine 27, 156–166 (2018).
2.	Gevaert, O., Villalobos, V., Sikic, B. I. & Plevritis, S. K. Identification of ovarian cancer driver genes by using module network integration of multi-omics data. Interface Focus 3, 20130013–20130013 (2013).
3.	Gevaert, O. MethylMix: an R package for identifying DNA methylation-driven genes. Bioinformatics 31, 1839–1841 (2015).

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


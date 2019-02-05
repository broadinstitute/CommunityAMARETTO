# CommunityAMARETTO

The goal of the Community-AMARETTO algorithm (Champion et al., EBioMedicine 2018) is to identify cell circuits and their drivers that are shared and distinct across biological systems. Specifically, Community-AMARETTO takes as input multiple regulatory networks inferred using the AMARETTO algorithm that are based on multi-omics and imaging data fusion. Next, Community-AMARETTO learns communities or subnetworks, in particular, regulatory modules comprising of cell circuits and their drivers, that are shared and distinct across multiple regulatory networks derived from multiple cohorts, diseases, or biological systems more generally, using the Girvan-Newman "edge betweenness community detection" algorithm (Girvan and Newman, Physical Review E. 2004).

REFERENCES

Champion M., Brennan K., Croonenborghs T., Gentles A.J., Pochet N., Gevaert O. (2018) Module analysis captures pancancer genetically and epigenetically deregulated cancer driver genes for smoking and antiviral response. EBioMedicine, 27:156-166. https://www.ebiomedicine.com/article/S2352-3964(17)30472-3/fulltext

GitHub:
AMARETTO: https://github.com/gevaertlab/AMARETTO
Community-AMARETTO: https://github.com/broadinstitute/CommunityAMARETTO

GenePattern: <under development>
AMARETTO: https://beta.genepattern.org/gp/pages/index.jsf?lsid=urn:lsid:broad.mit.edu:cancer.software.genepattern.module.analysis:00378:0.52
Community-AMARETTO: https://beta.genepattern.org/gp/pages/index.jsf?lsid=urn:lsid:broad.mit.edu:cancer.software.genepattern.module.analysis:00380:999999999
  
Docker: <under development>
AMARETTO: https://hub.docker.com/r/genepattern/docker-amaretto
Community-AMARETTO: https://hub.docker.com/r/genepattern/community-amaretto

## Installation

You can install the released version of CommunityAMARETTO from github:

``` r
library(devtools)
install_github("CommunityAMARETTO")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code
```


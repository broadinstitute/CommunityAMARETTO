#### How to use clinical data and analyse them
library("RMySQL")
library("stringr")

source("postprocessing/ClinicalVariables/ClinicalVariables.R")

# first connect isis genome (ssh -L 3307:localhost:3306 mchampion@isis-genome in a terminal)
CancerSite = c("LAML")
AllClinicalData <- DownloadClinicalData(CancerSite) # this is a big matrix with 740 clinical data for CancerSite
          # be careful COADREAD is seperated in "COAD" and "READ"

# here you clean the clinical data
ClinicalData <- CleanClinicalData(ClinicalData=AllClinicalData) # remove variables with missing values and clean some variables
# then choose which variables you want and correlate with gene expression, methylation,...
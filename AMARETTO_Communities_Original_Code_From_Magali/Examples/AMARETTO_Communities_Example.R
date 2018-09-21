#### Compare the modules, create and analyse the communities ####
source("package/AMARETTO_DevelVersion.R")
source("Pancancer_Communities/AMARETTO_Communities.R")
source("Other_codes/EnrichmentAnalysis/EA_createData.R")
source("Other_codes/EnrichmentAnalysis/EA_run.R")

library("plyr")  
library("igraph")

AllCancers <- c('BLCA','BRCA','COADREAD','GBM','HNSC','KIRC','LAML','LUAD','LUSC','OV','UCEC') # all the cancer types for your pancancer analysis

# notes: before running AMARETTO_Pancancer, you have to run AMARETTO on all cancers
# if you want to create your own file, use CreatePancancerData function in AMARETTO_Communities
load("Pancancer_Communities/AMARETTOPancancer11.rda")

# this function creates a p-value matrix measuring the association between all modules form all cancers
AMARETTOPancancer_results <- AMARETTO_Pancancer(AMARETTOPancancerData=AMARETTOPancancer11)
save(AMARETTOPancancer_results,file=".......Rdata")

CommunitiesResults <- AMARETTO_CreateCommunities(AMARETTOPancancer_results)
save(CommunitiesResults,file="......Rdata")

AnalysedResults <- AMARETTO_AnalyseCommunities(AMARETTOPancancer_results,CommunitiesResults,All.Communities=FALSE,CommNumber=c(13,15,18),AMARETTOPancancerData=AMARETTOPancancer11)     
    # if All.Communities = FALSE, you have to specify which Community you want to analyse (CommNumber)
    # this function analyse the communities: overlapping drivers and gene set enrichment analysis
save(AnalysedResults,file=".....Rdata")

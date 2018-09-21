### Code example for running gene set enrichment analysis
library("GSA");
library("limma");

# you should first set your working directory
Directory = "Postprocessing/EnrichmentAnalysis/" # path for the data
source("PostProcessing/EnrichmentAnalysis/EA_createData.R")
source("Postprocessing/EnrichmentAnalysis/EA_run.R")

# first create your dictionary of gene sets
DataBase <- CreateDataBase() # create the data base for GA (from MSigDB e.g.)
                           
### Enrichment analysis of a list of genes ###
GeneList <- c("CEACAM5","PIGR","CEACAM6","LCN2","SPRR3","PADI1","CXCL17","ST6GALNAC1","SCNN1A","PRSS3","TMC5","NCCRP1","CRABP2",
              "MUC1","C19orf21","GPRC5A","FUT2","ACPP","ANKRD56","CEACAM1","RASEF","MFSD4","POU2F3","RAB27B","LRG1","GPR37L1",
              "PRSS22","PRSS8","PLAC4","CLDN4")      
Results = EnrichmentAnalysis(Data=GeneList,DataBase=DataBase)
  # Results$ResultsThreshold gives you the significant GSEA results only: one line corresponds to one gene set,
  #           pvalueAdjust is the adjusted p-value, p-value, the non-adjusted p-value, 
  #           Intersection the number of overlapping genes (between your gene set and all other gene sets)
  #           GeneSetSize the number of genes in the gene set
  # Results$All.Results, the same but no threshold (all p-value results)
  # DataSize: the number of genes in your gene list

### Enrichment analysis on one module of AMARETTO ###
Cancer <- "BLCA" 
load("PancancerCommunities/AMARETTOPancancer11.rda")
ModuleNr <-3
AMARETTOresults <- AMARETTOPancancer11[[which(names(AMARETTOPancancer11)==Cancer)]]
ModuleMembers = AMARETTOresults$AllGenes[AMARETTOresults$ModuleMembership==ModuleNr]
Results = EnrichmentAnalysis(Data=ModuleMembers,DataBase=DataBase)

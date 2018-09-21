#### How to run AMARETTO
source("package/AMARETTO_DevelVersion.R")
source('package/TCGAdownload_DevelVersion.R') 
source('package/TCGAprocess_DevelVersion.R')

library("RCurl")
library("RMySQL")
library("stringr")
library("limma")
library("matrixStats")
library("doParallel")
library("ComplexHeatmap")
library("circlize")
library("impute")
library("glmnet")

TargetDirectory="/Users/mchampion/Desktop/PancancerAnalysis/Data/" # where to put the data 

# choose your cancer in the list below
AllCancers <- c('ACC','BLCA','BRCA','CESC','CHOL','COAD','COADREAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC',
                'LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ','SARC','STAD','TCGT','THCA','THYM','UCEC','UCS','UVM')
CancerSite = "BLCA"

# Downloading the data (MA, CNV and MET needed to run AMARETTO)
DataSetDirectories = Download_CancerSite(CancerSite,TargetDirectory,downloadData=FALSE,AMARETTO=TRUE) # if you want to download the data, set downloadData = TRUE

# Preprocessing all data
# if AMARETTO package is not installed, you have to load Batch Data 
# if UseMethylMix = TRUE you also need to dowload MethylData (only available for 25 of the cancer sites)
load("package/BatchData.rda")
load("package/Pancancer25_MethylStates.rda")
ProcessedData = Preprocess_CancerSite(CancerSite,DataSetDirectories,AMARETTO=TRUE)  
save(ProcessedData,file=paste0(TargetDirectory,"ProcessedData",CancerSite,".RData"))

# Or load direclty your data
load("package/ProcessedData_LAML.rda")
# /!\ AMARETTO uses methylation data from MethylMix, so use data from Pancancer25_MethylStates.rda
# /!\ before running AMARETTO (in case you use your own data sets) be sure that all genes in CNV matrix and MET matrix have expression levels (are in MA matrix)!

# initialization step
AMARETTOinit <- AMARETTO_Initialize(MA_matrix = ProcessedData$MA_TCGA,CNV_matrix = ProcessedData$CNV_TCGA,MET_matrix = ProcessedData$MET_TCGA,NrModules = 100,VarPercentage=75)  
save(AMARETTOinit,file=paste0(TargetDirectory,"AMARETTOinit",CancerSite,".RData"))

# main code
AMARETTOresults <- AMARETTO_Run(AMARETTOinit)
save(AMARETTOresults,file=paste0(TargetDirectory,"AMARETTOresults",CancerSite,".RData"))

# visualize modules
ModuleNr <- 54 # define the module you want to visualize
AMARETTO_VisualizeModule(AMARETTOinit,AMARETTOresults,CNV_matrix = ProcessedData$CNV_TCGA,MET_matrix = ProcessedData$MET_TCGA,ModuleNr=ModuleNr)



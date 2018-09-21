### code to download and preprocess data ###
library("RCurl")
library("RMySQL")
library("stringr")
library("limma")
library("impute")

# where the codes are (you need to be in the good folder, define your working directory!)
source('package/TCGAdownload_DevelVersion.R') 
source('package/TCGAprocess_DevelVersion.R')

TargetDirectory="/Users/mchampion/Desktop/PancancerAnalysis/Data/" # where to put the data 

########################################################
# Processing data manually
########################################################

CancerSite = "BLCA" # to be chosen from 'ACC','BLCA','BRCA','CESC','CHOL','COAD','COADREAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ','SARC','STAD','TCGT','THCA','THYM','UCEC','UCS','UVM'
  
# Downloading all the the data (mRNA, miRNA, CNV, MET and Mutation)
DataSetDirectories = Download_CancerSite(CancerSite,TargetDirectory,downloadData=FALSE,UseMethylMix=TRUE) # if you want to download the data, set downloadData = TRUE
  
# Preprocessing all data
# if AMARETTO package is not installed, you have to load Batch Data 
# if UseMethylMix = TRUE you also need to dowload MethylData (only available for 25 of the cancer sites)
# if you miRNA data, you finally need the csv file annotationMatureID 
# everything is in bitbucket
load("package/BatchData.rda")
load("package/Pancancer25_MethylStates.rda")
miRNA_annotation = read.csv2(paste0(getwd(),"/package/annotationMatureID.csv"),sep=",",header=TRUE)
ProcessedData = Preprocess_CancerSite(CancerSite,DataSetDirectories)  
    
# save the data
SaveFile = paste(TargetDirectory,'TCGA_',CancerSite,"_ProcessedData.RData",sep='')
save(file=SaveFile,ProcessedData)
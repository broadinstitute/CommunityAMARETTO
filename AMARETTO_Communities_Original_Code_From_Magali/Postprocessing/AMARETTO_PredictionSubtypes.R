#### predict molecular subtypes ###
library("AMARETTO")
library("ROCR")
library("cvTools")
AllCancers <- c("BLCA","BRCA","COADREAD","GBM","HNSC","KIRC","LAML","LUAD","LUSC","OV","UCEC")
TargetDirectory="/Users/mchampion/Desktop/PancancerAnalysis/Data/" # where to put the data 

AUCs <- list()
MSEs <- list()
Models <- list()
for (i in (1:length(AllCancers))){
  CancerSite <- AllCancers[i]
  
  ######## download and load the subtype data #########
  SubtypeData <- AMARETTO_LoadSubtypes(CancerSite,TargetDirectory,downloadData=FALSE)
  Subtype <- read.csv("/Users/mchampion/Desktop/PancancerAnalysis/AMARETTO/SubtypeBRCA.csv")
  SubtypeData <- Subtype[,3]
  names(SubtypeData) <- paste0(Subtype[,1],"-01")
  SubtypeData <- as.character(SubtypeData)
  names(SubtypeData) <- paste0(Subtype[,1],"-01")
  SubtypeData[which(SubtypeData=="Basal")] <- rep(1,length(which(SubtypeData=="Basal")))
  SubtypeData[which(SubtypeData=="LumA")] <- rep(2,length(which(SubtypeData=="LumA")))
  SubtypeData[which(SubtypeData=="LumB")] <- rep(3,length(which(SubtypeData=="LumB")))
  SubtypeData[which(SubtypeData=="Normal")] <- rep(4,length(which(SubtypeData=="Normal"))) 
  SubtypeData[which(SubtypeData=="Her2")] <- rep(5,length(which(SubtypeData=="Her2")))
  SubtypeData <- as.numeric(SubtypeData)
  names(SubtypeData) <- paste0(Subtype[,1],"-01")
  
  ######### Load the gene expression data module ########
  load(paste0("/Users/mchampion/Desktop/PancancerAnalysis/AMARETTO/initstep",CancerSite,".RData"))
  MAData <- AMARETTOinit$MA_matrix_Var
  MAmoduleMatrix <- c()
  for (j in (1:length(unique(AMARETTOinit$Clusters)))){
    MAmodule <- MAData[which(AMARETTOinit$Clusters==j),]
    MAmodule <- colMeans(MAmodule)
    MAmoduleMatrix <- rbind(MAmoduleMatrix,MAmodule)
  } 
  rownames(MAmoduleMatrix) <- paste0(paste0("M",c(1:100),"_"),CancerSite)
  SubtypeData <- SubtypeData[!is.na(SubtypeData[,which(colnames(SubtypeData)=="mRNAseq_cNMF")]),]
  SubtypeData <- SubtypeData[intersect(rownames(SubtypeData),colnames(MAmoduleMatrix)),]
  MAmoduleMatrix <- MAmoduleMatrix[,intersect(rownames(SubtypeData),colnames(MAmoduleMatrix))]
  
  ####### predicting all subtypes ########
  SubtypesResults <- AMARETTO_PredictSubtypes(ModuleMatrix,SubtypeData,nFoldCV=10)
  
  AUCs <- c(AUCs,list(colMeans(SubtypesResults$AUC)))
  MSEs <- c(MSEs,list(colMeans(SubtypesResults$MSE)))
  Models <- c(Models,list(SubtypesResults$Models))
}
names(AUCs)=names(MSEs)=names(Models) <- AllCancers

### Analyse the results ###
Modules <- list() # list of modules used to predict each subtype of each cancer
MeanModules <- list() # number of modules used to predict each subtype of each cancer
SubtypesCancers <- list() # overlap of modules used to predict subtype
for (i in (1:length(AllCancers))){
  CancerSite <- AllCancers[i]
  model <- Models[[i]]
  ModuleFold <- list()
  Mean <- c()
  Subtype <- list()
  for (s in (1:length(model))){ 
    ModuleModel <- list()
    Number <- c()
    modelfold <- model[[s]]
    for (j in (1:nFoldCV)){
      modelshort <- modelfold[[j]][-1]
      names(modelshort) <- paste0("M",c(1:100),"_",CancerSite)
      moduleModel <- names(which(abs(modelshort)>0))
      number <- length(names(which(abs(modelshort)>0)))
      ModuleModel <- c(ModuleModel,list(moduleModel))
      Number <- c(Number,number)
    } 
    subtype <- table(unlist(ModuleModel))      
    subtype <- subtype[order(subtype,decreasing=TRUE)]
    Subtype <- c(Subtype,list(subtype))   
    mean <- mean(Number)
    ModuleFold <- c(ModuleFold,list(ModuleModel))
    Mean <- c(Mean,mean)    
  }
  MeanModules <- c(MeanModules,list(Mean))
  names(ModuleFold) <- paste0("subtype",c(1:3)) 
  names(Subtype) <- paste0("subtype",c(1:3))  
  SubtypesCancers <- c(SubtypesCancers,list(Subtype))
  Modules <- c(Modules,list(ModuleFold))
}
names(MeanModules) <- AllCancers
names(Modules) <- AllCancers
names(SubtypesCancers)<- AllCancers

### Store the results, including the regulators
load("/Users/mchampion/Desktop/PancancerAnalysis/AMARETTO/Regulators_perCancer")
Reg_cancers <- list()
for (i in (1:length(AllCancers))){
  CancerSite <- AllCancers[i]
  Model <- Modules[[i]]
  Reg_subtype <- list()
  for (s in (1:length(Model))){
    model <- Model[[s]]
    modules_list <- unique(unlist(model))
    regexp <- "[[:digit:]]+"
    modules_list_nbrs <- as.numeric(str_extract(modules_list,regexp))
    reg <- RegPerCancer[[which(names(RegPerCancer)==AllCancers[i])]]
    reg <- reg[,modules_list_nbrs]
    I <- which(rowSums(reg)==0)
    if (length(I)>0){
      reg <- reg[-I,]
    }
    colnames(reg) <- str_replace_all(colnames(reg),"_","")
    colnames(reg) <- paste0(colnames(reg),"_",CancerSite)
    mult <- table(unlist(model))
    mult <- mult[colnames(reg)]
    for (j in (1:ncol(reg))){
      reg[,j] <- reg[,j]*mult[j]
    }
    Total <- rowSums(reg)
    reg <- cbind(reg,Total)
    reg <- reg[order(reg[,ncol(reg)],decreasing=TRUE),]
    Reg_subtype <- c(Reg_subtype,list(reg))
  }
  names(Reg_subtype) <- paste0("subtype",c(1:3))  
  Reg_cancers <- c(Reg_cancers,list(Reg_subtype))
}
names(Reg_cancers) <- AllCancers


AMARETTO_PredictSubtypes <- function(ModuleMatrix,SubtypeData,nFoldCV){
  AllModels <- list()
  AUC <- matrix(0,nrow=nFoldCV,ncol=length(unique(SubtypeData$mRNAseq_cNMF)))
  MSE <- matrix(0,nrow=nFoldCV,ncol=length(unique(SubtypeData$mRNAseq_cNMF)))
  #ROC <- list()
  for (s in (1:length(unique(SubtypeData$mRNAseq_cNMF)))){
    Subtype <- SubtypeData$mRNAseq_cNMF
    names(Subtype) <- rownames(SubtypeData)
    Subtype[-which(Subtype==s)] <- rep(0,length(Subtype[-which(Subtype==s)]))
    Subtype[which(Subtype>0)] <- rep(1,length(which(Subtype>0)))
    
    ####### cross-validation loop #######
    AssignGroups <- cvFolds(n=nrow(SubtypeData),K=nFoldCV)
    Models <- list()
    #PerfROC <- list()
    for (j in (1:nFoldCV)){
      TestSet <- AssignGroups$subsets[which(AssignGroups$which==j)]
      TrainSet <- c(1:nrow(SubtypeData))
      TrainSet <- TrainSet[-TestSet]
      TrainMatrix <- MAmoduleMatrix[,TrainSet]
      TestMatrix <- MAmoduleMatrix[,TestSet]
      TrainSubtype <- Subtype[colnames(TrainMatrix)]
      TestSubtype <- Subtype[colnames(TestMatrix)]
      
      ####### Adjust a linear model #######
      fit <- cv.glmnet(t(TrainMatrix),TrainSubtype, alpha = 1.0,family="binomial")
      bestNonZeroLambda <- fit$lambda.1se
      fitOpt <- glmnet(t(TrainMatrix),TrainSubtype,alpha = 1.0,lambda=bestNonZeroLambda,family="binomial")
      b_o<-coef(fitOpt,s = bestNonZeroLambda)
      Models <- c(Models,list(b_o))
      
      ######## Predict using this model ########
      Predictions = predict(fitOpt,type="class",t(TestMatrix)) 
      pred <- prediction(as.numeric(Predictions),TestSubtype)
      perf <- performance(pred,"tpr","fpr")
      auc <- performance(pred,"auc")
      AUC[j,s] <- unlist(auc@y.values)
      mse <- performance(pred,"rmse")
      MSE[j,s] <- unlist(mse@y.values)
#      PerfROC <- c(PerfROC,list(perf))      
    }
    AllModels <- c(AllModels,list(Models))
#    ROC <- c(ROC,list(PerfROC))
  }
  names(AllModels) = colnames(MSE)=colnames(AUC) <- paste0("Subtype",c(1:length(unique(SubtypeData$mRNAseq_cNMF))))
  #color = rainbow(nFoldCV)
  #for (s in (1:length(unique(SubtypeData$mRNAseq_cNMF)))){
    #for (j in (1:nFoldCV)){
      #perf <- ROC[[s]]
      #perf <- perf[[j]]
      #if (j==1){
        #plot(perf,col=color[j],main=paste0("ROC curve for subtype ",s," and cancer ",CancerSite))
      #} else {
        #plot(perf,col=color[j],add=TRUE)
      #}   
    #}
  #}
  list(Models=AllModels,MSE=MSE,AUC=AUC)
}

AMARETTO_LoadSubtypes <- function(CancerSite,TargetDirectory,downloadData){
  ######## download the subtype data #########
  command <- paste("mkdir ",TargetDirectory,sep="")
  if (!file.exists(TargetDirectory)) {
    #system(command)
    dir.create(file.path(TargetDirectory), showWarnings = FALSE)
  }
  gdacURL= "http://gdac.broadinstitute.org/runs/"
  TCGA_acronym_uppercase=toupper(CancerSite)
  fileType= "tar.gz"
  saveDir=TargetDirectory  
  
  dataType='analyses'
  if (CancerSite=="LAML"){
    dataFileTag='TB.Aggregate_Molecular_Subtype_Clusters.Level_4'
  } else {
    dataFileTag='TP.Aggregate_Molecular_Subtype_Clusters.Level_4'  
  }
  if (length(dataFileTag)==1) {
    Subtypedirectory=get_firehoseData(downloadData,saveDir,TCGA_acronym_uppercase,dataType,dataFileTag)    
  } else {
    Subtypedirectory=c()
    for (i in 1:length(dataFileTag)) {
      Subtypedirectory=c(MAdirectory,get_firehoseData(downloadData,saveDir,TCGA_acronym_uppercase,dataType,dataFileTag[i]))
    }        
  }
  if (CancerSite=="LAML") { 
    subtypestring='TB.mergedcluster.txt'
  } else {
    subtypestring='TP.mergedcluster.txt'
  }        
  Subtypefiles=dir(Subtypedirectory)
  MatchedFile=grep(subtypestring,Subtypefiles)        
  SubtypeData=read.csv(paste(Subtypedirectory,Subtypefiles[MatchedFile],sep=''),sep="\t",row.names=1,head=TRUE,na.strings=c("NA","null")) 
}
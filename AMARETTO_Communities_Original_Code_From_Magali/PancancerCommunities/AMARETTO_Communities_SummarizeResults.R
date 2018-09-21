### summarize the results obtained after running AMARETTO_Communities_Example

load("/Users/mchampion/Desktop/PancancerAnalysis/Communities/Communities.Rdata")
load("/Users/mchampion/Desktop/PancancerAnalysis/Communities/AnalysedResults.Rdata")
load("/Users/mchampion/Desktop/PancancerAnalysis/SurvivalAnalysis/ResultsClinical.Rdata")
load("/Users/mchampion/Desktop/PancancerAnalysis/SurvivalAnalysis/ResultsSurvival.Rdata")
Community <- substr(names(AnalysedResults),10,length(AnalysedResults))
CommSize <- Communities$CommInfo$CommSize[which(Communities$CommInfo$KeepCom==TRUE)]
CommInfoResults <- matrix(0,ncol=9,nrow=length(Community))
CommInfoResults[,1] <- CommSize
colnames(CommInfoResults) <- c("CommSize","BestEA","BestReg","BestClin","BestSurv","PercentageEA","PercentageReg","PercentageClin","PercentageSurv")
rownames(CommInfoResults) <- Community
for (i in (1:length(Community))){
  commnbr <- as.numeric(Community[i])
  EA <- AnalysedResults[[i]]$EA
  Reg <- AnalysedResults[[i]]$Regulators
  Clin <- ResultsClinical[[i]]$Info
  Surv <- ResultsSurvival[[i]]$Info
  if (nrow(EA)>0){
    CommInfoResults[i,2] <- EA[1,ncol(EA)]  
  }
  if (nrow(Reg)>0){
    CommInfoResults[i,3] <- Reg[1,ncol(Reg)]
  }  
  CommInfoResults[i,4] <- Clin[1,"OverlapMax"]
  CommInfoResults[i,8] <- Clin[1,"PercentageMax"]  
  CommInfoResults[i,5] <- Surv["OverlapMax"]
  CommInfoResults[i,9] <- Surv["PercentageMax"] 
}
CommInfoResults[,6] <- CommInfoResults[,"BestEA"]/CommInfoResults[,"CommSize"]
CommInfoResults[,7] <- CommInfoResults[,"BestReg"]/CommInfoResults[,"CommSize"]
Mean <- rowMeans(CommInfoResults[,c(6,7,8,9)])
CommInfoResults <- cbind(CommInfoResults,MeanResults=Mean)
CommInfoResults <- CommInfoResults[order(CommInfoResults[,"MeanResults"],decreasing=TRUE),]

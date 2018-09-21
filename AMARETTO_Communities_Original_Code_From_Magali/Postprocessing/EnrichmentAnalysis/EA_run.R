#### Enrichment analysis on one cluster ####
EnrichmentAnalysis <- function(Data,DataBase,Multi=FALSE,Universe=c(),Threshold=list(maxpvalue=0.05,mininter=1,maxgenesetsize=500)){
  # Data: our set of genes (on which we want to do EA)
  # DataBase: obtained by running function EA_createData (e.g. MSigDB or GeneSig)
  # Multi: if TRUE, multidim case, compute the universe before running the algorithm  
  #        by default, Multi is FALSE
  # Threshold: a list of three different thresholds - max p.value
  #                                                 - min intersection beween the DataBase and the cluster of genes 
  #                                                 - max number of genes in the DataBase
  
  library(doParallel)
  
  ##### 1st step: define the universe (If Multi=FALSE) #####
  # Note: not too large (to avoid small p.values) but must contain all genes of the DataBase and all interesting genes 
  if (Multi==FALSE){
      UniverseGenes <- union(unlist(DataBase),Data)
      UniverseSize <- length(UniverseGenes)
    } else{
      UniverseGenes = Universe
      UniverseSize = length(Universe)
    }
       
  ##### 2nd step: main loop #####
  n.cluster <- 4
  cl2 <- makeCluster(c(rep("localhost", n.cluster)), type = "SOCK")
  cluster = cl2
  
  registerDoParallel(n.cluster)
  getDoParWorkers()
  
  ResultsLoop = foreach(i=1:length(DataBase),.combine='rbind') %dopar% {
    Geneset <- DataBase[[i]]
    GenesetSize <- length(Geneset)
    
    Inter <- length(intersect(Geneset,Data)) # genes in common
    
    pvalue <- phyper(Inter-1,length(Data),UniverseSize-length(Data),GenesetSize,lower.tail=FALSE,log.p=FALSE)
    # when we compute only one p.value, it's unnecessary to compute the adjusted p.value 
    c(p.value=pvalue,Intersection=Inter,GeneSetSize=GenesetSize)
  }
  rownames(ResultsLoop) <- names(DataBase)
  stopCluster(cl2)
  
  ##### 3rd step: threshold #####
  pvalueAdjust <- p.adjust(ResultsLoop[,1],method="hochberg") # Adjusted p.value
  ResultsLoop = cbind(pvalueAdjust,ResultsLoop) 
  I1.threshold = which(ResultsLoop[,1]<Threshold[[1]])
  I2.threshold = which(ResultsLoop[,3]>Threshold[[2]])
  I3.threshold = which(ResultsLoop[,4]<Threshold[[3]])
  I.threshold = intersect(I1.threshold,I2.threshold)
  I.threshold = intersect(I.threshold,I3.threshold)
  ResultsThreshold <- ResultsLoop[I.threshold,]
  if (length(I.threshold)>1){
    ResultsOrd <- ResultsThreshold[order(ResultsThreshold[,1]),]   
  } else{
    ResultsOrd <- ResultsThreshold
  }
  list(ResultThreshold=ResultsOrd,All.results=ResultsLoop,DataSize=length(Data)) 
}

#### Enrichment analysis on several gene sets ####
#### this code is usued for Pancancer AMARETTO
EnrichmentAnalysis_multi <- function(Data,DataBase,Threshold){
    # Data: a list containing lists of genes
    # DataBase: obtained by running function EA_createData
    # Threshold: a list of three different thresholds - max p.value
    #                                                 - min intersection beween the DataBase and the cluster of genes 
    #                                                 - max number of genes in the DataBase
    
    ##### 1st step: define the universe #####
    # Note: not too large (to avoid small p.values) but must contain all genes of the DataBase and all interesting genes 
    UniverseGenes <- union(unlist(DataBase),unlist(Data))
    UniverseSize <- length(UniverseGenes)
    
    ##### 2nd step: sample the Data #####
    EnrichmentResults <- matrix(1,nrow=length(Data),ncol=length(DataBase))
    colnames(EnrichmentResults) <- names(DataBase)
    rownames(EnrichmentResults) <- names(Data)
    
    for (j in 1:length(Data)){
        Data_uni <- Data[[j]]
        Results <- EnrichmentAnalysis(Data=Data_uni,DataBase=DataBase,Multi=TRUE,Threshold=Threshold,Universe=UniverseGenes)
        I <- rownames(Results$ResultThreshold)
        if (length(I)>1){
            EnrichmentResults[j,I] <- Results$ResultThreshold[,1]
        } else{
            if (length(I)>0){
                EnrichmentResults[j,I] <- Results$ResultThreshold[1]
            } else {
            }
        }
    }
    ## Clean the data ##
    I <- c()
    for (i in (1:ncol(EnrichmentResults))){
        Test <- EnrichmentResults[,i]==rep(1,length(EnrichmentResults[,i]))
        if (length(which(Test==TRUE))==length(EnrichmentResults[,i])){
            I <- c(I,i)
        }
        else {
        }    
    }
    EnrichmentResults <- EnrichmentResults[,-I]
    EnrichmentResults
}
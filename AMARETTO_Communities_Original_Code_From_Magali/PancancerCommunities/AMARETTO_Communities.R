#### Run AMARETTO, compare the modules, create and analyse the communities #### 

############################################
##### Run AMARETTO Pancancer algorithm #####
############################################
AMARETTO_Pancancer <- function(AMARETTOPancancerData,InnerThreshold=c(pvalue=0.05,inter=5),Threshold=0.05){
    # AMARETTOPancancerData: list of AMARETTO results
    # Threshold: by default threshold=c(pvalue_min=0.05,inter_min=5)
    # Threshold.plot: threshold for plotting

    # 1st step: comparison between the modules #
    NumberCancer <- length(AMARETTOPancancerData) 
    
    GeneNames <- c()
    for (i in (1:length(AMARETTOPancancerData))){
      GeneNames <- c(GeneNames,AMARETTOPancancerData[[i]]["GeneNames"])
    }  
    names(GeneNames) <- names(AMARETTOPancancerData)
    AllGenes <- unique(unlist(GeneNames))  
    SizeAllGenes <- length(AllGenes)
    
    NrClusters <- c()
    for (i in (1:length(AMARETTOPancancerData))){
      NrClusters <- c(NrClusters,AMARETTOPancancerData[[i]]["N"])
    }
    names(NrClusters) <- names(AMARETTOPancancerData)
    N <- sum(unlist(N))
   
    Assign <- c()
    for (i in (1:length(AMARETTOPancancerData))){
      Assign <- c(Assign,AMARETTOPancancerData[[i]]["Assign"])
    }
    names(Assign) <- names(AMARETTOPancancerData)
    
    ResultsCompCancersp.value <- matrix(1,N,N)
    ResultsCompCancersinter <- matrix(0,N,N)
    Names <- c()
    for (i in (1:NumberCancer)){
        Cancer1 <- names(AMARETTOPancancerData)[i] 
        for (j in (1:NumberCancer)){
            if (j>i){
                Cancer2 <- names(AMARETTOPancancerData)[j]        
                CompCancer <- ModuleCompCancer(GeneNames,Assign,Cancers=names(AMARETTOPancancerData),Cancer1,Cancer2,NrClusters,SizeAllGenes)
                cat("\n","Module comparison between Cancer ", Cancer1," and Cancer ", Cancer2)            
                I1 = NrClusters[[i]]*(i-1)+1
                I2 = NrClusters[[i]]*i
                J1 = NrClusters[[j]]*(j-1)+1 
                J2 = NrClusters[[j]]*j
                ResultsCompCancersp.value[I1:I2,J1:J2] <- CompCancer$p.value
                ResultsCompCancersinter[I1:I2,J1:J2] <- CompCancer$inter
            }
            else {}
        }
        if (i<NumberCancer){
            Names <- c(Names,rownames(CompCancer$inter))
        }
        else {}
    }
    Names <- c(Names,colnames(CompCancer$inter))
    colnames(ResultsCompCancersinter) = colnames(ResultsCompCancersp.value) = rownames(ResultsCompCancersinter) = rownames(ResultsCompCancersp.value) = Names
    
    I.inter <- which(ResultsCompCancersinter<InnerThreshold['inter'])
    I.pvalue <- which(ResultsCompCancersp.value>InnerThreshold['pvalue'])
    I <- union(I.inter,I.pvalue)
    ResultsCompCancersp.value[I] <- 1
    ResultsCompCancersp.value[lower.tri(ResultsCompCancersp.value,diag=TRUE)] <- 0
    ResultsCompCancersp.value <- t(ResultsCompCancersp.value) + ResultsCompCancersp.value
    diag(ResultsCompCancersp.value) <- rep(1,nrow(ResultsCompCancersp.value))
    
    # 2nd step: plot the module network #
    ResultsPlot <- ResultsCompCancersp.value
    I <- which(ResultsCompCancersp.value>Threshold)
    ResultsPlot[I] <- 1
    ResultsPlot[lower.tri(ResultsPlot,diag=TRUE)] <- 1
    InitNetwork <- which(ResultsPlot < 1,arr.ind = TRUE)
    Network <- matrix(0,nrow=dim(InitNetwork)[1],ncol=4)
    for (i in 1:dim(InitNetwork)[1]){
      Network[i,1] <- rownames(ResultsPlot)[InitNetwork[i,1]]
      Network[i,2] <- colnames(ResultsPlot)[InitNetwork[i,2]]
      Network[i,3] <- ResultsPlot[InitNetwork[i,1],InitNetwork[i,2]]
    }
    Network[,4] <- -log(as.numeric(as.character(Network[,3])))
    
    colnames(Network) <- c('Node1','Node2',"p.value","score.pvalue")
    Network <- data.frame(Network)
    
    nodeNames <- union(Network[,1],Network[,2])
    nodeIDs <- c(1:length(nodeNames))  
    nodeAttributes <- cbind(nodeIDs,nodeNames)
    tmp <- data.frame(nodeIDs,nodeNames)
    colnames(tmp)[1] <- "nodeID1"
    colnames(tmp)[2] <- "Node1"
    EdgeMatrix <- join(tmp,Network,match="all",type="right",by=c("Node1"))
    tmp <- data.frame(nodeIDs,nodeNames)
    colnames(tmp)[1] <- "nodeID2"
    colnames(tmp)[2] <- "Node2"
    EdgeMatrix <- join(tmp,EdgeMatrix,match="all",type="right",by=c("Node2"))
    EdgeMatrix <- data.frame(EdgeMatrix$nodeID1,EdgeMatrix$nodeID2,
                                    EdgeMatrix$Node1,EdgeMatrix$Node2,EdgeMatrix$p.value,EdgeMatrix$score.pvalue)
    ModuleGraph <-  graph.data.frame(data.frame(EdgeMatrix),directed=FALSE,vertices=data.frame(nodeAttributes))
    
    graph.degrees <- degree(ModuleGraph)
    V(ModuleGraph)$size <- 2*sqrt(graph.degrees)
    Disease <- strsplit2(nodeAttributes[,2],"_")[,2]
    
    laytest <- layout.fruchterman.reingold(ModuleGraph)
    plot(ModuleGraph,layout=laytest,
         vertex.color=rainbow(length(AMARETTOPancancerData))[as.factor(Disease)],
         vertex.label=NA, 
         edge.width=as.numeric(as.character(EdgeMatrix$EdgeMatrix.score.pvalue))/10,
         main="Module network")
# to make the nodes transparent, vertex.color=adjustcolor(rainbow(length(Cancers))[as.factor(Disease)],alpha.f=.5)
    colours = rainbow(length(AMARETTOPancancerData))
    labels = names(AMARETTOPancancerData)
    legend("right",legend=names(AMARETTOPancancerData), col=colours, pch=19)    

ModuleMap <- list(p.value=ResultsCompCancersp.value,EdgeMatrix=EdgeMatrix,laytest=laytest)
return(ModuleMap)
}

####################################################
##### 2nd Step: create the communities Network #####
####################################################
AMARETTO_CreateCommunities <- function(ModuleMap,Threshold=FALSE,ratioCommSize=0.01, # ratio comm size, community network size
                                       MinCancer=2, # at least 2 cancers represented in each community
                                       ratioCancerSize=0.1, # the number of cancers min depends on the comm size
                                       ratioEdgesInOut=0.5){
    library("stringr")
    library("RColorBrewer")  
    library("plyr")  
    library("igraph")
    library("limma") 
    brewPal = c("Set3")
    
    ### 1st step: cluster the modules in community ###
    EdgeMatrix <- ModuleMap$EdgeMatrix
    if (Threshold>0){
      I <- which(as.numeric(as.character(EdgeMatrix$EdgeMatrix.p.value))<Threshold)
      EdgeMatrix <- EdgeMatrix[I,]
    } 
    nodeNames <- union(EdgeMatrix$EdgeMatrix.Node1,EdgeMatrix$EdgeMatrix.Node2)
    nodeIDs <- c(1:length(nodeNames))  
    nodeAttributes <- cbind(nodeIDs,nodeNames)
    tmp <- data.frame(nodeIDs,nodeNames)
    colnames(tmp)[1] <- "nodeID1"
    colnames(tmp)[2] <- "EdgeMatrix.Node1"
    EdgeMatrix <- join(tmp,EdgeMatrix[,-1],match="all",type="right",by=c("EdgeMatrix.Node1"))
    tmp <- data.frame(nodeIDs,nodeNames)
    colnames(tmp)[1] <- "nodeID2"
    colnames(tmp)[2] <- "EdgeMatrix.Node2"
    EdgeMatrix <- join(tmp,EdgeMatrix[,-3],match="all",type="right",by=c("EdgeMatrix.Node2"))
    EdgeMatrix <- data.frame(EdgeMatrix$nodeID1,EdgeMatrix$nodeID2,
                             EdgeMatrix$EdgeMatrix.Node1,EdgeMatrix$EdgeMatrix.Node2,EdgeMatrix$EdgeMatrix.p.value,EdgeMatrix$EdgeMatrix.score.pvalue)
    colnames(EdgeMatrix) <- c("nodeID1","nodeID2","Node1","Node2","p.value","score.pvalue") 
    ModuleGraph <-  graph.data.frame(data.frame(EdgeMatrix),directed=FALSE,vertices=data.frame(nodeAttributes))    
    graph.degrees <- degree(ModuleGraph)
    V(ModuleGraph)$size <- 2*sqrt(graph.degrees)
    
    comm <- edge.betweenness.community(ModuleGraph,directed=FALSE,merges=TRUE,modularity=TRUE,
                                       membership=TRUE,weights=as.numeric(as.character(EdgeMatrix$score.pvalue))) 
    message("There are ", length(unique(comm$membership))," different communities detected using weighted edges.")

    CommGraph <-  graph.data.frame(data.frame(EdgeMatrix),directed=FALSE,vertices=data.frame(nodeAttributes))
    names(comm$membership) <- V(CommGraph)$nodeNames
    
    ### 2nd step: some info about communities ###
    membership <- cbind(c(1:length(comm$membership)),comm$membership);
    colnames(membership) <- c("nodeID","community");
    numCommunitiesOrig <- length(unique(membership[,"community"]));
    membership[,"community"] <- as.factor(membership[,"community"]);   
    cancer <- strsplit2(rownames(membership),"_")[,2]
    membership <- cbind(membership,cancer)
    colnames(membership)[ncol(membership)] <- "studyNum";
    
    numEdgesInComm <- c()
    totalNumEdges <- c()    
    for(m in 1:nrow(membership)){ 
      commNum <- membership[m,"community"]
      Id <-  membership[m,"nodeID"]
      edgeMatrixIndices <- union(which(EdgeMatrix[,1]==Id), which(EdgeMatrix[,2]==Id))
      edgeMatrixShort <- EdgeMatrix[edgeMatrixIndices, ]
      edgeMatrixVector <- c(edgeMatrixShort[,1],edgeMatrixShort[,2])
      edgeMatrixVector <- edgeMatrixVector[-which(edgeMatrixVector ==Id)]
      
      totalNumEdges[m] <- nrow(edgeMatrixShort)
      membershipShort <- membership[match(edgeMatrixVector,membership[,"nodeID"]),]
      if (nrow(edgeMatrixShort)==1){
        numEdgesInComm[m] <- length(which(membershipShort["community"]==commNum)) 
      } else {
        numEdgesInComm[m] <- length(which(membershipShort[,"community"]==commNum))     
      }
    }
    
    fractEdgesInOut <- numEdgesInComm/totalNumEdges
    numEdgesNotInComm <- totalNumEdges - numEdgesInComm
    Perm <- sample(1:numCommunitiesOrig)
    Color <- rainbow(length(unique(comm$membership)))
    Color<-Color[Perm]
    Color<- Color[comm$membership]
    
    membership <- data.frame(membership,Color, cbind(totalNumEdges,numEdgesInComm,fractEdgesInOut,numEdgesNotInComm))
    
    commEdgeInfo <- data.frame()
    for(c in 1:numCommunitiesOrig){  
      #these will be double-counted as both nodes for each edge will be in the community
      numTotalEdgesInCommunity <- sum(membership[which(membership[,"community"]==c),"numEdgesInComm"])/2
      #these will not be double-counted.
      numTotalEdgesNotInCommunity <- sum(membership[which(membership[,"community"]==c),"numEdgesNotInComm"])
      
      fractEdgesInVsOut <- numTotalEdgesInCommunity/(numTotalEdgesNotInCommunity+numTotalEdgesInCommunity)
      
      numDatasetsPerCommunity <- length(unique(membership[which(membership[,"community"]==c), "studyNum"]))
      CommSize <- nrow(membership[which(membership[,"community"]==c),])
      fractDatasetsSize <- numDatasetsPerCommunity/CommSize
      tmp <- data.frame(numTotalEdgesInCommunity,numTotalEdgesNotInCommunity, fractEdgesInVsOut, numDatasetsPerCommunity,CommSize,fractDatasetsSize,stringsAsFactors=FALSE)
      commEdgeInfo <- rbind(commEdgeInfo,tmp)  
    }
    Color2 <- rainbow(length(unique(comm$membership)),alpha=0.3)
    Color2<-Color2[Perm]
    commEdgeInfo <- data.frame(commEdgeInfo,Color2)
    
    # rename the communities from the largest to the smallest #
    commEdgeInfo <- commEdgeInfo[order(commEdgeInfo$CommSize,decreasing=TRUE),] 
    NewComNumber <- c(1:nrow(commEdgeInfo))
    names(NewComNumber) <- rownames(commEdgeInfo)
    rownames(commEdgeInfo) <- NewComNumber
    membership$community <- match(membership$community,names(NewComNumber))
    
    ### Step 3: remove communities according to the thresholds ###    
    SmallComm<-which(commEdgeInfo[,"CommSize"]/nrow(nodeAttributes)<ratioCommSize)
    HomogenousComm <- which(commEdgeInfo[,"numDatasetsPerCommunity"]<MinCancer)
    HomogenousComm2 <- which(commEdgeInfo[,"fractDatasetsSize"]<ratioCancerSize)
    CommEx <- which(commEdgeInfo[,"fractEdgesInVsOut"]<ratioEdgesInOut)
    RemoveCommunities <- unique(c(SmallComm,HomogenousComm,HomogenousComm2,CommEx))
    message("There are ", length(RemoveCommunities)," communities to remove.")
    KeepComm <- c(1:numCommunitiesOrig)
    if (length(RemoveCommunities)>0){
      KeepComm <- KeepComm[-RemoveCommunities]  
    } 
    names(KeepComm) <- c(1:length(KeepComm))
    names(RemoveCommunities) <- c(length(KeepComm):numCommunitiesOrig)[-1]
    test <- c(1:nrow(commEdgeInfo))
    InfoRemove <- (test%in%KeepComm)
    commEdgeInfo <-cbind(commEdgeInfo,KeepCom=InfoRemove)      
    
    List <- list()
    for (i in (1:length(KeepComm))){
      Commun <- KeepComm[[i]]
      CommunMember <- as.numeric(as.vector(membership[which(membership[,"community"]==Commun),"nodeID"]))
      List <- c(List,list(CommunMember))
    }
    CommGraphSmall <-  graph.data.frame(data.frame(EdgeMatrix),directed=FALSE,vertices=data.frame(membership))
    graph.degrees <- degree(CommGraphSmall)
    V(CommGraphSmall)$size <- 2*sqrt(graph.degrees)
    
    plot(CommGraphSmall,layout=ModuleMap$laytest,
         # vertex.frame.color=rainbow(length(unique(comm$membership)))[comm$membership],
         vertex.color=as.character(membership$Color),
         vertex.label=NA, 
         #edge.color="grey",
         mark.groups = List,
         mark.col=as.character(commEdgeInfo$Color2[KeepComm]),
         mark.border=NA,
         edge.width=as.numeric(as.character(EdgeMatrix$score.pvalue))/10,
         main="Community network")
    
    Communities <- list(CommInfo=commEdgeInfo,NodeInfo=membership,EdgeMatrix=EdgeMatrix)
    return(Communities)
}

#############################################
##### 3rd Step: Analyse the communities #####
#############################################
AMARETTO_AnalyseCommunities <- function(AMARETTOPancancer_results,Communities,All.Communities,CommNumber,AMARETTOPancancerData){
    
    if (All.Communities==TRUE){
      CommNumber = which(Communities$CommInfo$KeepCom==TRUE)
    }
    
    ### 1st step: highlight the communities ###
    CommGraph <- graph.data.frame(data.frame(Communities$EdgeMatrix),directed=FALSE,vertices=data.frame(Communities$NodeInfo))
    graph.degrees <- degree(CommGraph)
    V(CommGraph)$size <- 2*sqrt(graph.degrees)
    
    List <- list()
    for (i in (1:length(CommNumber))){
      community=CommNumber[i]
      CommunMember <- as.numeric(as.vector(Communities$NodeInfo[which(Communities$NodeInfo[,"community"]==community),"nodeID"]))
      List <- c(List,list(CommunMember))
    }
    if (length(CommNumber)>1){
      title <- c(paste0("Zoom on Community ",CommNumber[1]),paste0(CommNumber[2:length(CommNumber)]))
      title <- paste0(title,collapse=",")
    } else {
      title <- paste0("Zoom on Community ",CommNumber)
    }
    I <- which(Communities$NodeInfo$community %in% CommNumber)
    ColorComm <- as.character(Communities$NodeInfo$Color)
    ColorComm[-I] <- rep("gray",length(ColorComm)-length(I))
    
    plot(CommGraph,layout=AMARETTOPancancer_results$laytest,
         # vertex.frame.color=rainbow(length(unique(comm$membership)))[comm$membership],
         vertex.color=ColorComm,
         vertex.label=NA, 
         #edge.color="grey",
         mark.groups = List,
         mark.col=as.character(Communities$CommInfo$Color2[CommNumber]),
         mark.border=NA,
         edge.width=as.numeric(as.character(Communities$EdgeMatrix$score.pvalue))/10,
         main=title)    
    
    ### 2nd step: Analyse the Communities ###
    DataBase <- CreateDataBase() # DataBase from GeneSetDB and MSigDB
    
    AnalysedResults <- list()
    Names <- c()
    for (i in (1:length(CommNumber))){
        Community <- CommNumber[i]
        CommData <- Communities$NodeInfo$community
        names(CommData) <- rownames(Communities$NodeInfo)
        cat("\n","Analysing Community", CommNumber[i])
        
        ## Find the regulators ##
        regulators <- FindRegulators(CommData=CommData,Community,AMARETTOPancancerData)
        
        ## Enrichment analysis ##
        resultEAComm <- PerformEA(CommData=CommData,Community,DataBase,AMARETTOPancancerData)
        
        results <- list(Regulators=regulators,EA=resultEAComm)
        AnalysedResults <- c(AnalysedResults,list(results))        
        Name <- paste0("Community",Community)
        Names <- c(Names,Name)
    }
    names(AnalysedResults) <- Names
    AnalysedResults
}

##################################
##### Other useful functions #####
##################################

#### Comparison between two cancers #####
ModuleCompCancer <- function(GeneNames,Assign,Cancers,Cancer1,Cancer2,NrClusters,SizeAllGenes){
    # Cancer 1 the first cancer
    # Cancer 2 the second cancer
    # AMARETTOPancancerData all the results from AMARETTO 
    library(doParallel)
    loc1 <- which(Cancers==Cancer1)
    loc2 <- which(Cancers==Cancer2)
    
    Resultsp.value <- matrix(1,nrow=NrClusters[[loc1]],ncol=NrClusters[[loc2]])
    Resultsinter <- matrix(0,nrow=NrClusters[[loc1]],ncol=NrClusters[[loc2]])
    for (i in (1:NrClusters[[loc1]])){
        NrModule1 <- i
        
        n.cluster <- 4
        cl2 <- makeCluster(c(rep("localhost", n.cluster)), type = "SOCK")
        cluster = cl2
        registerDoParallel(n.cluster)
        getDoParWorkers()
        
        ResultsLoop = foreach(j=1:NrClusters[[loc2]],.combine='rbind') %dopar% {
            NrModule2 <- j
            Results <- ModuleComp(GeneNames,Assign,Cancers,Cancer1,Cancer2,NrModule1,NrModule2,SizeAllGenes)
            c(Results)
        }
        ResultsLoop[,1] <- p.adjust(ResultsLoop[,1],method="hochberg")
        stopCluster(cl2)
        
        Resultsp.value[i,] <- ResultsLoop[,1]
        Resultsinter[i,] <- ResultsLoop[,2]
    }
    Max <- max(NrClusters[[loc1]],NrClusters[[loc2]])
    Name1 <- c()
    Name2 <- c()
    for (i in 1:Max){
        name1 <- paste0("M",i,"_",Cancer1)
        Name1 <- c(Name1,name1)
        name2 <- paste0("M",i,"_",Cancer2)
        Name2 <- c(Name2,name2)
    }
    rownames(Resultsp.value)=rownames(Resultsinter) <- Name1[1:NrClusters[[loc1]]] 
    colnames(Resultsp.value)=colnames(Resultsinter) <- Name2[1:NrClusters[[loc2]]]
    ResultsModuleComp <- list(p.value=Resultsp.value,inter=Resultsinter)
}

##### Comparison between two modules ####
ModuleComp <- function(GeneNames,Assign,Cancers,Cancer1,Cancer2,NrModule1,NrModule2,SizeAllGenes){
    ## Define the gene sets ##
    loc1 <- which(Cancers==Cancer1)
    loc2 <- which(Cancers==Cancer2)
    ModuleData1 = GeneNames[[loc1]]
    ModuleData1 <- ModuleData1[Assign[[loc1]]==NrModule1]
    ModuleSize1 <- length(ModuleData1)
    ModuleData2 = GeneNames[[loc2]]
    ModuleData2 <- ModuleData2[Assign[[loc2]]==NrModule2]
    ModuleSize2 <- length(ModuleData2)
    
    ## HyperGeoTest ##
    Inter <- length(intersect(ModuleData1,ModuleData2))
    pvalue <- phyper(Inter-1,ModuleSize2,SizeAllGenes-ModuleSize2,ModuleSize1,lower.tail=FALSE,log.p=FALSE)
    Results = c(p.value=pvalue,Intersection=Inter,Module1=ModuleSize1,Module2=ModuleSize2)
}

#### Analyse EA results ####
PerformEA <- function(CommData,Community,DataBase,AMARETTOPancancerData){
    ModuleName = names(CommData[which(CommData==Community)])
    GeneNames <- c()
    for (i in (1:length(AMARETTOPancancerData))){
      GeneNames <- c(GeneNames,AMARETTOPancancerData[[i]]["GeneNames"])
    }  
    names(GeneNames) <- names(AMARETTOPancancerData)
   
    Assign <- c()
    for (i in (1:length(AMARETTOPancancerData))){
      Assign <- c(Assign,AMARETTOPancancerData[[i]]["Assign"])
    }
    names(Assign) <- names(AMARETTOPancancerData)
    
    Data <- list()
    Name <- c()
    for (i in (1:length(ModuleName))){
        myData <- ModuleName[i]
        Cancer <- strsplit2(myData, "_")[,2] 
        ModuleNr <- strsplit2(myData, "_")[,1]
        ModuleNr <- as.numeric(str_replace_all(ModuleNr,"[[A-Z]]","")) 
        ModuleData = GeneNames[[Cancer]]
        ModuleData <- ModuleData[Assign[[Cancer]]==ModuleNr]
        Data <- c(Data,list(ModuleData))
        Name <- c(Name,myData)
    }
    names(Data)= Name 
    SizeComm <- length(ModuleName) 
    
    Threshold = list(maxpvalue=0.05,mininter=5,maxgenesetsize=500)
    resultsEAComm = EnrichmentAnalysis_multi(Data=Data,DataBase=DataBase,Threshold=Threshold)
    resultsEAComm <- t(resultsEAComm)
    
    ## Clean the results (only if there is something to clean!!) ##     
    if (nrow(resultsEAComm)>0){
      Mean <- c()
      Nr <- c()
      for (i in (1:nrow(resultsEAComm))){
          EAshort <- resultsEAComm[i,which(resultsEAComm[i,]<1)]
          mean <- mean(EAshort)
          Mean <- c(Mean,mean)
          nr <- length(EAshort)
          Nr <- c(Nr,nr)
      }
      resultsEAComm <- cbind(resultsEAComm,Mean,Nr)    
      resultsEAComm <- resultsEAComm[order(resultsEAComm[,ncol(resultsEAComm)],decreasing=TRUE),]
    }
    resultsEAComm
}

#### Find the regulators ####
FindRegulators <- function(CommData,Community,AMARETTOPancancerData){
    ModuleName = names(CommData[which(CommData==Community)])
    Reg <- list()
    Name <- c()
    for (i in (1:length(ModuleName))){
        myData <- ModuleName[i]
        Cancer <- strsplit2(myData, "_")[,2] 
        ModuleNr <- strsplit2(myData, "_")[,1]
        ModuleNr <- as.numeric(str_replace_all(ModuleNr,"[[A-Z]]","")) 
        res <- AMARETTOPancancerData[[Cancer]]["v"]
        I = which(abs(res$v[ModuleNr,])>0)
        Reg <- c(Reg,list(Regulator=names(I)))      
        Name <- c(Name,myData)       
    }
    names(Reg) = Name 
    
    ## Clean the results ##
    RegNames <- unique(unlist(Reg))
    Regulators <- matrix(0,ncol=length(ModuleName),nrow=length(RegNames))
    colnames(Regulators) <- Name
    rownames(Regulators) <- RegNames
    if (nrow(Regulators)>0){
      for (i in (1:length(Name))){
          Regulators[Reg[[Name[i]]],i] <- rep(1,length(Reg[[Name[i]]]))
       }
      Nr <- rowSums(apply(Regulators[,c(1:length(Name))], c(1, 2), as.numeric))
      Regulators <- cbind(Regulators,Nr)  
      Regulators <- Regulators[order(as.numeric(Regulators[,ncol(Regulators)]),decreasing=TRUE),]
    }
    Regulators <- data.frame(Regulators)
}

###### Create the AMARETTO results file #######
CreatePancancerData <- function(Cancers,ResultsDirectory){
  AMARETTOPancancerData <- list()
  for (i in (1:length(Cancers))){
    CancerSite <- Cancers[i]
    load(paste0(ResultsDirectory,"AMARETTOresults",CancerSite,".RData"))
    AMARETTOPancancerData <- c(AMARETTOPancancerData,list(AMARETTOresults))
  }
  names(AMARETTOPancancerData) <- Cancers
}
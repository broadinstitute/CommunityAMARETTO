#' cAMARETTO_ReadZipExport
#'
#' @param zipdirs A list of multiple AMARETTO results in a zipformat retrieved by the AMARETTO_ExportData function
#'
#' @return a list with AMARETTOinit and AMARETTOresults data objects from multiple runs
#'

cAMARETTO_ReadZipExport <- function(zipdirs){
  for (zipdir in zipdirs){
    unzip(zipdir)
  }
}

library(doParallel)
library(gtools)
library(tidyverse)

cAMARETTO_Results <- function(AMARETTOinit_all,AMARETTOresults_all,NrCores,output_dir){
  #test if names are matching
  if (all(names(AMARETTOinit_all) == names(AMARETTOresults_all))) {
    runnames<-names(AMARETTOresults_all)
    if (!length(unique(runnames)) == length(runnames)){
      stop("The run names are not unique. Give unique names.")
    }
  } else {
    stop("The names of the lists are not matching.")
  }
  
  dir.create(file.path(output_dir, "gmt_files"), recursive = FALSE, showWarnings = FALSE)
  
  # for each file a gmt for the modules
  for (run in runnames){
    gmt_file <- file.path(output_dir,"gmt_files", paste0(run, "_modules.gmt"))
    GmtFromModules(AMARETTOinit_all[[run]], AMARETTOresults_all[[run]], gmt_file, run)
  }
  # compare gmts pairwise between runs

  all_run_combinations <- as.data.frame(combinations(n=length(runnames), r=2, v=runnames, repeats.allowed=F))

  output_hgt_allcombinations <- apply(all_run_combinations, 1, function(x) {
    gmt_run1 <- file.path(output_dir, "gmt_files", paste0(x["V1"], "_modules.gmt"))
    gmt_run2 <- file.path(output_dir, "gmt_files", paste0(x["V2"], "_modules.gmt"))
    output_hgt_combination <- HyperGTestGeneEnrichment(gmt_run1, gmt_run2, NrCores)
    return(output_hgt_combination)
  })
  
  output_hgt_allcombinations <- do.call(rbind, output_hgt_allcombinations)
  output_hgt_allcombinations$padj <- p.adjust(output_hgt_allcombinations$p_value, method="BH")
  output_hgt_allcombinations <- output_hgt_allcombinations %>% 
                                  mutate(p_value=case_when(Geneset == Testset~NA_real_, TRUE~p_value))
  return(list(runnames=runnames, hgt_modules=output_hgt_allcombinations, NrCores=NrCores))
}

library(ComplexHeatmap)
library(tidyverse)
cAMARETTO_heatmap<-function(cAMARETTOresults,run1,run2){
  results_filtered<-cAMARETTOresults$hgt_modules %>% filter((grepl(run1,Testset)|grepl(run2,Testset))&(grepl(run1,Geneset)|grepl(run2,Geneset)))
  pvalue_matrix<-spread(results_filtered %>% select(Testset,Geneset,p_value),key=Geneset,value=p_value)
  pvalue_matrix<-column_to_rownames(pvalue_matrix,"Testset")
  pvalue_matrix<--log10(pvalue_matrix)
  Heatmap(pvalue_matrix, name = "Pvalues Compare Modules", column_title = "Regulator Genes\nExpression",show_column_names=TRUE,column_names_gp = gpar(fontsize = 8),row_names_gp = gpar(fontsize = 8),
                    column_title_gp = gpar(fontsize = 12, fontface = "bold"), col=colorRamp2(c(0, max(pvalue_matrix,na.rm = TRUE)), c("white", "darkred")))
}

library(tidyverse)
library(randomcoloR)
library(igraph)
cAMARETTO_ModuleNetwork<-function(cAMARETTOresults,pvalue,inter,color_list=NULL,edge_method=NULL){
  output_hgt_allcombinations_filtered<-cAMARETTOresults$hgt_modules %>% filter(padj<pvalue & n_Overlapping>=inter)
  node_information<-as.data.frame(unique(c(output_hgt_allcombinations_filtered$Geneset,output_hgt_allcombinations_filtered$Testset)))
  colnames(node_information)<-c("modulenames")
  node_information<-node_information %>% mutate(run=sub("_Module_.*$","",modulenames))
  module_network <- graph_from_data_frame(d=output_hgt_allcombinations_filtered,vertices=node_information,directed=FALSE)
  if (is.null(color_list)){
    color_list<-randomColor(length(cAMARETTOresults$runnames),luminosity="bright")
    names(color_list)<-cAMARETTOresults$runnames
  } else {
    cAMARETTOresults$runnames %in% names(color_list)
  }
  V(module_network)$color <- color_list[V(module_network)$run]
  V(module_network)$size<-2*sqrt(degree(module_network,mode="all"))
  if (!is.null(edge_method)){
    if (edge_method=="pvalue"){
      E(module_network)$width<- -(log10(E(module_network)$p_value))*0.2
    } else if (edge_method=="overlap"){
    E(module_network)$width<-E(module_network)$n_Overlapping/8
    }
  }
  layoutMN<-layout_with_fr(module_network)
  plot(module_network,vertex.frame.color=NA,layout=l,vertex.label=NA,main="Module network",edge.color="gray80")
  legendMN<-legend(x=-1.5,y=-1.1,legend=names(color_list), col=color_list, pch=19,bty="n")
  legendMN
  return(list(module_network=module_network,layoutMN=layoutMN,pvalue=pvalue,inter=inter,colMN=color_list))
}

library(tidyverse)
library(igraph)
cAMARETTO_IdentifyCom<-function(cAMARETTOnetworkM,color_list=NULL,ratioCommSize=0.01,MinCancer=2,ratioCancerSize=0.1,ratioEdgesInOut=0.5){
  
  comm <- edge.betweenness.community(cAMARETTOnetworkM$module_network,directed=FALSE, weights=NULL)
                                     #merges=TRUE,modularity=TRUE,
                                     #membership=TRUE,weights=output_hgt_allcombinations_filtered$padj)

  message("There are ", length(unique(comm$membership))," different communities detected using weighted edges.")

  names(comm$membership) <- V(cAMARETTOnetworkM$module_network)$name
  membership <- as.data.frame(cbind(c(1:length(comm$membership)),comm$membership))
  colnames(membership) <- c("nodeID","community")
  numCommunitiesOrig <- length(unique(membership[,"community"]))
  membership<-rownames_to_column(membership,"nodeName") %>% mutate(run=sub("_Module_.*$","",nodeName))
  
  Edges_Mnetwork<-as_data_frame(cAMARETTOnetworkM$module_network,what="edges")
  Nodes_Mnetwork<-as_data_frame(cAMARETTOnetworkM$module_network,what="vertices")
  
  for(m in 1:nrow(membership)){ 
    commNum <- membership[m,"community"]
    Id <-  membership[m,"nodeName"]
    edgeMatrixVector <- unlist(Edges_Mnetwork %>% filter(from==Id | to==Id) %>% select(from,to))
    edgeMatrixVector <- edgeMatrixVector[-which(edgeMatrixVector ==Id)]
    membership[m,"totalNumEdges"]<-length(edgeMatrixVector)
    membership[m,"numEdgesInComm"]<-nrow(membership[match(edgeMatrixVector,membership[,"nodeName"]),] %>% filter(community==commNum))
  }
  
  membership <- membership %>% mutate(fractEdgesInOut=numEdgesInComm/totalNumEdges,numEdgesNotInComm = totalNumEdges - numEdgesInComm)
  
  if (is.null(color_list)){
    color_list<-randomColor(length(numCommunitiesOrig),luminosity="bright")
    names(color_list)<-numCommunitiesOrig
    color_list<-color_list[comm$membership]
  } else {
    length(color_list)>=numCommunitiesOrig
  }

  commEdgeInfo <- membership %>% group_by(community) %>% summarise(numTotalEdgesInCommunity=sum(numEdgesInComm)/2,numTotalEdgesNotInCommunity=sum(numEdgesNotInComm),
                                                                   fractEdgesInVsOut=numTotalEdgesInCommunity/(numTotalEdgesNotInCommunity+numTotalEdgesInCommunity),
                                                                   numDatasetsPerCommunity=length(unique(run)),CommSize=n(),fractDatasetsSize=numDatasetsPerCommunity/CommSize,
                                                                   CommsizeFrac=CommSize/nrow(Nodes_Mnetwork))

  commEdgeInfo<-commEdgeInfo %>% arrange(-CommSize) %>% mutate(NewComNumber=row_number())
  membership<-left_join(membership,commEdgeInfo %>% select(community,NewComNumber))
  
  #Post Filter communities
   # ratio comm size, community network size
   # at least 2 cancers represented in each community
   # the number of cancers min depends on the comm size
  
  #filter
  KeepCommEdgeInfo <- commEdgeInfo %>% filter(CommsizeFrac>ratioCommSize & numDatasetsPerCommunity>MinCancer & fractDatasetsSize>ratioCancerSize & fractEdgesInVsOut>ratioEdgesInOut)
  message("There are ", nrow(commEdgeInfo) - nrow(KeepCommEdgeInfo)," communities to remove.")
  
  Nodes_Mnetwork<-left_join(Nodes_Mnetwork,membership%>% select(-run),by=c("name"="nodeName"))
  CommGraph <- graph.data.frame(Edges_Mnetwork,directed=FALSE,vertices=data.frame(Nodes_Mnetwork))
  graph.degrees <- igraph::degree(CommGraphSmall)
  V(CommGraphSmall)$size <- 2*sqrt(graph.degrees)
  
  community_list_df<-membership %>% filter(community %in% KeepCommEdgeInfo$community) %>%select(nodeName,community) %>% split(.$community)
  community_list<-lapply(community_list_df,function(x) unlist(x$nodeName))
  names(community_list)<-names(community_list_df)
  
  plot(CommGraph,layout=cAMARETTOnetworkM$layoutMN,
       vertex.color=as.character(Nodes_Mnetwork$color),
       vertex.label=NA,
       vertex.frame.color=NA,
       edge.color="gray80",
       mark.groups = community_list,
       #mark.col=as.character(commEdgeInfo$Color2[KeepComm]),
       mark.border=NA,
       main="Community network")

  legendMN<-legend(x=-1.5,y=-1.1,legend=names(cAMARETTOnetworkM$colMN), col=cAMARETTOnetworkM$colMN, pch=19,bty="n")
  legendMN
  return(list(CommGraph=CommGraph,community_list=community_list))
}

library(igraph)
cAMARETTO_GraphSelCom <- function(cAMARETTOnetworkC,com_selected){
  community_list=cAMARETTOnetworkC$community_list
  community_list_sel<-community_list[names(community_list) %in% com_selected]
  plot(cAMARETTOnetworkC$CommGraph,layout=cAMARETTOnetworkM$layoutMN,
       vertex.color=as.character(Nodes_Mnetwork$color),
       vertex.label=NA,
       vertex.frame.color=NA,
       edge.color="gray80",
       mark.groups = community_list_sel,
       #mark.col=as.character(commEdgeInfo$Color2[KeepComm]),
       mark.border=NA,
       main="Community network")
}

library(igraph)
cAMARETTO_HGTCom <- function(AMARETTOresults_all,cAMARETTOnetworkC,hyper_geo_reference="../../GetToKnowHowItWorks/H.C2CP.genesets_adapt.gmt",com_selected){
  
  for(community in cAMARETTOnetworkC$community_list){
    AMARETTOresults_all
    #retrieve target genes
    
    #retrieve regulator genes
  }
  
  #make community wise gmt
  regulators <- FindRegulators(CommData=CommData,Community,AMARETTOPancancerData)
    
  #enrichment analysis
  resultEAComm <- PerformEA(CommData=CommData,Community,DataBase,AMARETTOPancancerData)
    

}

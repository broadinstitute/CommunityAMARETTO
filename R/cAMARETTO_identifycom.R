#' @title cAMARETTO_IdentifyCom
#'
#' @param cAMARETTOnetworkM
#' @param color_list An optional list with colors
#' @param ratioCommSize
#' @param MinRuns
#' @param ratioRunSize
#' @param ratioEdgesInOut
#' @param edge_method Define edge weights based p-values or overlap between the gene sets.
#'
#' @return a list with the module network, layout for the network, used p-value, used overlap en colors
#' 
#' @import randomcoloR
#' @import tidyverse
#' @import igraph
#' @export
cAMARETTO_IdentifyCom <- function(cAMARETTOnetworkM, color_list=NULL, filterComm=TRUE, ratioCommSize=0.01, MinRuns=2, ratioRunSize=0.1, ratioEdgesInOut=0.5 ){
  
  comm <- edge.betweenness.community(cAMARETTOnetworkM$module_network, directed=FALSE, merges=TRUE, modularity=TRUE, membership=TRUE)
  
  message("There are ", length(unique(comm$membership)), " different communities detected using weighted edges.")
  
  names(comm$membership) <- V(cAMARETTOnetworkM$module_network)$name
  membership <- as.data.frame(cbind(c(1:length(comm$membership)), comm$membership))
  colnames(membership) <- c("nodeID", "community")
  numCommunitiesOrig <- length(unique(membership[, "community"]))
  membership<-rownames_to_column(membership, "nodeName") %>% mutate(run=sub("_Module_.*$", "", nodeName))
  
  Edges_Mnetwork <- igraph::as_data_frame(cAMARETTOnetworkM$module_network, what="edges")
  Nodes_Mnetwork <- igraph::as_data_frame(cAMARETTOnetworkM$module_network, what="vertices")
  
  for(m in 1:nrow(membership)){ 
    commNum <- membership[m, "community"]
    Id <-  membership[m, "nodeName"]
    edgeMatrixVector <- unlist(Edges_Mnetwork %>% filter(from==Id | to==Id) %>% select(from,to))
    edgeMatrixVector <- edgeMatrixVector[-which(edgeMatrixVector ==Id)]
    membership[m, "totalNumEdges"] <- length(edgeMatrixVector)
    membership[m, "numEdgesInComm"] <- nrow(membership[match(edgeMatrixVector, membership[,"nodeName"]), ] %>% filter(community==commNum))
  }
  
  membership <- membership %>% mutate(fractEdgesInOut=numEdgesInComm/totalNumEdges, numEdgesNotInComm = totalNumEdges - numEdgesInComm)

  commEdgeInfo <- membership %>% 
    group_by(community) %>% 
    summarise(numTotalEdgesInCommunity = sum(numEdgesInComm)/2, 
              numTotalEdgesNotInCommunity = sum(numEdgesNotInComm),
              fractEdgesInVsOut = numTotalEdgesInCommunity/(numTotalEdgesNotInCommunity+numTotalEdgesInCommunity),
              numDatasetsPerCommunity = length(unique(run)),
              CommSize = n(),
              fractDatasetsSize = numDatasetsPerCommunity/CommSize,
              CommsizeFrac = CommSize/nrow(Nodes_Mnetwork))
  
  commEdgeInfo <- commEdgeInfo %>% 
    arrange(-CommSize) %>% 
    mutate(NewComNumber=row_number())
  
  membership<-left_join(membership,commEdgeInfo %>% select(community, NewComNumber))
  
  #Post Filter communities
  # ratio comm size, community network size
  # at least 2 cancers represented in each community
  # the number of cancers min depends on the comm size
  
  if (filterComm ==TRUE) {
  KeepCommEdgeInfo <- commEdgeInfo %>% 
    filter(CommsizeFrac > ratioCommSize & numDatasetsPerCommunity > MinRuns & fractDatasetsSize > ratioRunSize & fractEdgesInVsOut > ratioEdgesInOut)
  
  message("There are ", nrow(commEdgeInfo) - nrow(KeepCommEdgeInfo)," communities to remove.")
  } else {
    KeepCommEdgeInfo <- commEdgeInfo
  }
  
  Nodes_Mnetwork <- left_join(Nodes_Mnetwork, membership %>% select(-run), by = c("name" = "nodeName"))
  CommGraph <- graph.data.frame(Edges_Mnetwork, directed=FALSE, vertices = data.frame(Nodes_Mnetwork))
  graph.degrees <- igraph::degree(CommGraph)
  V(CommGraph)$size <- 2*sqrt(graph.degrees)
  
  community_list_df <- membership %>% 
    filter(community %in% KeepCommEdgeInfo$community) %>%
    select(nodeName, community) %>% 
    split(.$community)
  
  community_list <- lapply(community_list_df, function(x) unlist(x$nodeName))
  names(community_list) <- names(community_list_df)
  
  
  if (is.null(color_list)){
    color_list <- randomColor(length(community_list), luminosity="light")
    names(color_list) <- names(community_list)
  } else {
    length(color_list) >= length(community_list)
  }
  
  plot(CommGraph,layout = cAMARETTOnetworkM$layoutMN,
       vertex.color = as.character(Nodes_Mnetwork$color),
       vertex.label = NA,
       vertex.frame.color = NA,
       edge.color = "gray80",
       mark.groups = community_list,
       mark.col = color_list,
       mark.border = NA,
       main = "Community network")
  
  legendMN <- legend(x = -1.5, y = -1.1+0.05*length(cAMARETTOnetworkM$colMN), legend = names(cAMARETTOnetworkM$colMN), col = cAMARETTOnetworkM$colMN, pch=19, bty="n",ncol=ceiling(length(cAMARETTOnetworkM$colMN)/5))
  legendMN
  
  legendCOM <- legend(x = 1.5, y = 1.5, legend=names(color_list), col=color_list, pch=19, bty="n",cex=max(0.9,1/(1+0.02*length(color_list))),ncol=ceiling(length(color_list)/15))
  legendCOM
  
  return(list(CommGraph=CommGraph, community_list=community_list, color_list=color_list))
}

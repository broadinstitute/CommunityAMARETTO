#' @title cAMARETTO_InformationTable
#'
#' @param cAMARETTOnetworkM The output of the Module Network function.
#' @param cAMARETTOnetworkC The output of the Identify Communities function.
#'
#' @return an information table about the communities
#'  containing the included nodes and overlapping genes per community.
#' 
#' @importFrom igraph as_data_frame degree E graph_from_data_frame
#'  layout_with_fr V graph.data.frame norm_coords edge.betweenness.community
#' @importFrom tibble add_row tibble column_to_rownames rownames_to_column
#' 
#' @examples 
#' try(
#' cAMARETTO_InformationTable(cAMARETTOnetworkM, cAMARETTOnetworkC)
#' )
#' @export
cAMARETTO_InformationTable <- function(cAMARETTOnetworkM, cAMARETTOnetworkC) {
  to<-from<-NULL
  Nodes_Cnetwork <- igraph::as_data_frame(cAMARETTOnetworkC$CommGraph,
                                          what="vertices")
  Edges_Cnetwork <- igraph::as_data_frame(cAMARETTOnetworkC$CommGraph,
                                          what="edges")
  Comm_Info <- tibble::tibble(community_numb = as.numeric(),
                              included_nodes = as.character(),
                              overlapping_genes = as.character())
  i=1
  for (comm in cAMARETTOnetworkC$community_list) {
    Edges_Cnetwork_comm <- Edges_Cnetwork %>% 
      filter(from %in% comm & to %in% comm)
    all_overlapping_genes <- unique(unlist(strsplit(
      Edges_Cnetwork_comm$Overlapping_genes,", ")))
    Comm_Info <- add_row(Comm_Info,
                         community_numb = as.numeric(
                           names(cAMARETTOnetworkC$community_list[i])),
                         included_nodes = paste(comm,collapse=", "),
                         overlapping_genes = paste(all_overlapping_genes,
                                                   collapse=", "))
    i = i + 1
  }
  return(Comm_Info)
}

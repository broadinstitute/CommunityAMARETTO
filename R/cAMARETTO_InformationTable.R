#' @title cAMARETTO_InformationTable
#'
#' @param cAMARETTOnetworkM
#' @param cAMARETTOnetworkC
#' @param ModuleNr
#'
#' @return a plot with one (or none) colored modules
#' 
#' @import igraph
#' @export
cAMARETTO_InformationTable <- function(cAMARETTOnetworkM, cAMARETTOnetworkC) {
  Nodes_Cnetwork <- igraph::as_data_frame(cAMARETTOnetworkC$CommGraph, what="vertices")
  Edges_Cnetwork <- igraph::as_data_frame(cAMARETTOnetworkC$CommGraph, what="edges")
  Comm_Info <- tibble(community_numb = as.numeric(), included_nodes = as.character(), overlapping_genes = as.character())
  i=1
  for (comm in cAMARETTOnetworkC$community_list) {
    Edges_Cnetwork_comm <- Edges_Cnetwork %>% 
      filter(from %in% comm & to %in% comm)
    all_overlapping_genes <- unique(unlist(strsplit(Edges_Cnetwork_comm$Overlapping_genes,", ")))
    Comm_Info <- add_row(Comm_Info, community_numb = as.numeric(names(cAMARETTOnetworkC$community_list[i])), included_nodes = paste(comm,collapse=", "), overlapping_genes = paste(all_overlapping_genes,collapse=", "))
    i = i + 1
  }
}

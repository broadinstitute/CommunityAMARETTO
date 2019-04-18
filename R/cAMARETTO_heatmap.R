#' @title cAMARETTO_heatmap
#'
#' Constructs a heatmap with the comparison between modules of two runs based on a hyper geometric test
#' 
#' @param cAMARETTOresults cAMARETTO_Run output
#' @param run1 name of run 1
#' @param run2 name of run 2
#'
#' @return None
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom tidyr spread
#' @importFrom tibble column_to_rownames
#' @importFrom grid gpar
#' @importFrom dplyr select filter
#' 
#' 
#' @examples 
#' cAMARETTO_heatmap(cAMARETTOresults, "LIHC", "GBM")
#' 
#' @export
cAMARETTO_heatmap<-function(cAMARETTOresults,run1,run2){
  results_filtered<-cAMARETTOresults$hgt_modules %>% dplyr::filter((grepl(run1,Geneset2)|grepl(run2,Geneset2))&(grepl(run1,Geneset1)|grepl(run2,Geneset1)))
  pvalue_matrix<-tidyr::spread(results_filtered %>% dplyr::select(Geneset2,Geneset1,p_value),key=Geneset1,value=p_value)
  pvalue_matrix<-tibble::column_to_rownames(pvalue_matrix,"Geneset2")
  pvalue_matrix<--log10(pvalue_matrix)
  Heatmap(pvalue_matrix, name = "Pvalues Compare Modules", column_title = "Regulator Genes\nExpression",show_column_names=TRUE,column_names_gp = gpar(fontsize = 8),row_names_gp = gpar(fontsize = 8),
          column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"), col=colorRamp2(c(0, max(pvalue_matrix,na.rm = TRUE)), c("white", "darkred")))
}

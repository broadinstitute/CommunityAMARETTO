#' @title cAMARETTO_heatmap
#'
#' @param cAMARETTOresults cAMARETTO_Run output
#' @param run1 name of run 1
#' @param run2 name of run 2
#'
#' @return a heatmap with the comparison between two runs
#' @import ComplexHeatmap
#' @import tidyverse
#' @import circlize
#' @export
cAMARETTO_heatmap<-function(cAMARETTOresults,run1,run2){
  results_filtered<-cAMARETTOresults$hgt_modules %>% filter((grepl(run1,Testset)|grepl(run2,Testset))&(grepl(run1,Geneset)|grepl(run2,Geneset)))
  pvalue_matrix<-spread(results_filtered %>% select(Testset,Geneset,p_value),key=Geneset,value=p_value)
  pvalue_matrix<-column_to_rownames(pvalue_matrix,"Testset")
  pvalue_matrix<--log10(pvalue_matrix)
  Heatmap(pvalue_matrix, name = "Pvalues Compare Modules", column_title = "Regulator Genes\nExpression",show_column_names=TRUE,column_names_gp = gpar(fontsize = 8),row_names_gp = gpar(fontsize = 8),
          column_title_gp = gpar(fontsize = 12, fontface = "bold"), col=colorRamp2(c(0, max(pvalue_matrix,na.rm = TRUE)), c("white", "darkred")))
}

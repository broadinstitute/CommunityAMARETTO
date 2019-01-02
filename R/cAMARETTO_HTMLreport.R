#' @title cAMARETTO_HTMLreport
#'
#' @param cAMARETTOresults
#' @param cAMARETTOnetworkM
#' @param cAMARETTOnetworkC
#' @param report_address
#'
#' @return a plot with one (or none) colored modules
#' 
#' @import igraph
#' @import DT
#' @import tidyverse
#' @import reshape2
#' @export
cAMARETTO_HTMLreport <- function(cAMARETTOresults,cAMARETTOnetworkM,cAMARETTOnetworkC,report_address="./"){
  
  #dataframe with modules per Run
  ComModulesLink <- stack(cAMARETTOnetworkC$community_list) %>% 
    rename(Community="ind", Module="values")
  
  #adding Modules that are not in Communities or Communities that are filtered out
  all_module_names <- unique(c(cAMARETTOresults$hgt_modules$Geneset,cAMARETTOresults$hgt_modules$Testset))
  all_module_names_fil<-all_module_names[!all_module_names %in% ComModulesLink$Module]
  
  Nodes_Mnetwork <- igraph::as_data_frame(cAMARETTOnetworkM$module_network, what="vertices")
  Module_no_Network <- all_module_names_fil[!all_module_names_fil %in% Nodes_Mnetwork$name]
  Module_no_Com <- all_module_names_fil[all_module_names_fil %in% Nodes_Mnetwork$name]
  
  ComModulesLink<-left_join(as.data.frame(all_module_names)%>% rename(Module="all_module_names"),ComModulesLink)
  
  ComModulesLink <- ComModulesLink %>% mutate(Community=ifelse(Module %in% Module_no_Network,"Not in Network",ifelse(Module %in% Module_no_Com, "Not in a Community",paste0("Community ", Community))))
  
  ComModulesLink <- ComModulesLink %>% 
    separate(Module, c("Run","ModuleNr"), "_Module_") %>% 
    mutate(ModuleNr = paste0("Module ", ModuleNr)) %>% 
    group_by(Community, Run) %>% 
    summarise(ModuleNrs=paste(ModuleNr, collapse = ", "))

  ComModulesLink <- dcast(ComModulesLink, Community~Run, fill=0)
  full_path <- normalizePath(report_address)
  rmarkdown::render(
    system.file("templates/TemplateIndexPage.Rmd", package = "CommunityAMARETTO"),
    output_dir = paste0(full_path, "/htmls/"),
    output_file = "index.html",
    params = list(
      cAMARETTOnetworkM = cAMARETTOnetworkM,
      cAMARETTOnetworkC = cAMARETTOnetworkC,
      ComModulesLink = ComModulesLink,
    ),
    quiet = TRUE
  )
}  
  
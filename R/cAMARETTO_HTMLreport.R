#' @title cAMARETTO_HTMLreport
#'
#' Creates a HTMLreport for the community AMARETTO results
#'
#' @param cAMARETTOresults The output of the Results function.
#' @param cAMARETTOnetworkM The output of the Module Network function.
#' @param cAMARETTOnetworkC The output of the Identify Communities function.
#' @param output_address The output repository for the HTML report.
#' @param HTMLsAMARETTOlist A list with AMARETTO reports to link with the Community AMARETTO report. If NULL, no links are added.
#' @param CopyAMARETTOReport Boolean to indicate if the AMARETTO reports needs to be copied in the AMARETTO report directory. In this way links are contained when moving the HTML directory.
#' @param hyper_geo_reference A reference gmt file to perform the Hyper Geometric Test.
#' @param NrCores Number of Cores to use during generation of the HTML report.
#' @param driverGSEA if TRUE, driver genes beside the target genes will also be included for hypergeometric test. 
#' @param PhenotypeTablesList List of Phenotype Association Tables for different AMARETTO runs.
#'
#' @return A set of HTMLs, giving caracteristics of the communities
#' 
#' @import igraph
#' @import DT
#' @importFrom rmarkdown render
#' @importFrom stringr str_order
#' @importFrom dplyr arrange group_by left_join mutate select summarise  rename  filter everything pull distinct
#' @importFrom tibble tibble rownames_to_column
#' @importFrom knitr knit_meta
#' @importFrom reshape2 dcast
#' @importFrom utils stack
#' @importFrom tidyr separate
#' @importFrom R.utils insert
#' @examples cAMARETTO_HTMLreport(cAMARETTOresults,cAMARETTOnetworkM, cAMARETTOnetworkC,HTMLsAMARETTOlist = HTMLsAMARETTOlist, hyper_geo_reference = gmtfile, output_address= "./")
#' 
#' @export
cAMARETTO_HTMLreport <- function(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC, PhenotypeTablesList = NULL,
                                 output_address="./", HTMLsAMARETTOlist=NULL, CopyAMARETTOReport = TRUE,
                                 hyper_geo_reference = NULL,
                                 driverGSEA=TRUE,NrCores=2){
  
  RunInfoList<-InitialCheckInputs(cAMARETTOresults,output_address,HTMLsAMARETTOlist,CopyAMARETTOReport,hyper_geo_reference)
  RunInfo<-RunInfoList$RunInfo
  RunInfo2<-RunInfoList$RunInfo2
  full_path<-RunInfoList$full_path
  HTMLsAMARETTOlist<-RunInfoList$HTMLsAMARETTOlist
  #==============================================================================================================
  # Extract main dataframes
  com_gene_df<-suppressWarnings(ComRunModGenInfo(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC))
  comm_info <-suppressWarnings(cAMARETTO_InformationTable(cAMARETTOnetworkM, cAMARETTOnetworkC))
  #==============================================================================================================
  Runs_AMARETTOs_info<-com_gene_df%>%dplyr::select(Run_Names,AMARETTOres)%>%dplyr::distinct()%>%dplyr::mutate(Run_Names = RunHyperLink(Run_Names,AMARETTOres,HTMLsAMARETTOlist,CopyAMARETTOReport))%>%filter(AMARETTOres==1)%>%select(Run_Names)
  #==============================================================================================================
  ComModulesLink<-CommunityModuleTableCreate(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC,HTMLsAMARETTOlist,CopyAMARETTOReport)
  #==============================================================================================================
  com_gene_df<-com_gene_df%>%dplyr::mutate(Color=sapply(as.numeric(Weights), function(x){
    if(is.na(x)){
      return("")
    }
    else if(x>0){
      return("darkred")
    }
    else if(x<0){
      return("darkblue")
    }
    else {
      return("darkgreen")
    }
  }))%>%dplyr::mutate(TypeColored=paste0('<font color=',Color,'>',Type,'</font>'))
  #==============================================================================================================
  #[stringr::str_order(Community, numeric = TRUE),]
  GeneComLink<-com_gene_df%>%filter(AMARETTOres==1)%>%dplyr::rename(GeneName = GeneNames)%>%
    dplyr::mutate(GeneName = paste0("<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=", GeneName, "\">", GeneName, "</a>"))%>%
    dplyr::select(c(GeneName,Community,TypeColored,Community_key,Community_type))%>%dplyr::rename(Type=TypeColored)%>%
    dplyr::mutate(Community = CommunityHyperLink(Community,Community_key,Community_type))%>%dplyr::arrange(GeneName)%>%select(-Community_key,-Community_type)%>%select(GeneName,everything())%>%distinct()
  #==============================================================================================================
  #GeneComLink<-GeneComLink[stringr::str_order(GeneComLink$Community, numeric = TRUE),]
  #==============================================================================================================
  #adding Community to driver genes table 
  Comm_Drivers<-com_gene_df%>%dplyr::filter(Type=="Driver")%>%
    dplyr::mutate(GeneNames=paste0("<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=", GeneNames, "\">", GeneNames, "</a>"))%>%
    dplyr::group_by(Community_key,Run_Names,Community_type,Community)%>%dplyr::summarise(Drivers=paste(unique(sort(GeneNames)),collapse = ", "))
  Comm_Drivers<-data.frame(Comm_Drivers)%>%dplyr::mutate(Community = CommunityHyperLink(Community,Community_key,Community_type))
  Comm_Drivers<-Comm_Drivers[stringr::str_order(Comm_Drivers$Community, numeric = TRUE),]
  #==============================================================================================================
  #HGT to test for gene set enrichment
   # avoid showing datatable size-related warnings.
  #==============================================================================================================
  # add phenotype table
  if (!is.null(PhenotypeTablesList)){
    phenotype_table_all<-CreatePhenotypeTable(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC, PhenotypeTablesList)
  }
  #==============================================================================================================
  # do hypergeometric test
  if(!is.null(hyper_geo_reference)){
    if (is.character(hyper_geo_reference)){
      all_hgt_output<-CreateHyperGeoTestAll(cAMARETTOresults,cAMARETTOnetworkM,cAMARETTOnetworkC,hyper_geo_reference,driverGSEA)
    }
    else if (is.data.frame(hyper_geo_reference)){
      all_hgt_output<-hyper_geo_reference
    }
    else{
      stop("hyper_geo_reference is not in the correct format!")
    }
  }

  #============================================================================================================== 
  options('DT.warn.size'=FALSE)
  buttons_list = list(list(extend ='csv'), list(extend ='excel'), list(extend = 'pdf', pageSize = 'A4', orientation = 'landscape'),list(extend ='print'), list(extend ='colvis'))
  columnDefs = list(list(className = 'dt-head-center', targets = "_all"),list(className = 'text-left',targets = "_all"))
  optionsList = list(deferRender=TRUE,
                     pageLength = 10,
                    lengthMenu = c(5, 10, 20, 50, 100),
                    keys = TRUE,
                    dom = "Blfrtip", 
                    buttons = buttons_list,
                    columnDefs = columnDefs,
                    paging = TRUE)
  #============================================================================================================== 
  # Community Pages : 
  for (ComNr in unique(com_gene_df$Community_key)){
    
    ModuleList<-com_gene_df%>%filter(Community_key==ComNr)%>%select(Run_Names,ModuleNr,AMARETTOres)%>%distinct()%>%mutate(ModuleNr=ModuleHyperLink(ModuleNr,Run_Names,AMARETTOres,HTMLsAMARETTOlist,CopyAMARETTOReport,page=2))%>%select(-AMARETTOres)
    
    DTML <- DT::datatable(ModuleList, 
                        class = "display",
                        filter = 'top',
                        extensions = c('Buttons','KeyTable'),
                        rownames = FALSE,
                        options = optionsList,
                        colnames = c("Data Set", "Module"),
                        escape = FALSE)
    
    #adding Gene-Module-Run tabel
    genelists_module<-com_gene_df%>%filter(AMARETTOres==1)%>%dplyr::filter(Community_key==ComNr)%>%dplyr::arrange(GeneNames)%>%dplyr::rename(Run=Run_Names)%>%dplyr::rename(ModuleName=ModuleNr)%>%dplyr::rename(Genes=GeneNames)%>%dplyr::select(-c(Type))%>%dplyr::rename(Type=TypeColored)
    genelists_module <- genelists_module%>%mutate(ModuleName = ModuleHyperLink(ModuleName,Run,AMARETTOres,HTMLsAMARETTOlist,CopyAMARETTOReport,page=2))%>%dplyr::mutate(Genes = paste0("<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=", Genes, "\">", Genes, "</a>"))%>%dplyr::select(c(Run,ModuleName,Genes,Type))

    DTGenes <- DT::datatable(genelists_module,
                             class = "display",
                             filter = 'top',
                             extensions = c('Buttons','KeyTable'),
                             rownames = FALSE,
                             options = optionsList,
                             colnames = c("Data Set", "Module", "Gene", "Gene Type"),
                             escape=FALSE)
    
    if(!is.null(hyper_geo_reference)) {
      DTGSEA<-create_hgt_datatable(all_hgt_output, com_table=TRUE, ComNr = ComNr)
  } else {
    DTGSEA <- "Genesets were not analysed as they were not provided."
  }
  # add phenotype table for each community page
  if (!is.null(PhenotypeTablesList)){
    phenotype_table_community<-phenotype_table_all%>%dplyr::filter(Community_key==ComNr)%>%select(Run_Names,everything())%>%dplyr::mutate(ModuleNr=ModuleHyperLink(ModuleNr,Run_Names,AMARETTOres,HTMLsAMARETTOlist,CopyAMARETTOReport,page=2))%>%select(-AMARETTOres)
    phenotype_table_community<-phenotype_table_community%>%dplyr::select(-Community_key,-Community,-Community_type)%>%arrange(q.value)
    DTPhC<-DT::datatable(phenotype_table_community,
                     class = "display",
                     filter = 'top',
                     extensions = c('Buttons','KeyTable'),
                     rownames = FALSE,
                     colnames=c("Data Set","Module","Phenotype","Statistics Test","P-value","FDR Q-value","Descriptive Statistics"),
                     options = optionsList, escape = FALSE)
  }
  else{ DTPhC = "Phenotype Statistical Analysis is not provided" }
  
  comm_name<-com_gene_df%>%filter(Community_key==ComNr)%>%select(Community,Community_key)%>%distinct()%>%pull(Community)
  if (grepl("Not in", comm_name)){
    ComTitle = comm_name
  }
  else{
    ComTitle = paste0("Community ",comm_name)
  }
  
  knitr::knit_meta(class=NULL, clean = TRUE)  # cleaning memory, avoiding memory to be overloaded
  rmarkdown::render(
      system.file("templates/community_page_template/TemplateCommunityPage.Rmd", package = "CommunityAMARETTO"),
      output_dir = paste0(full_path, "/communities"),
      output_file = paste0("Community_",ComNr,".html"),
      params = list(
        ComNr = ComNr,
        ComTitle = ComTitle,
        DTGSEA = DTGSEA,
        DTML = DTML,
        DTGenes = DTGenes,
        DTPhC = DTPhC,
        cAMARETTOnetworkM = cAMARETTOnetworkM,
        cAMARETTOnetworkC = cAMARETTOnetworkC
      ), quiet = TRUE)
  }
  file_remove<-suppressWarnings(suppressMessages(file.remove(paste0(full_path,"/communities/Community_",c(1:length(unique(com_gene_df$Community_key))),"_files"))))
  file_remove<-suppressWarnings(suppressMessages(file.remove(paste0(full_path,"index_files"))))
  #==============================================================================================================
  # index page : 
  DTRunInfo<-datatable(Runs_AMARETTOs_info,
            class = "display",
            filter = 'top',
            extensions = c('Buttons','KeyTable'),
            escape=FALSE,
            rownames = FALSE,
            options=optionsList,
            colnames = c("AMARETTO Report"))
  
  ComModulesLink<-ComModulesLink %>% dplyr::rename(Edges="numTotalEdgesInCommunity","Fraction Edges"="fractEdgesInVsOut","Fraction Size"="CommsizeFrac")
  DTComModulesLink<-datatable(ComModulesLink, 
            class = "display",
            filter = 'top',
            extensions = c('Buttons','KeyTable'),
            rownames = FALSE,
            options = optionsList,
            escape=FALSE) %>%
    formatSignif(c("Fraction Edges","Fraction Size"),2)

  DTGeneComLink<-datatable(GeneComLink, 
            class = "display",
            filter = 'top',
            extensions = c('Buttons','KeyTable'),
            rownames = FALSE,
            options = optionsList,
            colnames = c("Gene", "Community", "Gene Type"),
            escape=FALSE)

  if (!is.null(hyper_geo_reference)) {
    DTGSEAall <-create_hgt_datatable(all_hgt_output, com_table=FALSE)
    }
  else{
    DTGSEAall <- "Genesets were not analysed as they were not provided."
  }
  
  Comm_Drivers<-Comm_Drivers%>%select(-Community_key,-Community_type)%>%select(Community,everything())
  DTComDrivers<- DT::datatable(Comm_Drivers,
                          class = "display",
                          filter = 'top',
                          extensions = c('Buttons','KeyTable'),
                          rownames = FALSE,
                          options = optionsList,
                          colnames = c("Community", "Data Set", "Driver Genes"),
                          escape=FALSE)
  
  # add phenotype table for index page
  if (!is.null(PhenotypeTablesList)){
    phenotype_table_all<-phenotype_table_all%>%dplyr::mutate(ModuleNr=ModuleHyperLink(ModuleNr,Run_Names,AMARETTOres,HTMLsAMARETTOlist,CopyAMARETTOReport))%>%select(-AMARETTOres)
    phenotype_table_all<-phenotype_table_all%>%dplyr::mutate(Community = CommunityHyperLink(Community,Community_key,Community_type))%>%select(-Community_key,-Community_type)%>%select(Community,Run_Names,everything())%>%arrange(q.value,Community)
    DTPh<-DT::datatable(phenotype_table_all,
                     class = "display",
                     filter = 'top',
                     extensions = c('Buttons','KeyTable'),
                     rownames = FALSE,
                     options = optionsList,
                     colnames=c("Community","Data Set","Module","Phenotype","Statistics Test","P-value","FDR Q-value","Descriptive Statistics"),
                     escape=FALSE)
  }
  else{ DTPh = "Phenotype Statistical Analysis is not provided" }
  
  driversFreqTbl<-DriversSharedTbl(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC)%>%left_join(com_gene_df%>%select(Community_key,Community,Community_type)%>%distinct(),by="Community_key")%>%dplyr::mutate(Community = CommunityHyperLink(Community,Community_key,Community_type))%>%select(-Community_key,-Community_type)%>%select(Community,everything())
  DTdriverFreq<- DT::datatable(driversFreqTbl,
                               class = "display",
                               filter = 'top',
                               extensions = c('Buttons','KeyTable'),
                               rownames = FALSE,
                               options = optionsList, escape=FALSE)
  
  rmarkdown::render(
    system.file("templates/TemplateIndexPage.Rmd", package = "CommunityAMARETTO"),
    output_dir = full_path,
    output_file = "index.html",
    params = list(
      cAMARETTOnetworkM = cAMARETTOnetworkM,
      cAMARETTOnetworkC = cAMARETTOnetworkC,
      DTComModulesLink = DTComModulesLink,
      DTRunInfo =  DTRunInfo
    ), quiet = TRUE)
  
  rmarkdown::render(
    system.file("templates/TemplateIndexPage_RunsInfo.Rmd", package = "CommunityAMARETTO"),
    output_dir = full_path,
    output_file = "index_RunsInfo.html",
    params = list(
      DTRunInfo =  DTRunInfo
    ), quiet = TRUE)
  
  rmarkdown::render(
    system.file("templates/TemplateIndexPage_AllCommunities.Rmd", package = "CommunityAMARETTO"),
    output_dir = full_path,
    output_file = "index_AllCommunities.html",
    params = list(
      DTRunInfo =  DTRunInfo
    ), quiet = TRUE)
  
  rmarkdown::render(
    system.file("templates/TemplateIndexPage_Drivers.Rmd", package = "CommunityAMARETTO"),
    output_dir = full_path,
    output_file = "index_Drivers.html",
    params = list(
      DTdriverFreq = DTdriverFreq
    ), quiet = TRUE)
  
  rmarkdown::render(
    system.file("templates/TemplateIndexPage_Drivers2.Rmd", package = "CommunityAMARETTO"),
    output_dir = full_path,
    output_file = "index_Drivers2.html",
    params = list(
      DTComDrivers=DTComDrivers,
    ), quiet = TRUE)
  
  rmarkdown::render(
    system.file("templates/TemplateIndexPage_AllGenes.Rmd", package = "CommunityAMARETTO"),
    output_dir = full_path,
    output_file = "index_AllGenes.html",
    params = list(
      DTGeneComLink = DTGeneComLink
    ), quiet = TRUE)
  
  rmarkdown::render(
    system.file("templates/TemplateIndexPage_GenesetsEnrichment.Rmd", package = "CommunityAMARETTO"),
    output_dir = full_path,
    output_file = "index_GenesetsEnrichment.html",
    params = list(
      DTGSEAall = DTGSEAall
    ), quiet = TRUE)

  rmarkdown::render(
    system.file("templates/TemplateIndexPage_PhenoAssociation.Rmd", package = "CommunityAMARETTO"),
    output_dir = full_path,
    output_file = "index_PhenoAssociation.html",
    params = list(
      DTPh = DTPh
    ), quiet = TRUE)
  
}

#' @title HGTGeneEnrichmentList
#'
#' Calculates the p-values for unranked gene set enrichment based on two gmt files as input and the hyper geometric test.
#'
#' @param genelist The gmt file with reference gene set.
#' @param gmtfile The gmt file with gene sets to test. In our case, the gmt file of the modules.
#' @param NrCores Number of cores used for parallelization.
#' @param ref.numb.genes The total number of genes teste, standard equal to 45956 (MSIGDB standard).
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom stats p.adjust phyper
#' @keywords internal
#' @export
HGTGeneEnrichmentList <- function(genelist, gmtfile, NrCores, ref.numb.genes = 45956) {
    gmtset <- readGMT(gmtfile)  # the hallmarks_and_co2...
    ########################### Parallelizing :
    cluster <- parallel::makeCluster(c(rep("localhost", NrCores)), type = "SOCK")
    doParallel::registerDoParallel(cluster, cores = NrCores)
    resultHGT<-foreach::foreach(i = 1:length(gmtset), .combine = "rbind") %dopar% {
        l <- length(gmtset[[i]])
        # k <- sum(gmtset[[i]] %in% genelist)
        k <- length(intersect(gmtset[[i]],genelist))
        m <- ref.numb.genes
        n <- length(genelist)
        p1 <- stats::phyper(k - 1, l, m - l, n, lower.tail = FALSE)
          
        if (k > 0) {
          overlapping.genes <- intersect(gmtset[[i]],genelist)
          #overlapping.genes <- gmtset[[i]][gmtset[[i]] %in% genelist]
          overlapping.genes <- paste(overlapping.genes, collapse = ", ")
          c(Geneset = names(gmtset[i]),
            Testset = names(genelist),
            p_value = p1,
            n_Overlapping = k,
            n_RefGeneset = l,
            Overlapping_genes = overlapping.genes)
        }
    }
    
    parallel::stopCluster(cluster)
    resultHGT <- as.data.frame(resultHGT, stringsAsFactors = FALSE)
    resultHGT$p_value <- as.numeric(resultHGT$p_value)
    resultHGT$n_Overlapping <- as.numeric(resultHGT$n_Overlapping)
    resultHGT$n_RefGeneset<-as.numeric(resultHGT$n_RefGeneset)
    resultHGT[, "padj"] <- stats::p.adjust(resultHGT[, "p_value"], method = "BH")
    return(resultHGT)
}

#' Title ComRunModGenInfo
#'
#' @param cAMARETTOresults 
#' @param cAMARETTOnetworkM 
#' @param cAMARETTOnetworkC 
#'
#' @import igraph
#' @importFrom dplyr arrange rename left_join mutate
#' @importFrom utils stack 
#' @importFrom purrr map
#' @return a dataframe contaning all communities, runname, and modules relationships. 
#' @export
#'
#' @examples ComRunModGenInfo(cAMARETTOresults,cAMARETTOnetworkM,cAMARETTOnetworkC)
ComRunModGenInfo<-function(cAMARETTOresults,cAMARETTOnetworkM,cAMARETTOnetworkC){
  ComModulesLink <- utils::stack(cAMARETTOnetworkC$community_list) %>% dplyr::rename(Module="values", Community="ind")
  #all_module_names <- unique(c(cAMARETTOresults$hgt_modules$Geneset1,cAMARETTOresults$hgt_modules$Geneset2))
  all_module_names<-unique(cAMARETTOresults$all_genes_modules_df%>%mutate(Module=paste0(Run_Names,"|",ModuleNr))%>%pull(Module))
  modules_with_community<-ComModulesLink$Module
  suppressWarnings(ComModulesLink <- dplyr::left_join(as.data.frame(all_module_names) %>% dplyr::rename(Module="all_module_names"),ComModulesLink, by="Module"))
  ComModulesLink<-ComModulesLink%>%dplyr::mutate(ModuleNr=unlist(purrr::map(strsplit(Module,"\\|"),2)))%>% dplyr::mutate(Run_Names=unlist(purrr::map(strsplit(Module,"\\|"),1)))
  # so far we have all the modules with their community information, NA for modules with no community.
  Nodes_Mnetwork <- igraph::as_data_frame(cAMARETTOnetworkM$module_network, what="vertices")
  # Nodes_Mnetwork = all modules selected to be in the network
  Module_no_Network <-setdiff(all_module_names,Nodes_Mnetwork$name)
  Module_no_Com<-setdiff(Nodes_Mnetwork$name,modules_with_community)
  suppressMessages(ComModulesLink <-dplyr::left_join(cAMARETTOresults$all_genes_modules_df,ComModulesLink,by=c("Run_Names","ModuleNr")))
  
  ComModulesLink <- ComModulesLink %>%dplyr::arrange(as.numeric(Community))%>%dplyr::mutate(Community=ifelse(Module %in% Module_no_Network,paste0("Not in Network ",Run_Names),ifelse(Module %in% Module_no_Com, paste0("Not in Community ",Run_Names),Community)))

  ComModulesLink <- suppressWarnings(ComModulesLink %>%dplyr::left_join(data.frame(Community=unique(ComModulesLink$Community),Community_key=1:length(unique(ComModulesLink$Community))),by="Community"))
  ComModulesLink<-ComModulesLink%>%mutate(Community_type=sapply(Community, function(x){
    
    if (grepl("Not in Network", x)){
      return(1)
    }
    else if (grepl("Not in Community", x)){
      return(2)
    }
    else{
      return(0)
    }
  }))
  #ComModulesLink<-ComModulesLink%>%select(-Community) 
  return(ComModulesLink)
}

#' Title CreatePhenotypeTable
#'
#' @param cAMARETTOresults 
#' @param cAMARETTOnetworkM 
#' @param PhenotypeTablesList 
#' @param cAMARETTOnetworkC 
#'
#' @importFrom dplyr select distinct mutate mutate left_join select
#' @return results
#' @export
#'
#' @examples CreatePhenotypeTable(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC, PhenotypeTables)
CreatePhenotypeTable<-function(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC, PhenotypeTablesList){
  phenotype_table_all<-NULL
  CommunityRunModuleTable<-ComRunModGenInfo(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC)%>%dplyr::select(Run_Names,ModuleNr,Community,Community_key,Community_type,AMARETTOres)%>%dplyr::distinct()
  for (i in 1:length(PhenotypeTablesList)){
    if (is.null(PhenotypeTablesList[[i]])){
      next
    }
    phenotype_table<-PhenotypeTablesList[[i]]%>%dplyr::mutate(ModuleNr=gsub("Module ","Module_",ModuleNr))%>%dplyr::mutate(Run_Names=names(PhenotypeTablesList)[i])%>%dplyr::left_join(CommunityRunModuleTable,by=c("Run_Names","ModuleNr"))%>%dplyr::select(Community,Community_key,Community_type,AMARETTOres,Run_Names,ModuleNr,Phenotypes,Statistical_Test,p.value,q.value,Descriptive_Statistics)
    phenotype_table_all<-rbind(phenotype_table_all,phenotype_table)
  }
  return(phenotype_table_all)
}

#' Title CreateHyperGeoTestAll
#'
#' @param cAMARETTOresults 
#' @param cAMARETTOnetworkM 
#' @param cAMARETTOnetworkC 
#' @param hyper_geo_reference 
#' @param driverGSEA 
#' @param NrCores 
#'
#' @return Geneset Enrichment Analysis for the entire communities. 
#' @export
#'
#' @examples CreateHyperGeoTestAll(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC, hyper_geo_reference = './h.all.v6.2.symbols.gmt', driverGSEA=TRUE,NrCores=4)
CreateHyperGeoTestAll<-function(cAMARETTOresults,cAMARETTOnetworkM,cAMARETTOnetworkC,hyper_geo_reference,driverGSEA=TRUE,NrCores=4){
  print("Performing Geneset Enrichment Analysis ...")
  com_gene_df<-suppressWarnings(ComRunModGenInfo(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC))
  communities_all<-unique(com_gene_df$Community_key)
  all_hgt_output<-NULL
  if (is.null(cAMARETTOresults)){return(1)}
  for(i in 1:length(hyper_geo_reference)){
    for (ComNr in communities_all){
      print(ComNr)
      target_genes<-com_gene_df%>%dplyr::filter(Community_key==ComNr)%>%dplyr::filter(Type=="Target")%>%dplyr::arrange(GeneNames)%>%dplyr::pull(GeneNames)
      driver_genes<-com_gene_df%>%dplyr::filter(Community_key==ComNr)%>%dplyr::filter(Type=="Driver")%>%dplyr::arrange(GeneNames)%>%dplyr::pull(GeneNames)
      driver_genes_weights<-com_gene_df%>%dplyr::filter(Community_key==ComNr)%>%dplyr::filter(Type=="Driver")%>%dplyr::arrange(GeneNames)%>%dplyr::pull(Weights)
      if(driverGSEA){
        genelist<-unique(c(target_genes,driver_genes))
      }
      else{
        genelist<-unique(target_genes)
      }
      outputHGT <- HGTGeneEnrichmentList(genelist, hyper_geo_reference[i], NrCores = NrCores)
      if (nrow(outputHGT)>0){

        all_hgt_output<-rbind(all_hgt_output, outputHGT%>%mutate(Community_key=ComNr)) 
      }
      else{
        all_hgt_output<-rbind(all_hgt_output, NULL) 
      }
    }
    cat("The hyper geometric test results are calculated.\n")
  }
  print(head(all_hgt_output))
  print(dim(all_hgt_output))
  #Community_key
  utils::data(MsigdbMapping)
  MsigdbMapping<-MsigdbMapping%>%dplyr::mutate(url=paste0('<a href="http://software.broadinstitute.org/gsea/msigdb/cards/',geneset,'">',gsub("_"," ",geneset),'</a>'))
  all_hgt_output<-all_hgt_output%>%dplyr::left_join(MsigdbMapping,by=c("Geneset"="geneset"))%>%
    dplyr::mutate(description=ifelse(is.na(description),Geneset,description))%>%
    dplyr::mutate(Geneset=ifelse(is.na(url),Geneset,url))%>%dplyr::rename("Description"="description")%>%dplyr::select(-url)
  
  all_hgt_output <- all_hgt_output %>% dplyr::mutate(overlap_perc = n_Overlapping / n_RefGeneset) %>% dplyr::select(Community_key,Geneset, Description,n_RefGeneset, n_Overlapping, Overlapping_genes, overlap_perc, p_value, padj) %>% dplyr::arrange(padj)
  
  all_hgt_output<-all_hgt_output%>%left_join(com_gene_df%>%select(Community_key,Community,Community_type)%>%distinct(),by="Community_key")
  print("Geneset Enrichment Analysis is done!")
  return(all_hgt_output)
}

#' Title InitialCheckInputs
#'
#' @param cAMARETTOresults 
#' @param output_address 
#' @param HTMLsAMARETTOlist 
#' @param CopyAMARETTOReport 
#' @param hyper_geo_reference 
#'
#' @return RunInfo dataframe
#' @export
#'
#' @examples InitialCheckInputs(cAMARETTOresults,output_address="./",HTMLsAMARETTOlist,CopyAMARETTOReport=FALSE,hyper_geo_reference)
InitialCheckInputs<-function(cAMARETTOresults,output_address,HTMLsAMARETTOlist,CopyAMARETTOReport,hyper_geo_reference){
  if (!dir.exists(output_address)) {
    stop("Output directory is not existing.")
  } else {
    dir.create(file.path(output_address,"htmls"), showWarnings = FALSE)
    full_path <- file.path(normalizePath(output_address),"htmls")
    print(paste0("The output directory is: ",full_path))
  }
  if (!is.null(hyper_geo_reference)) {
    if(is.character(hyper_geo_reference)){
      if (!file.exists(hyper_geo_reference[1])) {
        stop("GMT for hyper geometric test is not existing.")
      }
    }
  }
  
  i=1
  if(is.null(HTMLsAMARETTOlist)==FALSE){
    if(!all(names(HTMLsAMARETTOlist) %in% cAMARETTOresults$runnames)==TRUE){
      stop("The RUN names don't match those of the cAMARETTOresults")
    }
    for(htmldir in HTMLsAMARETTOlist){
      if(!file.exists(htmldir)){
        stop(paste0("The AMARETTO ",names(HTMLsAMARETTOlist)[i] ," html directory is not existing."))
      }
      htmldir<-normalizePath(file.path(htmldir,"/AMARETTOhtmls/"))
      if (CopyAMARETTOReport==TRUE){
        dir.create(file.path(full_path,names(HTMLsAMARETTOlist)[i]),showWarnings = FALSE)
        file.copy(htmldir, file.path(full_path,names(HTMLsAMARETTOlist)[i]), recursive = TRUE)
        htmldir<-file.path(".",names(HTMLsAMARETTOlist)[i],"AMARETTOhtmls")      
      }
      HTMLsAMARETTOlist[i]<-htmldir
      i=i+1
    }
  }
  else{
    HTMLsAMARETTOlist=NULL
  }
  return(list(full_path=full_path,HTMLsAMARETTOlist=HTMLsAMARETTOlist))
}
  
#' Title CommunityModuleTableCreate
#'
#' @param cAMARETTOresults 
#' @param cAMARETTOnetworkM 
#' @param HTMLsAMARETTOlist 
#' @param CopyAMARETTOReport 
#' @param cAMARETTOnetworkC 
#'
#' @return Community To Module Dataframe
#' @export
#'
#' @examples CommunityModuleTableCreate (cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC, RunInfo)
CommunityModuleTableCreate<-function(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC,HTMLsAMARETTOlist,CopyAMARETTOReport){
  com_gene_df<-suppressWarnings(ComRunModGenInfo(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC))%>%dplyr::select(Community_key,Community,Community_type,AMARETTOres,Run_Names,ModuleNr)%>%distinct()
  ComModule<-com_gene_df %>% dplyr::mutate(ModuleLink=ModuleHyperLink(ModuleNr,Run_Names,AMARETTOres,HTMLsAMARETTOlist,CopyAMARETTOReport))
  ComModule<-ComModule%>%group_by(Community_key,Run_Names) %>% summarise(ModuleLinks=paste(ModuleLink,collapse = ", "))
  ComModule<-suppressMessages(reshape2::dcast(ComModule, Community_key~Run_Names, fill=0))
  ComModule<-suppressMessages(ComModule%>%left_join(com_gene_df%>%select(Community_key,Community,Community_type)%>%distinct()))
  ComModule<-ComModule%>%left_join(cAMARETTOnetworkC$commEdgeInfo%>%mutate(Community_key=Community)%>%select(Community_key,numTotalEdgesInCommunity,fractEdgesInVsOut,CommsizeFrac),by="Community_key")
  ComModule<-ComModule%>%dplyr::mutate(Community = CommunityHyperLink(Community,Community_key,Community_type) )
  ComModule<-ComModule%>%select(Community,sort(cAMARETTOresults$runnames),everything())%>%select(-Community_key,-Community_type)
  return(ComModule)
}

#' Title CommunityHyperLink
#'
#' @param Community CommunityHyperLink
#' @param Community_key 
#' @param Community_type 
#'
#' @return a hyperlink and presentable name for the communities used for datatables
#' @export
#'
#' @examples CommunityHyperLink(Community,Community_key,Community_type)
CommunityHyperLink<-function(Community,Community_key,Community_type){
  Comm_hyperLink<-paste0("<a href=\"./communities/",paste0("Community","_",Community_key),".html\">",paste0(ifelse((Community_type==0),"Community ",""),Community), "</a>")
  return(Comm_hyperLink)
}



#' Title DriversSharedTbl
#'
#' @param cAMARETTOresults 
#' @param cAMARETTOnetworkM 
#' @param cAMARETTOnetworkC 
#'
#' @return  Frequency of driver genes across different communities and datasets
#' @export
#'
#' @examples  DriversSharedTbl(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC)
DriversSharedTbl<-function(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC){
  com_gene_df<-suppressWarnings(ComRunModGenInfo(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC))
  ComDrivers<-com_gene_df%>%filter(Type=="Driver")%>%mutate(GeneNames = paste0("<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=", GeneNames, "\">", GeneNames, "</a>"))
  ComDrivers<-ComDrivers%>%group_by(Community_key,Run_Names)%>%summarise(GeneNames=paste(sort(unique(GeneNames)),collapse = ", "))
  ComDrivers<-suppressMessages(reshape2::dcast(ComDrivers, Community_key~Run_Names, fill=""))
  ComDrivers<-ComDrivers%>%left_join(ComDrivers%>%tidyr::unite("all_genes",-one_of("Community_key"),sep = ", "),by="Community_key")
  gene_freq_df<-NULL
  for (ComNr in unique(ComDrivers$Community_key)){
    gene_table<-table(strsplit(gsub("","",ComDrivers%>%filter(Community_key==ComNr)%>%pull(all_genes)),","))
    freq_tbl<-data.frame(GeneNames=names(gene_table),Freq=paste0("# Data Sets = ",as.numeric(gene_table)),Community_key=ComNr)
    gene_freq_df<-rbind(gene_freq_df,freq_tbl)
  }
  gene_freq_df<-gene_freq_df%>%group_by(Community_key,Freq)%>%summarise(GeneNames=paste(GeneNames,collapse = ", "))
  gene_freq_df<-suppressMessages(reshape2::dcast(gene_freq_df, Community_key~Freq, fill=""))
  # order columns
  columnnames<-c("Community_key",sort(colnames(gene_freq_df)[-1]))
  gene_freq_df<-gene_freq_df[,columnnames] 
  ComDrivers<-ComDrivers%>%left_join(gene_freq_df,by="Community_key")%>%select(-all_genes)
  return(ComDrivers)
}


#' Title ModuleHyperLink
#'
#' @param Module 
#' @param Run_Names 
#' @param AMARETTOres 
#' @param HTMLsAMARETTOlist 
#' @param CopyAMARETTOReport 
#' @param page 
#'
#' @return Hyperlinks for Modules
#' @export
#'
#' @examples ModuleHyperLink(Module,Run_Names,AMARETTOres,HTMLsAMARETTOlist,CopyAMARETTOReport,page=1)
ModuleHyperLink<-function(Module,Run_Names,AMARETTOres,HTMLsAMARETTOlist,CopyAMARETTOReport,page=1){
  # htmldir = "someLocalAdress/LIHC_Report_75/AMARETTOhtmls"
  # if CopyAMARETTOReport ==TRUE :
  #   htmldir = "./TCGA_LIHC/AMARETTOhtmls"
  if(is.null(HTMLsAMARETTOlist)==FALSE){
    RunInfo<-(data.frame(Run_Names=names(HTMLsAMARETTOlist),ModuleLink=HTMLsAMARETTOlist,stringsAsFactors = FALSE))
    if(CopyAMARETTOReport==TRUE & page == 2){
      RunInfo <- RunInfo %>% mutate(ModuleLink=gsub("^./","../",ModuleLink))
    }
    Mod_hyperLink<-data.frame(Run_Names=Run_Names,Module=Module,AMARETTOres=AMARETTOres,stringsAsFactors = FALSE)%>%dplyr::left_join(RunInfo,by="Run_Names")%>%mutate(mod_hyperlink=ifelse(AMARETTOres==1,paste0("<a href=",ModuleLink,"/modules/",gsub("Module_","module",Module),".html>",gsub("Module_","Module ",Module), "</a>"),Module))%>%dplyr::pull(mod_hyperlink)
  }
  else{
    Mod_hyperLink<-ifelse(AMARETTOres==1,gsub("Module_","Module ",Module),Module)
  }
  return(Mod_hyperLink)
}


#' Title RunHyperLink
#'
#' @param Run_Names 
#' @param AMARETTOres 
#' @param HTMLsAMARETTOlist 
#' @param CopyAMARETTOReport 
#' @param page 1 for index page and 2 for community page
#'
#' @return Hyperlinks for Run Names
#' @export
#'
#' @examples RunHyperLink(Run_Names, AMARETTOres, HTMLsAMARETTOlist, CopyAMARETTOReport, page=1)
#' 
RunHyperLink<-function(Run_Names,AMARETTOres,HTMLsAMARETTOlist,CopyAMARETTOReport,page=1){
  if(is.null(HTMLsAMARETTOlist)==FALSE){
    RunInfo<-(data.frame("Run_Names"=names(HTMLsAMARETTOlist),RunLink=HTMLsAMARETTOlist,stringsAsFactors = FALSE))
    if(CopyAMARETTOReport==TRUE & page == 2){
      RunInfo <- RunInfo %>% mutate(RunLink=gsub("^./","../",RunLink))
    }
    Run_hyperLink<-data.frame("Run_Names"=Run_Names,AMARETTOres=AMARETTOres,stringsAsFactors = FALSE)%>%dplyr::left_join(RunInfo,by="Run_Names")%>%mutate(run_hyperlink=ifelse(AMARETTOres==1,paste0("<a href=",RunLink,"/index.html>",Run_Names, "</a>"),Run_Names))%>%dplyr::pull(run_hyperlink)
  }
  else{
    Run_hyperLink<-Run_Names
  }
  return(Run_hyperLink)
}

#' Title cAMARETTO_Cytoscape
#'
#' @param cAMARETTOsList A list containing Community AMARETTO results with the following format : list(cAMARETTOresults = cAMARETTOresults,cAMARETTOnetworkM = cAMARETTOnetworkM, cAMARETTOnetworkC = cAMARETTOnetworkC)
#' @param communityReportURL Optional URL linked to the Community-AMARETTO HTML report index page. 
#' @param cytoscape_name A character, for naming the network which will be shown in cytoscape. 
#'
#' @return result
#' @import RCy3
#' @import igraph
#' @importFrom  purrr map
#' @import tidyverse
#' @export
#'
#' @examples  
#' cytoscape_name<-"cAMARETTO_Liver2DS"
#  cAMARETTOsList<-readRDS(file="./outputs/cAMARETTO_Liver2DS.rds")
#  communityReportURL<-"http://portals.broadinstitute.org/pochetlab/demo/cAMARETTO_Liver_2DS/" 
#' cAMARETTO_Cytoscape(cAMARETTOsList,communityReportURL = "",cytoscape_name="my_cytoscape")
cAMARETTO_Cytoscape<-function(cAMARETTOsList,communityReportURL = "",cytoscape_name="my_cytoscape"){

  graph<-cAMARETTOsList$cAMARETTOnetworkC$CommGraph
  runnames<-cAMARETTOsList$cAMARETTOresults$runnames
  runURLs<-data.frame(run=runnames)
  runURLs<-runURLs%>%dplyr::mutate(run_URL=paste0(communityReportURL,run))
  
  nodes_df<-igraph::as_data_frame(graph, what="vertices")
  nodes_df<-suppressWarnings(nodes_df%>%dplyr::left_join(runURLs,by="run"))
  nodes_df<-nodes_df%>%dplyr::mutate(Module_name=unlist(map(strsplit(name,"\\|"),2)))%>%dplyr::mutate(Module_name=gsub("Module_","module",Module_name))
  nodes_df<-nodes_df%>%dplyr::mutate(URL=ifelse(run %in% runnames,paste0('<a href="',run_URL,'/AMARETTOhtmls/modules/',Module_name,'.html">URL</a>'),''))
  nodes_df<-nodes_df%>%dplyr::select(-c("Module_name","run_URL","NewComNumber"))
  
  edges_df<-igraph::as_data_frame(graph, what="edges")
  edges_df<-suppressWarnings(edges_df%>%dplyr::mutate(from_Run=unlist(purrr::map(strsplit(from,"\\|"),1)))%>%dplyr::mutate(from_Module=unlist(map(strsplit(from,"\\|"),2)))%>%dplyr::mutate(from_Module=gsub("Module_","module",from_Module))%>%dplyr::left_join(runURLs,by=c("from_Run"="run"))%>%dplyr::rename(from_run_URL=run_URL)%>%dplyr::mutate(from_URL=ifelse(from_Run %in% runnames,paste0('<a href="',from_run_URL,'/AMARETTOhtmls/modules/',from_Module,'.html">URL</a>'),'')))
  edges_df<-suppressWarnings(edges_df%>%dplyr::mutate(to_Run=unlist(purrr::map(strsplit(to,"\\|"),1)))%>%dplyr::mutate(to_Module=unlist(map(strsplit(to,"\\|"),2)))%>%dplyr::mutate(to_Module=gsub("Module_","module",to_Module))%>%dplyr::left_join(runURLs,by=c("to_Run"="run"))%>%dplyr::rename(to_run_URL=run_URL)%>%dplyr::mutate(to_URL=ifelse(to_Run %in% runnames,paste0('<a href="',to_run_URL,'/AMARETTOhtmls/modules/',to_Module,'.html">URL</a>'),'')))
  edges_df<-edges_df%>%dplyr::select(-c("from_Run","from_Module","from_run_URL","to_Run","to_Module","to_run_URL"))%>%dplyr::rename(source_URL=from_URL,target_URL=to_URL)

  graph_new<-igraph::graph_from_data_frame(edges_df, directed=FALSE, vertices=nodes_df)
  try(RCy3::createNetworkFromIgraph(graph_new,cytoscape_name),silent=TRUE)
  return(graph_new)
}



#' Title create_hgt_datatable
#'
#' @param output_hgt 
#' @param com_table 
#' @param ComNr 
#'
#' @return DataTable
#'
#' @examples 
create_hgt_datatable<-function(output_hgt, com_table=FALSE, ComNr = 1){
  if (com_table){
    outputHGT<-all_hgt_output%>%filter(Community_key==ComNr)%>%select(-Community_key,-Community_type,-Community)
    if (nrow(outputHGT)>0){
      DTGSEA_colnames<-c("Gene Set Name","Gene Set Description","# Genes in Gene Set","# Genes in Overlap","Genes in Overlap","% Genes in overlap","P-value","FDR Q-value")
      outputHGT<-outputHGT%>%dplyr::arrange(padj)
      DTGSEA <- DT::datatable(outputHGT,
                              class = "display",
                              filter = 'top',
                              extensions = c('Buttons','KeyTable'),
                              rownames = FALSE,
                              options = optionsList,
                              colnames = DTGSEA_colnames , escape = FALSE) %>% 
        DT::formatSignif(c("p_value", "padj","overlap_perc"), 2) %>% 
        DT::formatStyle("overlap_perc", background = styleColorBar(c(0, 1), "lightblue"), backgroundSize = "98% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "center")%>%DT::formatStyle(columns = c(5), fontSize = '60%')
    } 
    else{
      DTGSEA <- "No significant overlap was identified in the geneset enrichment analysis."
    }
   return(DTGSEA)
  }
  else{
    all_hgt_output<-all_hgt_output%>%filter(n_Overlapping>2)%>%dplyr::mutate(Community=CommunityHyperLink(Community,Community_key,Community_type))%>%select(-Community_type,-Community_key)%>%select(Community,everything())%>% dplyr::arrange(padj,Community)
    DTGSEA_colnames<-c("Community","Gene Set Name","Gene Set Description","# Genes in Gene Set","# Genes in Overlap","Genes in Overlap","% Genes in overlap","P-value","FDR Q-value")
    all_hgt_output<-all_hgt_output %>% dplyr::mutate(Geneset = paste0("<a href=\"http://software.broadinstitute.org/gsea/msigdb/cards/", Geneset, "\">", gsub("_", " ", Geneset),"</a>"))
    DTGSEAall <- DT::datatable(all_hgt_output,
                               class = "display",
                               filter = 'top',
                               extensions = c('Buttons','KeyTable'),
                               rownames = FALSE,
                               options = optionsList,
                               colnames = DTGSEA_colnames, escape = FALSE)%>%
      DT::formatSignif(c("p_value", "padj","overlap_perc"), 2) %>% 
      DT::formatStyle("overlap_perc", background = styleColorBar(c(0, 1), "lightblue"), backgroundSize = "98% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "center")%>%DT::formatStyle(columns = c(6), fontSize = '60%')
    return(DTGSEAall)
  }
}

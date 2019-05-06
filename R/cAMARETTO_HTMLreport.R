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
#' @param hyper_geo_test_bool Boolean if Hyper Geometric Test needs to be performed.
#' @param hyper_geo_reference A reference gmt file to perform the Hyper Geometric Test.
#' @param MSIGDB Boolean if gmt is MSIGDB derived.
#' @param NrCores Number of Cores to use during generation of the HTML report.
#' @param driverGSEA if TRUE, driver genes beside the target genes will also be included for hypergeometric test. 
#' @param PhenotypeTablesList 
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
#' @examples
#' 
#' cAMARETTO_HTMLreport(cAMARETTOresults,cAMARETTOnetworkM, cAMARETTOnetworkC,HTMLsAMARETTOlist = HTMLsAMARETTOlist, hyper_geo_test_bool = TRUE, hyper_geo_reference = gmtfile, MSIGDB = TRUE, output_address= "./")
#' 
#' @export
cAMARETTO_HTMLreport <- function(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC, PhenotypeTablesList = NULL,
                                 output_address="./", HTMLsAMARETTOlist=NULL, CopyAMARETTOReport = TRUE,
                                 hyper_geo_test_bool = FALSE, hyper_geo_reference = NULL,
                                 MSIGDB = FALSE,driverGSEA=TRUE,
                                 NrCores=2){
  
  RunInfoList<-InitialCheckInputs(cAMARETTOresults,output_address,HTMLsAMARETTOlist,CopyAMARETTOReport,hyper_geo_test_bool,hyper_geo_reference)
  RunInfo<-RunInfoList$RunInfo
  RunInfo2<-RunInfoList$RunInfo2
  full_path<-RunInfoList$full_path
  HTMLsAMARETTOlist<-RunInfoList$HTMLsAMARETTOlist
  print("ali1")
  #==============================================================================================================
  # Extract main dataframes
  com_gene_df<-suppressWarnings(ComRunModGenInfo(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC))
  print("ali2")
  comm_info <-suppressWarnings(cAMARETTO_InformationTable(cAMARETTOnetworkM, cAMARETTOnetworkC))
  print("ali3")
  #==============================================================================================================
  Runs_AMARETTOs_info<-com_gene_df%>%dplyr::select(Run_Names,AMARETTOres)%>%dplyr::distinct()%>%dplyr::mutate(Run_Names = RunHyperLink(Run_Names,AMARETTOres,HTMLsAMARETTOlist,CopyAMARETTOReport))%>%filter(AMARETTOres==1)%>%select(Run_Names)
  #==============================================================================================================
  ComModulesLink<-CommunityModuleTableCreate(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC,HTMLsAMARETTOlist,CopyAMARETTOReport)
  #==============================================================================================================
  print("ali4")
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
  print("ali5")
  GeneComLink<-com_gene_df%>%filter(AMARETTOres==1)%>%dplyr::rename(GeneName = GeneNames)%>%
    dplyr::mutate(GeneName = paste0("<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=", GeneName, "\">", GeneName, "</a>"))%>%
    dplyr::select(c(GeneName,Community,TypeColored,Community_key,Community_type))%>%dplyr::rename(Type=TypeColored)%>%
    dplyr::mutate(Community = CommunityHyperLink(Community,Community_key,Community_type))%>%dplyr::arrange(GeneName)%>%select(-Community_key,-Community_type)%>%select(GeneName,everything())%>%distinct()
  #==============================================================================================================
  #GeneComLink<-GeneComLink[stringr::str_order(GeneComLink$Community, numeric = TRUE),]
  #==============================================================================================================
  #adding Community to driver genes table 
  print("ali6")
  Comm_Drivers<-com_gene_df%>%dplyr::filter(Type=="Driver")%>%
    dplyr::mutate(GeneNames=paste0("<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=", GeneNames, "\">", GeneNames, "</a>"))%>%
    dplyr::group_by(Community_key,Run_Names,Community_type,Community)%>%dplyr::summarise(Drivers=paste(unique(sort(GeneNames)),collapse = ", "))
  Comm_Drivers<-data.frame(Comm_Drivers)%>%dplyr::mutate(Community = CommunityHyperLink(Community,Community_key,Community_type))
  Comm_Drivers<-Comm_Drivers[stringr::str_order(Comm_Drivers$Community, numeric = TRUE),]
  print("ali7")
  #==============================================================================================================

  #HGT to test for gene set enrichment
   # avoid showing datatable size-related warnings.

  print("ali8")

  #==============================================================================================================
  # add phenotype table
  print("ali11")
  if (!is.null(PhenotypeTablesList)){
    phenotype_table_all<-CreatePhenotypeTable(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC, PhenotypeTablesList)
    print(head(phenotype_table_all))
  }
  #==============================================================================================================
  # do hypergeometric test
  if (hyper_geo_test_bool == TRUE) {
    all_hgt_output<-CreateHyperGeoTestAll(cAMARETTOresults,cAMARETTOnetworkM,cAMARETTOnetworkC,hyper_geo_reference,MSIGDB,driverGSEA)
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
    
    print(ComNr)

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
    
    print("maman joonam")
    
    if (hyper_geo_test_bool) {
        outputHGT<-all_hgt_output%>%filter(Community_key==ComNr)%>%select(-Community_key,-Community_type,-Community)
        if (nrow(outputHGT)>0){
          DTGSEA_colnames<-c("Gene Set Name","# Genes in Gene Set","# Genes in Overlap","Genes in Overlap","% Genes in overlap","P-value","FDR Q-value")
          if (MSIGDB == TRUE) {
            outputHGT<-outputHGT %>% dplyr::mutate(Geneset = paste0("<a href=\"http://software.broadinstitute.org/gsea/msigdb/cards/", Geneset, "\">", gsub("_", " ", Geneset),"</a>"))
            DTGSEA_colnames<-R.utils::insert(DTGSEA_colnames,2,"Gene Set Description")
          }
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
  } else {
    DTGSEA <- "Genesets were not analysed as they were not provided."
  }
  print("hadi")
  # add phenotype table for each community page
  if (!is.null(PhenotypeTablesList)){
    phenotype_table_community<-phenotype_table_all%>%dplyr::filter(Community_key==ComNr)%>%select(Run_Names,everything())%>%dplyr::mutate(ModuleNr=ModuleHyperLink(ModuleNr,Run_Names,AMARETTOres,HTMLsAMARETTOlist,CopyAMARETTOReport,page=2))%>%select(-AMARETTOres)
    DTPhC<-DT::datatable(phenotype_table_community%>%dplyr::select(-Community_key,-Community,-Community_type),
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
      system.file("templates/TemplateCommunityPage.Rmd", package = "CommunityAMARETTO"),
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
  #==============================================================================================================
  # index page : 
  DTRunInfo<-datatable(Runs_AMARETTOs_info,class = "display",
            filter = 'top',
            extensions = c('Buttons','KeyTable'),
            escape=FALSE,
            rownames = FALSE,
            options=optionsList,
            colnames = c("AMARETTO Report"))

  DTComModulesLink<-datatable(ComModulesLink %>% 
              dplyr::rename(Edges="numTotalEdgesInCommunity","Fraction Edges"="fractEdgesInVsOut","Fraction Size"="CommsizeFrac"), 
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

  if (hyper_geo_test_bool) {
    all_hgt_output<-all_hgt_output%>%filter(n_Overlapping>2)%>%dplyr::mutate(Community=CommunityHyperLink(Community,Community_key,Community_type))%>%select(-Community_type,-Community_key)%>%select(Community,everything())
    DTGSEA_colnames<-c("Community","Gene Set Name","# Genes in Gene Set","# Genes in Overlap","Genes in Overlap","% Genes in overlap","P-value","FDR Q-value")
    if (MSIGDB == TRUE) {
      DTGSEA_colnames<-R.utils::insert(DTGSEA_colnames,3,"Gene Set Description")
    }
    DTGSEAall <- DT::datatable(all_hgt_output %>% dplyr::mutate(Geneset = paste0("<a href=\"http://software.broadinstitute.org/gsea/msigdb/cards/", Geneset, "\">", gsub("_", " ", Geneset),"</a>")),
                        class = "display",
                        filter = 'top',
                        extensions = c('Buttons','KeyTable'),
                        rownames = FALSE,
                        options = optionsList,
                        colnames = DTGSEA_colnames, escape = FALSE)%>%
                        DT::formatSignif(c("p_value", "padj","overlap_perc"), 2) %>% 
                        DT::formatStyle("overlap_perc", background = styleColorBar(c(0, 1), "lightblue"), backgroundSize = "98% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "center")%>%DT::formatStyle(columns = c(6), fontSize = '60%')
    }
  else{
    DTGSEAall <- "Genesets were not analysed as they were not provided."
  }
  DTComDrivers<- DT::datatable(Comm_Drivers%>%select(-Community_key,-Community_type)%>%select(Community,everything()),
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
    DTPh<-DT::datatable(phenotype_table_all%>%dplyr::mutate(Community = CommunityHyperLink(Community,Community_key,Community_type))%>%select(-Community_key,-Community_type)%>%select(Community,Run_Names,everything()),
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
      DTComDrivers=DTComDrivers,
      DTGeneComLink = DTGeneComLink,
      DTRunInfo =  DTRunInfo,
      DTGSEAall = DTGSEAall,
      DTPh = DTPh,
      DTdriverFreq = DTdriverFreq,
    ), quiet = TRUE)
}

#' @title HGTGeneEnrichmentList
#'
#' Calculates the p-values for unranked gene set enrichment based on two gmt files as input and the hyper geometric test.
#'
#' @param genelist The gmt file with reference gene set.
#' @param gmtfile The gmt file with gene sets to test. In our case, the gmt file of the modules.
#' @param NrCores Number of cores used for parallelization.
#' @param ref.numb.genes The total number of genes teste, standard equal to 45 956 (MSIGDB standard).
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
        k <- sum(gmtset[[i]] %in% genelist)
        m <- ref.numb.genes
        n <- length(genelist)
        p1 <- stats::phyper(k - 1, l, m - l, n, lower.tail = FALSE)
          
        if (k > 0) {
          overlapping.genes <- gmtset[[i]][gmtset[[i]] %in% genelist]
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

#' @title GeneSetDescription
#'
#' @param filename The name of the gmt file.
#' @param MSIGDB TRUE or FALSE
#' @importFrom utils data
#' @return provides descriptions of MSIGDB Genesets
#' @keywords internal
#' @examples GeneSetDescription(filename,MSIGDB)
GeneSetDescription<-function(filename,MSIGDB){
  utils::data(MsigdbMapping)
  gmtLines<-strsplit(readLines(filename),"\t")
  gmtLines_description <- lapply(gmtLines, function(x) {
    c(x[[1]],x[[2]],length(x)-2)
  })
  gmtLines_description<-data.frame(matrix(unlist(gmtLines_description),byrow=T,ncol=3),stringsAsFactors=FALSE)
  rownames(gmtLines_description)<-NULL
  colnames(gmtLines_description)<-c("GeneSet","Description","NumberGenes")
  gmtLines_description$NumberGenes<-as.numeric(gmtLines_description$NumberGenes)
  if(MSIGDB){
    gmtLines_description$Description<-sapply(gmtLines_description$GeneSet, function(x) {
      index<-which(MsigdbMapping$geneset==x)
      ifelse(length(index)!=0, MsigdbMapping$description[index],gmtLines_description$Description[which(gmtLines_description$GeneSet==x)]) 
    })}
  return(gmtLines_description)
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
#' @param MSIGDB 
#' @param driverGSEA 
#' @param NrCores 
#'
#' @return Geneset Enrichment Analysis for the entire communities. 
#' @export
#'
#' @examples CreateHyperGeoTestAll(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC, hyper_geo_reference = './h.all.v6.2.symbols.gmt', MSIGDB=TRUE, driverGSEA=TRUE,NrCores=4)
CreateHyperGeoTestAll<-function(cAMARETTOresults,cAMARETTOnetworkM,cAMARETTOnetworkC,hyper_geo_reference,MSIGDB,driverGSEA,NrCores){
  print("Performing Geneset Enrichment Analysis ...")
  com_gene_df<-suppressWarnings(ComRunModGenInfo(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC))
  communities_all<-unique(com_gene_df$Community_key)
  all_hgt_output<-NULL
  GeneSetDescriptions <- GeneSetDescription(hyper_geo_reference,MSIGDB)
  for (ComNr in communities_all){
    target_genes<-com_gene_df%>%dplyr::filter(Community_key==ComNr)%>%dplyr::filter(Type=="Target")%>%dplyr::arrange(GeneNames)%>%dplyr::pull(GeneNames)
    driver_genes<-com_gene_df%>%dplyr::filter(Community_key==ComNr)%>%dplyr::filter(Type=="Driver")%>%dplyr::arrange(GeneNames)%>%dplyr::pull(GeneNames)
    driver_genes_weights<-com_gene_df%>%dplyr::filter(Community_key==ComNr)%>%dplyr::filter(Type=="Driver")%>%dplyr::arrange(GeneNames)%>%dplyr::pull(Weights)
    if(driverGSEA){
      genelist<-unique(c(target_genes,driver_genes))
    }
    else{
      genelist<-unique(target_genes)
    }
    outputHGT <- HGTGeneEnrichmentList(genelist, hyper_geo_reference, NrCores = NrCores)
    if (nrow(outputHGT)>0){
      outputHGT <- dplyr::left_join(outputHGT,GeneSetDescriptions, by = c(Geneset = "GeneSet")) %>%
        dplyr::mutate(overlap_perc = n_Overlapping / NumberGenes) %>% dplyr::select(Geneset, Description,n_RefGeneset, n_Overlapping, Overlapping_genes, overlap_perc, p_value, padj) %>% dplyr::arrange(padj)
      
      all_hgt_output<-rbind(all_hgt_output, outputHGT%>%mutate(Community_key=ComNr)) 
    }
    else{
      all_hgt_output<-rbind(all_hgt_output, NULL) 
    }
  }
  all_hgt_output<-all_hgt_output%>% dplyr::select(Community_key,everything())%>%left_join(com_gene_df%>%select(Community_key,Community,Community_type)%>%distinct(),by="Community_key")
  print("Geneset Enrichment Analysis is done!")
  return(all_hgt_output)
}


#' Title InitialCheckInputs
#'
#' @param cAMARETTOresults 
#' @param output_address 
#' @param HTMLsAMARETTOlist 
#' @param CopyAMARETTOReport 
#' @param hyper_geo_test_bool 
#' @param hyper_geo_reference 
#'
#' @return RunInfo dataframe
#' @export
#'
#' @examples InitialCheckInputs(cAMARETTOresults,output_address="./",HTMLsAMARETTOlist,CopyAMARETTOReport=FALSE,hyper_geo_test_bool = TRUE,hyper_geo_reference){
InitialCheckInputs<-function(cAMARETTOresults,output_address,HTMLsAMARETTOlist,CopyAMARETTOReport,hyper_geo_test_bool,hyper_geo_reference){
  if (!dir.exists(output_address)) {
    stop("Output directory is not existing.")
  } else {
    dir.create(file.path(output_address,"htmls"), showWarnings = FALSE)
    full_path <- file.path(normalizePath(output_address),"htmls")
    print(paste0("The output directory is: ",full_path))
  }
  if (hyper_geo_test_bool == TRUE) {
    if (!file.exists(hyper_geo_reference)) {
      stop("GMT for hyper geometric test is not existing.")
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
#' 
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
  ComDrivers$all_genes
  for (ComNr in unique(ComDrivers$Community_key)){
    gene_table<-table(strsplit(gsub("","",ComDrivers%>%filter(Community_key==ComNr)%>%pull(all_genes)),","))
    freq_tbl<-data.frame(GeneNames=names(gene_table),Freq=paste0("# Data Sets = ",as.numeric(gene_table)),Community_key=ComNr)
    gene_freq_df<-rbind(gene_freq_df,freq_tbl)
  }
  gene_freq_df<-gene_freq_df%>%group_by(Community_key,Freq)%>%summarise(GeneNames=paste(GeneNames,collapse = ", "))
  gene_freq_df<-suppressMessages(reshape2::dcast(gene_freq_df, Community_key~Freq, fill=""))
  ComDrivers<-ComDrivers%>%left_join(gene_freq_df,by="Community_key")%>%select(-all_genes)
  return(ComDrivers)
}


ModuleHyperLink<-function(Module,Run_Names,AMARETTOres,HTMLsAMARETTOlist,CopyAMARETTOReport,page=1){
  # htmldir = "someLocalAdress/LIHC_Report_75/AMARETTOhtmls"
  # if CopyAMARETTOReport ==TRUE :
  #   htmldir = "./LIHC_Report_75/AMARETTOhtmls"
  if(is.null(HTMLsAMARETTOlist)==FALSE){
    RunInfo<-rownames_to_column(as.data.frame(HTMLsAMARETTOlist),"Run_Names") %>% dplyr::rename(ModuleLink="HTMLsAMARETTOlist")
    if(CopyAMARETTOReport==TRUE & page == 2){
      RunInfo <- RunInfo %>% mutate(ModuleLink=sub("^./","../",ModuleLink))
    }
    Mod_hyperLink<-data.frame(Run_Names=Run_Names,Module=Module,AMARETTOres=AMARETTOres,stringsAsFactors = FALSE)%>%dplyr::left_join(RunInfo,by="Run_Names")%>%mutate(mod_hyperlink=ifelse(AMARETTOres==1,paste0("<a href=",ModuleLink,"/modules/",gsub("Module_","module",Module),".html>",gsub("Module_","Module ",Module), "</a>"),Module))%>%dplyr::pull(mod_hyperlink)
  }
  else{
    Mod_hyperLink<-ifelse(AMARETTOres==1,gsub("Module_","Module ",Module),Module)
  }
  return(Mod_hyperLink)
}


RunHyperLink<-function(Run_Names,AMARETTOres,HTMLsAMARETTOlist,CopyAMARETTOReport,page=1){
  if(is.null(HTMLsAMARETTOlist)==FALSE){
    RunInfo<-rownames_to_column(as.data.frame(HTMLsAMARETTOlist),"Run_Names") %>% dplyr::rename(RunLink="HTMLsAMARETTOlist")
    if(CopyAMARETTOReport==TRUE & page == 2){
      RunInfo <- RunInfo %>% mutate(RunLink=sub("^./","../",RunLink))
    }
    Run_hyperLink<-data.frame(Run_Names=Run_Names,AMARETTOres=AMARETTOres,stringsAsFactors = FALSE)%>%dplyr::left_join(RunInfo,by="Run_Names")%>%mutate(run_hyperlink=ifelse(AMARETTOres==1,paste0("<a href=",RunLink,"/index.html>",Run_Names, "</a>"),Run_Names))%>%dplyr::pull(run_hyperlink)
  }
  else{
    Run_hyperLink<-Run_Names
  }
  return(Run_hyperLink)
}

#' @title cAMARETTO_HTMLreport
#'
#' @param cAMARETTOresults
#' @param cAMARETTOnetworkM
#' @param cAMARETTOnetworkC
#' @param output_address
#' @param HTMLsAMARETTOlist
#' @param hyper_geo_test_bool
#' @param hyper_geo_refence
#' @param MSIGDB
#' @param GMTURL
#' @param NrCores
#'
#' @return A set of HTMLs, giving caracteristics of the communities
#' 
#' @import igraph
#' @import DT
#' @import tidyverse
#' @import reshape2
#' @import rmarkdown
#' @return
#' @export
cAMARETTO_HTMLreport <- function(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC,
                                 output_address="./", HTMLsAMARETTOlist=NULL,
                                 hyper_geo_test_bool = FALSE, hyper_geo_reference = NULL, MSIGDB = FALSE,GMTURL = FALSE,
                                 NrCores=2){
  
  #test parameters
  if (!dir.exists(output_address)) {
    stop("Output directory is not existing.")
  }
  
  if (hyper_geo_test_bool == TRUE) {
    if (!file.exists(hyper_geo_reference)) {
      stop("GMT for hyper geometric test is not existing.")
    }
  }
  
  i=1
  if(is.null(HTMLsAMARETTOlist)==FALSE){
    for(htmldir in HTMLsAMARETTOlist){
      if(!file.exists(htmldir)){
        stop("An AMARETTO html directory is not existing.")
      }
      htmldir<-normalizePath(htmldir,"/report_html/htmls/modules/module")
      HTMLsAMARETTOlist[i] <- htmldir
      i=i+1
    }
  }
  
  #dataframe with modules per Run
  ComModulesLink <- stack(cAMARETTOnetworkC$community_list) %>% 
    dplyr::rename(Module="values", Community="ind")
  
  #adding Modules that are not in Communities or Communities that are filtered out
  all_module_names <- unique(c(cAMARETTOresults$hgt_modules$Geneset,cAMARETTOresults$hgt_modules$Testset))
  suppressWarnings(ComModulesLink <- left_join(as.data.frame(all_module_names) %>% dplyr::rename(Module="all_module_names"),ComModulesLink, by="Module"))
  all_module_names<-all_module_names[!all_module_names %in% ComModulesLink$Module]
  
  Nodes_Mnetwork <- igraph::as_data_frame(cAMARETTOnetworkM$module_network, what="vertices")
  Module_no_Network <- all_module_names[!all_module_names %in% Nodes_Mnetwork$name]
  Module_no_Com <- all_module_names[all_module_names %in% Nodes_Mnetwork$name]
  
  ComModulesLink <- ComModulesLink %>% 
    mutate(Community=ifelse(Module %in% Module_no_Network,"Not in Network",ifelse(Module %in% Module_no_Com, "Not in a Community",paste0("<a href=\"./communities/Community_",Community,".html\">Community ",Community, "</a>")))) %>%
    separate(Module, c("Run","ModuleNr"), "_Module_") %>% 
    mutate(ModuleNr = paste0("Module ", ModuleNr)) %>% 
    group_by(Community, Run) %>% 
    summarise(ModuleNrs=paste(ModuleNr, collapse = ", "))

  ComModulesLink <- dcast(ComModulesLink, Community~Run, fill=0)
  full_path <- normalizePath(output_address)
  
  comm_info <- cAMARETTO_InformationTable(cAMARETTOnetworkM, cAMARETTOnetworkC)
  
  GeneComLink <- comm_info %>% 
    select(community_numb,overlapping_genes) %>% 
    mutate(overlapping_genes = strsplit(as.character(overlapping_genes), ", "), community_numb=paste0("<a href=\"./communities/Community_",community_numb,".html\">Community ",community_numb, "</a>")) %>% 
    unnest(overlapping_genes) %>% 
    dplyr::rename(Community="community_numb",GeneName="overlapping_genes") %>% 
    arrange(GeneName)
  
  rmarkdown::render(
    system.file("templates/TemplateIndexPage.Rmd", package = "CommunityAMARETTO"),
    output_dir = paste0(full_path, "/htmls/"),
    output_file = "index.html",
    params = list(
      cAMARETTOnetworkM = cAMARETTOnetworkM,
      cAMARETTOnetworkC = cAMARETTOnetworkC,
      ComModulesLink = ComModulesLink,
      GeneComLink = GeneComLink,
    ), quiet = TRUE)
  
  print("The index html is created.")
  
  #HGT to test for gene set enrichment
  if (hyper_geo_test_bool) {
    GeneSetDescriptions <- GeneSetDescription(hyper_geo_reference)
  }
  
  for (i in 1:nrow(comm_info)){
    community_info <- comm_info[i,]
    ComNr <- community_info$community_numb
    genelist <- unlist(strsplit(community_info$overlapping_genes,", "))
    ModuleList <- unlist(strsplit(community_info$included_nodes,", "))
    ModuleNames <- ModuleList
    
    if(is.null(HTMLsAMARETTOlist)==FALSE){
      for (i in 1:length(HTMLsAMARETTOlist)){
        ModuleList <- sub(names(HTMLsAMARETTOlist)[i], HTMLsAMARETTOlist[i], ModuleList)
      }
      ModuleList <- sub("_Module_","/report_html/htmls/modules/module",ModuleList)
      ModuleList <- paste0('<a href="file://',ModuleList,'.html"> ModuleLink </a>')
      names(ModuleList) <- ModuleNames
      ModuleList <- rownames_to_column(as.data.frame(ModuleList),"ModuleNames") %>% separate(ModuleNames,c("Run","ModuleNr"),"_Module_") %>% mutate(ModuleNr = sub("^","Module ",ModuleNr))
      DTML <- datatable(ModuleList, 
                        class = "display",
                        extensions = "Buttons",
                        rownames = FALSE,
                        options = list(pageLength = 10, dom = "Bfrtip", buttons = list(list(extend = 'csv',text = "Save CSV", title=paste0("ModulesCom",ComNr)))),
                        escape = FALSE)
    } else {
      ModuleList <- as.data.frame(ModuleList) %>% separate(ModuleList,c("Run","ModuleNr"),"_Module_") %>% mutate(ModuleNr = sub("^","Module ",ModuleNr))
      DTML <- datatable(ModuleList, 
                        class = "display",
                        extensions = "Buttons",
                        rownames = FALSE,
                        options = list(pageLength = 10, dom = "Bfrtip", buttons = list(list(extend = 'csv',text = "Save CSV", title=paste0("ModulesCom",ComNr)))),
                        escape = FALSE)
    }
    
    if (hyper_geo_test_bool) {
      outputHGT <- HGTGeneEnrichmentList(genelist, hyper_geo_reference, NrCores = NrCores)
      outputHGT <- left_join(outputHGT,GeneSetDescriptions, by = c(Geneset = "GeneSet")) %>%
        mutate(overlap_perc = n_Overlapping / NumberGenes) %>% dplyr::select(Geneset, Description, n_Overlapping, Overlapping_genes, overlap_perc, p_value, padj) %>% arrange(padj)
      
      if (MSIGDB == TRUE & GMTURL == FALSE) {
        DTGSEA <- datatable(outputHGT %>% mutate(Geneset = paste0("<a href=\"http://software.broadinstitute.org/gsea/msigdb/cards/", Geneset, "\">", gsub("_", " ", Geneset),"</a>")),
                          class = "display",
                          extensions = "Buttons",
                          rownames = FALSE,
                          options = list(pageLength = 10, dom = "Bfrtip", buttons = list(list(extend = 'csv',text = "Save CSV", title=paste0("HGTresults_Com",ComNr)))),
                          colnames = c("Gene Set Name", "Description", "# Genes in Overlap",  "Overlapping Genes", "Percent of GeneSet overlapping", "p-value", "FDR q-value"), escape = FALSE) %>% 
                          formatSignif(c("p_value", "padj","overlap_perc"), 2) %>% 
                          formatStyle("overlap_perc", background = styleColorBar(c(0, 1), "lightblue"), backgroundSize = "98% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "center")
      } else if (MSIGDB == TRUE & GMTURL == TRUE) {
        DTGSEA <- datatable(outputHGT %>% dplyr::select(-Description) %>% mutate(Geneset = paste0("<a href=\"http://software.broadinstitute.org/gsea/msigdb/cards/", Geneset, "\">", gsub("_", " ", Geneset),"</a>")),
                          class = "display",
                          extensions = "Buttons",
                          rownames = FALSE,
                          options = list(pageLength = 10, dom = "Bfrtip", buttons = list(list(extend = 'csv',text = "Save CSV", title=paste0("HGTresults_Com",ComNr)))),
                          colnames = c("Gene Set Name", "# Genes in Overlap",  "Overlapping Genes", "Percent of GeneSet overlapping", "p-value", "FDR q-value"), escape = FALSE) %>% 
                          formatSignif(c("p_value", "padj","overlap_perc"), 2) %>% 
                          formatStyle("overlap_perc", background = styleColorBar(c(0, 1), "lightblue"), backgroundSize = "98% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "center")
      } else if (MSIGDB == FALSE & GMTURL == TRUE) {  
        DTGSEA <- datatable(outputHGT %>% mutate(Geneset = paste0("<a href=\"", Description, "\">", gsub("_", " ", Geneset),"</a>")) %>% dplyr::select(-Description),
                          class = "display",
                          extensions = "Buttons",
                          rownames = FALSE,
                          options = list(pageLength = 10, dom = "Bfrtip", buttons = list(list(extend = 'csv',text = "Save CSV", title=paste0("HGTresults_Com",ComNr)))),
                          colnames = c("Gene Set Name", "# Genes in Overlap",  "Overlapping Genes", "Percent of GeneSet overlapping", "p-value", "FDR q-value"), escape = FALSE) %>% 
                          formatSignif(c("p_value", "padj","overlap_perc"), 2) %>% 
                          formatStyle("overlap_perc", background = styleColorBar(c(0, 1), "lightblue"), backgroundSize = "98% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "center")
      } else {
        DTGSEA <- datatable(outputHGT,
                          class = "display",
                          extensions = "Buttons",
                          rownames = FALSE,
                          options = list(pageLength = 10, dom = "Bfrtip", buttons = list(list(extend = 'csv',text = "Save CSV", title=paste0("HGTresults_Com",ComNr)))),
                          colnames = c("Gene Set Name", "# Genes in Overlap",  "Overlapping Genes", "Percent of GeneSet overlapping", "p-value", "FDR q-value"), escape = FALSE) %>% 
                          formatSignif(c("p_value", "padj","overlap_perc"), 2) %>% 
                          formatStyle("overlap_perc", background = styleColorBar(c(0, 1), "lightblue"), backgroundSize = "98% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "center")
      }
  } else {
    DTGSEA <- "Genesets were not analysed as they were not provided."
  }
    rmarkdown::render(
      system.file("templates/TemplateCommunityPage.Rmd", package = "CommunityAMARETTO"),
      output_dir = paste0(full_path, "/htmls/communities"),
      output_file = paste0("Community_",ComNr,".html"),
      params = list(
        ComNr = ComNr,
        DTGSEA = DTGSEA,
        DTML = DTML,
        cAMARETTOnetworkM = cAMARETTOnetworkM,
        cAMARETTOnetworkC = cAMARETTOnetworkC,
      ), quiet = TRUE)
  }
  
  print("The community htmls are created.")
}

#' @title HGTGeneEnrichmentList
#'
#' Calculates the p-values for unranked gene set enrichment based on two gmt files as input and the hyper geometric test.
#'
#' @param genelist The gmt file with reference gene set.
#' @param gmtfile The gmt file with gene sets to test. In our case, the gmt file of the modules.
#' @param NrCores Number of cores used for parallelization.
#' @param ref.numb.genes The total number of genes teste, standard equal to 45 956 (MSIGDB standard).
#'
#' @import doParallel
#' @keywords internal
#' @export
HGTGeneEnrichmentList <- function(genelist, gmtfile, NrCores, ref.numb.genes = 45956) {
    gmtset <- readGMT(gmtfile)  # the hallmarks_and_co2...
    
    ########################### Parallelizing :
    cluster <- makeCluster(c(rep("localhost", NrCores)), type = "SOCK")
    registerDoParallel(cluster, cores = NrCores)
    
    resultHGT<-foreach(i = 1:length(gmtset), .combine = "rbind") %dopar% {
        # print(i)
        l <- length(gmtset[[i]])
        k <- sum(gmtset[[i]] %in% genelist)
        m <- ref.numb.genes
        n <- length(genelist)
        p1 <- phyper(k - 1, l, m - l, n, lower.tail = FALSE)
          
        if (k > 0) {
          overlapping.genes <- gmtset[[i]][gmtset[[i]] %in% genelist]
          overlapping.genes <- paste(overlapping.genes, collapse = ", ")
          c(Geneset = names(gmtset[i]),
            Testset = names(genelist),
            p_value = p1,
            n_Overlapping = k,
            Overlapping_genes = overlapping.genes)
        }
    }
    
    stopCluster(cluster)
    resultHGT <- as.data.frame(resultHGT, stringsAsFactors = FALSE)
    resultHGT$p_value <- as.numeric(resultHGT$p_value)
    resultHGT$n_Overlapping <- as.numeric((resultHGT$n_Overlapping))
    resultHGT[, "padj"] <- p.adjust(resultHGT[, "p_value"], method = "BH")
    return(resultHGT)
}

#' @title GeneSetDescription
#' @param filename The name of the gmt file.
#'
#' @return
#' @keywords internal
#' @examples
#' @export
GeneSetDescription <- function(filename) {
  gmtLines <- strsplit(readLines(filename), "\t")
  gmtLines_description <- lapply(gmtLines, function(x) {c(x[[1]], x[[2]], length(x) - 2)})
  gmtLines_description <- data.frame(matrix(unlist(gmtLines_description), byrow = T, ncol = 3), stringsAsFactors = FALSE)
  rownames(gmtLines_description) <- NULL
  colnames(gmtLines_description) <- c("GeneSet", "Description", "NumberGenes")
  gmtLines_description$NumberGenes <- as.numeric(gmtLines_description$NumberGenes)
  return(gmtLines_description)
}

  
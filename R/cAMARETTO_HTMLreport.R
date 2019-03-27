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
#'
#' @return A set of HTMLs, giving caracteristics of the communities
#' 
#' @import igraph
#' @import DT
#' @import tidyverse
#' @import reshape2
#' @import rmarkdown
#' @examples
#' 
#' cAMARETTO_HTMLreport(cAMARETTOresults,cAMARETTOnetworkM, cAMARETTOnetworkC,HTMLsAMARETTOlist = HTMLsAMARETTOlist, hyper_geo_test_bool = TRUE, hyper_geo_reference = gmtfile, MSIGDB = TRUE, output_address= "./")
#' 
#' @export
cAMARETTO_HTMLreport <- function(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC,
                                 output_address="./", HTMLsAMARETTOlist=NULL, CopyAMARETTOReport = TRUE,
                                 hyper_geo_test_bool = FALSE, hyper_geo_reference = NULL,
                                 MSIGDB = FALSE,driverGSEA=TRUE,
                                 NrCores=2){
  
  #test parameters
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
    RunInfo<-rownames_to_column(as.data.frame(HTMLsAMARETTOlist),"Run") %>% rename(ModuleLink="HTMLsAMARETTOlist")
  }
  
  # Extract main dataframes
  com_gene_df<-suppressWarnings(ComRunModGenInfo(cAMARETTOresults,cAMARETTOnetworkC))
  comm_info <-  suppressWarnings(cAMARETTO_InformationTable(cAMARETTOnetworkM, cAMARETTOnetworkC))
  
  #dataframe with modules per Run
  ComModulesLink <- stack(cAMARETTOnetworkC$community_list) %>% 
    dplyr::rename(Module="values", Community="ind")
  
  all_module_names <- unique(c(cAMARETTOresults$hgt_modules$Geneset1,cAMARETTOresults$hgt_modules$Geneset2))
  suppressWarnings(ComModulesLink <- left_join(as.data.frame(all_module_names) %>% dplyr::rename(Module="all_module_names"),ComModulesLink, by="Module"))
  #adding Module Links
  if(is.null(HTMLsAMARETTOlist)==FALSE){
    suppressMessages(ComModulesLink<- left_join(ComModulesLink %>% mutate(Run=sub("\\|.*","",Module)), RunInfo))
    ComModulesLink <- ComModulesLink %>% mutate(ModuleLink=ifelse(Run %in% cAMARETTOresults$runnames,paste0(ModuleLink,"/modules/module",sub(".*_","",Module),".html"),NA))
  }
  #adding Modules that are not in Communities or Communities that are filtered out
  all_module_names <- ComModulesLink[is.na(ComModulesLink[,"Community"]),"Module"]
  Nodes_Mnetwork <- igraph::as_data_frame(cAMARETTOnetworkM$module_network, what="vertices")
  Module_no_Network <- all_module_names[!all_module_names %in% Nodes_Mnetwork$name]
  Module_no_Com <- all_module_names[all_module_names %in% Nodes_Mnetwork$name]
  
  ComModulesLink <- ComModulesLink %>%
    mutate(Community=ifelse(Module %in% Module_no_Network,"Not in Network",ifelse(Module %in% Module_no_Com, "Not in a Community",paste0("Community ",Community)))) %>%
    mutate(ModuleName = sub("^.*\\|","",Module)) %>% 
    mutate(ModuleName = sub("_"," ",ModuleName)) 
  if(is.null(HTMLsAMARETTOlist)==FALSE){
    ComModulesLink <- ComModulesLink %>% mutate(ModuleName = ifelse(Run %in% cAMARETTOresults$runnames, paste0('<a href="',ModuleLink,'">',ModuleName,'</a>'),ModuleName))
    RunInfo <- RunInfo %>% mutate(Run=paste0('<a href="',ModuleLink,'/index.html">',Run,'</a>'))
  } else {
    ComModulesLink <- ComModulesLink %>% mutate(Run=sub("\\|.*","",Module))
    RunInfo <- as.data.frame(cAMARETTOresults$runnames) 
    names(RunInfo) <- c("Run")
  }
  ComModulesLink <- ComModulesLink %>% 
    group_by(Community, Run) %>% 
    summarise(ModuleNames=paste(ModuleName, collapse = ", "))
  suppressMessages(ComModulesLink <- dcast(ComModulesLink, Community~Run, fill=0))
  suppressMessages(ComModulesLink <- left_join(ComModulesLink,cAMARETTOnetworkC$commEdgeInfo %>% dplyr::select(Community,numTotalEdgesInCommunity,fractEdgesInVsOut,CommsizeFrac) %>% mutate(Community=paste0("Community ",Community))))
  ComModulesLink <- ComModulesLink%>%mutate(Community = paste0("<a href=\"./communities/",sub(" ","_",Community),".html\">",Community, "</a>"))

  # adding genes to communities table
  
  com_gene_df<-com_gene_df%>%mutate(Color=sapply(as.numeric(Weights), function(x){
    if(is.na(x)){
      return("")
    }
    else if(x>0){
      return("darkred")
    }
    else {
      return("darkblue")
    }
  }))%>%mutate(TypeColored=paste0('<font color=',Color,'>',Type,'</font>'))
  #%>%mutate(Community=stringr::str_sort(Community, numeric = TRUE))
  GeneComLink<-com_gene_df%>%rename(GeneName = GeneNames)%>%
    dplyr::mutate(GeneName = paste0("<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=", GeneName, "\">", GeneName, "</a>"))%>%
    select(c(GeneName,TypeColored,Community))%>%rename(Type=TypeColored)%>%
    mutate(Community=paste0("<a href=\"./communities/",paste0("Community_",Community),".html\">",paste0("Community ",Community), "</a>"))
  
  #adding Community to driver genes table 
  Comm_Drivers<-com_gene_df%>%filter(Type=="Driver")%>%
    mutate(GeneNames=paste0("<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=", GeneNames, "\">", GeneNames, "</a>"))%>%
    group_by(Community,Run_Names)%>%summarise(Drivers=paste(unique(GeneNames),collapse = ", "))
  print("Hi Ali")
  Comm_Drivers<-data.frame(Comm_Drivers)%>%mutate(Community=paste0("<a href=\"./communities/",paste0("Community_",Community),".html\">",paste0("Community ",Community), "</a>"))
  print("Hi Omid")
  
  #HGT to test for gene set enrichment
  options('DT.warn.size'=FALSE) # avoid showing datatable size-related warnings.
  if (hyper_geo_test_bool) {
    all_hgt_output <- tibble("Community"=character(),"Geneset"=character(),"Description"=character(),"n_Overlapping"=numeric(),"Overlapping_genes"=character(),"overlap_perc"=numeric(),"p_value="=numeric(),"padj"=numeric())
    GeneSetDescriptions <- GeneSetDescription(hyper_geo_reference,MSIGDB)
  }
  if(is.null(HTMLsAMARETTOlist)==FALSE){
    RunInfo2<-rownames_to_column(as.data.frame(HTMLsAMARETTOlist),"Run") %>% rename(ModuleLink="HTMLsAMARETTOlist")
    if(CopyAMARETTOReport==TRUE){
      RunInfo2 <- RunInfo2 %>% mutate(ModuleLink=sub("^./","../",ModuleLink))
    }
  }
  
  for (i in 1:nrow(comm_info)){
    community_info <- comm_info[i,]
    ComNr <- community_info$community_numb
    target_genes<-com_gene_df%>%filter(Community==i)%>%filter(Type=="Target")%>%arrange(GeneNames)%>%pull(GeneNames)
    driver_genes<-com_gene_df%>%filter(Community==i)%>%filter(Type=="Driver")%>%arrange(GeneNames)%>%pull(GeneNames)
    driver_genes_weights<-com_gene_df%>%filter(Community==i)%>%filter(Type=="Driver")%>%arrange(GeneNames)%>%pull(Weights)
    
    ModuleList <- unlist(strsplit(community_info$included_nodes,", "))
    if(is.null(HTMLsAMARETTOlist)==FALSE){
      ModuleList <- as.data.frame(ModuleList) %>% separate(ModuleList,c("Run","ModuleName"),"\\|",extra = "merge")
      if (CopyAMARETTOReport==TRUE){
        suppressMessages(ModuleList <- left_join(ModuleList,RunInfo2) %>% mutate(ModuleName = ifelse(Run %in% cAMARETTOresults$runnames,paste0('<a href="',ModuleLink,'/modules/',sub("Module_","module",ModuleName),'.html">',sub("_"," ",ModuleName),'</a>'),sub("_"," ",ModuleName))))
      } else {
        suppressMessages(ModuleList <- left_join(ModuleList,RunInfo2) %>% mutate(ModuleName = ifelse(Run %in% cAMARETTOresults$runnames,paste0('<a href="',ModuleLink,'/modules/',sub("Module_","module",ModuleName),'.html">',sub("_"," ",ModuleName),'</a>'),sub("_"," ",ModuleName))))
      }  
      DTML <- datatable(ModuleList %>% select(-ModuleLink), 
                        class = "display",
                        extensions = "Buttons",
                        rownames = FALSE,
                        options = list(pageLength = 10, dom = "Bfrtip", buttons = list(list(extend = 'csv',text = "Save CSV", title=paste0("ModulesCom",ComNr)))),
                        escape = FALSE)
    } else {
      ModuleList <- as.data.frame(ModuleList) %>% separate(ModuleList,c("Run","ModuleName"),"\\|",extra = "merge") %>% mutate(ModuleName = sub("_"," ",ModuleName))
      DTML <- datatable(ModuleList, 
                        class = "display",
                        extensions = "Buttons",
                        rownames = FALSE,
                        options = list(pageLength = 10, dom = "Bfrtip", buttons = list(list(extend = 'csv',text = "Save CSV", title=paste0("ModulesCom",ComNr)))),
                        escape = FALSE)
    }
    
    if (hyper_geo_test_bool) {
      genelist<-ifelse(driverGSEA,unique(c(target_genes,driver_genes)),unique(target_genes))
      outputHGT <- HGTGeneEnrichmentList(genelist, hyper_geo_reference, NrCores = NrCores)
      if (nrow(outputHGT)>0){
        outputHGT <- left_join(outputHGT,GeneSetDescriptions, by = c(Geneset = "GeneSet")) %>%
          mutate(overlap_perc = n_Overlapping / NumberGenes) %>% dplyr::select(Geneset, Description, n_Overlapping, Overlapping_genes, overlap_perc, p_value, padj) %>% arrange(padj)
        all_hgt_output<-rbind(all_hgt_output, outputHGT %>% mutate(Community=paste0("Community ",ComNr)) %>% dplyr::select(Community,everything()))
        
        if (MSIGDB == TRUE) {
          DTGSEA <- datatable(outputHGT %>% mutate(Geneset = paste0("<a href=\"http://software.broadinstitute.org/gsea/msigdb/cards/", Geneset, "\">", gsub("_", " ", Geneset),"</a>")),
                            class = "display",
                            extensions = "Buttons",
                            rownames = FALSE,
                            options = list(pageLength = 10, dom = "Bfrtip", buttons = list(list(extend = 'csv',text = "Save CSV", title=paste0("HGTresults_Com",ComNr)))),
                            colnames = c("Gene Set Name", "Description", "# Genes in Overlap",  "Overlapping Genes", "Percent of GeneSet overlapping", "p-value", "FDR q-value"), escape = FALSE) %>% 
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
      } else{
        DTGSEA <- "No significant overlap was identified in the geneset enrichment analysis."
        }
  } else {
    DTGSEA <- "Genesets were not analysed as they were not provided."
  }
  
  #adding Gene-Module-Run tabel
  genelists_module<-com_gene_df%>%filter(Community==i)%>%arrange(GeneNames)%>%rename(Run=Run_Names)%>%rename(ModuleName=ModuleNr)%>%rename(Genes=GeneNames)%>%select(-c(Type))%>%rename(Type=TypeColored)

  if(is.null(HTMLsAMARETTOlist)==FALSE){
    genelists_module <- suppressMessages(left_join(genelists_module,RunInfo2) %>% 
      dplyr::mutate(ModuleName=paste0("Module ",ModuleName))%>%
      dplyr::mutate(ModuleName = ifelse(Run %in% cAMARETTOresults$runnames,paste0('<a href="',ModuleLink,'/modules/',sub("Module ","module",ModuleName),'.html">',ModuleName,'</a>'),ModuleName)) %>%
      dplyr::mutate(Genes = paste0("<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=", 
                                   Genes, "\">", Genes, "</a>"))%>%
      dplyr::select(-ModuleLink)%>%dplyr::select(c(Run,ModuleName,Genes,Type)))
  }
  DTGenes <- datatable(genelists_module,
                       class = "display",
                       extensions = "Buttons",
                       rownames = FALSE,
                       options = list(pageLength = 10, 
                                      dom = "Bfrtip", 
                                      buttons = list(list(extend = 'csv',text = "Save CSV", title="GeneModuleLink"))),
                       escape=FALSE)
  
  DTComDrivers<-datatable(Comm_Drivers,
                          class = "display",
                          extensions = "Buttons",
                          rownames = FALSE,
                          options = list(pageLength = 10, 
                                         dom = "Bfrtip", 
                                         buttons = list(list(extend = 'csv',text = "Save CSV", title="GeneModuleLink"))),
                          escape=FALSE)
  knitr::knit_meta(class=NULL, clean = TRUE)  # cleaning memory, avoiding memory to be overloaded
  rmarkdown::render(
      system.file("templates/TemplateCommunityPage.Rmd", package = "CommunityAMARETTO"),
      output_dir = paste0(full_path, "/communities"),
      output_file = paste0("Community_",ComNr,".html"),
      params = list(
        ComNr = ComNr,
        DTGSEA = DTGSEA,
        DTML = DTML,
        DTGenes = DTGenes,
        DTComDrivers=DTComDrivers,
        cAMARETTOnetworkM = cAMARETTOnetworkM,
        cAMARETTOnetworkC = cAMARETTOnetworkC
      ), quiet = TRUE)
  }
  if (hyper_geo_test_bool) {
    if (MSIGDB == TRUE) {
      DTGSEAall <- datatable(all_hgt_output %>% mutate(Geneset = paste0("<a href=\"http://software.broadinstitute.org/gsea/msigdb/cards/", Geneset, "\">", gsub("_", " ", Geneset),"</a>")),
                        class = "display",
                        extensions = "Buttons",
                        rownames = FALSE,
                        options = list(pageLength = 10, dom = "Bfrtip", buttons = list(list(extend = 'csv',text = "Save CSV", title=paste0("HGTresults_Com",ComNr)))),
                        colnames = c("Community","Gene Set Name", "Description", "# Genes in Overlap",  "Overlapping Genes", "Percent of GeneSet overlapping", "p-value", "FDR q-value"), escape = FALSE) %>% 
                        formatSignif(c("p_value", "padj","overlap_perc"), 2) %>% 
                        formatStyle("overlap_perc", background = styleColorBar(c(0, 1), "lightblue"), backgroundSize = "98% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "center")
    } else {
      DTGSEAall <- datatable(all_hgt_output,
                        class = "display",
                        extensions = "Buttons",
                        rownames = FALSE,
                        options = list(pageLength = 10, dom = "Bfrtip", buttons = list(list(extend = 'csv',text = "Save CSV", title=paste0("HGTresults_Com",ComNr)))),
                        colnames = c("Community","Gene Set Name", "# Genes in Overlap",  "Overlapping Genes", "Percent of GeneSet overlapping", "p-value", "FDR q-value"), escape = FALSE) %>% 
      formatSignif(c("p_value", "padj","overlap_perc"), 2) %>% 
      formatStyle("overlap_perc", background = styleColorBar(c(0, 1), "lightblue"), backgroundSize = "98% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "center")
    }
  } else{
    DTGSEAall <- "Genesets were not analysed as they were not provided."
  }
  rmarkdown::render(
    system.file("templates/TemplateIndexPage.Rmd", package = "CommunityAMARETTO"),
    output_dir = full_path,
    output_file = "index.html",
    params = list(
      cAMARETTOnetworkM = cAMARETTOnetworkM,
      cAMARETTOnetworkC = cAMARETTOnetworkC,
      ComModulesLink = ComModulesLink,
      GeneComLink = GeneComLink,
      RunInfo =  RunInfo %>% select(Run),
      DTGSEAall = DTGSEAall
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
#'
#' @param filename The name of the gmt file.
#' @param MSIGDB TRUE or FALSE
#'
#' @return
#' @keywords internal
#' @examples
#' @export
GeneSetDescription<-function(filename,MSIGDB){
  data(MsigdbMapping)
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
#' @param cAMARETTOnetworkC 
#'
#' @return a dataframe contaning all communities, runname, and modules relationships. 
#' @export
#'
#' @examples ComRunModGenInfo(cAMARETTOresults,cAMARETTOnetworkC)
ComRunModGenInfo<-function(cAMARETTOresults,cAMARETTOnetworkC){
  ComModulesLink <- stack(cAMARETTOnetworkC$community_list) %>% dplyr::rename(Module="values", Community="ind")
  all_module_names <- unique(c(cAMARETTOresults$hgt_modules$Geneset1,cAMARETTOresults$hgt_modules$Geneset2))
  suppressWarnings(ComModulesLink <- left_join(as.data.frame(all_module_names) %>% dplyr::rename(Module="all_module_names"),ComModulesLink, by="Module"))
  ComModulesLink<-ComModulesLink%>%mutate(ModuleNr=as.numeric(gsub("Module_","",unlist(map(strsplit(Module,"\\|"),2)))))%>% mutate(Run_Names=unlist(map(strsplit(Module,"\\|"),1)))
  all_module_names <- ComModulesLink[is.na(ComModulesLink[,"Community"]),"Module"]
  Nodes_Mnetwork <- igraph::as_data_frame(cAMARETTOnetworkM$module_network, what="vertices")
  Module_no_Network <- all_module_names[!all_module_names %in% Nodes_Mnetwork$name]
  Module_no_Com <- all_module_names[all_module_names %in% Nodes_Mnetwork$name]
  suppressMessages(ComModulesLink <-dplyr::left_join(cAMARETTOresults$all_genes_modules_df,ComModulesLink))
  ComModulesLink <- ComModulesLink %>%arrange(as.numeric(Community))%>%mutate(Community=ifelse(Module %in% Module_no_Network,"Not in Network",ifelse(Module %in% Module_no_Com, "Not in a Community",paste0("",Community))))
  return(ComModulesLink)
}






  

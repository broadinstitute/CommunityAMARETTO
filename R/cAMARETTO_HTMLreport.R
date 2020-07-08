#' @title cAMARETTO_HTMLreport
#' Creates a HTMLreport for the community AMARETTO results
#'
#' @param cAMARETTOresults The output of the Results function.
#' @param cAMARETTOnetworkM The output of the Module Network function.
#' @param cAMARETTOnetworkC The output of the Identify Communities function.
#' @param output_address The output repository for the HTML report.
#' @param HTMLsAMARETTOlist A list with AMARETTO reports to link with the 
#' Community AMARETTO report. If NULL, no links are added.
#' @param CopyAMARETTOReport Boolean to indicate if the AMARETTO reports
#' needs to be copied in the AMARETTO report directory.
#' In this way links are contained when moving the HTML directory.
#' @param hyper_geo_reference A reference gmt file to perform 
#' the Hyper Geometric Test.
#' @param NrCores Number of Cores to use during generation of the HTML report.
#' @param driverGSEA if TRUE, driver genes beside the target genes will also 
#' be included for hypergeometric test. 
#' @param hyper_geo_reference_gp Hypergeometric test table 
#' for genetic perturbation
#' @param hyper_geo_reference_cp Hypergeometric test table for
#'  chemical perturbation
#' @param PhenotypeTablesList List of Phenotype Association Tables 
#' for different AMARETTO runs.
#'
#' @return A set of HTMLs, giving caracteristics of the communities
#' @importFrom igraph as_data_frame degree E graph_from_data_frame
#'  layout_with_fr V graph.data.frame norm_coords edge.betweenness.community
#' @import DT
#' @import rmarkdown
#' @import utils
#' @importFrom stringr str_order
#' @importFrom dplyr arrange group_by left_join mutate select summarise  
#' rename  filter everything pull distinct mutate one_of pull summarise
#' @importFrom tibble add_row tibble column_to_rownames rownames_to_column
#' @importFrom knitr knit_meta 
#' @importFrom reshape2 dcast
#' @importFrom utils stack data
#' @importFrom tidyr separate unite
#' @importFrom R.utils insert
#' @examples 
#' try(
#' cAMARETTO_HTMLreport(cAMARETTOresults,
#'   cAMARETTOnetworkM,
#'   cAMARETTOnetworkC,
#'   HTMLsAMARETTOlist = HTMLsAMARETTOlist,
#'   hyper_geo_reference = gmtfile,
#'   output_address= "./")
#' )
#' @export
cAMARETTO_HTMLreport <- function(cAMARETTOresults=list(),
                                cAMARETTOnetworkM=list(),
                                cAMARETTOnetworkC=list(),
                                PhenotypeTablesList = NULL,
                                output_address ="./",
                                HTMLsAMARETTOlist = NULL,
                                CopyAMARETTOReport = TRUE,
                                hyper_geo_reference = NULL,
                                hyper_geo_reference_gp = NULL,
                                hyper_geo_reference_cp = NULL,
                                driverGSEA = TRUE,
                                NrCores=2){

    ##################################### Bioconductor Considerations :
    Run_Names <- AMARETTOres <- Weights <- Color<- Type<- GeneNames<-NULL
    GeneName <- Community <- TypeColored <- Community_key<-NULL
    Community_type<- ModuleNr<-NULL
    ModuleName <- Run <- Genes <- q.value<-NULL
    #####################################

    RunInfoList<-InitialCheckInputs(cAMARETTOresults=cAMARETTOresults,
                                    output_address=output_address,
                                    HTMLsAMARETTOlist=HTMLsAMARETTOlist,
                                    CopyAMARETTOReport=CopyAMARETTOReport,
                                hyper_geo_reference=hyper_geo_reference,
                                hyper_geo_reference_gp=hyper_geo_reference_gp,
                                hyper_geo_reference_cp=hyper_geo_reference_cp)
    RunInfo<-RunInfoList$RunInfo
    RunInfo2<-RunInfoList$RunInfo2
    full_path<-RunInfoList$full_path
    HTMLsAMARETTOlist<-RunInfoList$HTMLsAMARETTOlist
    #====================================================================
    # Extract main dataframes
    com_gene_df<-suppressWarnings(
        ComRunModGenInfo(cAMARETTOresults = cAMARETTOresults,
                        cAMARETTOnetworkM = cAMARETTOnetworkM,
                        cAMARETTOnetworkC = cAMARETTOnetworkC))
    comm_info <-suppressWarnings(
        cAMARETTO_InformationTable(cAMARETTOnetworkM = cAMARETTOnetworkM,
                        cAMARETTOnetworkC = cAMARETTOnetworkC))
    #====================================================================
    Runs_AMARETTOs_info<-com_gene_df%>%
        dplyr::select(Run_Names,AMARETTOres)%>%
        dplyr::distinct()%>%
        dplyr::mutate(Run_Names = RunHyperLink(Run_Names,AMARETTOres,
                    HTMLsAMARETTOlist,CopyAMARETTOReport))%>%
                        filter(AMARETTOres==1)%>%select(Run_Names)
    #====================================================================
    ComModulesLink<-CommunityModuleTableCreate(cAMARETTOresults = 
                                                cAMARETTOresults,
                                cAMARETTOnetworkM=cAMARETTOnetworkM,
                                cAMARETTOnetworkC=cAMARETTOnetworkC,
                                HTMLsAMARETTOlist=HTMLsAMARETTOlist,
                                CopyAMARETTOReport=CopyAMARETTOReport)
    #====================================================================
    com_gene_df<-com_gene_df%>%
        dplyr::mutate(Color=sapply(as.numeric(Weights), function(x){
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
    }))%>%dplyr::mutate(TypeColored=paste0('<font color=',
                                            Color,
                                            '>',
                                            Type,
                                            '</font>'))
    #========================================================
    #[stringr::str_order(Community, numeric = TRUE),]
    geneCardURL<-"<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene="
    GeneComLink<-com_gene_df%>%
        filter(AMARETTOres==1)%>%
        dplyr::rename(GeneName = GeneNames)%>%
        dplyr::mutate(GeneName = paste0(geneCardURL,
                                        GeneName,
                                        "\">",
                                        GeneName,
                                        "</a>"))%>%
        dplyr::select(c(GeneName,Community,TypeColored,
                        Community_key,Community_type))%>%
        dplyr::rename(Type=TypeColored)%>%
        dplyr::mutate(Community = CommunityHyperLink(Community,Community_key,
                                                        Community_type))%>%
        dplyr::arrange(GeneName)%>%
        select(-Community_key,-Community_type)%>%
        select(GeneName,everything())%>%
        distinct()
    #=======================================================================
    #GeneComLink<-GeneComLink[stringr::str_order(GeneComLink$Community,
    #numeric = TRUE),]
    #=======================================================================
    #adding Community to driver genes table 
    Comm_Drivers<-com_gene_df%>%
        dplyr::filter(Type=="Driver")%>%
        dplyr::mutate(GeneNames=paste0(geneCardURL,
                                        GeneNames,
                                        "\">",
                                        GeneNames,
                                        "</a>"))%>%
        dplyr::group_by(Community_key,Run_Names,Community_type,Community)%>%
        dplyr::summarise(Drivers=paste(unique(sort(GeneNames)),collapse = ", "))

    Comm_Drivers<-data.frame(Comm_Drivers)%>%
        dplyr::mutate(Community = CommunityHyperLink(Community,
                                                        Community_key,
                                                        Community_type))

    Comm_Drivers<-Comm_Drivers[stringr::str_order(Comm_Drivers$Community,
                                numeric = TRUE),]
    #===================================================================
    #HGT to test for gene set enrichment
    # avoid showing datatable size-related warnings.
    #===================================================================
    # add phenotype table
    if (!is.null(PhenotypeTablesList)){
        phenotype_table_all<-CreatePhenotypeTable(cAMARETTOresults,
                                                    cAMARETTOnetworkM,
                                                    cAMARETTOnetworkC,
                                                    PhenotypeTablesList)
    }
    #===================================================================
    # do hypergeometric test
    if(!is.null(hyper_geo_reference)){
        if (is.character(hyper_geo_reference)){
            all_hgt_output<-CreateHyperGeoTestAll(cAMARETTOresults,
                                                    cAMARETTOnetworkM,
                                                    cAMARETTOnetworkC,
                                                    hyper_geo_reference,
                                                    driverGSEA)
                                                }
        else if (is.data.frame(hyper_geo_reference)){
            all_hgt_output<-hyper_geo_reference
        }
    else{
        stop("hyper_geo_reference is not in the correct format!")
        }
    }
    #=====================================================================
    #=====================================================================
    # do hypergeometric test
    if(!is.null(hyper_geo_reference_gp)){
        if (is.character(hyper_geo_reference_gp)){
            all_hgt_output_gp<-CreateHyperGeoTestAll(cAMARETTOresults,
                                                    cAMARETTOnetworkM,
                                                    cAMARETTOnetworkC,
                                                    hyper_geo_reference_gp,
                                                    driverGSEA)
        }
        else if (is.data.frame(hyper_geo_reference_gp)){
            all_hgt_output_gp<-hyper_geo_reference_gp
        }
    else{
        stop("hyper_geo_reference_gp is not in the correct format!")
    }
    }
    #=====================================================================
    #=====================================================================
    # do hypergeometric test
    if(!is.null(hyper_geo_reference_cp)){
        if (is.character(hyper_geo_reference_cp)){
            all_hgt_output_cp<-CreateHyperGeoTestAll(cAMARETTOresults,
                                                    cAMARETTOnetworkM,
                                                    cAMARETTOnetworkC,
                                                    hyper_geo_reference_cp,
                                                    driverGSEA)
        }
        else if (is.data.frame(hyper_geo_reference_cp)){
        all_hgt_output_cp<-hyper_geo_reference_cp
        }
        else{
            stop("hyper_geo_reference is not in the correct format!")
        }
    }
    #===================================================================

    options('DT.warn.size'=FALSE)
    buttons_list = list(list(extend ='csv'),
                        list(extend ='excel'),
                        list(extend = 'pdf', pageSize = 'A4',
                        orientation = 'landscape'),
                        list(extend ='print'), list(extend ='colvis'))
    columnDefs = list(list(className = 'dt-head-center', targets = "_all"),
                    list(className = 'text-left',targets = "_all"))
    optionsList = list(deferRender=TRUE,
                        pageLength = 10,
                        lengthMenu = c(5, 10, 20, 50, 100),
                        keys = TRUE,
                        dom = "Blfrtip", 
                        buttons = buttons_list,
                        columnDefs = columnDefs,
                        paging = TRUE)
    #=========================================================================
    # Community Pages : 
    for (ComNr in unique(com_gene_df$Community_key)){
        ModuleList<-com_gene_df%>%
            filter(Community_key==ComNr)%>%
            select(Run_Names,ModuleNr,AMARETTOres)%>%
            distinct()%>%
            mutate(ModuleNr=ModuleHyperLink(ModuleNr,Run_Names,AMARETTOres,
                                    HTMLsAMARETTOlist,
                                    CopyAMARETTOReport,page=2))%>%
                                    select(-AMARETTOres)
    DTML <- DT::datatable(ModuleList, 
                            class = "display",
                            filter = 'top',
                            extensions = c('Buttons','KeyTable'),
                            rownames = FALSE,
                            options = optionsList,
                            colnames = c("Data Set", "Module"),
                            escape = FALSE)

    #adding Gene-Module-Run tabel
    genelists_module<-com_gene_df%>%
        filter(AMARETTOres==1)%>%
        dplyr::filter(Community_key==ComNr)%>%
        dplyr::arrange(GeneNames)%>%
        dplyr::rename(Run=Run_Names)%>%
        dplyr::rename(ModuleName=ModuleNr)%>%
        dplyr::rename(Genes=GeneNames)%>%
        dplyr::select(-c(Type))%>%
        dplyr::rename(Type=TypeColored)

    geneCardURL<-"<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene="
    genelists_module <- genelists_module%>%
        mutate(ModuleName = ModuleHyperLink(ModuleName,
                                            Run,
                                            AMARETTOres,
                                            HTMLsAMARETTOlist,
                                            CopyAMARETTOReport,
                                            page=2))%>%
                dplyr::mutate(Genes = paste0(geneCardURL,
                                            Genes,
                                            "\">",
                                            Genes,
                                            "</a>"))%>%
                dplyr::select(c(Run,ModuleName,Genes,Type))

    DTGenes <- DT::datatable(genelists_module,
                            class = "display",
                            filter = 'top',
                            extensions = c('Buttons','KeyTable'),
                            rownames = FALSE,
                            options = optionsList,
                            colnames = c("Data Set", "Module",
                                        "Gene", "Gene Type"),
                            escape=FALSE)
    
    if(!is.null(hyper_geo_reference)) {
        DTGSEA<-create_hgt_datatable(all_hgt_output,
                                    com_table=TRUE,
                                    ComNr = ComNr)
    } else {
        DTGSEA <- "Genesets were not analysed as they were not provided."
    }
    if(!is.null(hyper_geo_reference_gp)) {
        DTGSEA_gp<-create_hgt_datatable(all_hgt_output_gp,
                                        com_table=TRUE,
                                        ComNr = ComNr)
    } else {
        DTGSEA_gp <- "Genesets were not analysed as they were not provided."
    }
    if(!is.null(hyper_geo_reference_cp)) {
        DTGSEA_cp<-create_hgt_datatable(all_hgt_output_cp,
                                        com_table=TRUE,
                                        ComNr = ComNr)
    } else {
        DTGSEA_cp<- "Genesets were not analysed as they were not provided."
    }
    # add phenotype table for each community page
    if (!is.null(PhenotypeTablesList)){
        phenotype_table_community<-phenotype_table_all%>%
        dplyr::filter(Community_key==ComNr)%>%
        select(Run_Names,everything())%>%
        dplyr::mutate(ModuleNr=ModuleHyperLink(ModuleNr,
                                                Run_Names,
                                                AMARETTOres,
                                                HTMLsAMARETTOlist,
                                                CopyAMARETTOReport,
                                                page=2))%>%
        select(-AMARETTOres)

        phenotype_table_community<-phenotype_table_community%>%
            dplyr::select(-Community_key,-Community,-Community_type)%>%
            arrange(q.value)

        DTPhC<-DT::datatable(phenotype_table_community,
                            class = "display",
                            filter = 'top',
                            extensions = c('Buttons','KeyTable'),
                            rownames = FALSE,
                            colnames=c("Data Set","Module",
                                    "Phenotype","Statistics Test",
                                    "P-value","FDR Q-value",
                                    "Descriptive Statistics"),
                            options = optionsList,
                            escape = FALSE)
    }
    else{ DTPhC = "Phenotype Statistical Analysis is not provided" }
    
    comm_name<-com_gene_df%>%
        filter(Community_key==ComNr)%>%
        select(Community,Community_key)%>%
        distinct()%>%
        pull(Community)
    if (grepl("Not in", comm_name)){
        ComTitle = comm_name
    }
    else{
        ComTitle = paste0("Community ",comm_name)
    }
    # cleaning memory, avoiding memory to be overloaded
    knitr::knit_meta(class=NULL, clean = TRUE)  
    rmarkdown::render(
    system.file("templates/community_page_template/TemplateCommunityPage.Rmd",
        package = "CommunityAMARETTO"),
        output_dir = paste0(full_path, "/communities"),
        output_file = paste0("Community_",ComNr,".html"),
        params = list(
            ComNr = ComNr,
            ComTitle = ComTitle,
            DTGSEA = DTGSEA,
            DTGSEA_gp = DTGSEA_gp,
            DTGSEA_cp = DTGSEA_cp,
            DTML = DTML,
            DTGenes = DTGenes,
            DTPhC = DTPhC,
            cAMARETTOnetworkM = cAMARETTOnetworkM,
            cAMARETTOnetworkC = cAMARETTOnetworkC), quiet = TRUE)
    }
    file_remove<-suppressWarnings(
        suppressMessages(file.remove(paste0(full_path,
            "/communities/Community_",
            c(seq_len(length(unique(com_gene_df$Community_key)))),"_files"))))
    file_remove<-suppressWarnings(
        suppressMessages(file.remove(paste0(full_path,"index_files"))))
    #========================================================================
    # index page : 
    DTRunInfo<-datatable(Runs_AMARETTOs_info,
            class = "display",
            filter = 'top',
            extensions = c('Buttons','KeyTable'),
            escape=FALSE,
            rownames = FALSE,
            options=optionsList,
            colnames = c("AMARETTO Report"))
    
    ComModulesLink<-ComModulesLink %>%
        dplyr::rename(Edges="numTotalEdgesInCommunity",
                            "Fraction Edges"="fractEdgesInVsOut",
                            "Fraction Size"="CommsizeFrac")
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

    #=========================================
    if (!is.null(hyper_geo_reference)) {
        DTGSEAall <-create_hgt_datatable(all_hgt_output, com_table=FALSE)
    }
    else{
        DTGSEAall <- "Genesets were not analysed as they were not provided."
    }

    if (!is.null(hyper_geo_reference_gp)) {
        DTGSEAall_gp <-create_hgt_datatable(all_hgt_output_gp, com_table=FALSE)
    }
    else{
        DTGSEAall_gp <- "Genesets were not analysed as they were not provided."
    }

    if (!is.null(hyper_geo_reference_cp)) {
        DTGSEAall_cp <-create_hgt_datatable(all_hgt_output_cp, com_table=FALSE)
    }
    else{
        DTGSEAall_cp <- "Genesets were not analysed as they were not provided."
    }
    #=========================================
    Comm_Drivers<-Comm_Drivers%>%
        select(-Community_key,-Community_type)%>%
        select(Community,everything())
    DTComDrivers<- DT::datatable(Comm_Drivers,
                                class = "display",
                                filter = 'top',
                                extensions = c('Buttons','KeyTable'),
                                rownames = FALSE,
                                options = optionsList,
                                colnames = c("Community",
                                            "Data Set",
                                            "Driver Genes"),
                                escape = FALSE)

    # add phenotype table for index page
    if (!is.null(PhenotypeTablesList)){
        phenotype_table_all<-phenotype_table_all%>%
            dplyr::mutate(ModuleNr=ModuleHyperLink(ModuleNr,
                                                    Run_Names,
                                                    AMARETTOres,
                                                    HTMLsAMARETTOlist,
                                                    CopyAMARETTOReport))%>%
            select(-AMARETTOres)
        phenotype_table_all<-phenotype_table_all%>%
            dplyr::mutate(Community = CommunityHyperLink(Community,
                                                        Community_key,
                                                        Community_type))%>%
            select(-Community_key,-Community_type)%>%
            select(Community,Run_Names,everything())%>%
            arrange(q.value,Community)

    DTPh<-DT::datatable(phenotype_table_all,
                        class = "display",
                        filter = 'top',
                        extensions = c('Buttons','KeyTable'),
                        rownames = FALSE,
                        options = optionsList,
                        colnames=c("Community","Data Set","Module",
                                    "Phenotype","Statistics Test",
                                    "P-value","FDR Q-value",
                                    "Descriptive Statistics"),
                        escape = FALSE)
    }
    else{ DTPh = "Phenotype Statistical Analysis is not provided" }
    
    driversFreqTbl<-DriversSharedTbl(cAMARETTOresults,
                                    cAMARETTOnetworkM,
                                    cAMARETTOnetworkC)%>%
        left_join(com_gene_df%>%
                    select(Community_key,Community,Community_type)%>%
                    distinct(),by="Community_key")%>%
                    dplyr::mutate(Community = CommunityHyperLink(Community,
                                                Community_key,
                                                Community_type))%>%
        select(-Community_key,-Community_type)%>%
        select(Community,everything())

    DTdriverFreq<- DT::datatable(driversFreqTbl,
                                class = "display",
                                filter = 'top',
                                extensions = c('Buttons','KeyTable'),
                                rownames = FALSE,
                                options = optionsList,
                                escape = FALSE)

    rmarkdown::render(
        system.file("templates/TemplateIndexPage.Rmd",
        package = "CommunityAMARETTO"),
        output_dir = full_path,
        output_file = "index.html",
        params = list(
            cAMARETTOnetworkM = cAMARETTOnetworkM,
            cAMARETTOnetworkC = cAMARETTOnetworkC,
            DTComModulesLink = DTComModulesLink,
            DTRunInfo =  DTRunInfo),
        quiet = TRUE)

    rmarkdown::render(
        system.file("templates/TemplateIndexPage_RunsInfo.Rmd",
                    package = "CommunityAMARETTO"),
        output_dir = full_path,
        output_file = "index_RunsInfo.html",
        params = list(DTRunInfo =  DTRunInfo),
        quiet = TRUE)

    rmarkdown::render(
        system.file("templates/TemplateIndexPage_AllCommunities.Rmd",
                    package = "CommunityAMARETTO"),
        output_dir = full_path,
        output_file = "index_AllCommunities.html",
        params = list(DTRunInfo =  DTRunInfo), quiet = TRUE)

    rmarkdown::render(
        system.file("templates/TemplateIndexPage_Drivers.Rmd",
                    package = "CommunityAMARETTO"),
        output_dir = full_path,
        output_file = "index_Drivers.html",
        params = list(DTdriverFreq = DTdriverFreq),quiet = TRUE)

    rmarkdown::render(
        system.file("templates/TemplateIndexPage_Drivers2.Rmd",
                    package = "CommunityAMARETTO"),
        output_dir = full_path,
        output_file = "index_Drivers2.html",
        params = list(DTComDrivers=DTComDrivers), quiet = TRUE)
    
    rmarkdown::render(
    system.file("templates/TemplateIndexPage_AllGenes.Rmd",
                package = "CommunityAMARETTO"),
    output_dir = full_path,
    output_file = "index_AllGenes.html",
    params = list(DTGeneComLink = DTGeneComLink), quiet = TRUE)

    rmarkdown::render(
        system.file("templates/TemplateIndexPage_GenesetsEnrichment.Rmd",
                package = "CommunityAMARETTO"),
        output_dir = full_path,
        output_file = "index_GenesetsEnrichment.html",
        params = list(DTGSEAall = DTGSEAall), quiet = TRUE)

    rmarkdown::render(
        system.file("templates/TemplateIndexPage_GenesetsEnrichment_gp.Rmd",
                    package = "CommunityAMARETTO"),
        output_dir = full_path,
        output_file = "index_GenesetsEnrichment_gp.html",
        params = list(DTGSEAall_gp = DTGSEAall_gp), quiet = TRUE)

    rmarkdown::render(
        system.file("templates/TemplateIndexPage_GenesetsEnrichment_cp.Rmd",
        package = "CommunityAMARETTO"),
        output_dir = full_path,
        output_file = "index_GenesetsEnrichment_cp.html",
        params = list(DTGSEAall_cp = DTGSEAall_cp), quiet = TRUE)

    rmarkdown::render(
        system.file("templates/TemplateIndexPage_PhenoAssociation.Rmd",
        package = "CommunityAMARETTO"),
        output_dir = full_path,
        output_file = "index_PhenoAssociation.html",
        params = list(DTPh = DTPh), quiet = TRUE)
}

#' @title HGTGeneEnrichmentList
#'
#' Calculates the p-values for unranked gene set enrichment
#'  based on two gmt files as input and the hyper geometric test.
#'
#' @param genelist The gmt file with reference gene set.
#' @param gmtfile The gmt file with gene sets to test.
#'  In our case, the gmt file of the modules.
#' @param NrCores Number of cores used for parallelization.
#' @param ref.numb.genes The total number of genes teste,
#'  standard equal to 45956 (MSIGDB standard).
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom stats p.adjust phyper
#' @keywords internal
#' @return a list a dataframe for Hypergeometric test
#' @export
#'
#' @examples 
#' try(HGTGeneEnrichmentList(genelist,
#'   gmtfile,
#'   NrCores,
#'   ref.numb.genes = 45956))
#' 

HGTGeneEnrichmentList <- function(genelist, gmtfile, NrCores,
                                    ref.numb.genes = 45956) {

    i<-j<-NULL
    `%dopar%` <- foreach::`%dopar%`
    `%do%` <- foreach::`%do%`
    gmtset <- readGMT(gmtfile)  # the hallmarks_and_co2...
    ########################### Parallelizing :
    cluster <- parallel::makeCluster(c(rep("localhost", NrCores)),
                                    type = "SOCK")
    doParallel::registerDoParallel(cluster, cores = NrCores)
    resultHGT<-foreach::foreach(i = seq_len(length(gmtset)),
        .combine = "rbind") %dopar% {
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
    resultHGT[, "padj"] <- stats::p.adjust(resultHGT[, "p_value"],
                                            method = "BH")
    return(resultHGT)
}

#' @title ComRunModGenInfo
#'
#' @param cAMARETTOresults cAMARETTOresults object 
#' @param cAMARETTOnetworkM cAMARETTOnetworkM object
#' @param cAMARETTOnetworkC  cAMARETTOnetworkC object
#'
#' @importFrom igraph as_data_frame degree E graph_from_data_frame 
#' layout_with_fr V graph.data.frame norm_coords edge.betweenness.community
#' @importFrom dplyr arrange rename left_join mutate
#' @importFrom utils stack data
#' @importFrom purrr map
#' @return a dataframe contaning all communities, runname, and 
#' modules relationships. 
#' @export
#'
#' @examples 
#' try(
#' df<-ComRunModGenInfo(cAMARETTOresults,
#'  cAMARETTOnetworkM,
#'  cAMARETTOnetworkC)
#' )
#' 
ComRunModGenInfo<-function(cAMARETTOresults,
                            cAMARETTOnetworkM,
                            cAMARETTOnetworkC){

    ##################################### Bioconductor Considerations :
    Run_Names <- ModuleNr <- Module <- Community<-NULL
    #####################################
    
    ComModulesLink <- utils::stack(cAMARETTOnetworkC$community_list) %>%
        dplyr::rename(Module="values", Community="ind")
    #all_module_names <- unique(c(cAMARETTOresults$hgt_modules$Geneset1,
    #cAMARETTOresults$hgt_modules$Geneset2))
    all_module_names<-unique(cAMARETTOresults$all_genes_modules_df%>%
                        mutate(Module=paste0(Run_Names,"|",ModuleNr))%>%
                        pull(Module))
    modules_with_community<-ComModulesLink$Module
    suppressWarnings(
        ComModulesLink <- dplyr::left_join(as.data.frame(all_module_names) %>%
            dplyr::rename(Module="all_module_names"),ComModulesLink,
            by="Module"))
    ComModulesLink<-ComModulesLink%>%
        dplyr::mutate(ModuleNr=unlist(purrr::map(strsplit(Module,"\\|"),2)))%>%
        dplyr::mutate(Run_Names=unlist(purrr::map(strsplit(Module,"\\|"),1)))
    # so far we have all the modules with their community information,
    #NA for modules with no community.
    Nodes_Mnetwork <- igraph::as_data_frame(cAMARETTOnetworkM$module_network,
                            what="vertices")
    # Nodes_Mnetwork = all modules selected to be in the network
    Module_no_Network <-setdiff(all_module_names,Nodes_Mnetwork$name)
    Module_no_Com<-setdiff(Nodes_Mnetwork$name,modules_with_community)
    suppressMessages(
    ComModulesLink <-dplyr::left_join(cAMARETTOresults$all_genes_modules_df,
                                        ComModulesLink,
                                        by=c("Run_Names","ModuleNr")))
    ComModulesLink <- ComModulesLink %>%
        dplyr::arrange(as.numeric(Community))%>%
        dplyr::mutate(Community=ifelse(Module %in% Module_no_Network,
                                paste0("Not in Network ",Run_Names),
                                    ifelse(Module %in% Module_no_Com,
                                            paste0("Not in Community ",
                                                Run_Names),
                                    Community)))
    ComModulesLink <- suppressWarnings(ComModulesLink %>%
        dplyr::left_join(data.frame(Community=unique(ComModulesLink$Community),
            Community_key=seq_len(length(unique(ComModulesLink$Community)))),
            by="Community"))
    ComModulesLink<-ComModulesLink%>%
        mutate(Community_type=sapply(Community, function(x){
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

#' @title CreatePhenotypeTable
#'
#' @param cAMARETTOresults list of cAMARETTOresults results
#' @param cAMARETTOnetworkM cAMARETTOnetworkM object
#' @param PhenotypeTablesList a list of phenotypeassociation tables
#'  computed in AMARETTO aalysis
#' @param cAMARETTOnetworkC cAMARETTOnetworkC object
#'
#' @importFrom dplyr select distinct mutate mutate left_join select
#' @return results
#' @export
#'
#' @examples try(CreatePhenotypeTable(cAMARETTOresults,
#'    cAMARETTOnetworkM,
#'    cAMARETTOnetworkC,
#'    PhenotypeTables))
CreatePhenotypeTable<-function(cAMARETTOresults,
                                cAMARETTOnetworkM,
                                cAMARETTOnetworkC,
                                PhenotypeTablesList){
    Run_Names<-ModuleNr<-Community<-Community_key<-Community_type<-NULL
    AMARETTOres<-Phenotypes<-Statistical_Test<-NULL
    p.value<-q.value<-Descriptive_Statistics<-NULL
    phenotype_table_all<-NULL
    CommunityRunModuleTable<-ComRunModGenInfo(cAMARETTOresults,
                                                cAMARETTOnetworkM,
                                                cAMARETTOnetworkC)%>%
        dplyr::select(Run_Names,ModuleNr,Community,Community_key,
                    Community_type,AMARETTOres)%>%dplyr::distinct()
    for (i in seq_len(length(PhenotypeTablesList))){
        if (is.null(PhenotypeTablesList[[i]])){
            next
        }
    phenotype_table<-PhenotypeTablesList[[i]]%>%
        dplyr::mutate(ModuleNr=gsub("Module ","Module_",ModuleNr))%>%
            dplyr::mutate(Run_Names=names(PhenotypeTablesList)[i])%>%
            dplyr::left_join(CommunityRunModuleTable,
                    by=c("Run_Names","ModuleNr"))%>%
            dplyr::select(Community,Community_key,Community_type,
                        AMARETTOres,Run_Names,ModuleNr,
                        Phenotypes,Statistical_Test,p.value,
                        q.value,Descriptive_Statistics)
    phenotype_table_all<-rbind(phenotype_table_all,phenotype_table)
    }
    return(phenotype_table_all)
}

#' @title CreateHyperGeoTestAll
#'
#' @param cAMARETTOresults cAMARETTOresults object
#' @param cAMARETTOnetworkM cAMARETTOnetworkM object
#' @param cAMARETTOnetworkC cAMARETTOnetworkC  object
#' @param hyper_geo_reference a vector address for the gmt files.
#' @param driverGSEA TRUE for inclusion of driver genes
#'  in hypergeomertic test
#' @param NrCores Nr of core for parallel processing
#'
#' @importFrom foreach foreach %dopar% %do%
#' @return Geneset Enrichment Analysis for the entire communities. 
#' @export
#'
#' @examples try(
#' CreateHyperGeoTestAll(cAMARETTOresults,
#'  cAMARETTOnetworkM,
#'  cAMARETTOnetworkC,
#'  hyper_geo_reference = './h.all.v6.2.symbols.gmt',
#'  driverGSEA=TRUE,
#'  NrCores=4)
#' )
CreateHyperGeoTestAll<-function(cAMARETTOresults,cAMARETTOnetworkM,
                                cAMARETTOnetworkC,hyper_geo_reference,
                                driverGSEA=TRUE,NrCores=4){
    j<-Community_key<-Type<-GeneNames<-Weights<-NULL
    geneset<-description<-Geneset<-n_Overlapping<-NULL
    n_RefGeneset<-Description<-Overlapping_genes<-overlap_perc<-NULL
    p_value<-padj<-Community<-Community_type<-NULL

    print("Performing Geneset Enrichment Analysis ...")
    com_gene_df<-suppressWarnings(ComRunModGenInfo(cAMARETTOresults,
                                                    cAMARETTOnetworkM,
                                                    cAMARETTOnetworkC))
    communities_all<-unique(com_gene_df$Community_key)
    all_hgt_output<-NULL
    if (is.null(cAMARETTOresults)){return(1)}
    cluster <- parallel::makeCluster(c(rep("localhost", NrCores)),
                                    type = "SOCK")
    doParallel::registerDoParallel(cluster, cores = NrCores)
    for(i in seq_len(length(hyper_geo_reference))){
        hgt_output_gmt<-NULL
        hgt_output_gmt<-foreach::foreach(j =seq_len(length(communities_all)),
                                .combine = "rbind") %dopar% {
            #for (ComNr in communities_all){
            #print(ComNr)
            ComNr<-communities_all[j]
    
            target_genes<-com_gene_df%>%
                dplyr::filter(Community_key==ComNr)%>%
                dplyr::filter(Type=="Target")%>%
                dplyr::arrange(GeneNames)%>%
                dplyr::pull(GeneNames)
    
            driver_genes<-com_gene_df%>%
                dplyr::filter(Community_key==ComNr)%>%
                dplyr::filter(Type=="Driver")%>%
                dplyr::arrange(GeneNames)%>%
                dplyr::pull(GeneNames)

            driver_genes_weights<-com_gene_df%>%
                dplyr::filter(Community_key==ComNr)%>%
                dplyr::filter(Type=="Driver")%>%
                dplyr::arrange(GeneNames)%>%
                dplyr::pull(Weights)

            if(driverGSEA){
                genelist<-unique(c(target_genes,driver_genes))
            }
            else{
                genelist<-unique(target_genes)
            }
            outputHGT <- HGTGeneEnrichmentList(genelist,
                                hyper_geo_reference[i],
                                NrCores = NrCores)
            if (nrow(outputHGT)>0){
                return(outputHGT%>%mutate(Community_key=ComNr))
            }
            else{
                return(NULL) 
                }
        }
    all_hgt_output<-rbind(all_hgt_output,hgt_output_gmt)
    cat("The hyper geometric test with ",
            hyper_geo_reference[i],
            " gmt file is calculated.\n")
    }
    #Community_key
    utils::data(MsigdbMapping)
    msigdbURL<-'<a href="http://software.broadinstitute.org/gsea/msigdb/cards/'
    MsigdbMapping<-MsigdbMapping%>%
        dplyr::mutate(url=paste0(msigdbURL,
                            geneset,
                            '">',
                            gsub("_"," ",geneset),
                            '</a>'))
    all_hgt_output<-all_hgt_output%>%
        dplyr::left_join(MsigdbMapping,by=c("Geneset"="geneset"))%>%
        dplyr::mutate(description=ifelse(is.na(description),
                            Geneset,description))%>%
        dplyr::mutate(Geneset=ifelse(is.na(url),Geneset,url))%>%
        dplyr::rename("Description"="description")%>%
        dplyr::select(-url)

    all_hgt_output <- all_hgt_output %>%
        dplyr::mutate(overlap_perc = n_Overlapping / n_RefGeneset) %>%
        dplyr::select(Community_key,Geneset, Description,
                n_RefGeneset, n_Overlapping,
                Overlapping_genes, overlap_perc, p_value, padj) %>%
        dplyr::arrange(padj)

    all_hgt_output<-all_hgt_output%>%
        left_join(com_gene_df%>%
            select(Community_key,Community,Community_type)%>%
            distinct(),by="Community_key")
    print("Geneset Enrichment Analysis is done!")
    return(all_hgt_output)
}

#' @title InitialCheckInputs
#'
#' @param cAMARETTOresults cAMARETTOresults object
#' @param output_address output_address 
#' @param HTMLsAMARETTOlist A list with AMARETTO reports to link
#'  with the Community AMARETTO report. If NULL, no links are added.
#' @param CopyAMARETTOReport Boolean to indicate if the AMARETTO reports
#'  needs to be copied in the AMARETTO report directory.
#'   In this way links are contained when moving the HTML directory.
#' @param hyper_geo_reference comupted tables or addresses to gmt files
#'  for functional enrichment
#' @param hyper_geo_reference_gp comupted tables or addresses to gmt files
#'  for genetic perturbation
#' @param hyper_geo_reference_cp comupted tables or addresses to gmt files
#'  for chemical perturbation
#'
#' @return RunInfo dataframe
#'
#' @examples 
#' try(
#' InitialCheckInputs(cAMARETTOresults,
#'  output_address="./",
#'  HTMLsAMARETTOlist,
#'  CopyAMARETTOReport=FALSE,
#'  hyper_geo_reference)
#' )
InitialCheckInputs<-function(cAMARETTOresults,
                            output_address,
                            HTMLsAMARETTOlist,
                            CopyAMARETTOReport,
                            hyper_geo_reference,
                            hyper_geo_reference_gp,
                            hyper_geo_reference_cp){
    if (!dir.exists(output_address)) {
        stop("Output directory is not existing.")
    } else {
        dir.create(file.path(output_address,"cAMARETTOhtmls"),
                    showWarnings = FALSE)
        full_path <- file.path(normalizePath(output_address),"cAMARETTOhtmls")
        print(paste0("The output directory is: ",full_path))
    }
    if (!is.null(hyper_geo_reference)) {
        if(is.character(hyper_geo_reference)){
            if (!file.exists(hyper_geo_reference[1])) {
                stop("GMT for hyper geometric test is not existing.")
            }
        }
    }

    if (!is.null(hyper_geo_reference_gp)) {
        if(is.character(hyper_geo_reference_gp)){
            if (!file.exists(hyper_geo_reference_gp[1])) {
                stop("GMT for genetic perturbation hyper geometric test is
                    not existing.")
            }
        }
    }

    if (!is.null(hyper_geo_reference_cp)) {
        if(is.character(hyper_geo_reference_cp)){
            if (!file.exists(hyper_geo_reference_cp[1])) {
                stop("GMT for chemical perturbation hyper geometric test is
                not existing.")
            }
        }
    }
    if(!is.null(HTMLsAMARETTOlist)){
        if(!all(names(HTMLsAMARETTOlist) %in% cAMARETTOresults$runnames)){
            stop("The RUN names don't match those of the cAMARETTOresults")
        }
        i<-1
        for(htmldir in HTMLsAMARETTOlist){
            if(!file.exists(htmldir)){
                stop(paste0("The AMARETTO ",names(HTMLsAMARETTOlist)[i],
                    " html directory is not existing."))
            }
            htmldir<-normalizePath(file.path(htmldir,"/AMARETTOhtmls/"))
            if (CopyAMARETTOReport==TRUE){
                dir.create(file.path(full_path,names(HTMLsAMARETTOlist)[i]),
                    showWarnings = FALSE)
                file.copy(htmldir, file.path(full_path,
                        names(HTMLsAMARETTOlist)[i]),recursive = TRUE)
                htmldir<-file.path(".",names(HTMLsAMARETTOlist)[i],
                                    "AMARETTOhtmls")
            }
            HTMLsAMARETTOlist[i]<-htmldir
            i=i+1
        }
    }
    else{
        HTMLsAMARETTOlist<-NULL
    }
    return(list(full_path=full_path,HTMLsAMARETTOlist=HTMLsAMARETTOlist))
}

#' @title CommunityModuleTableCreate
#'
#' @param cAMARETTOresults cAMARETTOresults object
#' @param cAMARETTOnetworkM cAMARETTOnetworkM object
#' @param HTMLsAMARETTOlist List contating the address to index page for
#'  each data
#' @param CopyAMARETTOReport Boolean to indicate if the AMARETTO reports
#'  needs to be copied in the AMARETTO report directory.
#'   In this way links are contained when moving the HTML directory.
#' @param cAMARETTOnetworkC cAMARETTOnetworkC object
#'
#' @return Community To Module Dataframe
#' @export
#'
#' @examples
#'try( 
#'CommunityModuleTableCreate (cAMARETTOresults,
#' cAMARETTOnetworkM,
#' cAMARETTOnetworkC,
#' RunInfo)
#')
CommunityModuleTableCreate<-function(cAMARETTOresults,
                                    cAMARETTOnetworkM,
                                    cAMARETTOnetworkC,
                                    HTMLsAMARETTOlist,
                                    CopyAMARETTOReport){

    ##################################### Bioconductor Considerations :
    Community_key <- Community <- Community_type <-NULL
    AMARETTOres<- Run_Names<- ModuleNr<-NULL
    ModuleLink <- numTotalEdgesInCommunity <-NULL
    fractEdgesInVsOut <- CommsizeFrac<-NULL
    #####################################

    com_gene_df<-suppressWarnings(ComRunModGenInfo(cAMARETTOresults,
                                                    cAMARETTOnetworkM,
                                                    cAMARETTOnetworkC))%>%
        dplyr::select(Community_key,Community,Community_type,
                AMARETTOres,Run_Names,ModuleNr)%>%distinct()

    ComModule<-com_gene_df %>%
        dplyr::mutate(ModuleLink=ModuleHyperLink(ModuleNr,Run_Names,
                                                AMARETTOres,
                                                HTMLsAMARETTOlist,
                                                CopyAMARETTOReport))

    ComModule<-ComModule%>%
        group_by(Community_key,Run_Names) %>%
        summarise(ModuleLinks=paste(ModuleLink,collapse = ", "))

    ComModule<-suppressMessages(reshape2::dcast(ComModule,
                                        Community_key~Run_Names,
                                        fill=0))

    ComModule<-suppressMessages(
        ComModule%>%
            left_join(com_gene_df%>%
            select(Community_key,Community,Community_type)%>%
            distinct()))

    ComModule<-ComModule%>%
        left_join(cAMARETTOnetworkC$commEdgeInfo%>%
        mutate(Community_key=Community)%>%
        select(Community_key,numTotalEdgesInCommunity,
                fractEdgesInVsOut,CommsizeFrac),
                by="Community_key")

    ComModule<-ComModule%>%
        dplyr::mutate(Community = CommunityHyperLink(Community,
                                        Community_key,Community_type))
    ComModule<-ComModule%>%
        select(Community,sort(cAMARETTOresults$runnames),everything())%>%
        select(-Community_key,-Community_type)
    return(ComModule)
}

#' @title CommunityHyperLink
#'
#' @param Community CommunityHyperLink
#' @param Community_key Community_key 
#' @param Community_type Community_type
#'
#' @return a hyperlink and presentable name for 
#' the communities used for datatables
#' @export
#'
#' @examples try(CommunityHyperLink(Community,Community_key,Community_type))
CommunityHyperLink<-function(Community,Community_key,Community_type){
    Comm_hyperLink<-paste0("<a href=\"./communities/",
                        paste0("Community","_",Community_key),
                        ".html\">",paste0(ifelse((Community_type==0),
                        "Community ",""),Community),"</a>")
    return(Comm_hyperLink)
}

#' @title DriversSharedTbl
#'
#' @param cAMARETTOresults  cAMARETTOresults object
#' @param cAMARETTOnetworkM cAMARETTOnetworkM object
#' @param cAMARETTOnetworkC cAMARETTOnetworkC object
#'
#' @return  Frequency of driver genes across different communities and
#'  datasets
#' @export
#'
#' @examples
#' try(
#' DriversSharedTbl(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC)
#' )
DriversSharedTbl<-function(cAMARETTOresults,
                            cAMARETTOnetworkM,
                            cAMARETTOnetworkC){

    com_gene_df<-suppressWarnings(ComRunModGenInfo(cAMARETTOresults,
                                                    cAMARETTOnetworkM,
                                                    cAMARETTOnetworkC))

    Type<-GeneNames<-Community_key<-Run_Names<-NULL
    all_genes<-Freq<-NULL
    genecardURL<-"<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene="
    ComDrivers<-com_gene_df%>%filter(Type=="Driver")%>%
        mutate(GeneNames = paste0(genecardURL,
                                    GeneNames, "\">",
                                    GeneNames,
                                    "</a>"))
    ComDrivers<-ComDrivers%>%
        group_by(Community_key,Run_Names)%>%
        summarise(GeneNames=paste(sort(unique(GeneNames)),collapse = ", "))

    ComDrivers<-suppressMessages(reshape2::dcast(ComDrivers,
                            Community_key~Run_Names,fill=""))
    
    ComDrivers<-ComDrivers%>%
        left_join(ComDrivers%>%
            tidyr::unite("all_genes",-one_of("Community_key"),sep = ", "),
            by="Community_key")
    
    gene_freq_df<-NULL
    for (ComNr in unique(ComDrivers$Community_key)){
        gene_table<-table(strsplit(gsub("","",ComDrivers%>%
                                filter(Community_key==ComNr)%>%
                                    pull(all_genes)),","))
        freq_tbl<-data.frame(GeneNames=names(gene_table),
                            Freq=paste0("# Data Sets = ",
                                        as.numeric(gene_table)),
                            Community_key=ComNr)
        gene_freq_df<-rbind(gene_freq_df,freq_tbl)
        }
    gene_freq_df<-gene_freq_df%>%
        group_by(Community_key,Freq)%>%
        summarise(GeneNames=paste(GeneNames,collapse = ", "))
    
    gene_freq_df<-suppressMessages(reshape2::dcast(gene_freq_df,
                                                    Community_key~Freq,
                                                    fill=""))
    # order columns
    columnnames<-c("Community_key",sort(colnames(gene_freq_df)[-1]))
    gene_freq_df<-gene_freq_df[,columnnames] 
    ComDrivers<-ComDrivers%>%
        left_join(gene_freq_df,by="Community_key")%>%
        dplyr::select(-all_genes)
    return(ComDrivers)
}


#' @title ModuleHyperLink
#' Returns hyperlink to the html report of a given module and dataset.
#' @param Module Module Nr ex. 4
#' @param Run_Names Run_Names ex. "TCGA_LIHC"
#' @param AMARETTOres AMARETTO result object
#' @param HTMLsAMARETTOlist list containing the hyperlink address
#'  to html report of each dataset. 
#' @param CopyAMARETTOReport Boolean to indicate if the AMARETTO reports
#'  needs to be copied in the AMARETTO report directory.
#'   In this way links are contained when moving the HTML directory.
#' @param page 1 for index pages and 2 for community pages
#'
#' @return Hyperlinks for Modules
#' @export
#'
#' @examples 
#' try(ModuleHyperLink(Module,Run_Names,
#'    AMARETTOres,HTMLsAMARETTOlist,
#'    CopyAMARETTOReport,page=1))
ModuleHyperLink<-function(Module,
                            Run_Names,
                            AMARETTOres,
                            HTMLsAMARETTOlist,
                            CopyAMARETTOReport,
                            page=1){
    # htmldir = "someLocalAdress/LIHC_Report_75/AMARETTOhtmls"
    # if CopyAMARETTOReport ==TRUE :
    #   htmldir = "./TCGA_LIHC/AMARETTOhtmls"
    ModuleLink<-mod_hyperlink<-NULL

    if(is.null(HTMLsAMARETTOlist)==FALSE){
        RunInfo<-(data.frame(Run_Names=names(HTMLsAMARETTOlist),
                                ModuleLink=HTMLsAMARETTOlist,
                                stringsAsFactors = FALSE))
        if(CopyAMARETTOReport==TRUE & page == 2){
            RunInfo <- RunInfo %>% mutate(ModuleLink=gsub("^./",
                                                            "../",
                                                            ModuleLink))
            }
        Mod_hyperLink<-data.frame(Run_Names=Run_Names,
                                Module=Module,
                                AMARETTOres=AMARETTOres,
                                stringsAsFactors = FALSE)%>%
            dplyr::left_join(RunInfo,by="Run_Names")%>%
            mutate(mod_hyperlink=ifelse(AMARETTOres==1,
                                        paste0("<a href=",
                                                ModuleLink,
                                                "/modules/",
                                                gsub("Module_",
                                                    "module",
                                                    Module),
                                                ".html>",gsub("Module_",
                                                                "Module ",
                                                                Module),
                                            "</a>"),Module))%>%
            dplyr::pull(mod_hyperlink)
        }
    else{
        Mod_hyperLink<-ifelse(AMARETTOres==1,
                                gsub("Module_","Module ",Module),Module)
        }
    return(Mod_hyperLink)
}
#' @title RunHyperLink
#'
#' @param Run_Names a vector of Run names
#' @param AMARETTOres list of AMARETTOres object
#' @param HTMLsAMARETTOlist A list with AMARETTO reports to link
#'  with the Community AMARETTO report. If NULL, no links are added.
#' @param CopyAMARETTOReport Boolean to indicate if the AMARETTO reports
#'  needs to be copied in the AMARETTO report directory.
#'   In this way links are contained when moving the HTML directory.
#' @param page 1 for index page and 2 for community page
#'
#' @return Hyperlinks for Run Names
#' @export
#'
#' @examples try(RunHyperLink(Run_Names, AMARETTOres,
#'  HTMLsAMARETTOlist, CopyAMARETTOReport, page=1))
#' 
RunHyperLink<-function(Run_Names,AMARETTOres,
                        HTMLsAMARETTOlist,CopyAMARETTOReport,
                        page=1){
    RunLink<-run_hyperlink<-NULL
    if(is.null(HTMLsAMARETTOlist)==FALSE){
        RunInfo<-(data.frame("Run_Names"=names(HTMLsAMARETTOlist),
                            RunLink=HTMLsAMARETTOlist,
                            stringsAsFactors = FALSE))
        if(CopyAMARETTOReport==TRUE & page == 2){
            RunInfo <- RunInfo %>% mutate(RunLink=gsub("^./","../",RunLink))
            }
        Run_hyperLink<-data.frame("Run_Names"=Run_Names,
                                    AMARETTOres=AMARETTOres,
                                    stringsAsFactors = FALSE)%>%
            dplyr::left_join(RunInfo,by="Run_Names")%>%
            mutate(run_hyperlink=ifelse(AMARETTOres==1,
                                        paste0("<a href=",
                                                RunLink,"/index.html>",
                                                Run_Names, "</a>"),
                                        Run_Names))%>%
            dplyr::pull(run_hyperlink)
    }
    else{
        Run_hyperLink<-Run_Names
        }
    return(Run_hyperLink)
}
#' @title create_hgt_datatable
#'
#' @param output_hgt hyper geoetric test table
#' @param com_table TRUE if it is for community page, FALSE if index page.
#' @param ComNr community number
#' @import DT
#' @return DataTable
#'
#' @examples try(create_hgt_datatable(output_hgt, com_table=FALSE, ComNr = 1))
create_hgt_datatable<-function(output_hgt,
                                com_table=FALSE, ComNr = 1){
    Community_key<-Community_type<-Community<-padj<-n_Overlapping<-NULL
    #==========================================================================
    options('DT.warn.size'=FALSE)
    buttons_list = list(list(extend ='csv'),
                        list(extend ='excel'),
                        list(extend = 'pdf', pageSize = 'A4',
                            orientation = 'landscape'),
                        list(extend ='print'),
                        list(extend ='colvis'))
    columnDefs = list(list(className = 'dt-head-center', targets = "_all"),
                    list(className = 'text-left',targets = "_all"))
    optionsList = list(deferRender=TRUE,
                        pageLength = 10,
                        lengthMenu = c(5, 10, 20, 50, 100),
                        keys = TRUE,
                        dom = "Blfrtip", 
                        buttons = buttons_list,
                        columnDefs = columnDefs,
                        paging = TRUE)
    #=======================================================================
    if (com_table){
        outputHGT<-output_hgt%>%
        filter(Community_key==ComNr)%>%
        select(-Community_key,-Community_type,-Community)
    if (nrow(outputHGT)>0){
        DTGSEA_colnames<-c("Gene Set Name",
                            "Gene Set Description",
                            "# Genes in Gene Set",
                            "# Genes in Overlap",
                            "Genes in Overlap",
                            "% Genes in overlap",
                            "P-value",
                            "FDR Q-value")
    outputHGT<-outputHGT%>%dplyr::arrange(padj)
    DTGSEA <- DT::datatable(outputHGT,
                            class = "display",
                            filter = 'top',
                            extensions = c('Buttons','KeyTable'),
                            rownames = FALSE,
                            options = optionsList,
                            colnames = DTGSEA_colnames ,
                            escape = FALSE) %>% 
        DT::formatSignif(c("p_value",
                            "padj",
                            "overlap_perc"), 2) %>%
        DT::formatStyle("overlap_perc",
                        background = styleColorBar(c(0, 1),"lightblue"),
                        backgroundSize = "98% 88%",
                        backgroundRepeat = "no-repeat",
                        backgroundPosition = "center")%>%
        DT::formatStyle(columns = c(5), fontSize = '60%')
    }
        else{
            DTGSEA <- "No significant overlap was identified 
                        in the geneset enrichment analysis."
            }
    return(DTGSEA)
        }
    else{
        output_hgt<-output_hgt%>%filter(n_Overlapping>2)%>%
            dplyr::mutate(Community=CommunityHyperLink(Community,
                                                    Community_key,
                                                    Community_type))%>%
            select(-Community_type,-Community_key)%>%
            select(Community,everything())%>%
            dplyr::arrange(padj,Community)
        DTGSEA_colnames<-c("Community","Gene Set Name",
                        "Gene Set Description","# Genes in Gene Set",
                        "# Genes in Overlap","Genes in Overlap",
                        "% Genes in overlap","P-value","FDR Q-value")
    #output_hgt<-output_hgt %>% 
    # dplyr::mutate(Geneset = 
    #paste0("<a href=\"http://software.broadinstitute.org/gsea/msigdb/cards/",
    #                                Geneset,
    #                                "\">",
    #                                gsub("_", " ", Geneset),
    #                                "</a>"))
        DTGSEAall <- DT::datatable(output_hgt,
                                    class = "display",
                                    filter = 'top',
                                    extensions = c('Buttons','KeyTable'),
                                    rownames = FALSE,
                                    options = optionsList,
                                    colnames = DTGSEA_colnames,
                                    escape = FALSE)%>%
        DT::formatSignif(c("p_value", "padj","overlap_perc"), 2) %>% 
        DT::formatStyle("overlap_perc",
                            background = styleColorBar(c(0, 1),"lightblue"),
                            backgroundSize = "98% 88%",
                            backgroundRepeat = "no-repeat",
                            backgroundPosition = "center")%>%
        DT::formatStyle(columns = c(6), fontSize = '60%')
    return(DTGSEAall)
    }
}

#' Title drivers_communities_summary
#'
#' @param cAMARETTOresults cAMARETTOresults object
#' @param cAMARETTOnetworkM  cAMARETTOnetworkM object
#' @param cAMARETTOnetworkC  cAMARETTOnetworkC object
#'
#' @return a list containing 3 tables summarising drivers-comunities relationships.
#' @export
#'
#' @examples try(drivers_communities_summary(cAMARETTOresults,cAMARETTOnetworkM,cAMARETTOnetworkC))
drivers_communities_summary<-function(cAMARETTOresults,cAMARETTOnetworkM,cAMARETTOnetworkC){
    com_run_mod_gene<-CommunityAMARETTO::ComRunModGenInfo(cAMARETTOresults,
                                                          cAMARETTOnetworkM,
                                                          cAMARETTOnetworkC)
    tbl_freq<-com_run_mod_gene%>%filter(Type=="Driver")%>%
        group_by(GeneNames)%>%
        summarise(freq=length(ModuleNr))%>%
        arrange(-freq)
    high_genes<-tbl_freq$GeneNames

    tbl1<-com_run_mod_gene%>%
        filter(Type=="Driver")%>%
        mutate(ModuleNr=gsub("Module_","",ModuleNr))%>%
        select(GeneNames,Run_Names,ModuleNr,Type,Weights,Community)
    tbl2<-com_run_mod_gene%>%
        filter(Type=="Driver")%>%
        group_by(GeneNames,Run_Names)%>%
        summarise(freq=length(ModuleNr),modules=paste(ModuleNr,collapse = ","),communities=paste(Community,collapse = ","))%>%
        mutate(modules=gsub("Module_","",modules))%>%
        arrange(factor(GeneNames, levels = high_genes))
    tbl3<-com_run_mod_gene%>%
        filter(Type=="Driver")%>%
        mutate(modules=gsub("Module_","",ModuleNr))%>%
        mutate(driver_summary=paste0("[",Run_Names,", "," M",modules,", ","C",Community,", W=",signif(Weights,1),"]"))%>%
        select(GeneNames,driver_summary)%>%
        group_by(GeneNames)%>%
        summarise(freq=length(driver_summary),summary=paste(driver_summary,collapse = ","))%>%
        arrange(factor(GeneNames, levels = high_genes))
    return(list(tbl1=tbl1,
                tbl2=tbl2,
                tbl3=tbl3))
}


#' Title HubAuthority_scores_CommunityAMARETTO
#'
#' @param cAMARETTOresults cAMARETTOresults object
#' @param cAMARETTOnetworkM cAMARETTOnetworkM object
#' @param cAMARETTOnetworkC cAMARETTOnetworkC object
#'
#' @return a list containing Hub-score authority-score of genes and modules.
#' @importFrom igraph hub_score authority_score
#' @export
#' @examples try(Hubscore_CommunityAMARETTO(cAMARETTOresults,cAMARETTOnetworkM,cAMARETTOnetworkC))
HubAuthority_scores_CommunityAMARETTO<-function(cAMARETTOresults,cAMARETTOnetworkM,cAMARETTOnetworkC){
    ComRunModGenTbl<-CommunityAMARETTO::ComRunModGenInfo(cAMARETTOresults,
                                                         cAMARETTOnetworkM,
                                                         cAMARETTOnetworkC)
    from<-ComRunModGenTbl$GeneNames
    to<-paste0(ComRunModGenTbl$Run_Names,"&",ComRunModGenTbl$ModuleNr)
    weight<-ComRunModGenTbl$Weights
    df<-data.frame(from=from,to=to,weight=abs(weight))
    g <- graph_from_data_frame(df, directed=TRUE)
    ## calculate hub score
    hub_scores_vector<-hub_score(g, scale = TRUE)
    hub_scores_vector<-sort(hub_scores_vector$vector,decreasing = TRUE)
    hub_scores_df<-data.frame(genes=names(hub_scores_vector),
                              hub_score=hub_scores_vector,
                              rank=c(1:length(hub_scores_vector)))
    rownames(hub_scores_df)<-NULL
    ## calculate authority score
    authority_scores_vector<-authority_score(g, scale = TRUE)
    authority_scores_vector<-sort(authority_scores_vector$vector,decreasing = TRUE)
    authority_scores_df<-data.frame(genes=names(authority_scores_vector),
                                    hub_score=authority_scores_vector,
                                    rank=c(1:length(authority_scores_vector)))

    return(list("hub_score"=hub_scores_df,
                "authority_score"=authority_scores_df))
}


#' ########
#' 
#' # try(
#' #  #cAMARETTOsList<-readRDS (file="./outputs/cAMARETTO_Liver2DS.rds")
#' #  #communityReportURL<-c("http://portals.broadinstitute.org/
#' pochetlab/demo/cAMARETTO_Liver_2DS/")
#' #  cAMARETTO_Cytoscape(cAMARETTOsList,
#' communityReportURL = "",
#' cytoscape_name="cAMARETTO_Liver2DS")
#' # )
#' 
#' cAMARETTO_Cytoscape<-function(cAMARETTOsList,
#' communityReportURL = "",cytoscape_name="my_cytoscape"){
#' 
#'   graph<-cAMARETTOsList$cAMARETTOnetworkC$CommGraph
#'   runnames<-cAMARETTOsList$cAMARETTOresults$runnames
#'   runURLs<-data.frame(run=runnames)
#'   runURLs<-runURLs%>%dplyr::mutate(run_URL=paste0(communityReportURL,run))
#'   
#'   nodes_df<-igraph::as_data_frame(graph, what="vertices")
#'   nodes_df<-suppressWarnings(nodes_df%>%
#'                                dplyr::left_join(runURLs,by="run"))
#'   nodes_df<-nodes_df%>%
#'   dplyr::mutate(Module_name=unlist(map(strsplit(name,"\\|"),2)))%>%
#'     dplyr::mutate(Module_name=gsub("Module_","module",Module_name))
#'   nodes_df<-nodes_df%>%
#'     dplyr::mutate(URL=ifelse(run %in% runnames,
#'                              paste0('<a href="',run_URL,
#'                              '/AMARETTOhtmls/modules/',Module_name,
#'                              '.html">URL</a>'),
#'                              ''))
#'   nodes_df<-nodes_df%>%dplyr::select(-c("Module_name","run_URL",
#'   "NewComNumber"))
#'   
#'   edges_df<-igraph::as_data_frame(graph, what="edges")
#'   edges_df<-suppressWarnings(edges_df%>%
#'      dplyr::mutate(from_Run=unlist(purrr::map(strsplit(from,"\\|"),1)))%>%
#'      dplyr::mutate(from_Module=unlist(map(strsplit(from,"\\|"),2)))%>%
#'      dplyr::mutate(from_Module=gsub("Module_","module",from_Module))%>%
#'      dplyr::left_join(runURLs,by=c("from_Run"="run"))%>%
#'      dplyr::rename(from_run_URL=run_URL)%>%
#'      dplyr::mutate(from_URL=ifelse(from_Run %in% runnames,
#'                      paste0('<a href="',
#'                             from_run_URL,
#'                             '/AMARETTOhtmls/modules/',
#'                              from_Module,
#'                              '.html">URL</a>'),
#'                               '')))
#'   edges_df<-suppressWarnings(edges_df%>%
#'            dplyr::mutate(to_Run=unlist(
#'            purrr::map(strsplit(to,"\\|"),1)))%>%
#'            dplyr::mutate(to_Module=unlist(map(strsplit(to,"\\|"),2)))%>%
#'            dplyr::mutate(to_Module=gsub("Module_","module",to_Module))%>%
#'            dplyr::left_join(runURLs,by=c("to_Run"="run"))%>%
#'            dplyr::rename(to_run_URL=run_URL)%>%
#'            dplyr::mutate(to_URL=ifelse(to_Run %in% runnames,
#'                      paste0('<a href="',
#'                      to_run_URL,
#'                      '/AMARETTOhtmls/modules/',
#'                      to_Module,
#'                      '.html">URL</a>'),
#'                       '')))
#'   edges_df<-edges_df%>%
#'     dplyr::select(-c("from_Run","from_Module","from_run_URL","to_Run",
#'     "to_Module","to_run_URL"))%>%
#'     dplyr::rename(source_URL=from_URL,target_URL=to_URL)
#' 
#'   graph_new<-igraph::graph_from_data_frame(edges_df,
#'    directed=FALSE, vertices=nodes_df)
#'   try(RCy3::createNetworkFromIgraph(graph_new,cytoscape_name),silent=TRUE)
#'   return(graph_new)
#' }


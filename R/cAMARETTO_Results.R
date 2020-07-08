#' @title cAMARETTO_Results
#' 
#' To initiate the community AMARETTO this results functions performs
#'  hyper geometric tests between all modules and gene sets.
#'
#' @param AMARETTOresults_all A list of multiple
#'  AMARETTO_Run outputs. The names are run names.
#' @param NrCores Nr of Cores that can be used to
#'  calculated the results.
#' @param output_dir A directory that stores gmt files
#'  for the runs
#' @param gmt_filelist A list with gmt files that are added into
#'  the communities. The names of the list are used in the networks.
#'   NULL if no list is added.
#' @param drivers Boolean that defines if only targets or drivers
#'  and targets are used to calculate the HGT.
#'
#' @return a list with cAMARETTOresults data objects from multiple runs
#' @importFrom gtools combinations
#' @importFrom dplyr arrange group_by left_join mutate select summarise 
#' rename  filter everything pull distinct case_when
#' @importFrom stats p.adjust phyper
#' @importFrom devtools install_github
#' @examples 
#' try(
#' cAMARETTOresults <- cAMARETTO_Results(AMARETTOresults_all,
#'   gmt_filelist=list(ImmuneSignature = Cibersortgmt),
#'   NrCores = 4 ,
#'   output_dir = "./")
#' )
#' @export
cAMARETTO_Results <- function(AMARETTOresults_all,
                                NrCores=1,
                                output_dir="./",
                                gmt_filelist=NULL,
                                drivers = FALSE){

    #test if names are matching
    RunName1<-Geneset1<-RunName2<-Geneset2<-NULL
    runnames <- names(AMARETTOresults_all)
    if (!length(unique(runnames)) == length(runnames)){
        stop("The run names are not unique. Give unique names.")
    }
    dir.create(file.path(output_dir, "gmt_files"), recursive = FALSE,
                showWarnings = FALSE)
    # for each file a gmt for the modules
    create_gmt_filelist<-c()

    for (run in runnames){
        gmt_file <- file.path(output_dir,
                    "gmt_files",paste0(run, "_modules.gmt"))
        GmtFromModules(AMARETTOresults_all[[run]], gmt_file, run,
                        Drivers = drivers)
        create_gmt_filelist <- c(create_gmt_filelist,gmt_file)
    }
    names(create_gmt_filelist)<-runnames

    # add extra gmt files to compare with
    given_gmt_filelist <- c()
    if (!is.null(gmt_filelist)){
        for (gmt_file in gmt_filelist){
            gmt_file_path <- file.path(gmt_file)
            given_gmt_filelist <- c(given_gmt_filelist,gmt_file_path)
        }
    }
    names(given_gmt_filelist)<-names(gmt_filelist)
    all_gmt_files_list <- c(create_gmt_filelist,given_gmt_filelist)

    if (
        ! length(unique(names(all_gmt_files_list))) == 
        length(names(all_gmt_files_list))){
            stop("There is overlap between the gmt file names and run names")
    }

    if (length(all_gmt_files_list) < 2){
        stop("There are none or only one group given,
            community AMARETTO needs at least two groups.")
    }
    # compare gmts pairwise between runs
    all_run_combinations <- as.data.frame(gtools::combinations(
        n=length(all_gmt_files_list),
        r=2,
        v=names(all_gmt_files_list),
        repeats.allowed=FALSE))
    output_hgt_allcombinations <- apply(all_run_combinations, 1, function(x) {
        gmt_run1 <- all_gmt_files_list[x["V1"]]
        gmt_run2 <- all_gmt_files_list[x["V2"]]
        output_hgt_combination <- HyperGTestGeneEnrichment(gmt_run1,
                gmt_run2,
                runname1=as.character(x["V1"]),
                runname2=as.character(x["V2"]),
                NrCores=NrCores)
        return(output_hgt_combination)
        }
    )
    genelists<-lapply(all_gmt_files_list, function(x){
        genelist<-readGMT(x)
        names(genelist)<-gsub(" ","_",names(genelist))
        return(genelist)
        }
    )

    output_hgt_allcombinations <- do.call(rbind, output_hgt_allcombinations)
    output_hgt_allcombinations$padj <- stats::p.adjust(
    output_hgt_allcombinations$p_value,
        method="BH")
    output_hgt_allcombinations <- output_hgt_allcombinations %>% 
        dplyr::mutate(p_value=dplyr::case_when(
            Geneset1 == Geneset2~NA_real_, TRUE~p_value))
    output_hgt_allcombinations <- output_hgt_allcombinations %>%
        dplyr::mutate(Geneset1=ifelse(RunName1%in%names(given_gmt_filelist),
                            paste0(RunName1,"|",gsub(" ","_",Geneset1))
                            ,Geneset1),
                    Geneset2=ifelse(RunName2%in%names(given_gmt_filelist),
                            paste0(RunName2,"|",gsub(" ","_",Geneset2)),
                            Geneset2))

    # Extract relationship between genes and modules for all
    #AMARETTo files and the given gmt files.
    all_genes_modules_df<-Extract_Genes_Modules_All(AMARETTOresults_all,
                                        gmt_filelist)

    return(list(runnames=runnames,
                gmtnames=names(given_gmt_filelist),
                hgt_modules=output_hgt_allcombinations,
                genelists = genelists,
                all_genes_modules_df=all_genes_modules_df,
                NrCores=NrCores))
}

#' @title GmtFromModules
#'
#' @param AMARETTOresults A AMARETTO_Run output.
#' @param gmt_file A gmtfilename
#' @param run A runname
#' @param Drivers Add drivers to the gmt-file
#'
#' @return Creates a gmt file for a AMARETTO run
#' @importFrom utils write.table stack
#' @importFrom dplyr select arrange
#' @examples 
#' try(
#' GmtFromModules(AMARETTOresults,gmt_file,
#' run,Drivers=FALSE)
#' )
#' @export
GmtFromModules <- function(AMARETTOresults,gmt_file,run,Drivers=FALSE){

    Type<-GeneNames<-ModuleNr<-NULL
    ModuleMembership<-ExtractGenesInfo(AMARETTOresults,run)
    if (Drivers == FALSE){
        ModuleMembership<-dplyr::filter(ModuleMembership,Type=="Target")
    }
    ModuleMembership<-ModuleMembership%>%
        dplyr::select(GeneNames,ModuleNr)%>%
        dplyr::arrange(GeneNames)
    ModuleMembers_list<-split(ModuleMembership$GeneNames,
                                ModuleMembership$ModuleNr)
    names(ModuleMembers_list)<-paste0(run,"|Module_",names(ModuleMembers_list))
    utils::write.table(sapply(names(ModuleMembers_list),
                    function(x) paste(x,paste(ModuleMembers_list[[x]],
                    collapse="\t"),sep="\t")),
                    gmt_file,
                    quote = FALSE,
                    row.names = TRUE,
                    col.names = FALSE,
                    sep='\t')
}
#' @title readGMT
#'
#' @param filename A gmtfilename
#' @examples
#' try(
#' readGMT("file.gmt")
#' )
#' @return Reads a gmt file
#' @export
readGMT<-function(filename){
    gmtLines<-strsplit(readLines(filename),"\t")
    gmtLines_genes <- lapply(gmtLines, tail, -2)
    names(gmtLines_genes) <- sapply(gmtLines, head, 1)
    return(gmtLines_genes)
}

#' @title GmtFromModules
#'
#' @param gmtfile1 A gmtfilename that you want to compare
#' @param gmtfile2 A second gmtfile to compare with.
#' @param runname1 name of the first dataset.
#' @param runname2 name of the second dataset.
#' @param NrCores  Number of cores for parallel computing. 
#' @param ref.numb.genes The reference number of genes.
#' @importFrom doParallel registerDoParallel 
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom stats p.adjust phyper
#' 
#' @return Creates resultfile with p-values and padj when
#'  comparing two gmt files with a hyper geometric test.
#' @examples 
#' try(
#' HyperGTestGeneEnrichment(gmtfile1, gmtfile2,
#' runname1, runname2,NrCores, ref.numb.genes=45956)
#' )
#' @export
HyperGTestGeneEnrichment<-function(gmtfile1, gmtfile2,
                                    runname1, runname2,
                                    NrCores, ref.numb.genes=45956){
    i<-j<-NULL
    gmtfile1<-readGMT(gmtfile1) # our gmt_file_output_from Amaretto
    gmtfile2<-readGMT(gmtfile2)  # the hallmarks_and_co2...
    ###########################  Parallelizing :
    cluster <- parallel::makeCluster(c(rep("localhost", NrCores)),
                                    type = "SOCK")
    doParallel::registerDoParallel(cluster,cores=NrCores)
    resultloop<-foreach::foreach(j=seq_len(length(gmtfile2)),
                            .combine='rbind') %do% {
    #print(j)
    foreach::foreach(i=seq_len(length(gmtfile1)),.combine='rbind') %dopar% {
        #print(i)
        l<-length(gmtfile1[[i]])
        k<-sum(gmtfile1[[i]] %in% gmtfile2[[j]])
        m<-ref.numb.genes
        n<-length(gmtfile2[[j]])
        p1<-stats::phyper(k-1,l,m-l,n,lower.tail=FALSE)
        overlapping.genes<-gmtfile1[[i]][gmtfile1[[i]] %in% gmtfile2[[j]]]
        overlapping.genes<-paste(overlapping.genes,collapse = ', ')
        c(RunName1=runname1,
            RunName2=runname2,
            Geneset1=names(gmtfile1[i]),
            Geneset2=names(gmtfile2[j]),
            p_value=p1,n_Overlapping=k,
            Overlapping_genes=overlapping.genes)
    }
    }
    parallel::stopCluster(cluster)
    resultloop<-as.data.frame(resultloop,stringsAsFactors=FALSE)
    resultloop$p_value<-as.numeric(resultloop$p_value)
    resultloop$n_Overlapping<-as.numeric((resultloop$n_Overlapping))
    resultloop[,"padj"]<-stats::p.adjust(resultloop[,"p_value"],method='BH')
    return(resultloop)
}

#' Title  ExtractGenesInfo 
#'
#' @param AMARETTOresults AMARETTOreults object from an AMARETTO run.
#' @param run the run-name string , for example :"TCGA_LIHC",
#'  associated with the AMARETTOreults run. 
#'
#' @return returns a dataframe with columns of runname,
#'  module-number, genename, gene-types, and weights of the driver genes.
#' @export
#'
#' @examples 
#' try(
#' ExtractGenesInfo(AMARETTOresults,"TCGA-LIHC")
#' )
ExtractGenesInfo<-function(AMARETTOresults,run){
    ModuleMembership<-NULL
    for (ModuleNr in seq_len(AMARETTOresults$NrModules)){
        Targets<- names(
        AMARETTOresults$ModuleMembership[which(
        AMARETTOresults$ModuleMembership==ModuleNr),1])
        if (is.null(Targets)){next}
        Target_df<-data.frame(Run_Names=run,
                                ModuleNr=ModuleNr,
                                GeneNames=Targets,
                                Type="Target",
                                Weights=0,
                                stringsAsFactors = FALSE) 
        Drivers <- names(which(
            AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0))
        Weight_Driver<-AMARETTOresults$RegulatoryPrograms[ModuleNr,Drivers]
        Drivers_df<-data.frame(Run_Names=run,
                            ModuleNr=ModuleNr,
                            GeneNames=Drivers,
                            Type="Driver",
                            Weights=Weight_Driver,
                            stringsAsFactors = FALSE) 
    ModuleMembership <- rbind(ModuleMembership,Target_df,Drivers_df)
    }
    return(ModuleMembership)
}

#' Title Extract_Genes_Modules_All
#'
#' @param AMARETTOresults_all list of 
#' AMARETTOresults objects for different AMARETTO runs
#' @param gmt_filelist addresses to gmt files.
#'
#' @importFrom utils stack
#' @importFrom dplyr mutate rename select
#' @return a dataframe for all AMARETTO files and
#'  given gmt file with the following structure :  runname, ModuleNr, 
#'  genename, gene-types, and weights of the driver genes.
#' @export
#' @examples try(Extract_Genes_Modules_All(AMARETTOresults_all,gmt_filelist))
Extract_Genes_Modules_All<-function(AMARETTOresults_all,gmt_filelist){
    ModuleNr<-values<-Run_Names<-GeneNames<-NULL
    Type<-Weights<-AMARETTOres<-NULL
    all_genes_modules_df<-ind<-NULL
    for (run in names(AMARETTOresults_all)){
        all_genes_modules_df<-rbind(all_genes_modules_df,
                        ExtractGenesInfo(AMARETTOresults_all[[run]],run))
    }
    all_genes_modules_df<-all_genes_modules_df%>%
    dplyr::mutate(ModuleNr=paste0("Module_",ModuleNr))%>%
    dplyr::mutate(AMARETTOres=1)
    for (run in names(gmt_filelist)){
        gmt_genes_df<-utils::stack(readGMT(gmt_filelist[[run]]))%>%
            dplyr::mutate(Run_Names=run)%>%
            dplyr::rename(ModuleNr=ind)%>%
            dplyr::rename(GeneNames=values)%>%
            dplyr::mutate(Type="Target")%>%
            dplyr::mutate(Weights=0)%>%
            dplyr::mutate(AMARETTOres=0)%>%
            dplyr::select(Run_Names,ModuleNr,GeneNames,
                            Type,Weights,AMARETTOres)
    all_genes_modules_df<-rbind(all_genes_modules_df,gmt_genes_df)
    }
    return(all_genes_modules_df)
}
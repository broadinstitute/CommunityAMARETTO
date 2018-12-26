#' cAMARETTO_Results
#'
#' @param AMARETTOinit_all A list of multiple AMARETTO_Initialize outputs. The names are run names.
#' @param AMARETTOresults_all A list of multiple AMARETTO_Run outputs. The names are run names.
#' @param parallelparam BiocParallel BPPARAM
#' @param output_dir A directory that stores gmt files for the runs
#'
#' @return a list with AMARETTOinit and AMARETTOresults data objects from multiple runs
#' @import gtools
#' @import tidyverse

cAMARETTO_Results <- function(AMARETTOinit_all,AMARETTOresults_all,nrCores=1,output_dir="./",gmt_filelist=NULL){
  
  #test if names are matching
  if (all(names(AMARETTOinit_all) == names(AMARETTOresults_all))) {
    runnames <- names(AMARETTOresults_all)
    if (!length(unique(runnames)) == length(runnames)){
      stop("The run names are not unique. Give unique names.")
    }
  } else {
    stop("The names of the lists are not matching.")
  }
  
  dir.create(file.path(output_dir, "gmt_files"), recursive = FALSE, showWarnings = FALSE)
  
  # for each file a gmt for the modules
  create_gmt_filelist<-c()
  for (run in runnames){
    gmt_file <- file.path(output_dir,"gmt_files", paste0(run, "_modules.gmt"))
    GmtFromModules(AMARETTOinit_all[[run]], AMARETTOresults_all[[run]], gmt_file, run)
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
  
  if (! length(unique(names(all_gmt_files_list))) == length(names(all_gmt_files_list))){
    stop("There is overlap between the gmt file names and run names")
  }
  
  if (length(all_gmt_files_list) < 2){
    stop("There are none or only one group given, community AMARETTO needs at least two groups.")
  }
  
  # compare gmts pairwise between runs
  names_groups <- names(all_gmt_files_list)
  all_run_combinations <- as.data.frame(combinations(n=length(all_gmt_files_list), r=2, v=names_groups, repeats.allowed=F))
  
  output_hgt_allcombinations <- apply(all_run_combinations, 1, function(x) {
    gmt_run1 <- file.path(output_dir, "gmt_files", paste0(x["V1"], "_modules.gmt"))
    gmt_run2 <- file.path(output_dir, "gmt_files", paste0(x["V2"], "_modules.gmt"))
    output_hgt_combination <- HyperGTestGeneEnrichment(gmt_run1, gmt_run2)
    return(output_hgt_combination)
  })
  
  output_hgt_allcombinations <- do.call(rbind, output_hgt_allcombinations)
  output_hgt_allcombinations$padj <- p.adjust(output_hgt_allcombinations$p_value, method="BH")
  output_hgt_allcombinations <- output_hgt_allcombinations %>% 
                                    mutate(p_value=case_when(Geneset == Testset~NA_real_, TRUE~p_value))
  return(list(runnames=runnames, hgt_modules=output_hgt_allcombinations, nrCores))
}

#' GmtFromModules
#'
#' @param AMARETTOinit_all A list of multiple AMARETTO_Initialize outputs. The names are run names.
#' @param AMARETTOresults_all A list of multiple AMARETTO_Run outputs. The names are run names.
#' @param gmt_file A gmtfilename
#' @param run A runname
#'
#' @return Creates a gmt file for a AMARETTO run
#'
#' @import tidyverse

GmtFromModules <- function(AMARETTOinit,AMARETTOresults,gmt_file,run){
  
  ModuleMembership<-rownames_to_column(as.data.frame(AMARETTOresults$ModuleMembership),"GeneNames")
  NrModules<-AMARETTOresults$NrModules
  ModuleMembership<-ModuleMembership %>% arrange(GeneNames)
  
  ModuleMembers_list<-split(ModuleMembership$GeneNames,ModuleMembership$ModuleNr)
  names(ModuleMembers_list)<-paste0(run,"_Module_",names(ModuleMembers_list))
  
  write.table(sapply(names(ModuleMembers_list),function(x) paste(x,paste(ModuleMembers_list[[x]],collapse="\t"),sep="\t")),gmt_file,quote = FALSE,row.names = TRUE,col.names = FALSE,sep='\t')
}

#' readGMT
#'
#' @param filename A gmtfilename
#' @examples
#' readGMT("file.gmt")
#' @return Reads a gmt file

readGMT<-function(filename){
  gmtLines<-strsplit(readLines(filename),"\t")
  gmtLines_genes <- lapply(gmtLines, tail, -2)
  names(gmtLines_genes) <- sapply(gmtLines, head, 1)
  return(gmtLines_genes)
}

#' GmtFromModules
#'
#' @param gmtfile1 A gmtfilename that you want to compare
#' @param gmtfile2 A second gmtfile to compare with.
#' @param ref.numb.genes The reference number of genes.
#' 
#' @return Creates resultfile with p-values and padj when comparing two gmt files with a hyper geometric test.
#'
#' @import doParallel

HyperGTestGeneEnrichment<-function(gmtfile,testgmtfile,NrCores,ref.numb.genes=45956){
  
  test.gmt<-readGMT(testgmtfile) # our gmt_file_output_from Amaretto
  gmt.path<-readGMT(gmtfile)  # the hallmarks_and_co2...
  
  ###########################  Parallelizing :
  cluster <- makeCluster(c(rep("localhost", NrCores)), type = "SOCK")
  registerDoParallel(cluster,cores=NrCores)
  
  resultloop<-foreach(j=1:length(test.gmt), .combine='rbind') %do% {
    #print(j)
    foreach(i=1:length(gmt.path),.combine='rbind') %dopar% {
      #print(i)
      l<-length(gmt.path[[i]])
      k<-sum(gmt.path[[i]] %in% test.gmt[[j]])
      m<-ref.numb.genes
      n<-length(test.gmt[[j]])
      p1<-phyper(k-1,l,m-l,n,lower.tail=FALSE)
      
      overlapping.genes<-gmt.path[[i]][gmt.path[[i]] %in% test.gmt[[j]]]
      overlapping.genes<-paste(overlapping.genes,collapse = ', ')
      c(Geneset=names(gmt.path[i]),Testset=names(test.gmt[j]),p_value=p1,n_Overlapping=k,Overlapping_genes=overlapping.genes)
    }
  }
  
  stopCluster(cluster)
  resultloop<-as.data.frame(resultloop,stringsAsFactors=FALSE)
  resultloop$p_value<-as.numeric(resultloop$p_value)
  resultloop$n_Overlapping<-as.numeric((resultloop$n_Overlapping))
  resultloop[,"padj"]<-p.adjust(resultloop[,"p_value"],method='BH')
  return(resultloop)
}

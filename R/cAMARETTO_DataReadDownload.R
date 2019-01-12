#' @title cAMARETTO_Read
#'
#' @param AMARETTOdirectories 
#'
#' @return a list with the AMARETTOinit and AMARETTOresults
#' 
#' @import rlist
#' @export
cAMARETTO_Read<-function(AMARETTOdirectories){
  directories_to_unzip <- AMARETTOdirectories[grepl("zip$",AMARETTOdirectories)]
  # unzip directories if needed
  if(length(directories_to_unzip)>0){
    for(directory_to_unzip in directories_to_unzip){
      if(file.exists(directory_to_unzip)){
        unzip(directory_to_unzip)
      } else {
        stop(paste0("The ",directory_to_unzip," is not existing"))
      }
    }
  }
  
  AMARETTOdirectories <- lapply(AMARETTOdirectories, function(x) sub(".zip","",x))
  AMARETTOinit_all<-list()
  AMARETTOresults_all<-list()
  i = 1
  for (AMARETTOdirectory in AMARETTOdirectories){
    tmp = load(paste0(AMARETTOdirectory,"/amarettoInit.RData"))
    assign(paste0("AMARETTOinit_",names(AMARETTOdirectories)[i]),get(tmp))
    rm(tmp)
    AMARETTOinit_all<-list.append(AMARETTOinit_all,eval(parse(text=paste0("AMARETTOinit_",names(AMARETTOdirectories)[i]))))
    tmp = load(paste0(AMARETTOdirectory,"/amarettoResults.RData"))
    assign(paste0("AMARETTOresults_",names(AMARETTOdirectories)[i]),get(tmp))
    rm(tmp)
    AMARETTOresults_all<-list.append(AMARETTOresults_all,eval(parse(text=paste0("AMARETTOresults_",names(AMARETTOdirectories)[i]))))
    i = i + 1
  }
  names(AMARETTOinit_all)<-names(AMARETTOdirectories)
  names(AMARETTOresults_all)<-names(AMARETTOdirectories)
  return(list(AMARETTOinit_all = AMARETTOinit_all,AMARETTOresults_all = AMARETTOresults_all))
}

#' @title cAMARETTO_ExportResults
#'
#' @param cAMARETTOresults
#' @param cAMARETTOnetworkM
#' @param cAMARETTOnetworkC
#' @param output_address
#'
#' @return 
#' 
#' @export
cAMARETTO_ExportResults<-function(cAMARETTOresults,cAMARETTOnetworkM, cAMARETTOnetworkC, output_address="./"){

  if (!dir.exists(output_address)){
    stop("Output directory is not existing.")
  }
  
  #add a date stamp to the output directory
  output_dir<-paste0("cAMARETTOresults_",gsub("-|:","",gsub(" ","_",Sys.time())))
  dir.create(file.path(output_address,output_dir))

  #save rdata files for AMARETTO_Run and AMARETTO_Initialize output
  save(cAMARETTOresults, file=file.path(output_address,output_dir,"/cAMARETTOresults.RData"))
  save(cAMARETTOnetworkM, file=file.path(output_address,output_dir,"/cAMARETTOnetworkM.RData"))
  save(cAMARETTOnetworkC, file=file.path(output_address,output_dir,"/cAMARETTOnetworkC.RData"))

  zip(zipfile = file.path(output_address,output_dir),files=file.path(output_address,output_dir))
  unlink(output_dir,recursive = TRUE)
}

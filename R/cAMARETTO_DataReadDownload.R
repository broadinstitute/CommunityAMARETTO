#' @title cAMARETTO_Read
#' 
#' Reads the AMARETTOinit and AMARETTOresults from (zipped)
#'  AMARETTO directories. The names of the list are assigned as run names.
#'
#' @param unzipParentDirectory the address where the files are unzipped in
#' @param AMARETTOdirectories a list of AMARETTO directories 
#'
#' @return a list with the AMARETTOinit and AMARETTOresults
#' 
#' @importFrom  rlist list.append
#' @importFrom  stringr str_count
#' @importFrom  utils zip unzip
#' @examples
#' try(
#'  AMARETTO_all <- cAMARETTO_Read(list(LIHC="LIHC.zip",BLCA="BLCA.zip"))
#'  )
#' @export
cAMARETTO_Read<-function(AMARETTOdirectories,unzipParentDirectory=getwd()){
  
  directories_to_unzip <- AMARETTOdirectories[grepl("zip$",AMARETTOdirectories)]
  # unzip directories if needed
  if(length(directories_to_unzip)>0){
    for(i in 1:length(directories_to_unzip)){
      if(file.exists(directories_to_unzip[[i]])){
        if (
          stringr::str_count(directories_to_unzip[[i]],
          "/")>=str_count(directories_to_unzip[[i]],
              "\\\\")){
          #linux or windows addresses
          extdir=paste(sub("/$","",
                           unzipParentDirectory),
                       names(directories_to_unzip)[i],sep="/")
        }
        else{
          extdir=paste(sub("\\\\$","",unzipParentDirectory),
                       names(directories_to_unzip)[i],sep="\\")
        }
        unzip(directories_to_unzip[[i]],exdir=extdir,junkpaths = TRUE)
        AMARETTOdirectories[[i]]<-extdir
      } else {
        stop(paste0("The ",directories_to_unzip[[i]]," is not existing"))
      }
    }
  }
  
  AMARETTOdirectories <- lapply(AMARETTOdirectories,
                                function(x) sub(".zip","",x))
  AMARETTOinit_all<-list()
  AMARETTOresults_all<-list()
  i = 1
  for (AMARETTOdirectory in AMARETTOdirectories){
    tmp = load(paste0(AMARETTOdirectory,"/amarettoInit.RData"))
    assign(paste0("AMARETTOinit_",names(AMARETTOdirectories)[i]),get(tmp))
    rm(tmp)
    AMARETTOinit_all<-rlist::list.append(AMARETTOinit_all,
                                      eval(parse(text=paste0("AMARETTOinit_",
                                              names(AMARETTOdirectories)[i])),
                                              envir = environment()))
    tmp = load(paste0(AMARETTOdirectory,"/amarettoResults.RData"))
    assign(paste0("AMARETTOresults_",names(AMARETTOdirectories)[i]),get(tmp))
    rm(tmp)
    AMARETTOresults_all<-list.append(AMARETTOresults_all,eval(parse(
      text=paste0("AMARETTOresults_",
                           names(AMARETTOdirectories)[i])),
                           envir = environment()))
    i = i + 1
  }
  names(AMARETTOinit_all)<-names(AMARETTOdirectories)
  names(AMARETTOresults_all)<-names(AMARETTOdirectories)
  return(list(AMARETTOinit_all = AMARETTOinit_all,
              AMARETTOresults_all = AMARETTOresults_all))
}


#' Title cAMARETTO_HTML_Read
#'
#' @param unzipParentDirectory a directory address
#'  where the html report files are unzipped to
#' @param HTMLsAMARETTOZips a list of directories for each
#'  AMARETTO HTML.zip report, the output of AMARETTO.
#'
#' @return a named vector with directories of unzipped HTML reprorts.
#' @importFrom  utils zip unzip
#' @export
#' 
#' @examples 
#' try(
#' cAMARETTO_HTML_Read(list(TCGA_LIHC="TCGA_LIHC_report.zip",
#' TCGA_GBM="TCGA_GBM_report.zip"))
#' )

cAMARETTO_HTML_Read<-function(HTMLsAMARETTOZips,unzipParentDirectory=getwd()){
  HTMLsAMARETTOlist<-rep(NA,length(HTMLsAMARETTOZips))
  if (! is.null(HTMLsAMARETTOZips) ){
    for (i in 1:length(HTMLsAMARETTOZips)){
      if(file.exists(HTMLsAMARETTOZips[[i]])){
        extdir=sub("/$","",unzipParentDirectory)
        utils::unzip(HTMLsAMARETTOZips[[i]],exdir=extdir)
        print(paste(HTMLsAMARETTOZips[[i]],"is unzipped to",extdir))
        HTMLsAMARETTOlist[i]<-sub(".zip","",paste(extdir,
                                basename(HTMLsAMARETTOZips[[i]]),sep="/"))
      }
      else{
        stop(paste0("The ",HTMLsAMARETTOZips[[i]]," is not existing"))
      }
    }
    names(HTMLsAMARETTOlist)<-names(HTMLsAMARETTOZips)
  }
  print(HTMLsAMARETTOlist)
  return(HTMLsAMARETTOlist)
}


#' @title cAMARETTO_ExportResults
#'
#' @param cAMARETTOresults The output of the Results function.
#' @param cAMARETTOnetworkM The output of the Module Network function.
#' @param cAMARETTOnetworkC The output of the Identify Communities function.
#' @param output_address The address where the zipfile will be generated in. 
#' @importFrom  utils zip unzip
#' @return a zipfile containing cAmaretto results
#' @examples 
#' try(
#' cAMARETTO_ExportResults(cAMARETTOresults,cAMARETTOnetworkM,
#' cAMARETTOnetworkC, output_address="./")
#' )
#' @export
cAMARETTO_ExportResults<-function(cAMARETTOresults,cAMARETTOnetworkM,
                                  cAMARETTOnetworkC, output_address="./"){

  if (!dir.exists(output_address)){
    stop("Output directory is not existing.")
  }
  
  #add a date stamp to the output directory
  output_dir<-paste0("cAMARETTOresults_",gsub("-|:","",
                                              gsub(" ","_",Sys.time())))
  dir.create(file.path(output_address,output_dir))

  #save rdata files for AMARETTO_Run and AMARETTO_Initialize output
  save(cAMARETTOresults, file=file.path(output_address,
                       output_dir,"/cAMARETTOresults.RData"))
  save(cAMARETTOnetworkM, file=file.path(output_address,
                       output_dir,"/cAMARETTOnetworkM.RData"))
  save(cAMARETTOnetworkC, file=file.path(output_address,
                       output_dir,"/cAMARETTOnetworkC.RData"))

  utils::zip(zipfile = file.path(output_address,
                                 output_dir),
             files=file.path(output_address,output_dir))
  unlink(output_dir,recursive = TRUE)
}

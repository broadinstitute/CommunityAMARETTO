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

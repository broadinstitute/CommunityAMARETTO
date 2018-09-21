### Create the DataBase for gene enrichment analysis ###
CreateDataBase <- function(){

gmtFiles_symb <- list.files("PostProcessing/EnrichmentAnalysis",pattern=".gmt");
gmtNames_symb <- strsplit2(gmtFiles_symb,split="\\.")[,1];
MSigDB_lists_symb <- list();
MSigDB_lists_merged_symb <- list();

for(g in 1:length(gmtFiles_symb)){
  File = paste0("PostProcessing/EnrichmentAnalysis/",gmtFiles_symb[g])    
  tmp <- ReadGMT(File);
  MSigDB_lists_symb[[g]] <- tmp$genesets;
  names(MSigDB_lists_symb[[g]]) <- tmp$geneset.names;
  MSigDB_lists_merged_symb <- append(MSigDB_lists_merged_symb,MSigDB_lists_symb[[g]]);  
}

names(MSigDB_lists_symb) <- gmtNames_symb;

#REACTOME pathways are already included in some of the cancer pathway lists. so don't re-included the REACTOME here!
#many of these immun symbol lists are derived from GSE studies - interesting!
GSEA_base_MSigDB_lists_symb <- MSigDB_lists_symb

GSEA_base_MSigDB_lists_merged <- list();

for(e in 1:length(GSEA_base_MSigDB_lists_symb)){
  
  GSEA_base_MSigDB_lists_merged <- append(GSEA_base_MSigDB_lists_merged,GSEA_base_MSigDB_lists_symb[[e]]);
  if(any(duplicated(names(GSEA_base_MSigDB_lists_merged)))){
    cat("loop:",e,"\n")
    stop("error! Some pathway names are duplicated")
  }
}

return(GSEA_base_MSigDB_lists_merged)
}

### Read GMT files, based on GSA.read.gmt function (small changes) ###
ReadGMT <- function(filename){
    a = scan(filename, what = list("", ""), sep = "\t", quote = NULL, fill = T, flush = T, multi.line = F,quiet=TRUE)
    geneset.names = a[1][[1]]
    geneset.descriptions = a[2][[1]]
    dd = scan(filename, what = "", sep = "\t", quote = NULL,quiet=TRUE)
    nn = length(geneset.names)
    n = length(dd)
    ox = rep(NA, nn)
    ii = 1
    for (i in 1:nn) {
        #cat(i)
        while ((dd[ii] != geneset.names[i]) | (dd[ii + 1] != 
                                                   geneset.descriptions[i])) {
            ii = ii + 1
        }
        ox[i] = ii
        ii = ii + 1
    }
    genesets = vector("list", nn)
    for (i in 1:(nn - 1)) {
        #cat(i, fill = T)
        i1 = ox[i] + 2
        i2 = ox[i + 1] - 1
        geneset.descriptions[i] = dd[ox[i] + 1]
        genesets[[i]] = dd[i1:i2]
    }
    geneset.descriptions[nn] = dd[ox[nn] + 1]
    genesets[[nn]] = dd[(ox[nn] + 2):n]
    out = list(genesets = genesets, geneset.names = geneset.names, 
               geneset.descriptions = geneset.descriptions)
    class(out) = "GSA.genesets"
    return(out)
}
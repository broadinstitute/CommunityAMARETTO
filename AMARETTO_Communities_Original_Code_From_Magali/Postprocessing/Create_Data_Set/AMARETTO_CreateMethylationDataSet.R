
################################################################################################################
# Initializers
################################################################################################################

# loading scripts
RscriptsPath="/Users/ogevaert/Documents/WORK/Rscripts/" # this should be the only absolute path
source(paste(RscriptsPath,"TCGA_Pancancer_Scripts.R",sep=""))

# needed for the batch data
PANCANstring="Stanford/TCGApancancer/"
ProjectSTRsPANCAN=TCGA_GENERIC_InitializeProjects(PANCANstring,1)
DataDateString_Combined="data/2013-05-29/"
ResultsDateString_Combined="results/2013-05-29/"

################################################################################################################
# Creating data sets for each of the 12 cancers but using the original clustered data (so raw data)
################################################################################################################

MatchedComboCancers=c("COAD","LUAD","LUSC","READ","KIRC","UCEC","GBM")
for (i in 1:length(MatchedComboCancers)) {
    Cancer=MatchedComboCancers[i]
    cat("Loading data set",Cancer,"with 450k and 27k combined.\n")
    load(paste(ProjectSTRsPANCAN$ProjectRoot,DataDateString_Combined,"MET_",Cancer,"_450kand27k_Clustered.RData",sep=""))
    
    eval(parse(text=paste("ClusterResults=ClusterResults_",Cancer,"_450kand27k",sep="")))
    eval(parse(text=paste("rm(ClusterResults_",Cancer,"_450kand27k)",sep="")))
    MET_tumor=ClusterResults[[1]]
    MET_tumor=TCGA_GENERIC_CleanUpSampleNames(MET_tumor_450kand27k,12)
    MET_normal=ClusterResults[[2]]
    
    # saving in current project/working directory
    save(MET_tumor,MET_tumor, file=paste(Cancer,'_Clustered.rda',sep=''))
    
}
    
################################################################################################################
# Same but now only saving the MethylDriver genes. 
################################################################################################################

# Better to only store the results, the original stuff is way to bulky
AllCancers=c("BLCA_450kon27k","BRCA_450kand27k","COAD_450kand27k","GBM_450kand27k","HNSC_450kon27k","KIRC_450kand27k","LAML_27k","LUAD_450kand27k","LUSC_450kand27k","OV_27k","READ_450kand27k","UCEC_450kand27k","COADREAD_450kand27k")
AllCancers_Clean=c("BLCA","BRCA","COAD","GBM","HNSC","KIRC","LAML","LUAD","LUSC","OV","READ","UCEC","COADREAD")

MethylDrivers=list()
MethylStates=list()
for (i in 1:length(AllCancers)) {
    Cancer=AllCancers[i]
    cat("Loading MethylMix data for",AllCancers_Clean[i],"\n")
    DataFile=paste(ProjectSTRsPANCAN$ProjectRoot,ResultsDateString_Combined,"TCGA_Pancancer_",Cancer,"_BetaMixModel_MethylMixResults.RData",sep="")
    load(DataFile)
    
    # only saving the methylation drivers
    MethylDrivers[[i]]=TCGA_GENERIC_ExtractGeneNames(MethylMixResults$MethylationDrivers,1)        
    currentMethylStates=MethylMixResults$MethylationStates
    
    # saving the methylation states
    Genes=TCGA_GENERIC_ExtractGeneNames(rownames(currentMethylStates),0)
    rownames(currentMethylStates)=Genes
    MethylStates[[i]]=TCGA_GENERIC_MergeData(unique(Genes),currentMethylStates)
}

names(MethylDrivers)=AllCancers_Clean
names(MethylStates)=AllCancers_Clean

# Creating data sets
save(MethylDrivers,file=paste('Pancancer12_MethylDrivers.rda',sep=''))
save(MethylStates,file=paste('Pancancer12_MethylStates.rda',sep=''))




########################################################
# initializers
########################################################

RscriptsPath="/Users/ogevaert/Documents/WORK/Rscripts/" # this should be the only absolute path
source(paste(RscriptsPath,"TCGA_Pancancer_Scripts.R",sep=""))
source('AMARETTO.R')

PANCANstring="Stanford/TCGApancancer/"
ProjectSTRs=TCGA_GENERIC_InitializeProjects(PANCANstring,1)


########################################################
# Load LAML downloaded TCGA data
########################################################

LoadFile=paste(ProjectSTRs$ProjectRoot,'data/2014-12-07/AMARETTOdownload/TCGA_LAML_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData',sep='')
load(LoadFile)
ProcessedDataLAML=ProcessedData
MA_Data_Var=TCGA_GENERIC_GeneFiltering('MAD',ProcessedDataLAML$MA_TCGA,25)
Top1000Genes=rownames(MA_Data_Var)[1:1000]

ProcessedDataLAML$MA_TCGA=ProcessedDataLAML$MA_TCGA[Top1000Genes,]

class(ProcessedDataLAML$MA_TCGA) <- "matrix"

ProcessedDataLAML$MA_TCGA=as.single(ProcessedDataLAML$MA_TCGA)

SaveFile=paste(RscriptsPath,'amarettodevel/package/AMARETTO/data/ProcessedDataLAML_1000.rda',sep='')
save(ProcessedDataLAML, file=SaveFile)
pkgname <- "CommunityAMARETTO"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('CommunityAMARETTO')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("ComRunModGenInfo")
### * ComRunModGenInfo

flush(stderr()); flush(stdout())

### Name: ComRunModGenInfo
### Title: Title ComRunModGenInfo
### Aliases: ComRunModGenInfo

### ** Examples

## Not run: 
##D try(
##D df<-ComRunModGenInfo(cAMARETTOresults,cAMARETTOnetworkM,cAMARETTOnetworkC)
##D )
## End(Not run)



cleanEx()
nameEx("CommunityHyperLink")
### * CommunityHyperLink

flush(stderr()); flush(stdout())

### Name: CommunityHyperLink
### Title: Title CommunityHyperLink
### Aliases: CommunityHyperLink

### ** Examples

try(CommunityHyperLink(Community,Community_key,Community_type))



cleanEx()
nameEx("CommunityModuleTableCreate")
### * CommunityModuleTableCreate

flush(stderr()); flush(stdout())

### Name: CommunityModuleTableCreate
### Title: Title CommunityModuleTableCreate
### Aliases: CommunityModuleTableCreate

### ** Examples

try( 
CommunityModuleTableCreate (cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC, RunInfo)
)



cleanEx()
nameEx("CreateHyperGeoTestAll")
### * CreateHyperGeoTestAll

flush(stderr()); flush(stdout())

### Name: CreateHyperGeoTestAll
### Title: Title CreateHyperGeoTestAll
### Aliases: CreateHyperGeoTestAll

### ** Examples

try(
CreateHyperGeoTestAll(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC, hyper_geo_reference = './h.all.v6.2.symbols.gmt', driverGSEA=TRUE,NrCores=4)
)



cleanEx()
nameEx("CreatePhenotypeTable")
### * CreatePhenotypeTable

flush(stderr()); flush(stdout())

### Name: CreatePhenotypeTable
### Title: Title CreatePhenotypeTable
### Aliases: CreatePhenotypeTable

### ** Examples

try(CreatePhenotypeTable(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC, PhenotypeTables))



cleanEx()
nameEx("DriversSharedTbl")
### * DriversSharedTbl

flush(stderr()); flush(stdout())

### Name: DriversSharedTbl
### Title: Title DriversSharedTbl
### Aliases: DriversSharedTbl

### ** Examples

try(
DriversSharedTbl(cAMARETTOresults, cAMARETTOnetworkM, cAMARETTOnetworkC)
)



cleanEx()
nameEx("ExtractGenesInfo")
### * ExtractGenesInfo

flush(stderr()); flush(stdout())

### Name: ExtractGenesInfo
### Title: Title ExtractGenesInfo
### Aliases: ExtractGenesInfo

### ** Examples

try(ExtractGenesInfo(AMARETTOresults,"TCGA-LIHC"))



cleanEx()
nameEx("Extract_Genes_Modules_All")
### * Extract_Genes_Modules_All

flush(stderr()); flush(stdout())

### Name: Extract_Genes_Modules_All
### Title: Title Extract_Genes_Modules_All
### Aliases: Extract_Genes_Modules_All

### ** Examples

try(Extract_Genes_Modules_All(AMARETTOresults_all,gmt_filelist))



cleanEx()
nameEx("InitialCheckInputs")
### * InitialCheckInputs

flush(stderr()); flush(stdout())

### Name: InitialCheckInputs
### Title: Title InitialCheckInputs
### Aliases: InitialCheckInputs

### ** Examples

try(
InitialCheckInputs(cAMARETTOresults,output_address="./",HTMLsAMARETTOlist,CopyAMARETTOReport=FALSE,hyper_geo_reference)
)



cleanEx()
nameEx("ModuleHyperLink")
### * ModuleHyperLink

flush(stderr()); flush(stdout())

### Name: ModuleHyperLink
### Title: Title ModuleHyperLink
### Aliases: ModuleHyperLink

### ** Examples

try(ModuleHyperLink(Module,Run_Names,AMARETTOres,HTMLsAMARETTOlist,CopyAMARETTOReport,page=1))



cleanEx()
nameEx("RunHyperLink")
### * RunHyperLink

flush(stderr()); flush(stdout())

### Name: RunHyperLink
### Title: Title RunHyperLink
### Aliases: RunHyperLink

### ** Examples

try(RunHyperLink(Run_Names, AMARETTOres, HTMLsAMARETTOlist, CopyAMARETTOReport, page=1))




cleanEx()
nameEx("cAMARETTO_ColorOneModule")
### * cAMARETTO_ColorOneModule

flush(stderr()); flush(stdout())

### Name: cAMARETTO_ColorOneModule
### Title: cAMARETTO_ColorOneModule Creates a network with one (or none)
###   colored modules
### Aliases: cAMARETTO_ColorOneModule

### ** Examples

try(
cAMARETTO_ColorOneModule(cAMARETTOnetworkM, cAMARETTOnetworkC, 2)
)



cleanEx()
nameEx("cAMARETTO_Cytoscape")
### * cAMARETTO_Cytoscape

flush(stderr()); flush(stdout())

### Name: cAMARETTO_Cytoscape
### Title: Title cAMARETTO_Cytoscape
### Aliases: cAMARETTO_Cytoscape

### ** Examples

 
try(
cytoscape_name<-"cAMARETTO_Liver2DS"
cAMARETTOsList<-readRDS(file="./outputs/cAMARETTO_Liver2DS.rds")
communityReportURL<-"http://portals.broadinstitute.org/pochetlab/demo/cAMARETTO_Liver_2DS/" 
cAMARETTO_Cytoscape(cAMARETTOsList,communityReportURL = "",cytoscape_name="my_cytoscape")
)



cleanEx()
nameEx("cAMARETTO_HTML_Read")
### * cAMARETTO_HTML_Read

flush(stderr()); flush(stdout())

### Name: cAMARETTO_HTML_Read
### Title: Title cAMARETTO_HTML_Read
### Aliases: cAMARETTO_HTML_Read

### ** Examples

cAMARETTO_HTML_Read(list(TCGA_LIHC="TCGA_LIHC_report.zip",TCGA_GBM="TCGA_GBM_report.zip"))



cleanEx()
nameEx("cAMARETTO_HTMLreport")
### * cAMARETTO_HTMLreport

flush(stderr()); flush(stdout())

### Name: cAMARETTO_HTMLreport
### Title: cAMARETTO_HTMLreport Creates a HTMLreport for the community
###   AMARETTO results
### Aliases: cAMARETTO_HTMLreport

### ** Examples

## Not run: 
##D try(
##D cAMARETTO_HTMLreport(cAMARETTOresults,cAMARETTOnetworkM, cAMARETTOnetworkC,HTMLsAMARETTOlist = HTMLsAMARETTOlist, hyper_geo_reference = gmtfile, output_address= "./")
##D )
## End(Not run)



cleanEx()
nameEx("cAMARETTO_IdentifyCom")
### * cAMARETTO_IdentifyCom

flush(stderr()); flush(stdout())

### Name: cAMARETTO_IdentifyCom
### Title: cAMARETTO_IdentifyCom
### Aliases: cAMARETTO_IdentifyCom

### ** Examples

try(
cAMARETTOnetworkC<-cAMARETTO_IdentifyCom(cAMARETTOnetworkM,filterComm = FALSE)
)



cleanEx()
nameEx("cAMARETTO_InformationTable")
### * cAMARETTO_InformationTable

flush(stderr()); flush(stdout())

### Name: cAMARETTO_InformationTable
### Title: cAMARETTO_InformationTable
### Aliases: cAMARETTO_InformationTable

### ** Examples

try(
cAMARETTO_InformationTable(cAMARETTOnetworkM, cAMARETTOnetworkC)
)



cleanEx()
nameEx("cAMARETTO_ModuleNetwork")
### * cAMARETTO_ModuleNetwork

flush(stderr()); flush(stdout())

### Name: cAMARETTO_ModuleNetwork
### Title: cAMARETTO_ModuleNetwork Creates a module network.
### Aliases: cAMARETTO_ModuleNetwork

### ** Examples

try(
cAMARETTOnetworkM<-cAMARETTO_ModuleNetwork(cAMARETTOresults,0.10,5)
)



cleanEx()
nameEx("cAMARETTO_Read")
### * cAMARETTO_Read

flush(stderr()); flush(stdout())

### Name: cAMARETTO_Read
### Title: cAMARETTO_Read Reads the AMARETTOinit and AMARETTOresults from
###   (zipped) AMARETTO directories. The names of the list are assigned as
###   run names.
### Aliases: cAMARETTO_Read

### ** Examples

try(
AMARETTOdirectories <- list(LIHC="AMARETTOresults_20181102_142532.zip",BLCA="AMARETTOresults_20181102_142602.zip",GBM="AMARETTOresults_20181102_142636.zip")
AMARETTO_all <- cAMARETTO_Read(AMARETTOdirectories)
AMARETTOinit_all <- AMARETTO_all$AMARETTOinit_all
AMARETTOresults_all <- AMARETTO_all$AMARETTOresults_all
)



cleanEx()
nameEx("cAMARETTO_Results")
### * cAMARETTO_Results

flush(stderr()); flush(stdout())

### Name: cAMARETTO_Results
### Title: cAMARETTO_Results To initiate the community AMARETTO this
###   results functions performs hyper geometric tests between all modules
###   and gene sets.
### Aliases: cAMARETTO_Results

### ** Examples

try(
Cibersortgmt <- "ciberSort.gmt"
cAMARETTOresults <- cAMARETTO_Results(AMARETTOresults_all, gmt_filelist=list(ImmuneSignature = Cibersortgmt), NrCores = 4 , output_dir = "./")
)



cleanEx()
nameEx("cAMARETTO_heatmap")
### * cAMARETTO_heatmap

flush(stderr()); flush(stdout())

### Name: cAMARETTO_heatmap
### Title: cAMARETTO_heatmap Constructs a heatmap with the comparison
###   between modules of two runs based on a hyper geometric test
### Aliases: cAMARETTO_heatmap

### ** Examples

try(
cAMARETTO_heatmap(cAMARETTOresults, "LIHC", "GBM")
)



cleanEx()
nameEx("create_hgt_datatable")
### * create_hgt_datatable

flush(stderr()); flush(stdout())

### Name: create_hgt_datatable
### Title: Title create_hgt_datatable
### Aliases: create_hgt_datatable

### ** Examples

try(create_hgt_datatable(output_hgt, com_table=FALSE, ComNr = 1))



cleanEx()
nameEx("readGMT")
### * readGMT

flush(stderr()); flush(stdout())

### Name: readGMT
### Title: readGMT
### Aliases: readGMT

### ** Examples

readGMT("file.gmt")



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

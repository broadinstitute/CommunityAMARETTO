

# loading scripts
RscriptsPath="/Users/ogevaert/Documents/WORK/Rscripts/" # this should be the only absolute path
source(paste(RscriptsPath,"TCGA_Pancancer_Scripts.R",sep=""))

# needed for the batch data
PANCANstring="Stanford/TCGApancancer/"
ProjectSTRsPANCAN=TCGA_GENERIC_InitializeProjects(PANCANstring,1)

# loading the batch correction data
Filename=paste(ProjectSTRsPANCAN$ProjectRoot,"data/2013-05-06/BA_TCGA_Pancancer_Char15.txt",sep="")

# UPdate March 2015
Filename=paste(ProjectSTRsPANCAN$ProjectRoot,"data/2015-03-15/BA_TCGA_Pancancer_Char15_March2015.txt",sep="")
BatchData=read.table(Filename, header=T, sep='\t',comment.char='')
BatchData[,1]=gsub('\\.','-',BatchData[,1])
BatchData[,2]=gsub('\\.','-',BatchData[,2])

# Update June 2015
Filename=paste(ProjectSTRsPANCAN$ProjectRoot,"data/2015-06-10/BA_TCGA_Pancancer_Char15_June2015.txt",sep="")
BatchData=read.table(Filename, header=T, sep='\t',comment.char='')
BatchData[,1]=gsub('\\.','-',BatchData[,1])
BatchData[,2]=gsub('\\.','-',BatchData[,2])

# saving in current project/working directory
save(BatchData, file='BatchData.rda')
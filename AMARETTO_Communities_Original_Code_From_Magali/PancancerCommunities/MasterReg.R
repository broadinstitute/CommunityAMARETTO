library("RMySQL")
library("stringr")

### correlate GPX2 with smoking for more details
load("/Users/magali/Desktop/recherche/PancancerAnalysis/Communities/Communities.Rdata")
Community = 15 
ModuleCom <- rownames(Communities$NodeInfo[which(Communities$NodeInfo$community==Community),])  

# Select the modules with GPX2 as regulators
ModulesGPX <- c()
for (i in (1:length(ModuleCom))){
  CancerSite <- strsplit2(ModuleCom[i],"_")[,2]
  ModuleNr <- strsplit2(ModuleCom[i],"_")[,1]
  ModuleNr <- as.numeric(str_replace_all(ModuleNr,"M",""))
  load(paste0("/Users/magali/Desktop/recherche/PancancerAnalysis/AMARETTO/Results/results",CancerSite,".RData"))
  Reg <- colnames(AMARETTOresults$v)[which(abs(AMARETTOresults$v[ModuleNr,])>0)]
  if ("GPX2" %in% Reg){
    ModulesGPX <- c(ModulesGPX,ModuleCom[i])
  }
}

# load the clinical data
load("/Users/mchampion/Desktop/amarettodevel/Magali/SurvivalAnalysis/TCGA_barcode_index2.RData")
CodeCancer <- TCGA_barcode_index2
CodeCancer <- as.matrix(CodeCancer)
row.names(CodeCancer) <- CodeCancer[,1]

conn<- dbConnect(MySQL(), user="ogevaert", password="",host="127.0.0.1", dbname="StanfordData",port=3307)     
RawClinicalObs.All = dbGetQuery(conn, "SELECT * FROM TCGA_patient_03_09_2015 ;")     

rownames(RawClinicalObs.All) <- RawClinicalObs.All[,"row_names"]
RawClinicalObs.All <- RawClinicalObs.All[,-c(which(colnames(RawClinicalObs.All)=="row_names"),which(colnames(RawClinicalObs.All)=="bcr_patient_barcode"))]
OverlapSamples=intersect(CodeCancer[,1],rownames(RawClinicalObs.All))
RawClinicalObs.All <- RawClinicalObs.All[OverlapSamples,]
CodeCancer <- CodeCancer[OverlapSamples,]  
ClinicalVariables <- colnames(RawClinicalObs.All)

I <- grep(pattern="smok",ClinicalVariables)
Smoking <- ClinicalVariables[I]
Smoking<-c("tobacco_smoking_history","tobacco_smoking_age_started","stopped_smoking_year","year_of_tobacco_smoking_onset","number_pack_years_smoked") # after removing variables with only missing values
Smokingdata <- RawClinicalObs.All[,Smoking]
colnames(Smokingdata) <- c("Smoking_profile","Smoking_started_age","Smoking_stopped_year","Smoking_started_year","number_pack_years_smoked")

Impossible1<-intersect(which(Smokingdata$Smoking_stopped_year>0),which(Smokingdata$Smoking_profile=="Current smoker"))
Impossible2 <- intersect(which(Smokingdata$Smoking_stopped_year>0),which(Smokingdata$Smoking_profile=="Lifelong Non-smoker"))
Impossible3 <- intersect(which(Smokingdata$Smoking_started_year>0),which(Smokingdata$Smoking_profile=="Lifelong Non-smoker"))
Impossible4 <- intersect(which(Smokingdata$Smoking_started_age>0),which(Smokingdata$Smoking_profile=="Lifelong Non-smoker"))
Impossible5 <- intersect(which(Smokingdata$number_pack_years_smoked>0),which(Smokingdata$Smoking_profile=="Lifelong Non-smoker"))
Smokingdata[c(Impossible1,Impossible2,Impossible3,Impossible4,Impossible5),] <- matrix(NA,nrow=length(c(Impossible1,Impossible2,Impossible3,Impossible4,Impossible5)),ncol=5)
Smokingdata <- Smokingdata[-which(Smokingdata$Smoking_profile=="Current Reformed Smoker, Duration Not Specified"),]

SmokingprofileNum <- Smokingdata[,1] 
SmokingprofileNum[which(SmokingprofileNum=="Lifelong Non-smoker")] <- 0
SmokingprofileNum[which(SmokingprofileNum=="Current reformed smoker for < or = 15 years")] <- 2
SmokingprofileNum[which(SmokingprofileNum=="Current Reformed Smoker, Duration Not Specified")] <- NA
SmokingprofileNum[which(SmokingprofileNum=="Current reformed smoker for > 15 years")] <- 1
SmokingprofileNum[which(SmokingprofileNum=="Current smoker")] <- 3

NApatients <- which(rowSums(is.na(Smokingdata))==ncol(Smokingdata))
Smokingdata <- Smokingdata[-NApatients,]
# remove started and stopped, replace by number of years
NumberYears = as.numeric(Smokingdata$Smoking_stopped_year)-as.numeric(Smokingdata$Smoking_started_year)
Smokingdata <- cbind(Smokingdata,NumberYears)
SmokingdataNum <- cbind(Smoking_profile=SmokingprofileNum[-NApatients],Smokingdata[,-1])

# now correlation
#for (i in (1:length(ModulesGPX))){
#  CancerSite <- strsplit2(ModulesGPX[i],"_")[2]
#  ModuleNr <- strsplit2(ModulesGPX[i],"_")[1]
#  ModuleNr <- as.numeric(str_replace_all(ModuleNr,"M",""))
#  load(paste0("/Users/mchampion/Desktop/PancancerAnalysis/AMARETTO/results",CancerSite,".RData"))
#  load(paste0("/Users/mchampion/Desktop/PancancerAnalysis/AMARETTO/initstep",CancerSite,".RData"))
#  ModuleData <- AMARETTOinit$MA_matrix_Var[AMARETTOresults$GeneNames[which(AMARETTOresults$assign==ModuleNr)],]
#  colnames(ModuleData) <- str_replace_all(colnames(ModuleData),"-01","")
 
Cancers <- c("BLCA","HNSC","LUAD","LUSC") #,"UCEC") no smoking data for UCEC
Correlation <- matrix(0,ncol=2*ncol(Smokingdata),nrow=length(Cancers))
rownames(Correlation) <- Cancers
colnames(Correlation) <- rep(colnames(Smokingdata),each=2)

for (i in (1:length(Cancers))){  
  CancerSite = Cancers[i]
  load(paste0("/Users/mchampion/Desktop/PancancerAnalysis/AMARETTO/TCGA_",CancerSite,"_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData"))
  MA_matrix <- ProcessedData$MA_TCGA
  if (length(which(rownames(MA_matrix)=="GPX2"))>0){ # GPX2 present?
    GeneExpression <- MA_matrix["GPX2",]
    names(GeneExpression) <- str_replace_all(colnames(MA_matrix),"-01","")
    GeneExpression <- GeneExpression[intersect(names(GeneExpression),rownames(Smokingdata))]
  
    SmokingdataModule <- Smokingdata[intersect(names(GeneExpression),rownames(Smokingdata)),]
    I <- which(colSums(is.na(SmokingdataModule))==nrow(SmokingdataModule))
    if (length(I)>0){
      SmokingdataModule <- SmokingdataModule[,-I]
    }
    
    # j=1, smoking profile: qualitative
    test<- kruskal.test(SmokingdataModule[,1]~GeneExpression)
    Correlation[i,1:2] <- c(NA,test$p.value)
    
    for (j in (3:(ncol(Correlation)/2))){
      Var <- colnames(Correlation)[(j-1)*2+1]
      I <- which(colnames(SmokingdataModule)==Var)
      if (length(I)>0){
        correlation=cor(GeneExpression,as.numeric(as.character(SmokingdataModule[,I])),use="na.or.complete")
        test <- cor.test(GeneExpression,as.numeric(as.character(SmokingdataModule[,I])),method="pearson",exact=FALSE)
        Correlation[i,(2*(j-1)+1):(2*j)] <- c(correlation,test$p.value)
      }
    }
    
    CombinedData <- cbind(GeneExpression,SmokingdataModule[,1])
    CombinedData <- data.frame(CombinedData)
    CombinedData$GeneExpression <- as.numeric(as.character(CombinedData$GeneExpression))
    colnames(CombinedData) <- c("GPX2","SmokingProfile")
    p <- ggplot(CombinedData,aes(SmokingProfile,GPX2),main="test")
    print(p + geom_boxplot(aes(fill = SmokingProfile)) +ggtitle(paste0("Correlation of GPX2 expression for ",CancerSite," cancer with smoking profile"))+xlim("Lifelong Non-smoker","Current reformed smoker for > 15 years","Current reformed smoker for < or = 15 years","Current smoker")+labs(fill="Smoking Profile",x="",y="GPX2 expression"))
  
    CombinedData$SmokingProfile<-factor(CombinedData$SmokingProfile, levels = c("Current smoker","Current reformed smoker for < or = 15 years","Current reformed smoker for > 15 years","Lifelong Non-smoker"))
    mylogit<- glm(CombinedData$GPX2~CombinedData$SmokingProfile)
    summary(mylogit)
    wald.test(b = coef(mylogit), Sigma = vcov(mylogit), Terms = 2:4)
    
    library("coin")
    MolecularVarBinary <- CombinedData[-which(CombinedData$SmokingProfile %in% c("Current reformed smoker for < or = 15 years","Current reformed smoker for > 15 years","Current Reformed Smoker, Duration Not Specified")),]
    test<- wilcox_test(MolecularVarBinary$GPX2~MolecularVarBinary$SmokingProfile)
    Correlation[i,1] <- pvalue(test)
  }
}

# more details for HNSC
CancerSite="HNSC"
load(paste0("/Users/mchampion/Desktop/PancancerAnalysis/AMARETTO/initstep",CancerSite,".RData"))
GeneExpression <- AMARETTOinit$MA_matrix_Var["GPX2",]
names(GeneExpression) <- str_replace_all(colnames(AMARETTOinit$MA_matrix_Var),"-01","")
GeneExpression <- GeneExpression[intersect(names(GeneExpression),rownames(Smokingdata))]

SmokingdataModule <- Smokingdata[intersect(names(GeneExpression),rownames(Smokingdata)),]
I <- which(colSums(is.na(SmokingdataModule))==nrow(SmokingdataModule))
if (length(I)>0){
  SmokingdataModule <- SmokingdataModule[,-I]
}
library(ggplot2)
CombinedData <- cbind(GeneExpression,SmokingdataModule[,1])
CombinedData <- data.frame(CombinedData)
CombinedData$GeneExpression <- as.numeric(CombinedData$GeneExpression)
colnames(CombinedData) <- c("GPX2","SmokingProfile")
p <- ggplot(CombinedData,aes(SmokingProfile,GPX2),main="test")
print(p + geom_boxplot(aes(fill = SmokingProfile)) +ggtitle(paste0("Correlation of GPX2 expression for HNSC cancer with smoking profile"))+xlim("Lifelong Non-smoker","Current reformed smoker for > 15 years","Current reformed smoker for < or = 15 years","Current smoker"))


#### correlate OAS2 with immune response for more details
load("/Users/mchampion/Desktop/PancancerAnalysis/Communities/Communities.Rdata")
Community = 13 
ModuleCom <- rownames(Communities$NodeInfo[which(Communities$NodeInfo$community==Community),])  

# Select the modules with OAS2 as regulators
ModulesOAS2 <- c()
for (i in (1:length(ModuleCom))){
  CancerSite <- strsplit2(ModuleCom[i],"_")[2]
  ModuleNr <- strsplit2(ModuleCom[i],"_")[1]
  ModuleNr <- as.numeric(str_replace_all(ModuleNr,"M",""))
  load(paste0("/Users/mchampion/Desktop/PancancerAnalysis/AMARETTO/results",CancerSite,".RData"))
  Reg <- colnames(AMARETTOresults$v)[which(abs(AMARETTOresults$v[ModuleNr,])>0)]
  if ("OAS2" %in% Reg){
    ModulesOAS2 <- c(ModulesOAS2,ModuleCom[i])
  }
}

# load the clinical data
load("/Users/mchampion/Desktop/amarettodevel/Magali/SurvivalAnalysis/TCGA_barcode_index2.RData")
CodeCancer <- TCGA_barcode_index2
CodeCancer <- as.matrix(CodeCancer)
row.names(CodeCancer) <- CodeCancer[,1]

conn<- dbConnect(MySQL(), user="ogevaert", password="",host="127.0.0.1", dbname="StanfordData",port=3307)     
RawClinicalObs.All = dbGetQuery(conn, "SELECT * FROM TCGA_patient_03_09_2015 ;")     

rownames(RawClinicalObs.All) <- RawClinicalObs.All[,"row_names"]
RawClinicalObs.All <- RawClinicalObs.All[,-c(which(colnames(RawClinicalObs.All)=="row_names"),which(colnames(RawClinicalObs.All)=="bcr_patient_barcode"))]
OverlapSamples=intersect(CodeCancer[,1],rownames(RawClinicalObs.All))
RawClinicalObs.All <- RawClinicalObs.All[OverlapSamples,]
CodeCancer <- CodeCancer[OverlapSamples,]  
ClinicalVariables <- colnames(RawClinicalObs.All)

Immune <- c("lymphocytes_count","necrosis_percent")
Immunedata <- RawClinicalObs.All[,Immune]

NApatients <- which(rowSums(is.na(Immunedata))==ncol(Immunedata))
Immunedata <- Immunedata[-NApatients,]

# now correlation
Cancers <- c("BLCA","BRCA","GBM","HNSC","LUAD","LAML","LUSC","OV","UCEC") # immune data only for LAML
#Correlation <- matrix(0,ncol=2*ncol(Immunedata),nrow=length(Cancers))
#rownames(Correlation) <- Cancers
#colnames(Correlation) <- rep(colnames(Immunedata),each=2)

#for (i in (1:length(Cancers))){  
#  CancerSite = Cancers[i]
CancerSite = "LAML"
  load(paste0("/Users/mchampion/Desktop/PancancerAnalysis/AMARETTO/TCGA_",CancerSite,"_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData"))
  GeneExpression <- ProcessedData$MA_TCGA["OAS2",]
  names(GeneExpression) <- str_replace_all(colnames(ProcessedData$MA_TCGA),"-03","")
  GeneExpression <- GeneExpression[intersect(names(GeneExpression),rownames(Immunedata))]
  
  ImmunedataModule <- Immunedata[intersect(names(GeneExpression),rownames(Immunedata)),]
  I <- which(colSums(is.na(ImmunedataModule))==nrow(ImmunedataModule))
  if (length(I)>0){
    ImmunedataModule <- ImmunedataModule[,-I] # only data for lymphocytes_count
  }
  ImmunedataModule <- as.numeric(ImmunedataModule)
  
  #for (j in (1:length(colnames(SmokingdataModule)))){
    # plot(Gene,NstageNum,xlab="gene expression",ylab="N-stage")
    correlation=cor(GeneExpression,ImmunedataModule,use="na.or.complete")
    test <- cor.test(GeneExpression,ImmunedataModule,method="pearson",exact=FALSE)
    #Correlation[i,(2*(j-1)+1):(2*j)] <- c(correlation,test$p.value)
#  }
#}

# download data from gdac
source("/Users/mchampion/Desktop/amarettodevel/package/TCGAdownload_DevelVersion.R")
Cancers <- c("BLCA","BRCA","GBM","HNSC","LUAD","LUSC","OV","UCEC") # immune data only for LAML
Correlation <- matrix(0,ncol=2*11,nrow=length(Cancers))
rownames(Correlation) <- Cancers
colnames(Correlation) <- rep(c("lymphocyte_infiltration","monocyte_infiltration","neutrophil_infiltration","eosinophil_infiltration","granulocyte_infiltration","inflam_infiltration","necrosis","normal_cells","tumor_cells","stromal_cells","tumor_nuclei"),each=2)
for (i in (1:length(Cancers))){
  CancerSite <- Cancers[i]
  ClinData <- t(read.csv(paste0("/Users/mchampion/Desktop/PancancerAnalysis/Data/gdac_20151101/gdac.broadinstitute.org_",CancerSite,".Merge_Clinical.Level_1.2015110100.0.0/",CancerSite,".clin.merged.txt"),sep="\t",row.names=1,head=TRUE,na.strings=c("NA","null")))
  I <- grep("samples.sample.portions.portion.slides.slide.percent",colnames(ClinData))
  Variables <- colnames(ClinData)[I]
  Variables <- c("patient.bcr_patient_barcode",Variables)
  ClinData <- ClinData[,Variables]
  rownames(ClinData) <- toupper(as.character(ClinData[,1]))
  ClinData <- ClinData[,-1]
  ClinData <- data.frame(ClinData)
  
  NApatients <- which(rowSums(is.na(ClinData))==ncol(ClinData))
  if (length(NApatients)>0){
    ClinData <- ClinData[-NApatients,]
  }
  NAdata <- which(colSums(is.na(ClinData))==nrow(ClinData))
  if (length(NAdata)>0){
    ClinData <- ClinData[,-NAdata]
  }
  colnames(ClinData) <- strsplit2(colnames(ClinData),"percent_")[,2]
  
  load(paste0("/Users/mchampion/Desktop/PancancerAnalysis/AMARETTO/TCGA_",CancerSite,"_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData"))
  GeneExpression <- ProcessedData$MA_TCGA["OAS2",]
  names(GeneExpression) <- str_replace_all(colnames(ProcessedData$MA_TCGA),"-01","")
  if (CancerSite=="LAML"){
    names(GeneExpression) <- str_replace_all(colnames(ProcessedData$MA_TCGA),"-03","")
  }
  GeneExpression <- GeneExpression[intersect(names(GeneExpression),rownames(ClinData))]
  ClinData <- ClinData[intersect(names(GeneExpression),rownames(ClinData)),]
  
  for (j in (1:ncol(ClinData))){
    Var <- colnames(ClinData)[j]
    I <- which(colnames(Correlation)==Var)[1]
    correlation=cor(GeneExpression,as.numeric(as.character(ClinData[,j])),use="na.or.complete")
    test <- cor.test(GeneExpression,as.numeric(as.character(ClinData[,j])),method="pearson",exact=FALSE)
    Correlation[i,(I:(I+1))] <- c(correlation,test$p.value)
  }
}

# each cancer site
CancerSite="GBM"
load(paste0("/Users/mchampion/Desktop/PancancerAnalysis/AMARETTO/TCGA_",CancerSite,"_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData"))
GeneExpression <- ProcessedData$MA_TCGA["OAS2",]
names(GeneExpression) <- str_replace_all(colnames(ProcessedData$MA_TCGA),"-01","")
ClinData <- t(read.csv(paste0("/Users/mchampion/Desktop/PancancerAnalysis/Data/gdac_20151101/gdac.broadinstitute.org_",CancerSite,".Merge_Clinical.Level_1.2015110100.0.0/",CancerSite,".clin.merged.txt"),sep="\t",row.names=1,head=TRUE,na.strings=c("NA","null")))
Variables <- c("patient.bcr_patient_barcode","patient.samples.sample.portions.portion.slides.slide.percent_lymphocyte_infiltration")
ClinData <- ClinData[,Variables]
rownames(ClinData) <- toupper(as.character(ClinData[,1]))
ClinData <- ClinData[,-1]
NApatients <- which(is.na(ClinData)==TRUE)
if (length(NApatients)>0){
  ClinData <-ClinData[-NApatients]
}
ClinData <- ClinData[intersect(names(GeneExpression),names(ClinData))]
GeneExpression <- GeneExpression[intersect(names(GeneExpression),names(ClinData))]
tmp=50
LowLymphocyteInfiltration <- which(as.numeric(ClinData)<=tmp)
HighLymphocyteInfiltration <- which(as.numeric(ClinData)>tmp)
LymphocyteInfiltration <- rep(0,length(ClinData))
names(LymphocyteInfiltration) <- names(ClinData)
LymphocyteInfiltration[LowLymphocyteInfiltration] <- "Low"
LymphocyteInfiltration[HighLymphocyteInfiltration] <- "High"
library(ggplot2)
CombinedData <- cbind(GeneExpression,LymphocyteInfiltration)
CombinedData <- data.frame(CombinedData)
CombinedData$GeneExpression <- as.numeric(as.character(CombinedData$GeneExpression))
colnames(CombinedData) <- c("OAS2","LymphocyteInfiltration")
p <- ggplot(CombinedData,aes(LymphocyteInfiltration,OAS2),main="test")
print(p + geom_boxplot(aes(fill = LymphocyteInfiltration)) +ggtitle(paste0("Correlation of OAS2 expression for ",CancerSite," cancer with lymphocyte infiltration"))+xlim("Low","High"))

# same for BLCA
CancerSite <- "BLCA"
tmp <- 5 # need to be define 

# lymphocyte signature
load("/Users/mchampion/Desktop/amarettodevel/Magali/SurvivalAnalysis/TCGA_barcode_index2.RData")
CodeCancer <- TCGA_barcode_index2
CodeCancer <- as.matrix(CodeCancer)
row.names(CodeCancer) <- CodeCancer[,1]
CodeCancer[,2] <- str_replace_all(CodeCancer[,2],"COAD","COADREAD")
CodeCancer[,2] <- str_replace_all(CodeCancer[,2],"READ","COADREAD")
CodeCancer[,2] <- str_replace_all(CodeCancer[,2],"COADCOADREAD","COADREAD")
conn<- dbConnect(MySQL(), user="ogevaert", password="",host="127.0.0.1", dbname="StanfordData",port=3307)     
RawClinicalObs.All = dbGetQuery(conn, "SELECT * FROM TCGA_patient_03_09_2015 ;")     
rownames(RawClinicalObs.All) <- RawClinicalObs.All[,"row_names"]
RawClinicalObs.All <- RawClinicalObs.All[,-c(which(colnames(RawClinicalObs.All)=="row_names"),which(colnames(RawClinicalObs.All)=="bcr_patient_barcode"))]

Correlations <- matrix(0,ncol=5,nrow=length(AllCancers))
rownames(Correlations) <- AllCancers
colnames(Correlations) <- c("corr","p.value","corr_bin","p.value_bin","p.value2")
for (i in (1:length(AllCancers))){
  CancerSite <- AllCancers[i]
  load(paste0("/Users/mchampion/Desktop/PancancerAnalysis/Data/TCGA_",CancerSite,"_ProcessedData_MA_CNV.RData"))
  
  # Keep samples associated to this cancer (based on the barcode)
  I <- which(CodeCancer[,2]==CancerSite)
  Samples <- CodeCancer[I,1]
  OverlapSamples=intersect(Samples,rownames(RawClinicalObs.All))
  RawClinicalObs.All_cancer <- RawClinicalObs.All[OverlapSamples,]
  ClinicalVariables <- colnames(RawClinicalObs.All_cancer)
  
  # clean N-stages
  Nstage <- RawClinicalObs.All_cancer[,which(ClinicalVariables=="pathologic_n")]
  names(Nstage) <- rownames(RawClinicalObs.All_cancer)
  Nstage[which(Nstage=="NX")] <- NA
  I <- which(is.na(Nstage)==TRUE) # missing data
  if (length(I)<length(Nstage)){
  Nstage <- Nstage[-I]
  Nstage <- str_replace_all(Nstage,"[abc]","") 
  names(Nstage) <- rownames(RawClinicalObs.All_cancer)[-I]
  Nstage <- strsplit2(Nstage," ")[,1]
  Nstage <- str_replace_all(Nstage,"[abcmiol]","")
  names(Nstage) <- rownames(RawClinicalObs.All_cancer)[-I]
  NstageNum <- as.numeric(str_replace_all(Nstage,"N",""))
  names(NstageNum) <- rownames(RawClinicalObs.All_cancer)[-I]
  
  Nstage_bin <- NstageNum
  Nstage_bin[which(NstageNum>0)] <- rep("N+",length(which(NstageNum>0)))
  Nstage_bin[which(NstageNum==0)] <- rep("N0",length(which(NstageNum==0)))  
  Nstage_binNum <- Nstage_bin
  Nstage_binNum[which(Nstage_binNum=="N+")] <- rep(1,length(which(Nstage_binNum=="N+")))  
  Nstage_binNum[which(Nstage_binNum=="N0")] <- rep(0,length(which(Nstage_binNum=="N0")))  
  Nstage_binNum <- as.numeric(str_replace_all(Nstage_binNum,"N","")) 
  names(Nstage_binNum) <- rownames(RawClinicalObs.All_cancer)[-I]
  
  # now compute the correlations
  MAData <- ProcessedData$MA_TCGA
  colnames(MAData) <- str_replace_all(colnames(MAData),"-01","")
  colnames(MAData) <- str_replace_all(colnames(MAData),"-03","")
  OverlapSamples <- intersect(names(Nstage),colnames(MAData))
  Nstage <- Nstage[OverlapSamples]
  NstageNum <- NstageNum[OverlapSamples]
  Nstage_bin <- Nstage_bin[OverlapSamples]
  Nstage_binNum <- Nstage_binNum[OverlapSamples]
  MAData <- MAData[,OverlapSamples]
  
  Expression <- MAData["OAS2",]
  #plot(Expression,NstageNum,xlab="gene expression",ylab="N-stage")
  correlation=cor(Expression,NstageNum)
  test <- cor.test(Expression,NstageNum,method="pearson",exact=FALSE)
  Correlations[i,1:2] <- c(test$estimate,test$p.value)
  #plot(Expression,Nstage_binNum,xlab="gene expression",ylab="N-stage")
  correlation=cor(Expression,Nstage_binNum)
  test <- cor.test(Expression,Nstage_binNum,method="pearson",exact=FALSE)
  Correlations[i,3:4] <- c(test$estimate,test$p.value)
  
  # 3rd step plot the correlations using boxplots
  library(ggplot2)
  CombinedData <- cbind(Expression,Nstage)
  CombinedData <- data.frame(CombinedData)
    
  CombinedData$Expression <- as.numeric(CombinedData$Expression)
  colnames(CombinedData) <- c("OAS2","N_stage")
  p <- ggplot(CombinedData,aes(N_stage,OAS2),main="test")
  print(p + geom_boxplot(aes(fill = N_stage)) +ggtitle(paste0("Correlations of OAS2 with N_stage for ",CancerSite)))
  
  CombinedData <- cbind(Expression,Nstage_bin)
  CombinedData <- data.frame(CombinedData)
  CombinedData$Expression <- as.numeric(CombinedData$Expression)
  colnames(CombinedData) <- c("OAS2","N_stage")
  p <- ggplot(CombinedData,aes(N_stage,OAS2),main="test")
  print(p + geom_boxplot(aes(fill = N_stage)) +ggtitle(paste0("Correlations of OAS2 with N_stage for ",CancerSite)))

  library("coin")
  test<- wilcox_test(CombinedData$OAS2~CombinedData$N_stage)
  Correlations[i,5] <- pvalue(test)
  }
}

## correlate with gene signature
library("doBy")
CellData <- read.csv("/Users/magali/Desktop/recherche/PancancerAnalysis/TcellSignature/Tcelldata.csv",header=TRUE)
colnames(CellData) <- c("Gene","p.value","q.value","Ratio","Fold.change",'Attributes')
CellData[,"Attributes"] <- recodeVar(CellData[,"Attributes"], c("cold up vs hot", "cold down vs hot"), c("Up", "Down"))
CellData <- CellData[-which(CellData[,"Gene"]=="OAS2"),]
Results <- matrix(0,nrow=length(AllCancers),ncol=3)
rownames(Results) <- AllCancers
colnames(Results) <- c("p.value","correlation.estimate","sign")
for (i in (1:length(AllCancers))){
  CancerSite <- AllCancers[i]
  load(paste0("/Users/magali/Desktop/recherche/PancancerAnalysis/Data/TCGA_",CancerSite,"_ProcessedData_MA_CNV.Rdata"))
  MA_TCGA <- ProcessedData$MA_TCGA
  
  #### clean the data (overlap, scale, threshold) ####  
  OverlapGenes <- intersect(as.character(CellData[,"Gene"]),rownames(MA_TCGA))
  MA_TCGACold <- MA_TCGA[OverlapGenes,]
  CellData <- CellData[which(as.character(CellData[,"Gene"])%in%OverlapGenes),]      
  MA_TCGACold <- t(scale(t(MA_TCGACold))) 
  Upgenes <- as.character(CellData[which(CellData[,"Attributes"]=="Up"),"Gene"])
  Downgenes <- as.character(CellData[which(CellData[,"Attributes"]=="Down"),"Gene"])
  score <- colSums(MA_TCGACold[Upgenes,])-colSums(MA_TCGACold[Downgenes,])
  
  GeneExp <- MA_TCGA[which(rownames(MA_TCGA)=="OAS2"),which(colnames(MA_TCGA)%in% names(score))]
  GeneExp <- scale(GeneExp)
  names(GeneExp) <- names(GeneExp)
  plot(GeneExp,score,xlab="OAS2",ylab="score signature",main=paste0("Correlation between score signature and gene OAS2 for ",CancerSite))
  correlation=cor(GeneExp,score)
  reg1 <- lm(score~GeneExp)
  par(cex=.8)
  abline(reg1)
  
  test <- cor.test(GeneExp,score,method="spearman",exact=FALSE)
  Results[i,] <- c(p.value=test$p.value,correlation.estimate=test$estimate,sign=sign(test$estimate))
}

#### Kevin code
AllCancers <- c("BLCA","BRCA","COADREAD","GBM",'HNSC',"KIRC","LAML","LUAD","LUSC","OV","UCEC")
library(plyr)
library("doBy")
library(ggplot2)
library(coin)

#Get T cell marker gene data
#I used the T cell signature data that you have previously used, from the Spranger et al melanoma paper
CellData=read.csv("/Users/magali/Desktop/recherche/PancancerAnalysis/TcellSignature/TCellData.csv",header=TRUE)
colnames(CellData)=c("Gene","p.value","q.value","Ratio","Fold.change",'Attributes')
CellData[,"Attributes"] <- recodeVar(CellData[,"Attributes"], c("cold up vs hot", "cold down vs hot"), c("Up", "Down"))
CellDataTCell=CellData[CellData$Attributes=="Down","Gene"]
#These are the genes upregulated in T cell, i.e. T cell marker genes
CellDataNonTCell=CellData[CellData$Attributes=="Up","Gene"]
#These are the genes downregulated in T cell

correlation <- matrix(0,ncol=6,nrow=length(AllCancers))
colnames(correlation) = c("Corr","p-value","p-valueBin","signBin","CorrSign","p-valueSign")
rownames(correlation) <- AllCancers
for (i in (1:length(AllCancers))){
  CancerSite <- AllCancers[i]
  load(paste0("/Users/magali/Desktop/recherche/PancancerAnalysis/AMARETTO/Data/TCGA_",CancerSite,"_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData"))
  MAData <- ProcessedData$MA_TCGA
  #Subtract the mean for each gene from each sample
  MAcancerNorm=sweep(MAData, 1, rowMeans(MAData), "-")
  
  #Select the genes that are at least 2 fold above mean
  OverExpressedGene=list()
  #For each patient, select the genes that are >2, i.e. 2 fold higher than mean expression for that gene
  for(j in 1:ncol(MAcancerNorm)){
    OverExpressedGene[[j]]=MAcancerNorm[MAcancerNorm[,j]>=2,j]
  }
  names(OverExpressedGene)=colnames(MAcancerNorm)
  
  UnderExpressedGene=list()
  #For each patient, select the genes that are >2, i.e. 2 fold higher than mean expression for that gene
  for(j in 1:ncol(MAcancerNorm)){
    UnderExpressedGene[[j]]=MAcancerNorm[MAcancerNorm[,j]<=-2,j]
  }
  names(UnderExpressedGene)=colnames(MAcancerNorm)
  
  #Next, I do a hypergeometric test for each patient, to test for a signficiant overlap between the overexpressed genes for each patient, and the T cell marker genes 
  phyper_TILup_overexpressed=array(NA, c(length(OverExpressedGene),2))
  rownames(phyper_TILup_overexpressed)=names(OverExpressedGene)
  colnames(phyper_TILup_overexpressed)=c("P.value","Percentage_Overlap")
  
  for(j in 1:length(OverExpressedGene)){
    phyper_TILup_overexpressed[j,1]=phyper(length(intersect(names(OverExpressedGene[[j]]),CellDataTCell))-1,length(intersect(CellDataTCell, rownames(MAData))),length(setdiff(rownames(MAData),CellDataTCell)),length(names(OverExpressedGene[[j]])),lower.tail=F, log.p = FALSE)
    phyper_TILup_overexpressed[j,2]=length(intersect(names(OverExpressedGene[[j]]),CellDataTCell))/length(names(OverExpressedGene[[j]]))
  }
  phyper_TILup_overexpressed=as.data.frame(phyper_TILup_overexpressed)
  #FDR correct the p-values. I'm not sure this is strictly necessary
  phyper_TILup_overexpressed$Q.value=p.adjust(phyper_TILup_overexpressed$P.value, method="fdr")
  #Lastly, I just make a handy 1/2 categorical variable indicating whether each patient had a significant enrichment for T cell marker genes
  phyper_TILup_overexpressed$TILup_marker_gene_enrichment=rep(NA,nrow(phyper_TILup_overexpressed))
  phyper_TILup_overexpressed$TILup_marker_gene_enrichment[phyper_TILup_overexpressed$P.value>0.05]=1
  phyper_TILup_overexpressed$TILup_marker_gene_enrichment[phyper_TILup_overexpressed$P.value<0.05]=2
  
  # first linear correlations
  reg1 <- lm(phyper_TILup_overexpressed$Percentage_Overlap~MAData["OAS2",])
  par(cex=.8)
  plot(MAData["OAS2",], phyper_TILup_overexpressed$Percentage_Overlap, xlab="SMARCB1 expression",ylab="Gene signature",main=paste0("Correlation for ",CancerSite))
  abline(reg1)
  
  test <- cor.test(MAData["OAS2",],phyper_TILup_overexpressed$Percentage_Overlap,method="pearson",exact=FALSE)
  correlation[i,1:2] <- c(test$estimate,test$p.value)

  # 2nd correlation boxplots
  Combined <- cbind(OAS2=MAData["OAS2",],Tcell=phyper_TILup_overexpressed$TILup_marker_gene_enrichment)
  Combined[which(Combined[,2]==2),2] <- rep("Tcell",length(which(Combined[,2]==2)))
  Combined[which(Combined[,2]==1),2] <- rep("NoTcell",length(which(Combined[,2]==1)))
  Combined <- data.frame(Combined)
  Combined$OAS2 <- as.numeric(as.character(Combined$OAS2))
  p <- ggplot(Combined,aes(Tcell,OAS2),main="test")
  print(p + geom_boxplot(aes(fill = Tcell)) +ggtitle(paste0("Correlations of Tcell with OAS2 for cancer ",CancerSite))+
          labs(x = "T-cell status",y="OAS2 expression")+
        scale_x_discrete(breaks=c("NoTcell", "Tcell"),
                         labels=c("Non T-cell inflamed", "T-cell inflamed"))+
          scale_fill_discrete(name="T-cell status",
                                breaks=c("NoTcell", "Tcell"),
                                labels=c("Non T-cell inflamed", "T-cell inflamed")))
  
  test<- wilcox_test(Combined$OAS2~Combined$Tcell)
  correlation[i,3:4] <- c(pvalue(test),(statistic(test)>0))
  
  # another gene list
  TCellTranscripts=c("CD8A", "CCL2", "CCL3",
                     "CCL4", "CXCL9", "CXCL10", "ICOS", "GZMK", "IRF1", "HLA-DMA", "HLA-DMB", "HLA-DOA",
                     "HLA-DOB")
  Exp <- ProcessedData$MA_TCGA[which(rownames(MAData)%in%TCellTranscripts),]
  score <- colMeans(Exp)
  plot(score,MAData["OAS2",],xlab="Score signature",ylab="OAS2 expression",main=paste0("Correlation between OAS2 and Tcell signature for ",CancerSite))
  reg1 <- lm(MAData["OAS2",]~score)
  par(cex=.8)
  abline(reg1)
  
  test <- cor.test(MAData["OAS2",],score,method="pearson")
  correlation[i,5:6] <- c(test$estimate,test$p.value)  
} 

# correlate with PDL1 expression
AllCancers <- c("BLCA","BRCA","COADREAD","GBM",'HNSC',"LAML","LUAD","LUSC","OV","UCEC")
correlation <- matrix(0,nrow=length(AllCancers),ncol=2)
rownames(correlation) <- AllCancers
colnames(correlation) <- c("Correlation","p-value")
for (i in (1:length(AllCancers))){
  load(paste0("/Users/magali/Desktop/recherche/PancancerAnalysis/AMARETTO/Data/TCGA_",AllCancers[i],"_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData"))
  if (length(intersect(rownames(ProcessedData$MA_TCGA),"CD274"))>0){
    PDCD1 <- ProcessedData$MA_TCGA["CD274",]
    OAS2 <- ProcessedData$MA_TCGA["OAS2",] 
    reg1 <- lm(PDCD1~OAS2)
    par(cex=.8)
    plot(OAS2,PDCD1, xlab="OAS2 expression",ylab="CD274 expression",main=paste0("Correlation for ",AllCancers[i]))
    abline(reg1)
    test <- cor.test(OAS2,PDCD1,method="pearson",exact=FALSE)
    correlation[i,1:2] <- c(test$estimate,test$p.value)
  }
}
write.csv2(correlation,file="/Users/magali/Desktop/recherche/PancancerAnalysis/article/Images_all/correlationCD274.csv")

# virus data
Virus <- read.csv2("/Users/magali/Desktop/recherche/PancancerAnalysis/AMARETTO/Tang_Viruses_TCGA.txt",sep="\t")
Virus <- Virus[,c(1:10)]
Virus[,1] <- str_replace_all(Virus[,1],"COAD","COADREAD")
Virus[which(Virus[,1]=="READ"),1] <- rep("COADREAD",length(which(Virus[,1]=="READ")))
Correlation <- matrix(0,nrow=length(AllCancers),ncol=2)
colnames(Correlation) <- c("p-value","signe")
rownames(Correlation) <- AllCancers
Table <- matrix(0,nrow=length(AllCancers),ncol=2)
colnames(Table) <- c("No viral infection","Viral infection")
rownames(Table) <- AllCancers
for (i in (1:length(AllCancers))){
  CancerSite = AllCancers[i]
  Virus_cancer <- Virus[which(Virus[,1]==CancerSite),]
  Virus_cancer$Sample.barcode <- substr(Virus_cancer$Sample.barcode,0,15)
#  Virus_cancer <- Virus_cancer[which(substr(Virus_cancer$Sample.barcode, nchar(Virus_cancer$Sample.barcode)-1, nchar(Virus_cancer$Sample.barcode))=="01"),]
  load(paste0("/Users/magali/Desktop/recherche/PancancerAnalysis/AMARETTO/Data/TCGA_",CancerSite,"_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData"))
  Virus_cancer <- Virus_cancer[which(Virus_cancer$Sample.barcode %in% colnames(ProcessedData$MA_TCGA)),]  
  OAS2 <- ProcessedData$MA_TCGA["OAS2",intersect(Virus_cancer$Sample.barcode,colnames(ProcessedData$MA_TCGA))]
  Samples_virus <- rep("Viral infection",nrow(Virus_cancer))
  names(Samples_virus) <- Virus_cancer$Sample.barcode
  Samples_virus[which(Virus_cancer$Virus.description=="N/A")] <- rep("No viral infection",length(which(Virus_cancer$Virus.sequence.ID=="N/A")))
  if (length(table(Samples_virus))==1){
    Table[i,] <- c(table(Samples_virus),0)
  } else {
    Table[i,] <- table(Samples_virus)
  }
  
  if (length(which(Samples_virus=="Viral infection"))>0){
    if (length(Samples_virus)>length(OAS2)){
      I <- names(which(table(names(Samples_virus))>1))
      for (j in (1:length(I))){
        I2 <- which(names(Samples_virus)==I[j])
        test <- Samples_virus[I2]
        if (length(table(test)==1)){
          Samples_virus <- Samples_virus[-I2[-1]]
        } else {
          Samples_virus[I2] <- rep("Viral infection",length(I2))
          Samples_virus <- Samples_virus[-I2[-1]]
        }
      }
    }
    Combined <- cbind(OAS2=OAS2,Virus=Samples_virus)
    Combined <- data.frame(Combined)
    Combined$OAS2 <- as.numeric(as.character(Combined$OAS2))
    p <- ggplot(Combined,aes(Virus,OAS2),main="test")
    print(p + geom_boxplot(aes(fill = Virus)) +ggtitle(paste0("Correlations of viral infection with OAS2 for cancer ",CancerSite))
            + labs(x = "Viral infection",y="OAS2 expression")+
            scale_x_discrete(breaks=c("No viral infection", "Viral infection"),
                             labels=c("No viral infection", "Viral infection"))+
            scale_fill_discrete(name="Viral infection",
                                breaks=c("No viral infection", "Viral infection"),
                                labels=c("No viral infection", "Viral infection")))
  
    test<- wilcox_test(Combined$OAS2~Combined$Virus)
    Correlation[i,] <- c(pvalue(test),(statistic(test)>0))
  }
}

# xenobiotic gene signature
Xenob <- read.csv2("/Users/magali/Desktop/recherche/PancancerAnalysis/AMARETTO/Results/Xenobiotic.csv")
Xenob <- Xenob[,2]
Xenob_human <- read.csv2("/Users/magali/Desktop/recherche/PancancerAnalysis/AMARETTO/Results/Xenob_human.csv")
Xenob_human <- Xenob_human[,2]
AllCancers <- c("BLCA","BRCA","COADREAD","HNSC","LUAD","LUSC","OV","UCEC")
Results <- matrix(0,nrow=length(AllCancers),ncol=3)
for (i in (1:length(AllCancers))){
  CancerSite <- AllCancers[i]
  load(paste0("/Users/magali/Desktop/recherche/PancancerAnalysis/Data/TCGA_",CancerSite,"_ProcessedData_MA_CNV.Rdata"))
  MA_TCGA <- ProcessedData$MA_TCGA
  
  OverlapGenes <- intersect(Xenob,rownames(MA_TCGA))
  score <- colMeans(MA_TCGA[OverlapGenes,])
  
  GeneExp <- MA_TCGA[which(rownames(MA_TCGA)=="GPX2"),]
  plot(GeneExp,score,xlab="GPX2",ylab="score signature",main=paste0("Correlation between score signature and gene OAS2 for ",CancerSite))
  correlation=cor(GeneExp,score)
  reg1 <- lm(score~GeneExp)
  par(cex=.8)
  abline(reg1)
  
  test <- cor.test(GeneExp,score,method="spearman",exact=FALSE)
  Results[i,] <- c(p.value=test$p.value,correlation.estimate=test$estimate,sign=sign(test$estimate))
}
colnames(Results) <- c("p.value","correlation.estimate","sign")
rownames(Results) <- AllCancers

# more details on GPX2
Regs <- list()
# extract regulators
for (i in (1:length(ModulesGPX))){
  CancerSite <- strsplit2(ModulesGPX[i],"_")[,2]
  ModuleNr <- strsplit2(ModulesGPX[i],"_")[,1]
  ModuleNr <- as.numeric(str_replace_all(ModuleNr,"M",""))
  load(paste0("/Users/magali/Desktop/recherche/PancancerAnalysis/AMARETTO/Results/results",CancerSite,".RData"))
  Reg <- AMARETTOresults$GeneNames[which(AMARETTOresults$assign==ModuleNr)]
  Regs <- c(Regs,list(Reg))
}
names(Regs) <- ModulesGPX
# ModuleCom=c("M35_BLCA","M44_BRCA","M60_BRCA","M73_COADREAD","M24_HNSC","M90_HNSC","M15_LUSC","M36_LUSC") 
#perhaps we should add these modules which are not part of the community but are regulated by GPX2
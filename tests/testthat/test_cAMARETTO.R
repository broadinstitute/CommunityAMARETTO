# devtools::install_github("broadinstitute/CommunityAMARETTO",ref="develop")
# library(CommunityAMARETTO)
# context("CommunityAMARETTO output data objects testing")
# 
# data("TCGA_LIHC_data")
# data("CCLE_Liver_data")
# 
# 
# 
# AMARETTOresults_all<-list(TCGA_LIHC=TCGA_LIHC_data$AMARETTO_Results,
#                           CCLE_Liver=CCLE_Liver_data$AMARETTO_Results)
# 
# 
# cAMARETTOresults<-CommunityAMARETTO::cAMARETTO_Results(AMARETTOresults_all,
#                                                        gmt_filelist=NULL,
#                                                        NrCores = 4,
#                                                        output_dir = ".",
#                                                        drivers=TRUE)
# 
# cAMARETTOnetworkM<-cAMARETTO_ModuleNetwork(cAMARETTOresults,0.10,5)
# 
# cAMARETTOnetworkC<-cAMARETTO_IdentifyCom(cAMARETTOnetworkM,filterComm = FALSE)
# 
# 
# #Identify significantly connected subnetworks using
# #the Girvan-Newman algorithm
# cAMARETTOnetworkC<-cAMARETTO_IdentifyCom(cAMARETTOnetworkM,
#                                          filterComm = FALSE,
#                                          ratioCommSize = 0.01,
#                                          MinRuns = 2,
#                                          ratioRunSize = 0.1,
#                                          ratioEdgesInOut = 0.5,
#                                          plot_network = FALSE)
# 
# 
# 


# 
# testthat::test_that("Checking cAMARETTOresults object if it is
# in decent shape",{
#   expect_equal(length(cAMARETTOresults$runnames)==
#   length(AMARETTOinit_all),TRUE)
# 
# })
# 
# testthat::test_that("Checking cAMARETTOnetworkM object to
# be in decent shape",{
#   expect_equal(dim(cAMARETTOnetworkM$layoutMN)[1],40)
#   expect_equal(dim(cAMARETTOnetworkM$layoutMN)[2],2)
# 
# })
# 
# testthat::test_that("Checking cAMARETTOnetworkC object
# to be in decent shape",{
#   expect_equal(dim(cAMARETTOnetworkC$commEdgeInfo)[1],21)
#   expect_equal(dim(cAMARETTOnetworkC$commEdgeInfo)[2],9)
# })

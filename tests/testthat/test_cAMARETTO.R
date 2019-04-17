devtools::install_github("broadinstitute/CommunityAMARETTO",ref="develop")
library(CommunityAMARETTO)
context("CommunityAMARETTO output data objects testing")


data("AMARETTOinit_all")
data("AMARETTOresults_all")

cAMARETTOresults<-cAMARETTO_Results(AMARETTOinit_all =AMARETTOinit_all ,
                                    AMARETTOresults_all = AMARETTOresults_all,
                                    gmt_filelist = NULL,
                                    NrCores = 4,
                                    drivers = TRUE)

cAMARETTOnetworkM<-cAMARETTO_ModuleNetwork(cAMARETTOresults,
                                           pvalue = 0.05,
                                           inter = 5,
                                           plot_network = FALSE)

#Identify significantly connected subnetworks using the Girvan-Newman algorithm
cAMARETTOnetworkC<-cAMARETTO_IdentifyCom(cAMARETTOnetworkM,
                                         filterComm = FALSE,
                                         ratioCommSize = 0.01,
                                         MinRuns = 2,
                                         ratioRunSize = 0.1,
                                         ratioEdgesInOut = 0.5,
                                         plot_network = FALSE)






testthat::test_that("Checking cAMARETTOresults object if it is in decent shape",{
  expect_equal(length(cAMARETTOresults$runnames)==length(AMARETTOinit_all),TRUE)
  
})

testthat::test_that("Checking cAMARETTOnetworkM object to be in decent shape",{
  expect_equal(dim(cAMARETTOnetworkM$layoutMN)[1],40)
  expect_equal(dim(cAMARETTOnetworkM$layoutMN)[2],2)
  
})

testthat::test_that("Checking cAMARETTOnetworkC object to be in decent shape",{
  expect_equal(dim(cAMARETTOnetworkC$commEdgeInfo)[1],21)
  expect_equal(dim(cAMARETTOnetworkC$commEdgeInfo)[2],9)
})

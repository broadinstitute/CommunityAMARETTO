library(devtools)
#install.packages("/Users/mohsennabian/Documents/GitHub/CommunityAMARETTO",repos = NULL, type = "source")
install_github("broadinstitute/CommunityAMARETTO",
                         ref='develop',
                         force = 'TRUE',
                         upgrade='never')
library(CommunityAMARETTO)
context("CommunityAMARETTO output data objects testing")
TCGA_LIHC_data<-CommunityAMARETTO::TCGA_LIHC_data
scHCV_data<-CommunityAMARETTO::scHCV_data
scHBV_data<-CommunityAMARETTO::scHBV_data
#CCLE_Liver_data<-CommunityAMARETTO::CCLE_Liver_data


AMARETTOresults_all<-list(TCGA_LIHC=TCGA_LIHC_data$AMARETTO_Results,
                          scHCV=scHCV_data$AMARETTO_Results,
                          scHBV=scHBV_data$AMARETTO_Results)


cAMARETTOresults<-CommunityAMARETTO::cAMARETTO_Results(AMARETTOresults_all,
                                                       gmt_filelist=NULL,
                                                       NrCores = 4,
                                                       output_dir = ".",
                                                       drivers=TRUE)

cAMARETTOnetworkM<-CommunityAMARETTO::cAMARETTO_ModuleNetwork(cAMARETTOresults,
                                           pvalue=0.10,
                                           inter = 5,
                                           plot_network = FALSE)

cAMARETTOnetworkC<-CommunityAMARETTO::cAMARETTO_IdentifyCom(cAMARETTOnetworkM,
                                                            filterComm = FALSE,
                                                            plot_network = FALSE)


unlink("gmt_files", recursive = TRUE)

test_that("Check if results are calculated", {
  expect_equal(is.null(cAMARETTOresults), FALSE)
  expect_equal(is.null(cAMARETTOnetworkM), FALSE)
  expect_equal(is.null(cAMARETTOnetworkC), FALSE)
  expect_equal(is.null(cAMARETTOresults$hgt_modules), FALSE)
  expect_equal(is.null(cAMARETTOresults$genelists), FALSE)
})
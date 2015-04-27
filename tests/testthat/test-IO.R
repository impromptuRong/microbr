################################################################################
# Test IO and general access method
################################################################################
library("microbr")
library("testthat")
################################################################################
test_that("Initialization method", {
  data("oral")
  data("oral_raw")
  ob <- physet(otu_table = oral_raw$rawdata, sample_data = oral_raw$metainfo, 
               tax_table = oral_raw$taxonomy, phy_tree = oral_raw$rawtree)
  expect_that(ob, is_equivalent_to(oral))
  expect_that(slotNames(oral), 
              is_equivalent_to(c("otu_table", "sample_data", "tax_table", 
                                 "phy_tree", "edge_len", "edge_com", 
                                 "edge_mat", "seqdep")))
})
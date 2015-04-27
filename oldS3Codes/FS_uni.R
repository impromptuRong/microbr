############################################################
#' Feature Selection wapper for phyloseq class
#' Detach phyloseq class to otu-table and grouping information
#' Call FS_basic.R

################   Basic Feature Selection   ###############
FS_uni <- function(physeq, samtype, name, method=c("StaTest", "SPLSDA", "ISA", "RF", "ENet", "SVMRFE"), kernel="linear"){
  data <- data.frame(otu_table(physeq))
  sample_data <- data.frame(sample_data(physeq))
  # tree <- phy_tree(physeq)
  # taxonomy <- data.frame(tax_table(physeq))
  
  ID <- sample.names(physeq)
  feature <- species.names(physeq)
  group <- sample_data[,samtype]
  FS_basic(data, group, ID, name, kernel=kernel, method=method)
}

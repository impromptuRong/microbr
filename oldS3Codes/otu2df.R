######################################################
#' Modified Version of phyloseq:::otu2df
#######   Modified Version of phyloseq:::otu2df   #######
otu2df <- function(physeq, keepOnlyTheseTaxa=NULL, keepOnlyTheseSample=NULL, threshold=0, normalize=TRUE){
  #taxMap=rank.names(physeq), samMap=sample.variables(physeq), 
  #######  Define Order  #######
  speorder <- species.names(physeq)
  samorder <- sample.names(physeq)
  if(taxa_are_rows(physeq)){	physeq <- t(physeq)	}
  dfdata <- data.frame(otu_table(physeq))[samorder, speorder]
  # Normalization Part
  if(is.numeric(normalize)){
    dfdata <- dfdata/normalize
  } else if(normalize){
    dfdata <- dfdata/rowSums(dfdata)
  } else {	dfdata <- dfdata	}
  samples <- data.frame(sample_data(physeq))[samorder, ]
  species <- data.frame(tax_table(physeq))[speorder, ]
  #######  Organize Tables  #######
  if(!is.null(keepOnlyTheseTaxa)){
    dfdata <- dfdata[,speorder %in% keepOnlyTheseTaxa]
    species <- subset(species, speorder %in% keepOnlyTheseTaxa)
    species <- droplevels(species)
  }
  if(!is.null(keepOnlyTheseSample)){
    dfdata <- dfdata[samorder %in% keepOnlyTheseSample,]
    samples <- subset(samples, samorder %in% keepOnlyTheseSample)
    samples <- droplevels(samples)
  }
  #######  Reshape Abundance Table  #######
  tmp.Abu <- reshape(dfdata, direction="long", varying=colnames(dfdata), times=colnames(dfdata), v.names="Abu", timevar="TaxaGroup", idvar="Sample", ids=rownames(dfdata))
  tmp.sam <- samples[rep(1:nrow(samples),times=nrow(species)),]
  tmp.spe <- species[rep(1:nrow(species), each=nrow(samples)),]
  df <- data.frame(tmp.Abu, tmp.sam, tmp.spe)
  #######  Apply Restriction to df  ########
  df <- df[df$Abu>=threshold,]
  return(df)
}

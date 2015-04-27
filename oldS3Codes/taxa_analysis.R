######################################################
#' Basic Taxa Abundance Analysis and alpha diversity
#' Rarefaction analysis with random shuffling
#' 
#########   Taxanomy Analysis Fundemental   #########
taxa_analysis <- function(physeq, samtype="Sample", spetype="TaxaGroup", color=NULL, keepOnlyTheseTaxa=NULL, keepOnlyTheseSample=NULL, 
                          facet_formula=.~TaxaGroup, OTUpoints=FALSE, labelOTUs=FALSE, threshold=0, normalize=TRUE, trans=trans, rarefy=FALSE, permu=100){
  Abu_df <- otu2df(physeq, keepOnlyTheseTaxa, keepOnlyTheseSample, threshold=threshold, normalize=TRUE)
  
  #######   boxplot   #######
  if(OTUpoints){
    pbox_hor <- ggplot(Abu_df, aes_string(x=spetype, y="Abu", fill=samtype)) + geom_boxplot(outlier.colour="NA") + 
      labs(y="Relative Abundance", x=spetype) + facet_grid(facet_formula, scales="free", space="free", drop=TRUE) + 
      geom_point(aes_string(color=samtype), position="jitter", size=1.5, alpha=I(1/2)) + 
      scale_fill_manual(values=color) + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1))
  } else {
    pbox_hor <- ggplot(Abu_df, aes_string(x=spetype, y="Abu", fill=samtype)) + geom_boxplot(outlier.size=1.5) + 
      labs(y="Relative Abundance", x=spetype) + facet_grid(facet_formula, scales="free", space="free", drop=TRUE) + 
      scale_fill_manual(values=color) + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1))
  }
  #####   Get facet formula   #####
  ft <- asOneFormula(facet_formula, as.formula(paste("~",samtype,"+",spetype,sep="")))
  var_form <- rev(strsplit(toString(ft),"\\, ")[[1]])[1]
  summary_formula <- as.formula(paste("Abu",var_form,sep="~"))
  Dsummary <- summaryBy(summary_formula, Abu_df, FUN=c(mean, sd, length))
  
  #####   Bar Plot with Error bar   #####
  pbar_hor <- ggplot(Dsummary, aes_string(x=spetype, y="Abu.mean", fill=samtype)) + geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymin=Abu.mean-Abu.sd/sqrt(Abu.length), ymax=Abu.mean+Abu.sd/sqrt(Abu.length)), width=.25, position=position_dodge(.9)) + 
    theme_bw() + facet_grid(facet_formula, scales="free", space="free", drop=TRUE) + theme(axis.text.x=element_text(angle=90, hjust=1)) + 
    scale_fill_manual(values=color) + labs(y="Relative Abundance", x=spetype)
  if(OTUpoints){
    pbar_hor <- pbar_hor + geom_point(data=Abu_df, aes_string(x=spetype, y="Abu", color=samtype), position="jitter", size=1.5, alpha=I(1/2))
  }
  
  #####   Estimate Richness   #####
  if(!is.null(keepOnlyTheseTaxa)){	physeq <- prune_species(keepOnlyTheseTaxa, physeq)	}
  if(!is.null(keepOnlyTheseSample)){	physeq <- prune_samples(keepOnlyTheseSample, physeq)}
  sumfun <- function(x, ...){c(mean=mean(x, ...), var=var(x, ...), sd=sd(x, ...), length=length(x))}
  # unlockBinding("estimate_richness", as.environment("package:phyloseq"))
  # assignInNamespace("estimate_richness", my.estimate_richness, ns="phyloseq", envir="package:phyloseq")
  # assign("estimate_richness", my.estimate_richness, "package:phyloseq")
  # rich.plot <- plot_richness_estimates(physeq, samtype, color=samtype) + scale_color_manual(values=color) + theme_bw()
  result <- rarefy_analysis(physeq, group=samtype, color=color, rarefy=rarefy, permu=100)
  
  #####   Summary Result   #####
  result[["tax.boxplot"]] <- pbox_hor
  result[["tax.barplot"]] <- pbar_hor
  return(result)
}


##############    Rarefaction Analysis   #################
rarefy_analysis <- function(physeq, group="group", color=NULL, rarefy=FALSE, permu=100, ...){
  if(class(physeq)=="phyloseq"){
    otu_table <- otu_table(physeq)@.Data
    if(otu_table(physeq)@taxa_are_rows){otu_table <- t(otu_table)}
    if(class(group)=="character"){group <- data.frame(sample_data(physeq))[, group]}
  } else {
    otu_table <- as.matrix(physeq)
  }
  if(class(otu_table)!="matrix" || class(group)=="character"){stop("Input Data is Invalid.\n")}
  
  #####   Estimate alpha Diversity, Richness and Eveness   #####
  eco <- my.estimate_richness(otu_table, split=TRUE)
  df_eco <- eco[,c(1,2,4,6,7,8)]
  eco_test <- para_t.test(df_eco, group, p.adj="none", ...)
  df_eco <- reshape(df_eco, direction="long", varying=colnames(df_eco), v.names="value", timevar="index", times=colnames(df_eco), idvar="sample", ids=rownames(df_eco))
  df_eco$group <- rep(group, times=6)
  df_eco_summary <- summaryBy(.~group+index, data=df_eco, FUN=function(x, ...){c(mean=mean(x, ...), sd=sd(x, ...), length=length(x))}, na.rm=TRUE)
  pbox_rich <- ggplot(df_eco[grep("^S", df_eco$index),], aes(x=group, y=value, fill=group)) + geom_boxplot(outlier.size=1.5) + 
    labs(y="Richness", x="") + facet_grid(.~index, scales="free", space="free", drop=TRUE) + 
    scale_fill_manual(values=color) + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1))
  pbar_rich <- ggplot(df_eco_summary[grep("^S", df_eco_summary$index),], aes(x=group, y=value.mean, fill=group)) + geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymin=value.mean-value.sd/sqrt(value.length), ymax=value.mean+value.sd/sqrt(value.length)), width=.25, position=position_dodge(.9)) + 
    labs(y="Richness", x="") + facet_grid(.~index, scales="free", space="free", drop=TRUE) +
    scale_fill_manual(values=color) + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1))
  pbox_alpha <- ggplot(df_eco[-grep("^S", df_eco$index),], aes(x=group, y=value, fill=group)) + geom_boxplot(outlier.size=1.5) + 
    labs(y="Diversity", x="") + facet_grid(.~index, scales="free", space="free", drop=TRUE) + 
    scale_fill_manual(values=color) + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1))
  pbar_alpha <- ggplot(df_eco_summary[-grep("^S", df_eco_summary$index),], aes(x=group, y=value.mean, fill=group)) + geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymin=value.mean-value.sd/sqrt(value.length), ymax=value.mean+value.sd/sqrt(value.length)), width=.25, position=position_dodge(.9)) + 
    labs(y="Diversity", x="") + facet_grid(.~index, scales="free", space="free", drop=TRUE) + 
    scale_fill_manual(values=color) + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1))
  
  result <- list()
  result[["tax.diversity"]] <- eco
  result[["tax.diver.test"]] <- eco_test
  result[["tax.richness"]][["boxplot"]] <- pbox_rich
  result[["tax.richness"]][["barplot"]] <- pbar_rich
  result[["tax.alphaindex"]][["boxplot"]] <- pbox_alpha
  result[["tax.alphaindex"]][["barplot"]] <- pbar_alpha
  if(!rarefy){return(result)}
  
  #####   bootstrap Rarefaction Curve   #####
  sample_permute <- function(seq_list, dep_list, sample=NULL, permu=100){
    dep_eco <- lapply(dep_list, function(x){replicate(permu, unlist(my.estimate_richness(t(as.matrix(table(sample(seq_list, x)))))), simplify=TRUE)})
    dep_mean <- sapply(dep_eco, function(x){apply(x[c(1,2,4,6,7,8),],1,mean)})
    rownames(dep_mean) <- paste(rownames(dep_mean), "mean", sep=".")
    dep_sd <- sapply(dep_eco, function(x){apply(x[c(1,2,4,6,7,8),],1,sd)/sqrt(permu)})
    rownames(dep_sd) <- paste(rownames(dep_sd), "se", sep=".")
    result <- data.frame(sample=sample, seq_dep=dep_list, t(rbind(dep_mean, dep_sd)))
    return(result)
  }
  
  #####    Plot rarefaction curve for each group    #####
  rarefy_group <- lapply(levels(group), function(k){
    seqlist <- rep(colnames(otu_table), colSums(otu_table[group==k,]));
    deplist <- c(2^unique(floor(log(c(1:length(seqlist)),2))), length(seqlist));
    sample_permute(seqlist, deplist, sample=k, permu=permu);
  })
  df_rarefy_group <- do.call(rbind, rarefy_group)
  df_rarefy_group <- data.frame(sample=rep(df_rarefy_group$sample,6), seq_dep=rep(df_rarefy_group$seq_dep,6), mean=unlist(df_rarefy_group[,3:8]), se=unlist(df_rarefy_group[,9:14]), 
                                index=factor(rep(c("S.obs","S.chao1","S.ACE","shannon","simpson","pielou"), each=nrow(df_rarefy_group)), c("S.obs","S.chao1","S.ACE","shannon","simpson","pielou")))
  
  rarefy_cat <- ggplot(df_rarefy_group, aes(x=seq_dep, y=mean, group=sample, colour=sample)) + geom_line() + geom_point() + 
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.25, stat="identity", position=position_dodge(.9)) + labs(y="Index", x="") + 
    scale_colour_manual(values=color) + facet_wrap(~index, scales="free", drop=TRUE) + theme_bw()
  
  #####    Plot rarefaction curve for each sample    #####
  rarefy_sample <- lapply(rownames(otu_table), function(k){
    seqlist <- rep(colnames(otu_table), otu_table[k,]);
    deplist <- c(2^unique(floor(log(c(1:length(seqlist)),2))), length(seqlist));
    sample_permute(seqlist, deplist, sample=k, permu=permu);
  })
  df_rarefy_sample <- do.call(rbind, rarefy_sample)
  df_rarefy_sample$group <- rep(group, sapply(rarefy_sample, nrow))
  
  rarefy_sam.shannon <- ggplot(df_rarefy_sample, aes(x=seq_dep, y=shannon.mean, group=sample, colour=group)) + geom_line() + geom_point() + 
    geom_errorbar(aes(ymin=shannon.mean-shannon.se, ymax=shannon.mean+shannon.se), width=.25, position=position_dodge(.9)) + labs(y="shannon", x="") + 
    scale_colour_manual(values=color) + facet_wrap(~sample, nrow=ceiling(sqrt(nrow(otu_table))), scales="free", drop=TRUE) + theme_bw()
  rarefy_sam.simpson <- ggplot(df_rarefy_sample, aes(x=seq_dep, y=simpson.mean, group=sample, colour=group)) + geom_line() + geom_point() + 
    geom_errorbar(aes(ymin=simpson.mean-simpson.se, ymax=simpson.mean+simpson.se), width=.25, position=position_dodge(.9)) + labs(y="simpson", x="") + 
    scale_colour_manual(values=color) + facet_wrap(~sample, nrow=ceiling(sqrt(nrow(otu_table))), scales="free", drop=TRUE) + theme_bw()
  rarefy_sam.pielou <- ggplot(df_rarefy_sample, aes(x=seq_dep, y=pielou.mean, group=sample, colour=group)) + geom_line() + geom_point() + 
    geom_errorbar(aes(ymin=pielou.mean-pielou.se, ymax=pielou.mean+pielou.se), width=.25, position=position_dodge(.9)) + labs(y="pielou", x="") + 
    scale_colour_manual(values=color) + facet_wrap(~sample, nrow=ceiling(sqrt(nrow(otu_table))), scales="free", drop=TRUE) + theme_bw()
  rarefy_sam.S.obs <- ggplot(df_rarefy_sample, aes(x=seq_dep, y=S.obs.mean, group=sample, colour=group)) + geom_line() + geom_point() + 
    geom_errorbar(aes(ymin=S.obs.mean-S.obs.se, ymax=S.obs.mean+S.obs.se), width=.25, position=position_dodge(.9)) + labs(y="Rarefy", x="") + 
    scale_colour_manual(values=color) + facet_wrap(~sample, nrow=ceiling(sqrt(nrow(otu_table))), scales="free", drop=TRUE) + theme_bw()
  rarefy_sam.S.chao1 <- ggplot(df_rarefy_sample, aes(x=seq_dep, y=S.chao1.mean, group=sample, colour=group)) + geom_line() + geom_point() + 
    geom_errorbar(aes(ymin=S.chao1.mean-S.chao1.se, ymax=S.chao1.mean+S.chao1.se), width=.25, position=position_dodge(.9)) + labs(y="S.Chao1", x="") + 
    scale_colour_manual(values=color) + facet_wrap(~sample, nrow=ceiling(sqrt(nrow(otu_table))), scales="free", drop=TRUE) + theme_bw()
  rarefy_sam.S.ACE <- ggplot(df_rarefy_sample, aes(x=seq_dep, y=S.ACE.mean, group=sample, colour=group)) + geom_line() + geom_point() + 
    geom_errorbar(aes(ymin=S.ACE.mean-S.ACE.se, ymax=S.ACE.mean+S.ACE.se), width=.25, position=position_dodge(.9)) + labs(y="S.ACE", x="") + 
    scale_colour_manual(values=color) + facet_wrap(~sample, nrow=ceiling(sqrt(nrow(otu_table))), scales="free", drop=TRUE) + theme_bw()
  
  result[["tax.rarefy"]][["group"]] <- rarefy_cat
  result[["tax.rarefy"]][["sample"]][["shannon"]] <- rarefy_sam.shannon
  result[["tax.rarefy"]][["sample"]][["simpson"]] <- rarefy_sam.simpson
  result[["tax.rarefy"]][["sample"]][["pielou"]] <- rarefy_sam.pielou
  result[["tax.rarefy"]][["sample"]][["S.obs"]] <- rarefy_sam.S.obs
  result[["tax.rarefy"]][["sample"]][["S.chao1"]] <- rarefy_sam.S.chao1
  result[["tax.rarefy"]][["sample"]][["S.ACE"]] <- rarefy_sam.S.ACE
  return(result)
}


my.estimate_richness <- function(physeq, split=TRUE){
  # Check input Data and flipping/splitting data
  if(class(physeq)=="phyloseq"){
    otu_table <- otu_table(physeq)@.Data
    if(otu_table(physeq)@taxa_are_rows){otu_table <- t(otu_table)}
  } else{
    otu_table <- as.matrix(physeq)
  }
  if(!split){otu_table <- as.matrix(colSums(otu_table))}
  if(ncol(otu_table)==1){otu_table <- t(otu_table)}
  
  # Check for Input, singletons, and fraction Data.
  # These metrics only really meaningful if singletons are included.
  if(class(otu_table)!="matrix"){stop("Input Data is Invalid.\n")}
  if(sum(ceiling(otu_table)) > sum(otu_table)){
    warning("Ceiling Data were used.\n",
            "Results of Richness Estimation are probably unreliable.\n")
  }
  if(!any(otu_table==1)){
    warning("The experiment object you have provided does not have any singletons. This is highly suspicious.\n",
            "Results of richness estimates are probably unreliable, or wrong.\n",
            "if you have already trimmed low-abundance taxa from the data.\n",
            "It is recommended that you find the un-trimmed data and retry.")
  }
  
  # Calculate standard alpha diversity, Richness and Eveness
  shannon <- diversity(otu_table, index="shannon", MARGIN=1, base=exp(1))
  simpson <- diversity(otu_table, index="simpson", MARGIN=1, base=exp(1))
  richness <- estimateR(ceiling(otu_table))
  pielou <- ifelse(specnumber(otu_table)>1, shannon/log(specnumber(otu_table)), 0)
  eco <- data.frame(t(rbind(richness, shannon, simpson, pielou)))
  return(eco)
}

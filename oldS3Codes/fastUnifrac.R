
fastUnifrac <- function(physet, tree=NULL, seq.depth=NULL, parallel=FALSE, method="Edge") {
  if (class(physet) == "physet") {
    otu_table <- otu_table(physet)
    if (is.null(tree))
      tree <- phy_tree(physet)
    seq.depth <- seqdep(physet)
  } else {
    otu_table <- Matrix(as.matrix(physet))
  } 
  
  if (class(tree) != "phylo" || !is.rooted(tree) || is.null(tree$edge.length))
    stop("Invalid phylogenetic tree: not a phylo class, not a rooted tree, no branch lengths!")
  if (is.null(tree$node.label))
    tree$node.label <- paste("Branch", 1:tree$Nnode, sep="_")
  if (class(otu_table) != "Matrix")
    stop("Invalid otu table!")
  if (is.null(seq.depth)) seq.depth <- rowSums(otu_table)
  otu_table <- otu_table/seq.depth
  
  common_taxa <- intersect(colnames(otu_table), tree$tip.label)
  if (length(common_taxa) < length(tree$tip.label)) {
    warning("Not all taxa in tree are included in otu table")
    tree <- drop.tip(tree, setdiff(tree$tip.label, common_taxa))
  }
  if (length(common_taxa) < length(colnames(otu_taxa)))
    warning("Not all taxa in otu table are included in tree")
  otu_table <- otu_table[, tree$tip.label]
  
  Nsample <- nrow(otu_table)
  Ntaxa <- length(tree$tip.label)
  Nnode <- tree$Nnode
  Nedge <- nrow(tree$edge)
  edge_length <- tree$edge.length
  # if((Ntaxa + Nnode - 1)!=Nedge) stop("Tree structure error")

  extend_par <- function(edge, tree=tree, Ntaxa=Ntaxa) {
    int.nodes <- tree$edge[edge, 2]
    if (int.nodes <= Ntaxa)
      return(int.nodes)
    sons <- c()
    repeat {
      int.nodes <- tree$edge[which(tree$edge[, 1] %in% int.nodes), 2]
      if (length(int.nodes) == 0)
        break
      sons <- c(sons, int.nodes)
    }
    sons <- sons[sons <= Ntaxa]
    return(sons)
  }
  extend <- function(edge) extend_par(edge, tree, Ntaxa)
  
  ###   Non-parallel and foreach   ###
  if (parallel) {
    edge_list <- try(foreach(edge = 1:Nedge) %dopar% extend(edge))
    if (isTRUE(all.equal(class(edge_list), "try-error"))) {
      sfInit(parallel=TRUE, cpus=8)
      # sfClusterSetupRNG(type="SPRNG", seed=12345)
      edge_list <- sfLapply(1:Nedge, extend_par, tree=tree, Ntaxa=Ntaxa)
    }
  } else {
    system.time(edge_list <- lapply(1:Nedge, extend))
  }
  edge_name <- c(tree$tip.label, tree$node.label)[tree$edge[, 2]]
  if (is.rooted(tree)) {
    edge_list[[Nedge+1]] <- c(1:Nedge)
    edge_name <- c(edge_name, tree$node.label[which(table(c(tree$edge))==2) - Ntaxa])
    edge_length <- c(edge_length, 0)
  }
  
  edge_bool <- Matrix(sapply(edge_list, function(x){1:Ntaxa %in% x})*1);
  rownames(edge_bool) <- tree$tip.label
  colnames(edge_bool) <- edge_name;
  edge_matrix <- otu_table %*% edge_bool;
  
  ###   Organize Return List   ###
  method = c("Edge","dPCoA","non_w","w.non_nor","w.nor") %in% method;
  if (!sum(method)) {
    warning("Input Methods are Invalid!\nReturn edge_list, edge_matrix, edge_bool Only.\nProper Methods are: Edge; dPCoA; non_w; w.non_nor; w.nor");
    method = c(1, 0, 0, 0, 0)
  }
  
  result <- list();
  if(method[1]){
    result[["edge_list"]] <- edge_list
    result[["edge_bool"]] <- edge_bool
    result[["edge_matrix"]] <- edge_matrix
  }
  if(method[2]){
    D <- matrix(0, Nsample, Nsample, dimnames=list(rownames(otu_table), rownames(otu_table)))
    for(i in 1:Nsample){
      for(j in 1:Nsample){
        D[i,j] <- (edge_matrix[i,]-edge_matrix[j,])^2 %*% edge_length
      }
    }
    result[["dist"]][["dPCoA"]] = as.dist(D);
  }
  if(method[3]){
    D <- matrix(0, Nsample, Nsample, dimnames=list(rownames(otu_table), rownames(otu_table)))
    for(i in 1:Nsample){
      for(j in 1:Nsample){
        D[i,j] <- (abs((edge_matrix[i,]>0)-(edge_matrix[j,]>0)) %*% edge_length)/(abs((edge_matrix[i,]>0)|(edge_matrix[j,]>0)) %*% edge_length)
      }
    }
    result[["dist"]][["non_w"]] = as.dist(D);
  }
  if(method[4]){
    D <- matrix(0, Nsample, Nsample, dimnames=list(rownames(otu_table), rownames(otu_table)))
    for(i in 1:Nsample){
      for(j in 1:Nsample){
        D[i,j] <- abs(edge_matrix[i,]-edge_matrix[j,]) %*% edge_length
      }
    }
    result[["dist"]][["w.non_nor"]] = as.dist(D);
  }
  if(method[5]){
    D <- matrix(0, Nsample, Nsample, dimnames=list(rownames(otu_table), rownames(otu_table)))
    leaf_depth <- setNames(node.depth.edgelength(tree), c(tree$tip.label, tree$node.label))[colnames(otu_table)]
    for(i in 1:Nsample){
      for(j in 1:Nsample){
        D[i,j] <- (abs(edge_matrix[i,]-edge_matrix[j,]) %*% edge_length)/((otu_table[i,]+otu_table[j,]) %*% leaf_depth)
      }
    }
    result[["dist"]][["w.nor"]] = as.dist(D);
  }
  return(result)
}

################################################################################
# Physet Ploting Options
# Plot ordination
# Plot rarefaction
################################################################################
#' @title Physet Ploting Options
#' @description This function is used to generate ploting options for objects. 
#' 
#' @param ... The sample and taxa information if available. \cr 
#' Current supported plotting options: \cr
#' \itemize{
#' \item{opts}{ Update current \code{opts} with new parameters. }
#' \item{info}{ Build up \code{opts} with new sample and taxa info. }
#' \item{color}{ The colors used for plotting. }
#' \item{colorVar}{ The variable used to specify colors. }
#' \item{fillVar}{ The variable used to specify fill colors for shapes. }
#' \item{shape}{ The shapes used for plotting. }
#' \item{shapeVar}{ The variable used to specify shapes. }
#' \item{labelVar}{ The variable used to specify labels. }
#' \item{axes.2d}{ Axes ID as \dQuote{x} and \dQuote{y} for 2d plot. }
#' \item{axes.3d}{ Axes ID as for 3d plot. }
#' \item{keepOnly}{ A list with attributes \code{sample} and \code{taxa}. This 
#' parameter decides which sample and taxa are showned in figures. }
#' \item{heatmapscale}{ The scale used to generate heatmap color. }
#' \item{heatmaptrans}{ The pre-treatment method to data in heatmap. }
#' \item{outputName}{ The prefix name used for output tables and figures. }
#' \item{...}{ Unused parameters for now. Possibly added to make the function 
#' compatiblewith old source code for igraph and 3d plotting. }
#' }
#' @return A list includes all necessary plotting parameters for plot method 
#' in \code{microbr}. 
#' 
#' @seealso \code{link{plot.ordination}}
#' @examples
#' data(oral)
#' colorVar <- list(sample = "Group", taxa = "Domain")
#' shapeVar <- list(sample = "Periodontitis")
#' fillVar <- colorVar
#' opts <- phyplotOptions(info = oral, colorVar = colorVar, shapeVar = shapeVar, fillVar = fillVar)
#' ord <- ordinate(oral, "NMDS", method = "unifrac.w.un")
#' res <- plot(ord, opts)
#' res$biplot
#' res$splitplot
#' 
#' @rdname phyplotOptions
#' @export
#' 
# Unused: max.dist = 0.3, normalize = TRUE, flip = FALSE, distname = NULL, 
phyplotOptions <- function(...) {
  # inputs <- c(mget(ls()), list(...))
  inputs <- list(...)
  opts <- inputs$opts
  inputs$opts <- NULL
  if (is.null(opts))
    opts <- list(color = .phycolors,  
                 colorVar = list(sample = NULL, taxa = NULL), 
                 fill = .phycolors, 
                 fillVar = list(sample = NULL, taxa = NULL), 
                 shape = .physhapes, 
                 shapeVar = list(sample = NULL, taxa = NULL), 
                 labelVar = list(sample = "Sam_ID", taxa = "Tax_ID"), 
                 axes.2d = c(1, 2), 
                 axes.3d = NULL, 
                 keepOnly = list(sample = NULL, taxa = NULL), 
                 heatmapscale = c("white", "steelblue"), 
                 heatmaptrans = scales::log_trans(4), 
                 outputName = NULL)
  
  # Update each slots based on inputs
  info <- inputs$info
  inputs$info <- NULL
  for (slotN in names(inputs))
    opts[[slotN]] <- inputs[[slotN]]
  
  # Update plotting information
  if (!is.null(info)) {
    if (inherits(info, "physet"))
      info <- list(sample = sample_data(info), taxa = tax_table(info))
    sinfo <- data.frame(id.type = rep("sample", nrow(info$sample)), 
                        row.names=rownames(info$sample),stringsAsFactors=FALSE)
    tinfo <- data.frame(id.type = rep("taxa", nrow(info$taxa)), 
                        row.names = rownames(info$taxa),stringsAsFactors=FALSE)
    # update colorVar
    if (is.null(opts$colorVar$sample)){
      sinfo$colorVar <- "sample"
    } else {
      sinfo$colorVar <- info$sample[[opts$colorVar$sample]]
    }
    if (is.null(opts$colorVar$taxa)) {
      tinfo$colorVar <- "taxa"
    } else {
      tinfo$colorVar <- info$taxa[[opts$colorVar$taxa]]
    }
    # update fillVar
    if (is.null(opts$fillVar$sample)){
      sinfo$fillVar <- "sample"
    } else {
      sinfo$fillVar <- info$sample[[opts$fillVar$sample]]
    }
    if (is.null(opts$fillVar$taxa)) {
      tinfo$fillVar <- "taxa"
    } else {
      tinfo$fillVar <- info$taxa[[opts$fillVar$taxa]]
    }
    # update shapeVar
    if (is.null(opts$shapeVar$sample)) {
      sinfo$shapeVar <- "sample"
    } else {
      sinfo$shapeVar <- info$sample[[opts$shapeVar$sample]]
    }
    if (is.null(opts$shapeVar$taxa)) {
      tinfo$shapeVar <- "taxa"
    } else {
      tinfo$shapeVar <- info$taxa[[opts$shapeVar$taxa]]
    }
    # update labelVar
    if (is.null(opts$labelVar$sample)) {
      sinfo$labelVar <- ""
    } else {
      sinfo$labelVar <- info$sample[[opts$labelVar$sample]]
    }
    if (is.null(opts$labelVar$taxa)) {
      tinfo$labelVar <- ""
    } else {
      tinfo$labelVar <- info$taxa[[opts$labelVar$taxa]]
    }
    sinfo[is.na(sinfo)] <- "NA"
    tinfo[is.na(tinfo)] <- "NA"
    opts$info <- list(sample = sinfo, taxa = tinfo)
  }
  
  ## Update color, fill and shape
  if (!is.null(opts$info)) {
    if (is(opts$color, "function"))
      opts$color <- opts$color(list(opts$info$sample$colorVar, opts$info$taxa$colorVar))
    if (is(opts$fill, "function"))
      opts$fill <- opts$fill(list(opts$info$sample$fillVar, opts$info$taxa$fillVar))
    if (is(opts$shape, "function"))
      opts$shape <- opts$shape(list(opts$info$sample$shapeVar, opts$info$taxa$shapeVar))
  }
  return(opts)
}
################################################################################
#' @title Plotting method for Ordination Object
#' @description The function used to plot ordination class. 
#' 
#' @param x A ordination object generated by \code{\link{ordinate}}.
#' @param opts The plotting options generated by \code{\link{phyplotOptions}}. 
#' @param info A object contains sample and taxa information. This parameter is 
#' used to update \code{opts} with . \code{info} be a list 
#' object with name \dQuote{sample} and \dQuote{taxa}. Or any object inherits
#' \code{\link{physet}} class.
#' @param ... Other parameters passed to ordination and phydist function.  
#' @return A list of objects including: biplot, splitplot(sample and taxa in 2 
#' different panel with same scale), eigen plot, heatmap and plot data frame. 
#' (Note: 3D plot are outputed inside the function because it cannot be saved 
#' as a valid object.) 
#' 
#' @seealso 
#' \code{\link[phyloseq]{plot_ordination}}, 
#' \code{\link[phyloseq]{plot_phyloseq}}
#' 
#' @examples
#' data(oral)
#' colorVar <- list(sample = "Group", taxa = "Domain")
#' shapeVar <- list(sample = "Periodontitis")
#' fillVar <- colorVar
#' opts <- phyplotOptions(colorVar = colorVar, shapeVar = shapeVar, fillVar = fillVar)
#' ord <- ordinate(oral, "CAP", . ~ Group, method = "unifrac.w.un")
#' res <- plot(ord, opts, oral)
#' res$biplot
#' res$splitplot
#' res$eigenplot
#' res$plotDF
#' 
#' @rdname plot-ordination
#' @import ggplot2
#' @export
plot.ordination <- function(x, opts, info, ...) {
  ## Left Join sample info and taxa info into ordination scores
  if (!missing(info))
    opts <- phyplotOptions(opts = opts, info = info)
  scoreS <- x$coord$sample
  scoreT <- x$coord$taxa
  if (!is.null(opts$keepOnly$sample))
    scoreS <- scoreS[opts$keepOnly$sample, ]
  if (!is.null(opts$keepOnly$taxa))
    scoreT <- scoreT[opts$keepOnly$taxa, ]
  DF <- merge(rbind(scoreS, scoreT), rbind(opts$info$sample, opts$info$taxa), 
              by = "row.names", all.x = TRUE, sort = FALSE)
  rownames(DF) <- DF$Row.names
  DF <- DF[,-1]
  color <- opts$color
  shape <- opts$shape
  
  #####   Split Plot   #####
  id.type <- NULL
  sp <- ggplot(data = DF, 
               mapping = aes_string(x = names(DF)[1], y = names(DF)[2], 
                                    color = "colorVar", shape = "shapeVar", 
                                    size = "id.type", fill = "fillVar", 
                                    na.rm = TRUE)) + theme_bw() + 
    geom_point() + scale_x_continuous() + scale_y_continuous() + 
    geom_text(aes_string(label = "labelVar"), size = 3, vjust = 2) + 
    scale_color_manual(name = opts$colorVar$sample, values = color) + 
    scale_shape_manual(name = opts$shapeVar$sample, values = shape) + 
    scale_fill_manual(values = color, guide = "none") + 
    scale_size_manual(values = c(sample = 5, taxa = 2), guide = "none") + 
    facet_wrap( ~ id.type, nrow = 1)
  # theme(title="geom_point", plot.title=theme_text(size=40, vjust=1.5))
  
  #####  Biplot with Arrow  #####
  bp <- ggplot() + theme_bw() + scale_x_continuous() + scale_y_continuous() + 
    geom_point(data = subset(DF, id.type == "sample"), 
               mapping = aes_string(x = names(DF)[1], y = names(DF)[2], 
                                    color = "colorVar", shape = "shapeVar", 
                                    size = "id.type", fill = "fillVar", 
                                    na.rm = TRUE)) + 
    geom_segment(data = subset(DF, id.type == "taxa"), 
                 mapping = aes_string(x = 0, xend = names(DF)[1], 
                                      y = 0, yend = names(DF)[2], 
                                      colour = "colorVar", na.rm = TRUE), 
                 arrow = grid::arrow(length = grid::unit(0.2, "cm")), 
                                     size = 1, alpha=0.75) + 
    geom_text(data = DF, 
              mapping = aes_string(x = names(DF)[1], y = names(DF)[2], 
                                   label = "labelVar", na.rm = TRUE), 
              size = 3, vjust = 2) + 
    scale_color_manual(name = opts$colorVar$sample, values = color) + 
    scale_shape_manual(name = opts$shapeVar$sample, values = shape) + 
    scale_fill_manual(values = color, guide = "none") + 
    scale_size_manual(values = c(sample = 5, taxa = 1), guide = "none") + 
    geom_hline(aes(x = 0), size=.2) + geom_vline(aes(y = 0), size=.2)
  
  #####   Stress/gof, Eigen Plot   #####
  eigenplot <- NULL
  #######   Plot Heatmap, seldom used   #######
  # heatmap <- my.plot_heatmap(physet, PCoA, label="ID", species.label="Taxa")
  heatmap <- NULL
  
  return(list(biplot = bp, splitplot = sp, eigenplot = eigenplot, 
              plotDF = DF, opts = opts))
  
  #######   3D Plot, seldom used   #######
  #   if(!is.null(td.axes)){
  #     plotdata <- PCoA$vectors[,td.axes]
  #     if(!is.null(keepOnlyTheseSample)){
  #       plotdata <- plotdata[keepOnlyTheseSample,]
  #       group <- factor(as.vector(group[keepOnlyTheseSample]))
  #     }  
  #     cloud3d(plotdata, labels=group, filename=paste(name, "dist", distname, "pcoa3d","wrl",sep="."), type="vrml", pointstyle="s", 
  #             metalabels=rownames(plotdata), cols=color, scalefac=4, autoscale="independent", lab.axis=paste("PC",td.axes,sep=""), 
  #             col.axis="darkblue", showaxis=TRUE, col.lab="darkgray", col.bg="white", cex.lab=1, htmlout=NULL, hwidth=1200, hheight=800, showlegend=TRUE)
  #     plotdata <- NMDS.3d$points[,td.axes]
  #     if(!is.null(keepOnlyTheseSample)){
  #       plotdata <- plotdata[keepOnlyTheseSample,]
  #       group <- factor(as.vector(group[keepOnlyTheseSample]))
  #     }
  #     cloud3d(plotdata, labels=group, filename=paste(name,"dist",distname,"nmds3d","wrl",sep="."), type="vrml", pointstyle="s", 
  #             metalabels=rownames(plotdata), cols=color, scalefac=4, autoscale="independent", lab.axis=paste("MDS",td.axes,sep=""), 
  #             col.axis="darkblue", showaxis=TRUE, col.lab="darkgray", col.bg="white", cex.lab=1, htmlout=NULL, hwidth=1200, hheight=800, showlegend=TRUE)
  #   }
  
  #######   IGraph Network, seldom used, add max.dist->phyplotOptions  #######
  #   igraph <- list()
  #   if(!is.null(max.dist)){
  #     for(i in 1:length(max.dist)){
  #       dist <- max.dist[i]
  #       ig <- make_network(physet, type="samples", distance=dmat, max.dist=dist, keep.isolates=TRUE)
  #       if(length(ig[[3]])){
  #         igraph[[paste("igraph", round(dist,2), sep="_")]] <- plot_network(ig, physet, color=samtype, shape=shape, line_weight=0.4, label=NULL)
  #       }
  #     }
  #   }
}
################################################################################
#' @title Plotting method for Rarefaction Object
#' @description The function used to plot rarefaction curves. 
#' 
#' @param x A ordination object generated by \code{\link{ordinate}}. 
#' @param opts The plotting options generated by \code{\link{phyplotOptions}}. 
#' @param info A object conains sample and taxa information. This can be a list 
#' object with name \dQuote{sample} and \dQuote{taxa}. Or any object inherits
#' \code{\link{physet}} class.
#' @param ... Other parameters passed to ordination and phydist function. 
#' @return A list of objects including: biplot, splitplot(sample and taxa in 2 
#' different panel with same scale), eigen plot, heatmap and plot data frame. 
#' 
#' @seealso
#' \code{\link[phyloseq]{plot_richness}}, 
#' \code{\link[phyloseq]{plot_phyloseq}}
#' 
#' @examples
#' data(oral)
#' colorVar <- list(sample = "Group", taxa = "Domain")
#' shapeVar <- list(sample = "Periodontitis")
#' fillVar <- colorVar
#' opts <- phyplotOptions(colorVar = colorVar, shapeVar = shapeVar, fillVar = fillVar)
#' raf <- rarefaction(oral, . ~ Group, rarefy = TRUE, permutation = 100)
#' res <- plot(raf, opts, oral)
#' res$tax.richness$boxplot
#' res$tax.richness$barplot
#' res$tax.alphaindex$boxplot
#' res$tax.alphaindex$barplot
#' res$tax.rarefy$group
#' res$tax.rarefy$sample
#' 
#' @rdname plot-rarefy
#' @import ggplot2
#' @export
plot.rarefy <- function(x, opts, info, ...) {
  ## Left Join sample info and taxa info into diversity scores. 
  if (!missing(info))
    opts <- phyplotOptions(opts = opts, info = info)
  color <- opts$color
  shape <- opts$shape
  print(color)
  print(shape)
  
  res <- list()
  ## Diversity and Richness plot
  eco <- x$alpha.diver[, c(1,2,4,6,7,8)]
  RN <- names(eco)
  if (!is.null(opts$keepOnly$sample))
    eco <- eco[opts$keepOnly$sample, ]
  eco <- reshape(eco, direction = "long", varying = colnames(eco), 
                 v.names = "value", timevar = "index", times = colnames(eco), 
                 idvar = "Sam_ID", ids = rownames(eco))
  plotDF <- merge(eco, opts$info$sample, by.x ="Sam_ID", by.y = "row.names", 
                  all.x = TRUE, sort = FALSE)
  plotSB <- doBy::summaryBy(. ~ colorVar+index, plotDF, FUN=function(p, ...) {
    c(mean = mean(p, ...), sd = sd(p, ...), len = length(p))}, na.rm = TRUE)
  plotSB$value.se <- plotSB$value.sd/sqrt(plotSB$value.len)
  plotSB$lower <- plotSB$value.mean - plotSB$value.se
  plotSB$upper <- plotSB$value.mean + plotSB$value.se
  
  ## box plot
  index <- NULL
  pbox_rich <- ggplot(data = plotDF[grep("^S", plotDF$index), ], 
                      aes_string(x="colorVar", y="value", fill="colorVar")) + 
    geom_boxplot(outlier.size = 1.5) + scale_fill_manual(values = color) + 
    labs(y = "Richness", x = "") + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    facet_grid(. ~ index, scales = "free", space = "free", drop = TRUE)
  pbox_alpha <- ggplot(data = plotDF[-grep("^S", plotDF$index), ], 
                       aes_string(x="colorVar", y="value", fill="colorVar")) + 
    geom_boxplot(outlier.size=1.5) + scale_fill_manual(values = color) + 
    labs(y = "Diversity", x = "") + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    facet_grid(. ~ index, scales = "free", space = "free", drop = TRUE)
  ## bar plot
  pbar_rich <- ggplot(data = plotSB[grep("^S", plotSB$index), ], 
         aes_string(x="colorVar", y="value.mean", fill="colorVar")) + 
    geom_bar(stat = "identity", position = "dodge", colour="black") + 
    geom_errorbar(aes_string(ymin = "lower", ymax = "upper"), colour = "black", 
                  width = .25, position = position_dodge(.9)) + 
    scale_fill_manual(values = color) + 
    labs(y = "Richness", x = "") + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    facet_grid(. ~ index, scales = "free", space = "free", drop = TRUE)
  pbar_alpha <- ggplot(data = plotSB[-grep("^S", plotSB$index), ], 
                       aes_string(x="colorVar", y="value.mean", fill="colorVar")) + 
    geom_bar(stat = "identity", position = "dodge", colour = "black") + 
    geom_errorbar(aes_string(ymin = "lower", ymax = "upper"), colour = "black", 
                  width = .25, position = position_dodge(.9)) + 
    scale_fill_manual(values = color) + 
    labs(y = "Diversity", x = "") + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    facet_grid(. ~ index, scales = "free", space = "free", drop = TRUE)
  
  res[["tax.richness"]][["boxplot"]] <- pbox_rich
  res[["tax.richness"]][["barplot"]] <- pbar_rich
  res[["tax.alphaindex"]][["boxplot"]] <- pbox_alpha
  res[["tax.alphaindex"]][["barplot"]] <- pbar_alpha
  if (is.null(x$rarefy)) 
    return(res)
  
  ## Rarefaction curves by Groups
  RN <- names(x$alpha.diver[, c(1,2,4,6,7,8)])
  plotDF <- x$rarefy$grouping
  plotDF <- data.frame(sample = rep(plotDF$sample, 6), 
                       seq_dep = rep(plotDF$seq_dep, 6), 
                       mean = unlist(plotDF[, 3:8]), 
                       se = unlist(plotDF[, 9:14]), 
                       index = factor(rep(RN, each=nrow(plotDF)), RN))
  plotDF$lower <- plotDF$mean - plotDF$se
  plotDF$upper <- plotDF$mean + plotDF$se
  pR_G <- ggplot(plotDF, aes_string(x = "seq_dep", y = "mean", 
                                    colour = "sample"), na.rm = TRUE) + 
    geom_line() + geom_point() + 
    geom_errorbar(aes_string(ymin = "lower", ymax = "upper"), width = .25, 
                  stat = "identity", position = position_dodge(.9)) + 
    scale_colour_manual(values = color) + 
    labs(y = "Index", x = "") + theme_bw() + 
    facet_wrap( ~ index, scales = "free", drop = TRUE)

  ## Rarefaction curves by Samples
  pR_S <- list()
  plotDF <- x$rarefy$sample
  NM <- length(unique((plotDF$sample)))
  plotDF <- merge(x$rarefy$sample, opts$info$sample, by.x = "sample", by.y = 0, 
                  all.x = TRUE, sort = FALSE, suffixes = c(".R", ".M"))
  for (k in RN) {
    mean <- paste(k, "mean", sep = ".")
    se <- paste(k, "se", sep = ".")
    upper <- paste(k, "upper", sep = ".")
    lower <- paste(k, "lower", sep = ".")
    plotDF[[lower]] <- plotDF[[mean]] - plotDF[[se]]
    plotDF[[upper]] <- plotDF[[mean]] + plotDF[[se]]
    pR_S[[k]] <- ggplot(plotDF, aes_string(x = "seq_dep", y = mean, 
                                           colour = "colorVar"), na.rm = TRUE) + 
      geom_line() + geom_point() + 
      geom_errorbar(aes_string(ymin = lower, ymax = upper), width = .25, 
                    stat = "identity", position = position_dodge(.9)) + 
      scale_colour_manual(values = color) + 
      labs(y = k, x = "") + theme_bw() + 
      facet_wrap( ~ sample, nrow = ceiling(sqrt(NM)), 
                  scales = "free", drop = TRUE)
  }
  res[["tax.rarefy"]][["group"]] <- pR_G
  res[["tax.rarefy"]][["sample"]] <- pR_S
  return(res)
}
################################################################################


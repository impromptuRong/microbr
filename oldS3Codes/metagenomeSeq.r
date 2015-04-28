######################################################
#' MetagenomeSeq Statistical Analysis based on Nature Method
#' Lots of compatibility problem for now. 
#' Not good result on real data, needs further evaluation
library(metagenomeSeq)

########   Read in Data and create MRecperiment Object  ########
##   Loading count data  ##
dataDirectory <- system.file("extdata", package = "metagenomeSeq")
lung = load_meta(file.path(dataDirectory, "CHK_NAME.otus.count.csv"))
dim(lung$counts)
##    Loading Taxonomy   ##
taxa = read.delim(file.path(dataDirectory, "CHK_otus.taxonomy.csv"), stringsAsFactors = F)[, 2]
otu = read.delim(file.path(dataDirectory, "CHK_otus.taxonomy.csv"), stringsAsFactors = F)[, 1]
##    Loading metadata   ##
clin = load_phenoData(file.path(dataDirectory, "CHK_clinical.csv"), tran = TRUE)
ord = match(colnames(lung$counts), rownames(clin))
clin = clin[ord, ]
head(clin[1:2, ])
##  MRexperiment object  ##
phenotypeData = as(clin, "AnnotatedDataFrame")
phenotypeData


########   Sample MRexperiment Object   ########
data(lungData)
lungData
OTUdata = as(lung$taxa, "AnnotatedDataFrame")
varLabels(OTUdata) = "taxa"
OTUdata

obj = newMRexperiment(lung$counts,phenoData=phenotypeData,featureData=OTUdata)
# Links to a paper providing further details can be included optionally.
# experimentData(obj) = annotate::pmid2MIAME("21680950")
obj

########   Normlization factors and Normilzation   ########
##    Calculating normalization factors    ##
data(lungData)
p = cumNormStat(lungData)
# To calculate the scaling factors we simply run cumNorm
lungData = cumNorm(lungData, p = p)

##    Exporting data    ##
mat = MRcounts(lungData, norm = TRUE, log = TRUE)[1:5, 1:5]
exportMat(mat, output = file.path(dataDirectory, "tmp.tsv"))
exportStats(lungData[, 1:5], output = file.path(dataDirectory, "tmp.tsv"))
head(read.csv(file = file.path(dataDirectory, "tmp.tsv"), sep = "nt"))

########   Statistical testing   ########
##    Fit Zero-inflated Gaussian mixture model    ##
controls = grep("Extraction.Control", pData(lungData)$SampleType)
lungTrim = lungData[, -controls]
sparseFeatures = which(rowSums(MRcounts(lungTrim) > 0) < 10)
lungTrim = lungTrim[-sparseFeatures, ]
lungp = cumNormStat(lungTrim, pFlag = TRUE, main = "Trimmed lung data")
lungTrim = cumNorm(lungTrim, p = lungp)

smokingStatus = pData(lungTrim)$SmokingStatus
bodySite = pData(lungTrim)$SampleType
mod = model.matrix(~smokingStatus + bodySite)
settings = zigControl(maxit = 10, verbose = TRUE)
fit = fitZig(obj = lungTrim, mod = mod, control = settings)

##    Exporting fits    ##
taxa = sapply(strsplit(as.character(fData(lungTrim)$taxa), split = ";"), function(i){i[length(i)]})
head(MRcoefs(fit, taxa = taxa, coef = 2))

##    Testing presence/absence    ##
data(mouseData)
classes = pData(mouseData)$diet
res = MRfisher(mouseData[1:5, ], cl = classes)
# Warning - the p-value is calculating 1 despite a high odd's
# ratio.
head(res)


########   Visualization of features   ########
##    Structural overview    ##
data(mouseData)
trials = pData(mouseData)$diet
heatmapColColors = brewer.pal(12, "Set3")[as.integer(factor(trials))]
heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)
# plotMRheatmap
plotMRheatmap(obj = mouseData, n = 200, cexRow = 0.4, cexCol = 0.4, trace = "none", col = heatmapCols, ColSideColors = heatmapColColors)
# plotCorr
plotCorr(obj = mouseData, n = 200, cexRow = 0.25, cexCol = 0.25, trace = "none", dendrogram = "none", col = heatmapCols)

data(mouseData)
cl = factor(pData(mouseData)$diet)
# plotOrd
plotOrd(mouseData, tran = TRUE, usePCA = FALSE, useDist = TRUE, bg = cl, pch = 21, xlab = "1st coordinate", ylab = "2nd coordinate")
# plotRare
res = plotRare(mouseData, cl = cl, ret = TRUE, pch = 21, bg = cl)

tmp = lapply(levels(cl), function(lv) lm(res[, "ident"] ~ res[, "libSize"] - 1, subset = cl == lv))
for (i in 1:length(levels(cl))){
    abline(tmp[[i]], col = i)
}
legend("topleft", c("Diet 1", "Diet 2"), text.col = c(1, 2), box.col = NA)


##    Feature specific    ##
head(MRtable(fit, coef = 2, taxa = 1:length(fData(lungTrim)$taxa)))
patients = sapply(strsplit(rownames(pData(lungTrim)), split = "_"), function(i){i[3]})
pData(lungTrim)$patients = patients
classIndex = list(smoker = which(pData(lungTrim)$SmokingStatus == "Smoker"))
classIndex$nonsmoker = which(pData(lungTrim)$SmokingStatus == "NonSmoker")
otu = 779
# plotOTU
plotOTU(lungTrim, otu = otu, classIndex, main = "Neisseria meningitidis")
# Now multiple OTUs annotated similarly
x = fData(lungTrim)$taxa[otu]
otulist = grep(x, fData(lungTrim)$taxa)
# plotGenus
plotGenus(lungTrim, otulist, classIndex, labs = FALSE, main = "Neisseria meningitidis")
lablist <- c("S", "NS")
axis(1, at = seq(1, 6, by = 1), labels = rep(lablist, times = 3))

citation("metagenomeSeq")
sessionInfo()






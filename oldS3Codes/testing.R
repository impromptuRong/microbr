################################################################################
# cap.default, cap.dist, cap.physet
x <- as.matrix(otu_table(oral))/seqdep(oral)
ds <- sample_data(oral)
df <- data.frame(ds[,c("Group", "Age")], x)
dis <- phydist(x)

# cap.default
phycap(x, . ~ Group + Age, data = ds, method = "bray", comm = NULL)
phycap(dis, . ~ Group, data = ds, comm = x)
phycap(x, . ~ Group + Age, data = ds, method = "bray")
phycap(dis, . ~ Group + Age, data = df)
phycap(df, . ~ Group + Age)
phycap(x, . ~ 1, data = ds, method = "bray")
phycap(dis, . ~ 1, comm = x)
phycap(x, . ~ sample_data(oral)$Group, method = "bray")


x <- as.matrix(edge_mat(oral))/seqdep(oral)
bw <- edge_len(oral)
cap(x, . ~ Group + Age, data = ds, method = "unifrac.w.un", dfun = phydist, bw = bw)

phycap(oral, otu_table ~ Group)
phycap(oral, . ~ Age)
ds$Group2 <- ds$Group
phycap(oral, . ~ Group2, data = ds)
phycap(oral, edge_mat ~ Group2, data = ds)
phycap(oral, otu_table ~ Group, comm = x)
phycap(oral, . ~ 1, data = ds, comm = x, method = "euclidean")

# cca.default, cca.physet, rda.default, rda.physet
x <- as.matrix(otu_table(oral))
ds <- sample_data(oral)
df <- data.frame(ds[, c("Group", "Age")], x)
dis <- phydist(x)

phycca(x, . ~ sample_data(oral)$Group)
phycca(x, . ~ Group, data = ds)
phycca(df, . ~ Group + Age)
phycca(x, . ~ 1)
phycca(oral, otu_table ~ Group)
phycca(oral, . ~ Age)
ds$Group2 <- ds$Group
phycca(oral, . ~ Group2, data = ds)
phycca(oral, edge_mat ~ Group2, data = ds)

# mds.default, mds.dist, mds.physet
phymds(x)
phymds(x, . ~ 1, comm = edge_mat(oral))
phymds(x, . ~ 1, method = "bray", comm = edge_mat(oral))
phymds(edge_mat(oral), method = "unifrac.uw", dfun = phydist, bw = edge_len(oral))
phymds(dis)
phymds(dis, comm = edge_mat(oral), expand = FALSE)
phymds(dis, . ~ 1, method = "bray", comm = edge_mat(oral))
phymds(oral)
phymds(oral, method = "unifrac.w.un")
phymds(oral, edge_mat ~ 1, method = "unifrac.uw")
phymds(oral, . ~ 1, method = "bray", comm = edge_mat(oral))


ordinate(dis, "NMDS", formula = . ~ abc, comm = edge_mat(oral), expand = FALSE)
ordinate(dis, "NMDS", comm = edge_mat(oral), expand = FALSE)
ordinate(dis, "CAP", . ~ Group, data = ds, comm = x)
ordinate(dis, "NMDS", data = ds, comm = edge_mat(oral), expand = FALSE, k = 3)
ordinate(oral, "CCA", otu_table ~ Group)
ordinate(oral, "PCoA", method = "unifrac.uw")
ordinate(oral, "CAP", edge_mat ~ Group, distance = "unifrac.w.un")



data(oral)
colorVar <- list(sample = "Group", taxa = "Domain")
shapeVar <- list(sample = "Periodontitis")
fillVar <- colorVar
opts <- phyplotOptions(colorVar = colorVar, shapeVar = shapeVar, fillVar = fillVar)
ord <- ordinate(oral, "NMDS", method = "unifrac.w.un")
res <- plot(ord, oral, opts)
res$biplot
res$splitplot
ordinate(oral, "CCA", formula=otu_table ~ Group, methods = "unifrac.w.un")

data(oral)
opts <- phyplotOptions(colorVar = list(sample = "Group", taxa = "Domain"), 
                       shapeVar = list(sample = "Periodontitis"), 
                       fillVar = list(sample = "Group", taxa = "Domain"))
res <- ordinate_analysis(oral, . ~ Group, "unifrac.w.un", opts)

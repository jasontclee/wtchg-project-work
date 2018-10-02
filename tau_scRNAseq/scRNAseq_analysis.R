# The Seurat tutorial is provided at "https://satijalab.org/seurat/immune_alignment.html"
# We'll start by loading the packages needed for this analysis

library(Seurat)

# Setting up the working directory and loading the h5 files

setwd("/home/jlee631/todd_labwork/tau_singlecell/5_tau_cell_metrics")

KO <- Read10X_h5("KO_filtered_gene_bc_matrices_h5.h5", ensg.names = FALSE)
H1 <- Read10X_h5("H1_filtered_gene_bc_matrices_h5.h5", ensg.names = FALSE)

# Conversion of the h5 files to the Seurat format for use in the Seurat pipeline
# First, we create the Seurat object from the h5 file

KO <- CreateSeuratObject(raw.data = KO, min.cells = 3, min.genes = 0, project = "TauscRNAseqKO")
KO@meta.data$cond <- "KO"

# We will use nGene, nUMI and also percentage of mitochondrial genes present as a way to exclude doublets and dying cells

mito_genes_KO <- grep(pattern = "^mt-", x = rownames(x = KO@data), value = TRUE)
percent_mito_KO <- Matrix::colSums(KO@raw.data[mito_genes_KO, ])/Matrix::colSums(KO@raw.data)
KO <- AddMetaData(object = KO, metadata = percent_mito_KO, col.name = "percent.mito")

# Violin Plot allows us to visualize what the distribution is before filtering

KO_pre <- VlnPlot(object = KO, features.plot = c("nGene", "nUMI", "percent.mito"))
KO_pre_l <- length(KO@meta.data$nGene)

# Now we do this for the H1 dataset

H1 <- CreateSeuratObject(raw.data = H1, min.cells = 3, min.genes = 0, project = "TauscRNAseqH1")
H1@meta.data$cond <- "H1"

# We will use nGene, nUMI and also percentage of mitochondrial genes present as a way to exclude doublets and dying cells

mito_genes_H1 <- grep(pattern = "^mt-", x = rownames(x = H1@data), value = TRUE)
percent_mito_H1 <- Matrix::colSums(H1@raw.data[mito_genes_H1, ])/Matrix::colSums(H1@raw.data)
H1 <- AddMetaData(object = H1, metadata = percent_mito_H1, col.name = "percent.mito")

# Violin Plot allows us to visualize what the distribution is before filtering

H1_pre <- VlnPlot(object = H1, features.plot = c("nGene", "nUMI", "percent.mito"))
H1_pre_l <- length(H1@meta.data$nGene)

# Filtering steps

KO <- FilterCells(KO, subset.names = c("nGene", "nUMI", "percent.mito"), low.thresholds = c(1200, -Inf, -Inf), high.thresholds = c(5500, 60000, 0.025))
KO_post <- VlnPlot(object = KO, features.plot = c("nGene", "nUMI", "percent.mito"))
KO_post_l <- length(KO@meta.data$nGene)

H1 <- FilterCells(H1, subset.names = c("nGene", "nUMI", "percent.mito"), low.thresholds = c(1200, -Inf, -Inf), high.thresholds = c(5500, 60000, 0.025))
H1_post <- VlnPlot(object = H1, features.plot = c("nGene", "nUMI", "percent.mito"))
H1_post_l <- length(H1@meta.data$nGene)

# Data normalization after filtering steps

KO <- NormalizeData(KO)
KO <- ScaleData (KO, display.progress = F)

H1 <- NormalizeData(H1)
H1 <- ScaleData (H1, display.progress = F)

# Variable genes in the set are identified

KO <- FindVariableGenes(KO, do.plot = F)
H1 <- FindVariableGenes(H1, do.plot = F)
g.1 <- head(rownames(KO@hvg.info),1000)
g.2 <- head(rownames(H1@hvg.info),1000)
genes.use <- unique(c(g.1,g.2))
genes.use <- intersect(genes.use, rownames(KO@scale.data))
genes.use <- intersect(genes.use, rownames(H1@scale.data))

# A CCA is performed, similar to PCA but with both groups

combined <- RunCCA(KO, H1, add.cell.id1 = "KO", add.cell.id2 = "H1", genes.use = genes.use, num.cc = 30)

p1 <- DimPlot(object = combined, reduction.use = "cca", group.by = "cond")
p2 <- VlnPlot(object = combined, features.plot = "CC1", group.by = "cond")
p3 <- plot_grid(p1,p2)

png(filename = "CCA_DimVln_Plot.png")
p3
dev.off()

# This graph shows us a plot of how many dimensions we will need to use for the analysis

PrintDim(object = combined, reduction.type = "cca", dims.print = 1:2, genes.print = 10)

p1 <- MetageneBicorPlot(combined, grouping.var = "cond", dims.eval = 1:30)

png(filename = "tSNE_dim.png")
p1
dev.off()

p3 <- DimHeatmap(object = combined, reduction.type = "cca", cells.use = 500, dim.use = 1:9, do.balanced = TRUE)

png(filename = "DimHeatmap.png")
p3
dev.off()

combined <- AlignSubspace(combined, reduction.type = "cca", grouping.var = "cond", dims.align = 1:15)

p1 <- VlnPlot(object = combined, features.plot = "ACC1", group.by = "cond")
p2 <- VlnPlot(object = combined, features.plot = "ACC2", group.by = "cond")
p3 <- plot_grid(p1, p2)

png(filename = "Acc_Vln_Plot.png")
p3
dev.off()

combined <- RunTSNE(combined, reduction.use = "cca.aligned", dims.use = 1:15, do.fast = T)
combined <- FindClusters(combined, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:15)

# Plots the tSNE graphs that show you the clusters for this dataset

p1 <- TSNEPlot(combined, do.return = T, pt.size = 0.5, group.by = "cond")
p2 <- TSNEPlot(combined, do.label = T, do.return = T, pt.size = 0.5)
p3 <- plot_grid(p1, p2)

png(filename = "tSNE_clusters_dualpanel.png")
p3
dev.off()

png(filename = "tSNE_clusters_pre.png")
p2
dev.off()

# Looping through the 10 different cluster groups that we obtained and finding the top genes that are most significant
# Already wrote the csv files, only delete the hashtags if you need to re-write them to the drive

for (i in 0:9) {
  v <- c("group_0","group_1","group_2","group_3","group_4","group_5","group_6","group_7","group_8","group_9")
  a <- FindConservedMarkers(combined, ident.1 = i, grouping.var = "cond", print.bar = FALSE)
  write.csv(a, file = v[i+1])
}

# Can use this command line to search for the genes of interest plotted against the cluster that it falls within
# The first 8 genes are the top significant hits for the group, anything else that follows is just things that come to mind

FeaturePlot(object = combined, features.plot = c("Nupr1", "Gcg", "Ins1", "Sst","S100a4","Plvap","Ppy","Stmn1", "Ctrb1"), min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)

# Labelling the graph with the new identities

new.ident <- c("beta1", "alpha", "beta2", "delta", "beta3", "endothelial", "pp", "acinar", "proliferating.undefined","leukocyte")
for (i in 0:9) {
  combined <- RenameIdent(object = combined, old.ident.name = i, new.ident.name = new.ident[i+1])
}

# Plotting the tSNE again, this time with the new labels

p3 <- TSNEPlot(combined, do.label = T, pt.size = 0.5)

png(filename = "tSNE_clusters_post.png")
p3
dev.off()

# Split plot graph to show expression levels between groups, change the variables if you change the name of the cluster groups

combined@ident <- factor(combined@ident, levels = (c("beta1", "alpha", "beta2", "delta", "beta3", "endothelial", "pp", "acinar", "proliferating.undefined","leukocyte")))
markers.to.plot <- c("Ins1", "Ins2", "Gcg", "Sst", "MAPT", "Nupr1", "Ttr", "Cenpf", "Plvap", "Stmn1", "Cd74")
sdp <- SplitDotPlotGG(combined, genes.plot = rev(markers.to.plot), cols.use = c("blue", "red"), x.lab.rot = T, plot.legend = T, dot.scale = 8, do.return = T, grouping.var = "cond")

# Next part is massive
# This part doesn't need to be changed, it's just formatting the plot

LabelPoint <- function(plot, genes, exp.mat, adj.x.t = 0, adj.y.t = 0, adj.x.s = 0, 
                       adj.y.s = 0, text.size = 2.5, segment.size = 0.1) {
  for (i in genes) {
    x1 <- exp.mat[i, 1]
    y1 <- exp.mat[i, 2]
    plot <- plot + annotate("text", x = x1 + adj.x.t, y = y1 + adj.y.t, 
                            label = i, size = text.size)
    plot <- plot + annotate("segment", x = x1 + adj.x.s, xend = x1, y = y1 + 
                              adj.y.s, yend = y1, size = segment.size)
  }
  return(plot)
}

LabelUR <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.r.t = 0.15, adj.u.s = 0.05, 
                    adj.r.s = 0.05, ...) {
  return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = adj.r.t, 
                    adj.y.s = adj.u.s, adj.x.s = adj.r.s, ...))
}

LabelUL <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05, 
                    adj.l.s = 0.05, ...) {
  return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = -adj.l.t, 
                    adj.y.s = adj.u.s, adj.x.s = -adj.l.s, ...))
}

# This part we are changing

beta1 <- SubsetData(combined, ident.use = "beta1", subset.raw = T)
beta1 <- SetAllIdent(beta1, id = "cond")
avgbeta1 <- log1p(AverageExpression(beta1, show.progress = FALSE))
avgbeta1$gene <- rownames(avgbeta1)

beta2 <- SubsetData(combined, ident.use = "beta2", subset.raw = T)
beta2 <- SetAllIdent(beta2, id = "cond")
avgbeta2 <- log1p(AverageExpression(beta2, show.progress = FALSE))
avgbeta2$gene <- rownames(avgbeta2)

genes.to.label1 = c("Gcg", "Ppy", "Sst", "Pyy", "Nupr1")
genes.to.label2 = c("Ins1", "Ins2")
genes.to.label3 = c("Ttr", "Fau")

p1 <- ggplot(avgbeta1, aes(H1, KO)) + geom_point() + ggtitle("Beta 1 DE")
p1 <- LabelUR(p1, genes = c(genes.to.label1, genes.to.label2), avgbeta1, 
              adj.u.t = 0.3, adj.u.s = 0.23)
p1 <- LabelUL(p1, genes = genes.to.label3, avgbeta1, adj.u.t = 0.5, adj.u.s = 0.4, 
              adj.l.t = 0.25, adj.l.s = 0.25)

p2 <- ggplot(avgbeta2, aes(H1, KO)) + geom_point() + ggtitle("Beta 2 DE")
p2 <- LabelUR(p2, genes = c(genes.to.label1, genes.to.label3), avgbeta2, 
              adj.u.t = 0.3, adj.u.s = 0.23)
p2 <- LabelUL(p2, genes = genes.to.label2, avgbeta2, adj.u.t = 0.5, adj.u.s = 0.4, 
              adj.l.t = 0.25, adj.l.s = 0.25)
plot_grid(p1, p2)

# Finding the differentially expressed genes
# ONLY RUN THIS ONCE CAUSE OTHERWISE YOU'LL JUST KEEP TACKING THE LABEL AT THE END :)

combined@meta.data$celltype.cond <- paste0(combined@ident, "_", combined@meta.data$cond)
combined <- StashIdent(combined, save.name = "celltype")
combined <- SetAllIdent(combined, id = "celltype.cond")

# Run this to find the differentially expressed genes

beta1_response <- FindMarkers(combined, ident.1 = "beta1_H1", ident.2 = "beta1_KO", print.bar = FALSE)
write.csv(beta1_response, file = "beta1_DEG.csv")

alpha_response <- FindMarkers(combined, ident.1 = "alpha_H1", ident.2 = "alpha_KO", print.bar = FALSE)
write.csv(alpha_response, file = "alpha_DEG.csv")

beta2_response <- FindMarkers(combined, ident.1 = "beta2_H1", ident.2 = "beta2_KO", print.bar = FALSE)
write.csv(beta2_response, file = "beta2_DEG.csv")

delta_response <- FindMarkers(combined, ident.1 = "delta_H1", ident.2 = "delta_KO", print.bar = FALSE)
write.csv(delta_response, file = "delta_DEG.csv")






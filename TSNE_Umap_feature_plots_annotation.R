library(tidyverse)
library(Seurat)
library(SeuratData)
library(ggplot2)

#load data
cluster_annotation <- readr::read_csv("annotated_clusters.csv")
seurat_object <- readr::read_rds("seurat_object.rds")

sample_name = "NSCLC - scRNA-seq (GSE148071)"

#######################
#run tSNE
seurat_object <- RunTSNE(seurat_object, reduction = "pca", dims = 1:20, perplexity=1000)

plot_name = paste(sample_name,"_tSNE_celltype.png",sep="")

p <- DimPlot(seurat_object, reduction = "tsne", label = TRUE, pt.size = 3)
p + labs(title = "TSNE, dims:20, perplexity:1000")

ggsave(plot_name, width=12, height=10, dpi="retina", path="plots")
dev.off()

#UMAPS

#run UMAP
seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:20, n.components = 2, n.neighbors = 30,
                 n.epochs = 300, min.dist = 1, learning.rate = 5, spread = 2)

p <- DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 2)
p + labs(title = "TSNE Lung: n=30,ep=300,dist=1,rate=5, spread=2")

#annotate UMAP clusters
celltype <- cluster_annotation$`Cell Type`
new.cluster.ids <- celltype
names(new.cluster.ids) <- levels(seurat_object)
seurat_object <- RenameIdents(seurat_object, new.cluster.ids)
seurat_object$Celltype <- seurat_object@active.ident

#run UMAP again after annotation
seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:20, n.components = 2, n.neighbors = 30,
                 n.epochs = 300, min.dist = 1, learning.rate = 5, spread = 2)

#plot UMAP by sample source (orig.identity)
plot_name = paste(sample_name,"_UMAP_group.png",sep="")

p <- DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 2, group.by = "orig.ident")
p + labs(title = "UMAP: n=30,ep=300,dist=1,rate=5, spread=2, group")
ggsave(plot_name, width=12, height=10, dpi="retina", path="plots")
dev.off()

#plot UMAP by cluster id
plot_name = paste(sample_name,"_UMAP_clusters.png",sep="")

p <- DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 2, group.by = "seurat_clusters")
p + labs(title = "UMAP: n=30,ep=300,dist=1,rate=5, spread=2, cluster")
ggsave(plot_name, width=12, height=10, dpi="retina", path="plots")
dev.off()

#plot UMAP by Cell Type
plot_name = paste(sample_name,"_UMAP_celltype.png",sep="")

p <- DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 2)
p + labs(title = "UMAP: n=30,ep=300,dist=1,rate=5, spread=2, celltype")
ggsave(plot_name, width=12, height=10, dpi="retina", path="plots")
dev.off()

#######################
#feature plots
plot_name = paste(sample_name, "_feature_plot_umap.png", sep="")

features <- c("AKR1C1", "PERP", "BIRC5", "AKR1C1", "MS4A6A", "CTSL", "TRBC2", "JCHAIN",
              "SERPINB13", "SCEL", "PLK1", "SCGB3A1", "DCN", "LPL", "VWF", "RGS5")


p <- FeaturePlot(seurat_object, features = features)
p + labs(width=1200, height = 1000) & theme(plot.title = element_text(size=10))
#p + labs(title = "Feature Plot Mouse Lung: 33 vs 25 (UMAP)")
ggsave(plot_name, width=12, height=10, dpi="retina", path="plots", device="png")
dev.off()

#ridge plot
plot_name = paste(sample_name, "_feature_plot_ridge.png", sep="")

p <- RidgePlot(seurat_object, features=features, ncol=2)
p + labs(width=1200, height = 1000) & theme(plot.title = element_text(size=10))
#p + labs(title = "Feature Plot Mouse Lung: 33 vs 25 (Ridge)")

ggsave(plot_name, width=12, height=30, dpi="retina", path="plots")
dev.off()

#violin plot
plot_name = paste(sample_name, "_feature_plot_violin.png", sep="")

p <- VlnPlot(seurat_object,features=features)
p + labs(width=1200, height = 1000) & theme(plot.title = element_text(size=10))
#p + labs(title = "Violin Plot Mouse Lung: 33 vs 25 (Violin)")

ggsave(plot_name, width=12, height=22, dpi="retina", path="plots")
dev.off()

#######################
#heatmaps

library(gplots)

#heatmap for cluster by group
plot_name = paste(sample_name, "_heatmap_group.png", sep="")
p <- DoHeatmap(subset(seurat_object, downsample = 100), features = features, size = 3, group.by = "orig.ident")
p + labs(title = plot_name)
ggsave(plot_name, width=10, height=7, dpi="retina", path="plots")
dev.off()

#heatmap for cluster by cell type (with legend)
plot_name = paste(sample_name, "_heatmap_celltype.png", sep="")
p <- DoHeatmap(subset(seurat_object, downsample = 100), features = features, size = 3, group.by='Celltype')
p + labs(title = plot_name)
ggsave(plot_name, width=10, height=7, dpi="retina", path="plots")
dev.off()

#heatmap for cluster by cell type (no legend)
plot_name = paste(sample_name, "_heatmap1.png", sep="")
p <- DoHeatmap(seurat_object, features = VariableFeatures(seurat_object)[1:100], cells = 1:500, size = 4,
          angle = 90)
p + labs(title = plot_name)
ggsave(plot_name, width=10, height=16, dpi="retina", path="plots")
dev.off()

#Complex heatmap for cluster by cell type and cluster
library(Scillus)

plot_heatmap(dataset = seurat_object,
             n = 2,
             markers = VariableFeatures(seurat_object)[1:100],
             sort_var = c("Celltype","Group"),
             anno_var = c("Celltype", "Sample"),
             anno_colors = list("Set2","Paired"),
             hm_limit = c(-1,0, 1),
             hm_colors = c("purple","black","yellow"))


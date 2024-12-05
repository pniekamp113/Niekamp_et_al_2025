#load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(tibble)
library(RColorBrewer)
library(monocle3)
library(SeuratWrappers)

setwd() #set working directory

# Subsetting of CD8a clusters only and re-clustering
# Select CD8 clusters
Idents(combined) <- "tree.ident"
clusters_to_subset <- c("7", "8", "9", "11")
cells_to_keep <- which(Idents(combined) %in% clusters_to_subset)
Cd8a_subset <- subset(combined, cells = cells_to_keep)

#Re-clusters CD8 cluster
Cd8a_subset <- NormalizeData(Cd8a_subset, layer = "counts.Gene Expression")
Cd8a_subset <- FindVariableFeatures(Cd8a_subset, layer = "counts.Gene Expression", nfeatures = 2000)
Cd8a_subset <- ScaleData(Cd8a_subset)
Cd8a_subset <- RunPCA(Cd8a_subset)
Cd8a_subset <- RunUMAP(Cd8a_subset,dims = 1:40)
Cd8a_subset <- FindNeighbors(Cd8a_subset, dims = 1:40)
Cd8a_subset <- FindClusters(Cd8a_subset, resolution = 0.4)
Cd8a_subset <- BuildClusterTree(Cd8a_subset,reorder.numeric = T,reorder = T)

#check for successful selection of Cd8a+ cells
p1 <- DimPlot(Cd8a_subset, label = TRUE)
p2 <- DimPlot(Cd8a_subset, group.by = "HTO_maxID")
p3 <- FeaturePlot(Cd8a_subset, "Cd8a")
cowplot::plot_grid(p1, p2,p3)

#Finalize CD8a by excluding C5
Idents(Cd8a_subset) <- "tree.ident"
clusters_to_subset <- c("1", "2", "3", "4", "6", "7")
cells_to_keep <- which(Idents(Cd8a_subset) %in% clusters_to_subset)
Cd8a_final <- subset(Cd8a_subset, cells = cells_to_keep)

Cd8a_subset <- NormalizeData(Cd8a_final, layer = "counts.Gene Expression")
Cd8a_subset <- FindVariableFeatures(Cd8a_subset, layer = "counts.Gene Expression", nfeatures = 2000)
Cd8a_subset <- ScaleData(Cd8a_subset)
Cd8a_subset <- RunPCA(Cd8a_subset)
Cd8a_subset <- RunUMAP(Cd8a_subset,dims = 1:40)
Cd8a_subset <- FindNeighbors(Cd8a_subset, dims = 1:40)
Cd8a_subset <- FindClusters(Cd8a_subset, resolution = 0.4)
Cd8a_subset <- BuildClusterTree(Cd8a_subset,reorder.numeric = T,reorder = T)

#plot new clusters
p1 <- DimPlot(Cd8a_subset, label = TRUE)
p2 <- DimPlot(Cd8a_subset, group.by = "HTO_maxID")
p3 <- FeaturePlot(Cd8a_subset, "Cd8a")
cowplot::plot_grid(p1, p2,p3)

#Assign labels and re-name clusters
hto_labels <- c(
  "TotalSeq-C0301-HTO1" = "WT",
  "TotalSeq-C0302-HTO2" = "RARa-KO",
  "TotalSeq-C0303-HTO3" = "RARa-TG"
)
labels <- hto_labels[Cd8a_subset$HTO_maxID]
names(labels) <- colnames(Cd8a_subset)
Cd8a_subset <- AddMetaData(Cd8a_subset, metadata = labels, col.name = "labels")

# Rename clusters
new_cluster_names <- c(
 "1" = "Teff3",
 "2" = "Early Teff",
 "3" = "Naive-Like",
 "4" = "Teff2",
 "5" = "Teff1",
 "6" = "Proliferating Teff"
)

# Re-order the factor levels
levels(Cd8a_subset) <- c("Proliferating Teff", "Teff1", "Teff3", "Teff2", "Early Teff", "Naive-Like")


#Figure XA+B
p1 <- DimPlot(Cd8a_subset, label = TRUE) + NoLegend()
p2 <- DimPlot(Cd8a_subset, group.by = "labels") + NoLegend()
cowplot::plot_grid(p1, p2)

saveRDS(Cd8a_subset, file "Cd8a_subset.rds")

#Heatmap for signature gene comparison (for Fig XC)
gene_list <- read.csv("signature_genes.csv")
genes <- gene_list$gene
genes_in_data <- intersect(genes, rownames(Cd8a_subset))
all_markers <- FindAllMarkers(Cd8a_subset, 
                              only.pos = FALSE, 
                              logfc.threshold = 0, 
                              min.pct = 0,
                              return.thresh = Inf)

filtered_results <- all_markers[all_markers$gene %in% genes_in_data, ]

data_subset <- filtered_results %>%
  select(gene, avg_log2FC, cluster)

data_pivoted <- data_subset %>%
  pivot_wider(names_from = cluster, values_from = avg_log2FC)

write.csv(data_pivoted, file = "DEGs_for_heatmap_sorted.csv")


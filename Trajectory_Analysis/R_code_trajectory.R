# Trajectory code using Monocle3 for manuscript

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


# Load Cd8a subset seurat object
Cd8a_subset <- readRDS("Cd8a_subset.rds")

#set unused layers to NULL for SeuratWrapper function
Cd8a_subset[["RNA"]]@layers$`counts.Antibody Capture` <- NULL
Cd8a_subset[["RNA"]]@layers$`data.Antibody Capture` <- NULL

# run SeuratWrappers::as.cell_data_set
cds <- SeuratWrappers::as.cell_data_set(Cd8a_subset)

# run Monocle3 to predict trajectory
cds <- cluster_cells(cds)
plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by = "partition")
cds <- learn_graph(cds, use_partition = FALSE)

# identify root cluster (Naive-like = cluster 3)
root_cells <- which(cds$tree.ident == 3)
root_cells_ids <- colnames(cds)[root_cells]
cds <- order_cells(cds, root_cells = root_cells_ids)

# Figure SX
p1 <- DimPlot(Cd8a_subset) + NoLegend()
p2 <- plot_cells(cds, color_cells_by = "pseudotime", label_branch_points = FALSE, label_leaves = FALSE, show_trajectory_graph = FALSE, cell_size = 1) + NoLegend()
cowplot::plot_grid(p1, p2)

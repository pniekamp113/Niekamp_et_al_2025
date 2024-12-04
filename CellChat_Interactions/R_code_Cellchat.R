# Cellchat code for manuscript

#Load required libraries
library(Seurat)
library(patchwork)
library(CellChat)


# Load seurat object
seurat_obj <- readRDS("scRNA-seq-object.rds")


# Subset without apoptotic cell clusters
clusters_to_subset <- c("C1 CD4", "C2 Naive-Like", "C3 Early Teff", "C4 Teff 3", "C5 Teff 1", "C6 Teff 2", "C7 Proliferating Teff", "C8 Myeloid", "C9 MC38 cells")
cells_to_keep <- which(Idents(seurat_obj) %in% clusters_to_subset)
for_cell_chat <- subset(seurat_obj, cells = cells_to_keep)

# Create cellchat object
cellchat <- createCellChat(object = for_cell_chat)
CellChatDB <- CellChatDB.mouse 
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use

# Identify Genes for CellChat
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Project expression to PPI network to impute expression
cellchat <- projectData(cellchat, PPI.mouse)

# Compute Interaction Network  
cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)


# convert cellchat to DF and look up all possible pathways
df.cellchat <- subsetCommunication(cellchat)
#write as CSV file
write.csv(df.cellchat, "cellchat_LRpairs.csv")

# compute centrality
cellchat <- netAnalysis_computeCentrality(cellchat)

#save CellChat object for later use.
saveRDS(cellchat, file = "cellchat_results.rds")
cellchat <- readRDS("cellchat_results.rds")

#select source and target cells / can be removed if all cell clusters are vizualized 
source <- c("C1 CD4", "C2 Naive-Like", "C3 Early Teff", "C4 Teff 3", "C5 Teff 1", "C6 Teff 2", "C7 Proliferating Teff", "C8 Myeloid", "C9 MC38 cells")
target <- c("C1 CD4", "C2 Naive-Like", "C3 Early Teff", "C4 Teff 3", "C5 Teff 1", "C6 Teff 2", "C7 Proliferating Teff", "C8 Myeloid", "C9 MC38 cells")

# Results Figure SX
netVisual_individual(cellchat,  signaling = "CCL", pairLR.use = "CCL5_CCR1",  sources.use = source, targets.use = target, vertex.receiver = vertex.receiver, vertex.label.cex = 1.5)
netVisual_individual(cellchat,  signaling = "CCL", pairLR.use = "CCL3_CCR1",  sources.use = source, targets.use = target, vertex.receiver = vertex.receiver, vertex.label.cex = 1.5)
netVisual_individual(cellchat,  signaling = "CCL", pairLR.use = "CCL4_CCR5",  sources.use = source, targets.use = target, vertex.receiver = vertex.receiver, vertex.label.cex = 1.5)
netVisual_individual(cellchat,  signaling = "CCL", pairLR.use = "CCL3_CCR5",  sources.use = source, targets.use = target, vertex.receiver = vertex.receiver, vertex.label.cex = 1.5)
netVisual_individual(cellchat,  signaling = "CCL", pairLR.use = "CCL5_CCR5",  sources.use = source, targets.use = target, vertex.receiver = vertex.receiver, vertex.label.cex = 1.5)
netVisual_individual(cellchat,  signaling = "IFN-II", pairLR.use = "IFNG_IFNGR1_IFNGR2",  sources.use = source, targets.use = target, vertex.receiver = vertex.receiver,vertex.label.cex = 1.5)
netVisual_individual(cellchat,  signaling = "CXCL", pairLR.use = "CXCL16_CXCR6",  sources.use = source, targets.use = target, vertex.receiver = vertex.receiver, vertex.label.cex = 1.5)
netVisual_individual(cellchat,  signaling = "PD-L1", pairLR.use = "CD274_PDCD1",  sources.use = source, targets.use = target, vertex.receiver = vertex.receiver, vertex.label.cex = 1.5)
netVisual_individual(cellchat,  signaling = "PDL2", pairLR.use = "PDCD1LG2_PDCD1",  sources.use = source, targets.use = target, vertex.receiver = vertex.receiver, vertex.label.cex = 1.5)


# Results Figure SX 
pathways <- c("CCL", "CXCL", "PD-L1", "PDL2", "IL2", "TNF", "IFN-II")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", signaling = pathways, color.heatmap = "Purples", height = 3, width = 5, font.size = 8)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", signaling = pathways, color.heatmap = "Purples", height = 3, width = 5, font.size = 8)


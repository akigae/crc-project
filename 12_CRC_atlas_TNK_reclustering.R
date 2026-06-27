library(Seurat)
library(schard)
library(anndata)

# Load h5ad as Seurat----
snhx = schard::h5ad2seurat("~/projects/CRC_xenium/CAF_projext/validation_2/atlas/atlas_t.h5ad")

snhx <- subset(snhx, sample_type == "primary tumor")
snhx <- subset(snhx, microsatellite_status %in% c("MSI", "MSI-H", "MSS"))

# Re-clustering----
snhx <- NormalizeData(snhx)
snhx <- FindVariableFeatures(snhx)
snhx <- ScaleData(snhx)

snhx <- FindNeighbors(snhx, dims = 1:10, reduction = "XscANVI_")
snhx <- FindClusters(snhx, resolution = c(0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))

snhx <- RunUMAP(snhx, dims = 1:10, reduction = "XscANVI_")

# DEG analysis
Idents(snhx) <- "SCT_snn_res.0.8"
markers <- FindAllMarkers(snhx, only.pos = T, min.pct = 0.25)

# pick up only top10 significant genes (avg_log2FC > 0.6)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.6) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

top10$cluster <- factor(top10$cluster, levels = c(unique(top10$cluster)))
top10 <- top10[order(top10$cluster), ]

snhx_filtered <- ScaleData(snhx, features = top10$gene)
p <- DoHeatmap(snhx_filtered, features = top10$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()

# check the heatmap to check whether the clusters are well separated at this resolution
ggplot2::ggsave(filename = "snhx_t_0.8.png", plot = p, width = 20, height = 10)


VlnPlot(snhx, c("CD3D", "CD3E", "CD4", "CD40LG", "CD8A", "CD8B", "TCF7", "CCR7", "SELL", 
                "FOXP3", "CTLA4", "IL2RA", "IL17A", "CCL20", "KLRB1",　"MAML2", "KLF12", "SIK3",
                "ITGAE", "GZMA", "GZMB", "GZMK", "EOMES", "CD27", "GNLY", "ZNF683", "FCRL6",
                "TRDC", "TRDV1", "TRGC1", "MKI67", "TOP2A", "CCNB1",
                "FCGR3A", "TYROBP", "S1PR5","IL4I1", "KIT", "LST1"), 
        group.by = "RNA_snn_res.0.8", stack = T, flip = T, pt.size = 0, add.noise = T)


snhx$t_sub <- snhx$RNA_snn_res.0.8
snhx$t_sub <- gsub("^12$", "γδ T cell", snhx$t_sub)
snhx$t_sub <- gsub("^13$", "NK cell", snhx$t_sub)
snhx$t_sub <- gsub("^16$", "ILC", snhx$t_sub)
snhx$t_sub <- gsub("^11$", "Proliferating_T", snhx$t_sub)

snhx$t_sub <- gsub("^20$", "CD4T_Naive", snhx$t_sub)
snhx$t_sub <- gsub("^19$", "CD4T_Naive", snhx$t_sub)
snhx$t_sub <- gsub("^4$", "CD4T_Naive", snhx$t_sub)
snhx$t_sub <- gsub("^3$", "CD4T_Naive", snhx$t_sub)
snhx$t_sub <- gsub("^6$", "CD4T_Naive", snhx$t_sub)
snhx$t_sub <- gsub("^18$", "Treg", snhx$t_sub)
snhx$t_sub <- gsub("^7$", "Treg", snhx$t_sub)
snhx$t_sub <- gsub("^5$", "Treg", snhx$t_sub)
snhx$t_sub <- gsub("^1$", "CD4T_Th17", snhx$t_sub)

snhx$t_sub <- gsub("^17$", "CD8T_Naive", snhx$t_sub)
snhx$t_sub <- gsub("^14$", "CD8T_Trm", snhx$t_sub)
snhx$t_sub <- gsub("^8$", "CD8T_Trm", snhx$t_sub)
snhx$t_sub <- gsub("^0$", "CD8T_GZMK", snhx$t_sub)
snhx$t_sub <- gsub("^9$", "CD8T_GZMK", snhx$t_sub)
snhx$t_sub <- gsub("^2$", "CD8T_GZMK", snhx$t_sub)

snhx$t_sub <- gsub("^10$", "Activated_T", snhx$t_sub)
snhx$t_sub <- gsub("^15$", "Activated_T", snhx$t_sub)

# Subset NK cell & γδ T cell----
NK_sub <- subset(snhx, t_sub %in% c("NK cell", "γδ T cell"))

NK_sub <- NormalizeData(NK_sub)
NK_sub <- FindVariableFeatures(NK_sub)
NK_sub <- ScaleData(NK_sub)

NK_sub <- FindNeighbors(NK_sub, dims = 1:10, reduction = "XscANVI_")
NK_sub <- FindClusters(NK_sub, resolution = c(0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))

NK_sub <- RunUMAP(NK_sub, dims = 1:10, reduction = "XscANVI_")

NK_markers <- FindAllMarkers(NK_sub, only.pos = T, min.pct = 0.25)

VlnPlot(NK_sub, c("CD3D", "CD3E", "CD8A", "CD8B", 
                  "TRDC", "GNLY", "ZNF683", "TRDV1", "TRDV2", "FCGR3A"), 
        group.by = "RNA_snn_res.0.5", stack = T, flip = T, pt.size = 0, add.noise = T)

NK_sub$t_sub <- NK_sub$RNA_snn_res.0.5
NK_sub$t_sub <- gsub("^0$", "NK cell", NK_sub$t_sub)
NK_sub$t_sub <- gsub("^8$", "NK cell", NK_sub$t_sub)
NK_sub$t_sub <- gsub("^10$", "NK cell", NK_sub$t_sub)
NK_sub$t_sub <- gsub("^5$", "NK cell", NK_sub$t_sub)
NK_sub$t_sub <- gsub("^7$", "NK cell", NK_sub$t_sub)
NK_sub$t_sub <- gsub("^1$", "NK cell", NK_sub$t_sub)

NK_sub$t_sub <- gsub("^6$", "CD8T_Temra", NK_sub$t_sub)
NK_sub$t_sub <- gsub("^11$", "CD8T_Temra", NK_sub$t_sub)

NK_sub$t_sub <- gsub("^3$", "γδ T cell", NK_sub$t_sub)
NK_sub$t_sub <- gsub("^2$", "γδ T cell", NK_sub$t_sub)
NK_sub$t_sub <- gsub("^4$", "γδ T cell", NK_sub$t_sub)
NK_sub$t_sub <- gsub("^9$", "γδ T cell", NK_sub$t_sub)

DimPlot(NK_sub, reduction = "umap", label = T, group.by = "t_sub")

# Merge annotation data----
meta_1 <- snhx@meta.data["t_sub"]
meta_2 <- NK_sub@meta.data["t_sub"]

unique(meta_1$t_sub)
meta_1 <- subset(meta_1, t_sub %in% c("Activated_T", "CD8T_GZMK", "CD4T_Naive", "CD8T_Trm", "Treg", "Proliferating_T", 
                                      "CD4T_Th17", "ILC", "CD8T_Naive", "CD4T_Tex"))

merged_meta <- rbind(meta_1, meta_2)

snhx <- AddMetaData(snhx, merged_meta, col.name = "t_sub_2")
DimPlot(snhx, reduction = "umap", group.by = "t_sub_2")

# subset activated T----
act_sub <- subset(snhx, t_sub_2 == "Activated_T")

act_sub <- FindNeighbors(act_sub, dims = 1:10, reduction = "XscANVI_")
act_sub <- FindClusters(act_sub, resolution = c(0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))

act_sub <- RunUMAP(act_sub, dims = 1:10, reduction = "XscANVI_")
DimPlot(act_sub, reduction = "umap", group.by = "RNA_snn_res.0.2", label = T)
FeaturePlot(act_sub, features = "pct_counts_mito", reduction = "umap")

Idents(act_sub) <- "RNA_snn_res.0.2"
markers_2 <- FindAllMarkers(act_sub, min.pct = 0.25, only.pos = T)

VlnPlot(act_sub, c("CD3D", "CD3E", "CD4", "CD40LG", "CD8A", "CD8B", "TCF7", "CCR7", "SELL", 
                   "FOXP3", "CTLA4", "IL2RA", "IL17A", "CCL20", "KLRB1",　"MAML2", "KLF12", "SIK3",
                   "ITGAE", "GZMA", "GZMB", "GZMK", "EOMES", "CD27", "GNLY", "ZNF683", "FCRL6",
                   "TRDC", "TRDV1", "TRGC1", "MKI67", "TOP2A", "CCNB1",
                   "FCGR3A", "TYROBP", "S1PR5","IL4I1", "KIT", "LST1"), 
        group.by = "RNA_snn_res.0.3", stack = T, flip = T, pt.size = 0, add.noise = T)

act_sub$t_sub_2 <- act_sub$RNA_snn_res.0.2

act_sub$t_sub_2 <- gsub("^0$", "CD4T_MAML2", act_sub$t_sub_2)
act_sub$t_sub_2 <- gsub("^1$", "CD8T_Trm", act_sub$t_sub_2)
act_sub$t_sub_2 <- gsub("^2$", "CD4T_MAML2", act_sub$t_sub_2)
act_sub$t_sub_2 <- gsub("^3$", "CD4T_MAML2", act_sub$t_sub_2)
act_sub$t_sub_2 <- gsub("^4$", "CD4T_MAML2", act_sub$t_sub_2)

DimPlot(act_sub, reduction = "umap", group.by = "t_sub_2", label = T)

# merge data_2----
meta_1 <- snhx@meta.data["t_sub_2"]
meta_2 <- act_sub@meta.data["t_sub_2"]

unique(meta_1$t_sub_2)
meta_1 <- subset(meta_1, t_sub_2 %in% c("CD8T_GZMK", "CD4T_Naive", "CD8T_Trm", "Treg", "Proliferating_T", 
                                       "CD4T_Th17", "ILC", "CD8T_Naive", "γδ T cell" , "NK cell", "CD8T_Temra"))

merged_meta <- rbind(meta_1, meta_2)

snhx <- AddMetaData(snhx, merged_meta, col.name = "t_sub_3")

DimPlot(snhx, reduction = "umap", group.by = "t_sub_3")

# Re-make umap----
snhx <- NormalizeData(snhx)
snhx <- FindVariableFeatures(snhx)
snhx <- ScaleData(snhx)

snhx <- FindNeighbors(snhx, dims = 1:10, reduction = "XscANVI_")
snhx <- RunUMAP(snhx, dims = 1:10, reduction = "XscANVI_")

DimPlot(snhx, reduction = "umap", group.by = "t_sub_3")

saveRDS(snhx, "~/projects/CAF_projext/validation_2/atlas/snhx_anno.rds")

write.csv(snhx@meta.data, 
          file="~/projects/CAF_projext/validation_2/atlas/snhx_meta.csv",
          row.names=TRUE)

# marker plot----
snhx$t_sub_3 <- factor(snhx$t_sub_3, levels = c("ILC", "NK cell", "Proliferating_T", "γδ T cell", 
                                                "CD8T_Temra", "CD8T_GZMK", "CD8T_Trm", "CD8T_Naive", 
                                                "CD4T_MAML2", "CD4T_Th17", "Treg", "CD4T_Naive"))

genes <- c("CD3D", "CD3E", "CD4", "CD40LG", "CD8A", "CD8B", "TCF7", "CCR7", "SELL", 
           "FOXP3", "CTLA4", "IL2RA", "IL17A", "CCL20", "KLRB1",　"MAML2", "KLF12", "SIK3",
           "ITGAE", "GZMA", "GZMB", "GZMK", "EOMES", "CD27", "GNLY", "ZNF683", "FCRL6",
           "TRDC", "TRDV1", "TRGC1", "MKI67", "TOP2A", "CCNB1",
           "FCGR3A", "TYROBP", "S1PR5","IL4I1", "KIT", "LST1")

p1 <- DotPlot(snhx, features = genes, cols="RdBu", group.by = "t_sub_3", dot.scale = 10) + RotatedAxis()

ggsave("~/projects/CAF_projext/validation_2/atlas/atlas_t_dotplot.tiff", 
       plot = p1, units = "in", width = 14.6, height = 4.6, dpi = 300)

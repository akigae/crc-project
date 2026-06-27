library(Seurat)
library(ggplot2)
library(dplyr)

# Level 2 annotation
merged_obj <- readRDS("~/projects/CRC_xenium/xenium_data/xenium_obj.rds")

#### TNK ####
tnk.obj <- subset(merged_obj, Cell_Cluster_level1 == "TNK")

tnk.obj <- RunPCA(tnk.obj, npcs = 50, features = rownames(tnk.obj))
tnk.obj <- FindNeighbors(tnk.obj, reduction = "pca", dims = 1:30)
tnk.obj <- FindClusters(tnk.obj, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6))
tnk.obj <- RunUMAP(tnk.obj, reduction = "pca", dims = 1:30, reduction.name = "umap.tnk")

VlnPlot(tnk.obj, group.by = "SCT_snn_res.0.6", c("CD8A", "CD8B", "CD4"), 
        stack = T, flip = T, pt.size = 0, add.noise = T)

VlnPlot(tnk.obj, group.by = "SCT_snn_res.0.6", c("DCN", "FN1", "LUM", 'ACTA2',
                                                 "FCGR2A", "CD163", "S100A9", "CD68", "CD86", "MARCO", 
                                                 "EPCAM", "CEACAM6", "CEACAM1", "EGFR", "MET",
                                                 "IGHG1", "IGHGP", "IGHG3", "KIT", "CPA3", "MS4A1", 
                                                 "FCMR", "TNFRSF13C", "MPEG1", "GZMB",
                                                 "S100B", "NCAM1", "SPARC", "CD2", "CD3E", "IL2RB",
                                                 "PLVAP", "NOTCH3", "RGS5", "SOCS3", "STAT6", "SOCS3"), 
        stack = T, flip = T, pt.size = 0, add.noise = T)

Idents(tnk.obj) <- "SCT_snn_res.0.6"
markers <- FindAllMarkers(tnk.obj, min.pct = 0.25, only.pos = T)

FeaturePlot(tnk.obj, features = c("CD4", "CD8A", "TCF7", "CCR7", "SELL",
                                  "GZMK", "CXCR6", "ITGAE", "GZMA", "ENTPD1",
                                  "ZNF683", "GZMB"), reduction = "umap.tnk")

VlnPlot(tnk.obj, group.by = "SCT_snn_res.0.6", c("CD8A", "CD8B", "CD4", "CD40LG", "ITGAE", "GZMA", "GZMB",  
                                                 "GZMK", "EOMES", "CCL5", "TCF7", "CCR7", "SELL", 
                                                 "GNLY", "ZNF683", "KLRC1", 
                                                 "CXCL13", "FOXP3", "CTLA4", "TIGIT", "HAVCR2"), 
        stack = T, flip = T, pt.size = 0, add.noise = T)

tnk.obj$tnk_sub <- tnk.obj$SCT_snn_res.0.6

# annotation
tnk.obj$tnk_sub <- gsub("^3$", "CD8T_Trm", tnk.obj$tnk_sub)
tnk.obj$tnk_sub <- gsub("^11$", "CD8T_Trm", tnk.obj$tnk_sub)
tnk.obj$tnk_sub <- gsub("^13$", "CD8T_Trm", tnk.obj$tnk_sub)
tnk.obj$tnk_sub <- gsub("^7$", "CD8T_Trm", tnk.obj$tnk_sub)

tnk.obj$tnk_sub <- gsub("^14$", "CD8T_Naive", tnk.obj$tnk_sub)

tnk.obj$tnk_sub <- gsub("^1$", "CD8T_Eff_GNLY", tnk.obj$tnk_sub)

tnk.obj$tnk_sub <- gsub("^4$", "CD8T_Tem_GZMK", tnk.obj$tnk_sub)

tnk.obj$tnk_sub <- gsub("^9$", "CD4T", tnk.obj$tnk_sub)
tnk.obj$tnk_sub <- gsub("^0$", "CD4T", tnk.obj$tnk_sub)
tnk.obj$tnk_sub <- gsub("^2$", "CD4T", tnk.obj$tnk_sub)
tnk.obj$tnk_sub <- gsub("^6$", "CD4T", tnk.obj$tnk_sub)
tnk.obj$tnk_sub <- gsub("^5$", "CD4T", tnk.obj$tnk_sub)
tnk.obj$tnk_sub <- gsub("^8$", "CD4T", tnk.obj$tnk_sub)
tnk.obj$tnk_sub <- gsub("^10$", "CD4T", tnk.obj$tnk_sub)
tnk.obj$tnk_sub <- gsub("^12$", "CD4T", tnk.obj$tnk_sub)

DimPlot(tnk.obj, label = T, reduction = "umap.tnk", group.by = "tnk_sub")

# CD4T----
CD4.obj <- subset(tnk.obj, tnk_sub == "CD4T")
CD4.obj <- RunPCA(CD4.obj, npcs = 30, features = rownames(CD4.obj))

CD4.obj <- FindNeighbors(CD4.obj, reduction = "pca", dims = 1:30)
CD4.obj <- FindClusters(CD4.obj, resolution = 0.2)
CD4.obj <- RunUMAP(CD4.obj, reduction = "pca", dims = 1:30, reduction.name = "umap.tnk")

Idents(CD4.obj) <- "SCT_snn_res.0.2"
cd4_marker <- FindAllMarkers(CD4.obj, only.pos = T)

VlnPlot(CD4.obj, group.by = "SCT_snn_res.0.2", c("FOXP3", "CTLA4", "TIGIT",
                                                 "TCF7", "CCR7", "SELL", "IL7R", "CXCL13"), 
        stack = T, flip = T, pt.size = 0, add.noise = T)

CD4.obj$CD4_sub <- CD4.obj$SCT_snn_res.0.2

CD4.obj$CD4_sub <- gsub("^1$", "Treg", CD4.obj$CD4_sub)
CD4.obj$CD4_sub <- gsub("^4$", "Treg", CD4.obj$CD4_sub)

CD4.obj$CD4_sub <- gsub("^0$", "CD4T_Naive", CD4.obj$CD4_sub)
CD4.obj$CD4_sub <- gsub("^2$", "CD4T_Naive", CD4.obj$CD4_sub)

CD4.obj$CD4_sub <- gsub("^3$", "CD4T_CXCL13", CD4.obj$CD4_sub)

DimPlot(CD4.obj, label = T, reduction = "umap.tnk", group.by = "CD4_sub")

# merge annotation
meta_1 <- CD4.obj@meta.data["CD4_sub"]
meta_2 <- tnk.obj@meta.data["tnk_sub"]

meta_2 <- subset(meta_2, tnk_sub %in% c("CD8T_Eff_GNLY", "CD8T_Tem_GZMK", "CD8T_Trm", "CD8T_Naive"))
colnames(meta_1)[1] <- "tnk_sub"
mixed_df <- rbind(meta_1, meta_2)

tnk.obj <- AddMetaData(tnk.obj, mixed_df, col.name = "tnk_sub")
DimPlot(tnk.obj, label = T, reduction = "umap.tnk", group.by = "tnk_sub")

#### Myeloid ####
mye.obj <- subset(merged_obj, Cell_Cluster_level1 == "Myeloid")

mye.obj <- RunPCA(mye.obj, npcs = 50, features = rownames(mye.obj))

mye.obj <- FindNeighbors(mye.obj, reduction = "pca", dims = 1:30, k.param = 20)
mye.obj <- FindClusters(mye.obj, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6))
mye.obj <- RunUMAP(mye.obj, reduction = "pca", dims = 1:30, reduction.name = "umap.myeloid")

DimPlot(mye.obj, label = T, reduction = "umap.myeloid", group.by = "SCT_snn_res.0.2")

# DEG analysis
Idents(mye.obj) <- "SCT_snn_res.0.3"
markers <- FindAllMarkers(mye.obj, only.pos = TRUE, min.pct = c(0.25))

# pick up only top10 significant genes (avg_log2FC > 0.6)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.6) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

top10$cluster <- factor(top10$cluster, levels = c(unique(top10$cluster)))
top10 <- top10[order(top10$cluster), ]

top10

merged_obj_filtered <- ScaleData(mye.obj, features = top10$gene)
p <- DoHeatmap(merged_obj_filtered, features = top10$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()

# check the heatmap to check whether the clusters are well separated at this resolution
ggplot2::ggsave(filename = "mye_filtered_0.3.png", plot = p, width = 20, height = 10)

VlnPlot(mye.obj, group.by = "SCT_snn_res.0.2", c("FCGR2A", "ITGAX",  "CD163", "CD68", "S100A9", "CXCL5", 
                                                 "CD1C", "S100B", "IRF8", "CCR7", "CD80", "CD86", "S100A9"), 
        stack = T, flip = T, pt.size = 0, add.noise = T)

# annotation
mye.obj$SCT_snn_res.0.2 <- factor(mye.obj$SCT_snn_res.0.2, 
                                  levels = c("0", "1", "2", "3", "4", "5", "6", "7"))

mye.obj$sub_fin <- mye.obj$SCT_snn_res.0.2

mye.obj$sub_fin <- gsub("^0$", "C1QB_mac", mye.obj$sub_fin)
mye.obj$sub_fin <- gsub("^2$", "C1QB_mac", mye.obj$sub_fin)
mye.obj$sub_fin <- gsub("^3$", "C1QB_mac", mye.obj$sub_fin)
mye.obj$sub_fin <- gsub("^6$", "C1QB_mac", mye.obj$sub_fin)

mye.obj$sub_fin <- gsub("^1$", "Neutrophil", mye.obj$sub_fin)
mye.obj$sub_fin <- gsub("^4$", "CXCL5_mac", mye.obj$sub_fin)
mye.obj$sub_fin <- gsub("^5$", "Dendritic cell", mye.obj$sub_fin)
mye.obj$sub_fin <- gsub("^7$", "CCR7_mac", mye.obj$sub_fin)

#### Plasma ####
plasma.obj <- subset(merged_obj, Cell_Cluster_level1 == "Plasma")
plasma.obj <- RunPCA(plasma.obj, npcs = 50, features = rownames(plasma.obj))

### use PCA
# Determine percent of variation associated with each PC(https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html)
pct <- plasma.obj[["pca"]]@stdev / sum(plasma.obj[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs

plasma.obj <- FindNeighbors(plasma.obj, reduction = "pca", dims = 1:15)
plasma.obj <- FindClusters(plasma.obj, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
plasma.obj <- RunUMAP(plasma.obj, reduction = "pca", dims = 1:15, reduction.name = "umap.plasma")

DimPlot(plasma.obj, label = T, reduction = "umap.plasma", group.by = "SCT_snn_res.0.3")

# DEG analysis
Idents(plasma.obj) <- "SCT_snn_res.0.4"
markers_plasma <- FindAllMarkers(plasma.obj, only.pos = TRUE, min.pct = c(0.25))

VlnPlot(plasma.obj, group.by = "SCT_snn_res.0.4", c("JCHAIN", "IGHM", "MZB1", "XBP1",
                                                    "IGHG1", "IGHG2", "IGHG3", "IGHG4", 
                                                    "DCN", "LUM", "FN1"), 
        stack = T, flip = T, pt.size = 0, add.noise = T)

# annotation
plasma.obj$SCT_snn_res.0.4 <- factor(plasma.obj$SCT_snn_res.0.4, 
                                  levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"))

plasma.obj$sub_fin <- plasma.obj$SCT_snn_res.0.4

plasma.obj$sub_fin <- gsub("^0$", "IgG_Plasma", plasma.obj$sub_fin)
plasma.obj$sub_fin <- gsub("^2$", "IgG_Plasma", plasma.obj$sub_fin)
plasma.obj$sub_fin <- gsub("^4$", "IgG_Plasma", plasma.obj$sub_fin)
plasma.obj$sub_fin <- gsub("^6$", "IgG_Plasma", plasma.obj$sub_fin)
plasma.obj$sub_fin <- gsub("^7$", "IgG_Plasma", plasma.obj$sub_fin)
plasma.obj$sub_fin <- gsub("^10$", "IgG_Plasma", plasma.obj$sub_fin)
plasma.obj$sub_fin <- gsub("^11$", "IgG_Plasma", plasma.obj$sub_fin)

plasma.obj$sub_fin <- gsub("^1$", "IgG_Plasma_Stromal", plasma.obj$sub_fin)
plasma.obj$sub_fin <- gsub("^9$", "IgM_Plasma", plasma.obj$sub_fin)
plasma.obj$sub_fin <- gsub("^3$", "Plasma", plasma.obj$sub_fin)
plasma.obj$sub_fin <- gsub("^5$", "Plasma", plasma.obj$sub_fin)
plasma.obj$sub_fin <- gsub("^8$", "Plasma", plasma.obj$sub_fin)

#### Epithelial ####
epi.obj <- subset(merged_obj, Cell_Cluster_level1 == "Epithelial")

epi.obj <- RunPCA(epi.obj, npcs = 50, features = rownames(epi.obj))
epi.obj <- FindNeighbors(epi.obj, reduction = "pca", dims = 1:30)
epi.obj <- FindClusters(epi.obj, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))

epi.obj <- RunUMAP(epi.obj, reduction = "pca", dims = 1:30, reduction.name = "umap")

DimPlot(epi.obj, label = T, reduction = "umap", group.by = "SCT_snn_res.0.3")

Idents(epi.obj) <- "SCT_snn_res.0.2"
markers_epi <- FindAllMarkers(epi.obj, only.pos = TRUE, min.pct = c(0.25))

VlnPlot(epi.obj, group.by = "SCT_snn_res.0.2", c("IFITM3", "CCL28"), 
        stack = T, flip = T, pt.size = 0, add.noise = T)

epi.obj$SCT_snn_res.0.2 <- factor(epi.obj$SCT_snn_res.0.2, 
                                  levels = c("0", "1", "2", "3", "4", "5", "6", "7", 
                                             "8", "9", "10", "11", "12"))

epi.obj$sub_fin <- epi.obj$SCT_snn_res.0.2

epi.obj$sub_fin <- gsub("^5$", "CCL28_Epithelial", epi.obj$sub_fin)
epi.obj$sub_fin <- gsub("^9$", "CCL28_Epithelial", epi.obj$sub_fin)

epi.obj$sub_fin <- gsub("^0$", "IFITM3_Epithelial", epi.obj$sub_fin)
epi.obj$sub_fin <- gsub("^1$", "IFITM3_Epithelial", epi.obj$sub_fin)
epi.obj$sub_fin <- gsub("^2$", "IFITM3_Epithelial", epi.obj$sub_fin)
epi.obj$sub_fin <- gsub("^3$", "IFITM3_Epithelial", epi.obj$sub_fin)
epi.obj$sub_fin <- gsub("^4$", "IFITM3_Epithelial", epi.obj$sub_fin)
epi.obj$sub_fin <- gsub("^6$", "IFITM3_Epithelial", epi.obj$sub_fin)
epi.obj$sub_fin <- gsub("^7$", "IFITM3_Epithelial", epi.obj$sub_fin)
epi.obj$sub_fin <- gsub("^8$", "IFITM3_Epithelial", epi.obj$sub_fin)
epi.obj$sub_fin <- gsub("^11$", "IFITM3_Epithelial", epi.obj$sub_fin)
epi.obj$sub_fin <- gsub("^12$", "IFITM3_Epithelial", epi.obj$sub_fin)

#### Fibroblast ####
fib.obj <- subset(merged_obj, Cell_Cluster_level1 == "Fibroblast")

fib.obj <- RunPCA(fib.obj, npcs = 30, features = rownames(fib.obj))

fib.obj <- FindNeighbors(fib.obj, reduction = "pca", dims = 1:30)
fib.obj <- FindClusters(fib.obj, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6))
fib.obj <- RunUMAP(fib.obj, reduction = "pca", dims = 1:30, reduction.name = "umap.fibro")

VlnPlot(fib.obj, group.by = "SCT_snn_res.0.2", c("DCN", "FN1", "LUM", 'ACTA2',
                                                 "FCGR2A", "CD163", "S100A9", "CD68", "CD86", "MARCO", 
                                                 "EPCAM", "CEACAM6", "CEACAM1", "EGFR", "MET",
                                                 "IGHG1", "IGHGP", "JCHAIN", "KIT", "CPA3", "MS4A1", 
                                                 "FCMR", "TNFRSF13C", "MPEG1", "GZMB",
                                                 "S100B", "NCAM1", "SPARC", "CD2", "CD3E", "IL2RB",
                                                 "PLVAP", "NOTCH3", "RGS5", "SOCS3", "STAT6", "SOCS3"), 
        stack = T, flip = T, pt.size = 0, add.noise = T)

DimPlot(fib.obj, label = T, reduction = "umap.fibro", group.by = "SCT_snn_res.0.2")

Idents(fib.obj) <- "SCT_snn_res.0.2"
markers_fib <- FindAllMarkers(fib.obj, only.pos = TRUE, min.pct = c(0.25))

VlnPlot(fib.obj, group.by = "SCT_snn_res.0.2", c("DCN", "FN1", "LUM", "ACTA2",
                                                 "MGP", "CXCL12", "CXCL14", "SDC1", "CCL11", "SPARCL1"), 
        stack = T, flip = T, pt.size = 0, add.noise = T)

fib.obj$sub_fin <- fib.obj$SCT_snn_res.0.2

fib.obj$sub_fin <- gsub("^6$", "CCL11_high_Fib", fib.obj$sub_fin)
fib.obj$sub_fin <- gsub("^8$", "CCL21_Fib", fib.obj$sub_fin)
fib.obj$sub_fin <- gsub("^9$", "CCL21_Fib", fib.obj$sub_fin)
fib.obj$sub_fin <- gsub("^10$", "CCL21_Fib", fib.obj$sub_fin)
fib.obj$sub_fin <- gsub("^1$", "CXCL14_high_Fib", fib.obj$sub_fin)
fib.obj$sub_fin <- gsub("^3$", "IL6_Fib", fib.obj$sub_fin)
fib.obj$sub_fin <- gsub("^4$", "MGP_high_Fib", fib.obj$sub_fin)
fib.obj$sub_fin <- gsub("^0$", "SDC1_high_Fib", fib.obj$sub_fin)
fib.obj$sub_fin <- gsub("^2$", "SPARCL1_high_Fib", fib.obj$sub_fin)

fib.obj$sub_fin <- gsub("^7$", "Enteric_glial", fib.obj$sub_fin)

DimPlot(fib.obj, label = T, reduction = "umap.fibro", group.by = "sub_fin")


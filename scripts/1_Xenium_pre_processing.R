library(Seurat)
library(ggplot2)
library(dplyr)
library(sctransform)

#### 1. Create a Seurat object for slide_1 ####

# load Xenium data
xenium.obj <- LoadXenium("~/projects/CRC_xenium/xenium_data/xenium_output_folder/MSI", 
                         fov = "fov", segmentations="cell", flip.xy=TRUE)

# load cell IDs and coordinates for the selected region
tissue1_cells <- read.csv("~/projects/CRC_xenium/xenium_data/divide_sample/slides_1/S25_08502_border_cells_stats.csv", skip = 2, header = TRUE)

# Get a separate seurat object for the selected cells
tissue_1_obj <- subset(xenium.obj, cells = tissue1_cells$Cell.ID)

ImageFeaturePlot(tissue_1_obj, fov = "fov",
                 features = c("nCount_Xenium"),
                 max.cutoff="q95")+
  scale_y_reverse()

# add sample ID----
# S25_08502
tissue1_cells$sample_id <- "S25_08502_tumor"
S25_08502_tumor <- tissue1_cells
rownames(S25_08502_tumor) <- S25_08502_tumor$Cell.ID
head(S25_08502_tumor)

## Repeat the above steps for each sample to generate per-sample result objects ##

# merge sample objects
merge_df <- rbind(S25_08502_border, S25_08502_tumor, 
                  S25_04910_border, S25_04910_tumor, 
                  S24_23828_border, S24_23828_tumor, 
                  S24_22134_border, S24_22134_tumor, 
                  S24_21659_border, S24_21659_tumor)

xenium.obj <- AddMetaData(xenium.obj, merge_df["sample_id"], col.name = "sample_id")

saveRDS(xenium.obj, "~/projects/CRC_xenium/xenium_data/slide_1.rds")


#### 2. Create a Seurat object for slide_2 ####

# load Xenium data
xenium.obj.2 <- LoadXenium("~/projects/CRC_xenium/xenium_data/xenium_output_folder/MSS", 
                           fov = "fov", segmentations="cell", flip.xy=TRUE)

# load cell IDs and coordinates for the selected region
tissue1_cells <- read.csv("~/projects/CRC_xenium/xenium_data/divide_sample/slides_2/S25_17181_border_cells_stats.csv", skip = 2, header = TRUE)

# Get a separate seurat object for the selected cells
tissue_1_obj <- subset(xenium.obj.2, cells = tissue1_cells$Cell.ID)
tissue_1_obj

ImageFeaturePlot(tissue_1_obj, fov = "fov",
                 features = c("nCount_Xenium"),
                 max.cutoff="q95")+
  scale_y_reverse()

# add sample ID----
# S25_17181
tissue1_cells$sample_id <- "S25_17181_border"
S25_17181_border <- tissue1_cells
rownames(S25_17181_border) <- S25_17181_border$Cell.ID
head(S25_17181_border)

## Repeat the above steps for each sample to generate per-sample result objects ##

# merge sample objects
merge_df <- rbind(S25_14704_border, S25_14704_tumor, 
                  S25_15675_border, S25_15675_tumor, 
                  S25_17181_border, S25_17181_tumor, 
                  S25_16883_border, S25_16883_tumor, 
                  S25_16689_border, S25_16689_tumor)

xenium.obj.2 <- AddMetaData(xenium.obj.2, merge_df["sample_id"], col.name = "sample_id")
table(xenium.obj.2@meta.data["sample_id"])

saveRDS(xenium.obj.2, "~/projects/CRC_xenium/xenium_data/slide_2.rds")

xenium.obj@meta.data$slide_id <- "Slide_1"
xenium.obj.2@meta.data$slide_id <- "Slide_2"

hist(xenium.obj$nCount_Xenium, xlim=c(0,200), breaks = 5000)
abline(v = 20, col="red")

hist(xenium.obj.2$nCount_Xenium, xlim=c(0,200), breaks = 5000)
abline(v = 20, col="red")


#### 3. QC ####
options(future.globals.maxSize = 20 * 1024^3)
merged_obj <- merge(xenium.obj, xenium.obj.2)

saveRDS(merged_obj, "~/projects/CRC_xenium/xenium_data/before_QC.rds")

VlnPlot(merged_obj, features = c("nFeature_Xenium", "nCount_Xenium"), 
        ncol = 2, pt.size = 0, group.by = "slide_id")

## filtering cells based on detected transcripts per cell
## Upper limit cutoff - 98th percentile
## Lower limit cutoff - 20 determined by the previous plot

thres <- quantile(merged_obj$nCount_Xenium, c(0.98))
hist(merged_obj$nCount_Xenium, xlim=c(0,600), breaks = 5000)
abline(v = 20, col="red")
abline(v = thres[1], col="red")

merged_obj <- subset(merged_obj, subset = nCount_Xenium >= 20 & nCount_Xenium <= thres[1])

VlnPlot(merged_obj, features = c("nFeature_Xenium", "nCount_Xenium"), 
        ncol = 2, pt.size = 0, group.by = "slide_id")

saveRDS(merged_obj, "~/projects/CRC_xenium/xenium_data/after_qc.rds")

merged_obj <- readRDS("~/projects/CRC_xenium/xenium_data/after_qc.rds")
hist(merged_obj$nCount_Xenium, xlim=c(0,200), breaks = 5000)

merged_obj <- JoinLayers(merged_obj)

#### 4. Pre-processing ####
options(future.globals.maxSize = 20 * 1024^3)
DefaultAssay(merged_obj) <- "Xenium"

merged_obj <- SCTransform(merged_obj, assay = "Xenium")
merged_obj <- RunPCA(merged_obj, npcs = 30, features = rownames(merged_obj))
merged_obj <- RunUMAP(merged_obj, dims = 1:30)

merged_obj <- FindNeighbors(merged_obj, reduction = "pca", dims = 1:30)
merged_obj <- FindClusters(merged_obj, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))

DimPlot(merged_obj, group.by = "SCT_snn_res.0.2", label = T)

VlnPlot(merged_obj, group.by = "SCT_snn_res.0.2", c("DCN", "FN1", "LUM", 'ACTA2',
                                                    "FCGR2A", "CD163", "S100A9", "CD68", "CD86", "MARCO", 
                                                    "EPCAM", "CEACAM6", "CEACAM1", "CEACAM3", "EGFR", "MET",
                                                    "IGHG1", "IGHGP", "IGHG3", "KIT", "CPA3", "MS4A1", 
                                                    "FCMR", "TNFRSF13C", "MPEG1", "GZMB",
                                                    "S100B", "NCAM1", "SPARC", "CD2", "CD3E", "IL2RB",
                                                    "PLVAP", "NOTCH3", "RGS5", "SOCS3", "STAT6"), 
        stack = T, flip = T, pt.size = 0, add.noise = T)

# DEG analysis
Idents(merged_obj) <- "SCT_snn_res.0.2"
vsALL <- FindAllMarkers(merged_obj, min.pct = 0.25, only.pos = T)

# pick up top10 significant genes (avg_log2FC > 0.6)
vsALL %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.6) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

top10$cluster <- factor(top10$cluster, levels = c(unique(top10$cluster)))
top10 <- top10[order(top10$cluster), ]

top10

merged_obj_filtered <- ScaleData(merged_obj, features = top10$gene)
p <- DoHeatmap(merged_obj_filtered, features = top10$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()

# check the heatmap to check whether the clusters are well separated at this resolution
ggplot2::ggsave(filename = "~/projects/CRC_xenium/xenium_data/merged_obj_filtered_0.2.png", plot = p, width = 20, height = 10)

# annotation----
merged_obj$SCT_snn_res.0.2 <- factor(merged_obj$SCT_snn_res.0.2, 
                                     levels = c("0", "1", "2", "3", "4", "5", 
                                                "6", "7", "8", "9", "10", "11", "12", 
                                                "13", "14", "15", "16"))

level1.cluster.ids <- c("Epithelial", "Epithelial", "Fibroblast", "Plasma", "T/NK", "Myeloid", "Endothelial/Mural", 
                        "Epithelial", "B", "Epithelial", "Epithelial", "Myeloid", "Epithelial", "Epithelial", 
                        "Fibroblast", "Mast", "Fibroblast")

Idents(merged_obj) <- "SCT_snn_res.0.2"

names(level1.cluster.ids) <- levels(merged_obj)

merged_obj <- RenameIdents(merged_obj, level1.cluster.ids)
merged_obj[["level1.cluster.ids"]] <- Idents(object = merged_obj)

merged_obj <- AddMetaData(merged_obj, merged_obj$level1.cluster.ids, 
                          col.name = "Cell_Cluster_level1")

DimPlot(merged_obj, label = T, label.size = 4, reduction = "umap", 
        group.by = "Cell_Cluster_level1") + NoLegend()

VlnPlot(merged_obj, group.by = "Cell_Cluster_level1", c("EPCAM", "CEACAM6", "CEACAM1", "CEACAM3",
                                                        "FCGR2A", "CSF1R", "MMP12", "CD68", "CD86", "MARCO", 
                                                        "CD2", "CD3E", "IL7R", 
                                                        "IGHG1", "IGHGP", "IGHG3", 
                                                        "DCN", "FN1", "LUM", 
                                                        "PLVAP", "NOTCH3", "RGS5", 
                                                        "MS4A1", "FCMR", "TNFRSF13C", 
                                                        "KIT", "CPA3", "MS4A2"), 
        stack = T, flip = T, pt.size = 0, add.noise = T)

saveRDS(merged_obj, "~/projects/CRC_xenium/xenium_data/after_level1_annotation.rds")

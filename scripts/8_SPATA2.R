# Load required packages
library(SPATA2)
library(SPATAData)
library(Seurat)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)

# Load data----
xenium_obj <- readRDS("~/projects/CRC_xenium/xenium_data/xenium_obj.rds")

MSI_obj <- subset(xenium_obj, slide_id == "Slide_1")
MSS_obj <- subset(xenium_obj, slide_id == "Slide_2")

# Make SPATA object----
spata_object <- initiateSpataObjectXenium(
  sample_name = "MSS",
  directory_xenium = "~/projects/CRC_xenium/xenium_data/xenium_output_folder/MSS/")
spata_object

# Retrieve the cells that remained after QC
cells_keep <- Cells(MSS_obj)
cells_keep <- sub("_2$", "", cells_keep)
head(cells_keep)

spata_object1 <- subsetSpataObject(
  object = spata_object,
  barcodes = cells_keep,
  spatial_proc = TRUE,
  opt = "keep"
)

# Extract the SCT matrix from the Seurat object
sct_mat <- GetAssayData(MSS_obj, assay = "SCT", slot = "data")
colnames(sct_mat) <- sub("_2$", "", colnames(sct_mat))

# Add SCT matrix to SPATA object
spata_object1 <- addProcessedMatrix(
  object = spata_object1,
  proc_mtr = sct_mat,
  mtr_name = "SCT"
)

meta_df <- MSS_obj@meta.data
meta_df$barcodes <- sub("_2$", "", rownames(meta_df))
rownames(meta_df) <- meta_df$barcodes
meta_df$barcodes <- NULL

meta_df <- 
  tibble::rownames_to_column(meta_df, var = "barcodes") %>% 
  tibble::as_tibble()

# Add metadata to SPATA object
spata_object1 <- 
  addFeatures(
    object = spata_object1,
    feature_df = meta_df, 
    feature_names = c("region_annotation", "region_exclude"), 
    overwrite = TRUE
  )

saveRDS(spata_object1, "~/projects/CRC_xenium/xenium_data/SPATA/MSS_SPATA.obj")

## Repeat the above procedure to generate MSI_SPATA.obj ##

# Add signature score----
spata_object1 <- readRDS("~/projects/CRC_xenium/xenium_data/SPATA/MSS_SPATA.obj")

signatures <- list(
  CD8T_activation = c("GNLY", "NKG7", "GZMB", "PRF1", "GZMA", "GZMK", "LAMP1"),
  chemotaxis = c("CCL2", "CCL3", "CCL4", "CCL5", "CCL7", 
                 "CCL8", "CCL11", "CCL13", "CCL14", "CCL15", 
                 "CCL16", "CCL19", "CCL20", "CCL21", "CCL26", 
                 "CCL28", "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", 
                 "CXCL9", "CXCL10", "CXCL11", "CXCL12", "CXCL13", "CXCL14", "CXCL16")
)

spata_object1 <- addSignature(spata_object1, class = "UD", name = "CHEMOTAXIS", molecules = signatures$chemotaxis)
spata_object1 <- addSignature(spata_object1, class = "UD", name = "CYTOTOXIC_ACTIVATION", molecules = signatures$CD8T_activation)

# Create separate SPATA objects for each sample----
sample_coord <- read_csv("~/projects/CRC_xenium/SPATA/area/S25-16689_border_cells_stats.csv", skip = 2)

# Get after QC coordinates
meta <- getMetaDf(spata_object1)
cells_keep <- meta$barcodes[meta$sample_id %in% c("S25_16689_border")]
cells_keep　<- intersect(cells_keep, sample_coord$`Cell ID`)

S25_16689 <- subsetSpataObject(
  object = spata_object1,
  barcodes = cells_keep,
  opt = "keep",
  spatial_proc = TRUE
)

S25_16689 <- activateMatrix(
  object = S25_16689,
  mtr_name = "SCT"
)

# Put trajectory----
S25_16689 <- createSpatialTrajectories(object = S25_16689)

p <- plotSurface(
  object = S25_16689,
  color_by = "region_annotation",
  outline = F, pt_clrsp = "BuPu", pt_size = 0.7
)

p +
  scale_color_manual(
    values = c(
      Normal = "#7199B8",
      Tumor  = "#E5989B"
    )
  ) +
  geom_path(
    data = S25_16689@spatial@trajectories[["tumor_to_interface"]]@segment,
    aes(x = x_orig, y = y_orig),
    color = "black",
    linewidth = 1,
    arrow = arrow(
      type = "closed",
      length = unit(0.25, "cm")
    )
  )


# Chemotaxis score----
# Line plot
plotStsLineplot(
  S25_16689,
  id = "tumor_to_interface",
  variables = "UD_CHEMOTAXIS"
)

# Surface plot
plotSurface(S25_16689, color_by = "UD_CHEMOTAXIS", 
            outline = T, pt_clrsp = "BuPu", pt_size = 0.6)

# Heatmap
plotStsHeatmap(
  object = S25_16689,
  variables = signatures$chemotaxis,
  clrsp = "Reds 3",
  id = "tumor_to_interface",
  arrange_rows = "maxima"
)


# Cytotoxic score----
immune_sub <- subsetSpataObject(
  object = S25_16689,
  barcodes = S25_16689@meta_obs$barcodes[S25_16689@meta_obs$Major_cell_type %in% c("TNK")]
)

# Line plot
plotStsLineplot(
  immune_sub,
  id = "tumor_to_interface",
  variables = "UD_CYTOTOXIC_ACTIVATION"
)

# Heatmap
plotStsHeatmap(
  object = immune_sub,
  variables = signatures$CD8T_activation,
  clrsp = "Reds 3",
  id = "tumor_to_interface",
  arrange_rows = "maxima"
)

# Surface plot
df <- getCoordsDf(
  object = S25_16689,
  variables = c("UD_CYTOTOXIC_ACTIVATION", "Major_cell_type")
)

df$UD_CYTOTOXIC_ACTIVATION_immune_only <- ifelse(
  df$Major_cell_type == "TNK",
  df$UD_CYTOTOXIC_ACTIVATION,
  0
)

S25_16689 <- addFeatures(
  object = S25_16689,
  feature_df = df[, c("barcodes", "UD_CYTOTOXIC_ACTIVATION_immune_only")],
  overwrite = TRUE
)

plotSurface(
  S25_16689,
  color_by = "UD_CYTOTOXIC_ACTIVATION_immune_only",
  outline = TRUE,
  pt_size = 0.6
)

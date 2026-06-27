library(Seurat)
library(dplyr)

merged_obj <- readRDS("~/projects/CRC_xenium/xenium_data/xenium_obj.rds")

### making input data for CN tool----
meta <- merged_obj@meta.data
meta$cell_id <- rownames(meta)

# Get coordinates
img1 <- merged_obj@images[["fov"]]
img2 <- merged_obj@images[["fov.2"]]

fov_1 <- GetTissueCoordinates(img1, which = "centroids")
fov_2 <- GetTissueCoordinates(img2, which = "centroids")

coord_df = rbind(fov_1, fov_2)
rownames(coord_df) <- coord_df$cell

# Rename columns to match the CN tool input format
coord_df$XMin <- coord_df$x
coord_df$XMax <- coord_df$x
coord_df$YMin <- coord_df$y
coord_df$YMax <- coord_df$y
colnames(coord_df)[3] <- "cell_id"

meta$Class <- meta$sample_id
meta <- left_join(meta, coord_df, by = "cell_id")

meta <- meta[, c("cell_id", "sample_id", "XMin", "XMax", "YMin", "YMax", "sub_type", "Phenotype", "slide_id")]

# Rename columns to match the CN tool input format
colnames(meta)[2] <- "Class"
colnames(meta)[7] <- "Allsubtypes"

# Save as a CSV file
write.csv(meta, "~/projects/CRC_xenium/xenium_data/CN_resluts/xenium_AllSubtypes.csv", row.names = FALSE)

# Identify CN----
# Run the previously established Python code
# https://www.nature.com/articles/s43018-024-00824-y

# After identifying CNs----
gc_csd_CN <- read_csv(file = "~/projects/CRC_xenium/xenium_data/CN_resluts/ALL_cells_r=50_k=20_CN=10.csv")
gc_csd_CN_anno <- gc_csd_CN %>% mutate(All_CN=paste0('CN',neighborhood10))
gc_csd_CN_anno <- as.data.frame(gc_csd_CN_anno)
rownames(gc_csd_CN_anno) <- gc_csd_CN_anno$cell_id

merged_obj <- AddMetaData(merged_obj, gc_csd_CN_anno["All_CN"])

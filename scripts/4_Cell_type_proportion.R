library(Seurat)
library(tidyverse)
library(pheatmap)

# Load data
xenium_obj <- readRDS("~/projects/CRC_xenium/xenium_data/xenium_obj.rds")
meta_data <- xenium_obj@meta.data

# -------------------------------------------------------------------------
# 1. Major cell type
# -------------------------------------------------------------------------

# Calculate cell-type proportions for each sample
sample_prop <- meta_data %>%
  group_by(sample_id, slide_id, Phenotype, Major_cell_type) %>%
  tally() %>%
  group_by(sample_id) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# Calculate the mean proportion for each of the four groups
group_mean_matrix <- sample_prop %>%
  mutate(group_id = paste(slide_id, Phenotype, sep = "_")) %>%
  group_by(group_id, Major_cell_type) %>%
  summarise(mean_proportion = mean(proportion), .groups = 'drop') %>%
  pivot_wider(names_from = group_id, values_from = mean_proportion, values_fill = 0) %>%
  remove_rownames() %>% 
  column_to_rownames("Major_cell_type") %>%
  as.matrix()

# Create column annotations for the heatmap
annotation_col <- data.frame(
  row.names = c("MSI_Interface", "MSS_Interface", "MSI_Tumor", "MSS_Tumor"),
  Status    = c("MSI", "MSS", "MSI", "MSS"),
  Region    = c("Interface", "Interface", "Tumor", "Tumor")
)

new_ann_colors = list(
  Status = c(MSI = "#E27E7E", MSS = "#76B8B2"),
  Region = c(Interface = "#D2C5B3", Tumor = "#8A9A9F")
)

# Generate a heatmap using pheatmap with Z-score scaling
mypalette <- colorRampPalette(c("#6BAED6", "white", "#FB6A4A"))(100)

t <- pheatmap(
  mat = group_mean_matrix,
  scale = "row",                    
  annotation_col = annotation_col,     
  cluster_cols = FALSE,             
  cluster_rows = FALSE,       
  color = mypalette,
  annotation_colors = new_ann_colors,
  main = "Mean Cell Type Proportions (Row Z-score)",
  fontsize = 10,
  gaps_col = c(2, 4)
)

t

# -------------------------------------------------------------------------
# 1. Sub type
# -------------------------------------------------------------------------

# Calculate cell-type proportions for each sample
sample_prop <- meta_data %>%
  group_by(sample_id, slide_id, Phenotype, sub_type) %>%
  tally() %>%
  group_by(sample_id) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# Calculate the mean proportion for each of the four groups
group_mean_matrix <- sample_prop %>%
  mutate(group_id = paste(slide_id, Phenotype, sep = "_")) %>%
  group_by(group_id, sub_type) %>%
  summarise(mean_proportion = mean(proportion), .groups = 'drop') %>%
  pivot_wider(names_from = group_id, values_from = mean_proportion, values_fill = 0) %>%
  remove_rownames() %>% 
  column_to_rownames("sub_type") %>%
  as.matrix()

# Create column annotations for the heatmap
annotation_col <- data.frame(
  row.names = c("MSI_Interface", "MSS_Interface", "MSI_Tumor", "MSS_Tumor"),
  Status    = c("MSI", "MSS", "MSI", "MSS"),
  Region    = c("Interface", "Interface", "Tumor", "Tumor")
)

new_ann_colors = list(
  Status = c(MSI = "#E27E7E", MSS = "#76B8B2"),
  Region = c(Interface = "#D2C5B3", Tumor = "#8A9A9F")
)

group_mean_matrix_lv2 <- group_mean_matrix[c("CD4T_CXCL13", "CD4T_Naive", "CD8T_Eff_GNLY", "CD8T_Eff_GZMK", 
                                             "CD8T_Naive", "CD8T_Trm", "Treg"), ]

# Generate a heatmap using pheatmap with Z-score scaling
mypalette <- colorRampPalette(c("#6BAED6", "white", "#FB6A4A"))(100)

t <- pheatmap(
  mat = group_mean_matrix_lv2,
  scale = "row",                    
  annotation_col = annotation_col,     
  cluster_cols = FALSE,             
  cluster_rows = FALSE,       
  color = mypalette,
  annotation_colors = new_ann_colors,
  main = "Mean Cell Type Proportions (Row Z-score)",
  fontsize = 10,
  gaps_col = c(2, 4)
)

t


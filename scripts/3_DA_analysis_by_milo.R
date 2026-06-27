library(miloR)
library(Seurat)
library(scater)
library(patchwork)
library(SingleCellExperiment)
library(dplyr)
library(ggbeeswarm)
library(ggplot2)

# load data
merged_obj <- readRDS("~/projects/CRC_xenium/xenium_data/xenium_obj.rds")

# Make Milo object
sce <- as.SingleCellExperiment(merged_obj, assay = "SCT")
sce <- Milo(sce)

traj_milo <- miloR::buildGraph(sce, k = 10, d = 25, reduced.dim = "PCA")
traj_milo <- makeNhoods(traj_milo, k = 10, d = 25, prop = 0.1, refined = TRUE)
plotNhoodSizeHist(traj_milo)

traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), samples="sample_id")

# vs Phenotype
colData(traj_milo)$Sample <-colData(traj_milo)$sample_id
traj_design <- data.frame(colData(traj_milo))[,c("Sample", "Phenotype", "status")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$Sample
traj_design <- traj_design[colnames(nhoodCounts(traj_milo)), , drop=FALSE]

da_results <- testNhoods(traj_milo, design = ~ Phenotype, 
                         design.df = traj_design, 
                         fdr.weighting = "graph-overlap", norm.method = c("TMM"))

da_results_tmm_lv1 <- annotateNhoods(traj_milo, da_results, 
                                     coldata_col = "Major_cell_type")

da_results_tmm_lv1$celltype <- ifelse(da_results_tmm_lv1$Major_cell_type_fraction < 0.7, "Mixed", da_results_tmm_lv1$Major_cell_type)
plotDAbeeswarm(da_results_tmm_lv1, group.by = "celltype")

# Visualizing
alpha <- 0.1
set.seed(12345)

color_positive <- "#7FD1B9"
color_negative <- "#FF8E9E"

da_results_tmm_lv1 <- subset(da_results_tmm_lv1, celltype %in% c("Epithelial", "Fibroblast", "Endothelial/Mural", "TNK", 
                                                                 "B", "Plasma", "Myeloid", "Mast", "Enteric Glial"))

celltype_order <- c("Epithelial", "Fibroblast", "Endothelial/Mural", "TNK", 
                    "B", "Plasma", "Myeloid", "Mast", "Enteric Glial")

pl1 <- da_results_tmm_lv1 %>%
  group_by(celltype) %>%
  mutate(mean_logFC_all = median(logFC)) %>%
  ungroup() %>%
  mutate(
    nhood_anno = factor(celltype, levels = rev(celltype_order))
  ) %>%
  mutate(direction = case_when(
    SpatialFDR < alpha & logFC > 0 ~ "Up",
    SpatialFDR < alpha & logFC <= 0 ~ "Down",
    TRUE ~ "Not_Significant"
  )) %>%
  ggplot(aes(nhood_anno, logFC)) +
  geom_quasirandom(data = . %>% filter(direction == "Not_Significant"), 
                   size = 1, color = "grey90", alpha = 0.4, width = 0.3) +
  geom_quasirandom(data = . %>% filter(direction != "Not_Significant"), 
                   aes(color = direction), size = 1, alpha = 0.8, width = 0.3) +
  geom_point(data = . %>% distinct(nhood_anno, mean_logFC_all), 
             aes(x = nhood_anno, y = mean_logFC_all), color = "black", size = 2) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) + 
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  scale_color_manual(values = c("Up" = color_positive, "Down" = color_negative)) +
  theme_bw(base_size = 18) +
  labs(x = "", y = "log-Fold Change") +
  scale_y_continuous(breaks = seq(-6, 6, by = 2)) +
  guides(color = "none")   + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

pl1

library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(ggalluvial)
library(ggbreak)

atlas_obj <- readRDS("~/datas/atlas_only_ICI.rds")
atlas_sub <- subset(atlas_obj, sample_type %in% c("primary tumor"))
atlas_sub <- subset(atlas_sub, atlas_cell_type_coarse %in% c("Plasma cell"))

# Umap
atlas_sub <- FindNeighbors(atlas_sub, reduction = "X_scANVI")
atlas_sub <- RunUMAP(atlas_sub, dims = 1:10, reduction = "X_scANVI")
p1 = DimPlot(atlas_sub, group.by = "atlas_cell_type_fine", 
             reduction = "umap") + coord_fixed()
p1

ggsave("ICI_plasma_umap.png", plot = p1, 
       dpi = 400, 
       width = 9, height = 3.5, units = "in")

# Isotype annotation----
### MAX CHAIN approach
ix <- which(rownames(atlas_obj@assays[["originalexp"]]@counts) %in% c("IGHD", "IGHE", "IGHM", "IGHG1", "IGHG2", 
                                                                      "IGHG3", "IGHG4", "IGHA1", "IGHA2"))
x <- as.data.frame(t(as.matrix(atlas_obj@assays[["originalexp"]]@counts[ix,])))
colnames(x) <- rownames(atlas_obj@assays[["originalexp"]])[ix]
x$max.chain <- colnames(x)[max.col(x[,1:9])]
x$max.chain[x$max.chain %in% c("IGHA1","IGHA2")] <- 'IGHA'
atlas_obj$final_isotype_R <- x$max.chain

ISO_MAIN <- c("IGHD","IGHE","IGHM","IGHG1","IGHG2","IGHG3","IGHG4")
A_GENES  <- c("IGHA1","IGHA2")

cts <- LayerData(atlas_obj, assay = "originalexp", layer = "counts")
main_present <- intersect(ISO_MAIN, rownames(cts))
a_present    <- intersect(A_GENES,  rownames(cts))

# cells x genes（main）
mat_main <- t(as.matrix(cts[main_present, , drop=FALSE]))

# Sum IGHA1 and IGHA2
A_sum <- if (length(a_present) > 0) {
  Matrix::colSums(cts[a_present, , drop=FALSE])
} else {
  rep(0, ncol(cts))
}

cts[a_present, , drop=FALSE]
A_sum[1:5]

mat <- cbind(mat_main, IGHA = as.numeric(A_sum))

# MAX gene/class per cell
max_chain <- colnames(mat)[max.col(mat, ties.method="first")]

ntie <- apply(mat, 1, function(v) sum(v == max(v)))
max_chain[ntie > 1] <- "unknown"
atlas_obj$final_isotype_Amerge <- max_chain

new_col <- paste(atlas_obj$atlas_cell_type_coarse, atlas_obj$final_isotype_Amerge, sep="_")

# Put final isotype annotation to metadata
atlas_obj$iso_celltype <- new_col
table(atlas_obj$iso_celltype)

saveRDS(atlas_obj@meta.data, "~/datas/validation/atlas_b_isotype.rds")

# Proportion bar plot----
all_meta <- atlas_obj@meta.data
filter_meta <- subset(all_meta, sample_type %in% c("primary tumor"))
filter_meta <- subset(filter_meta, atlas_cell_type_coarse %in% c("Plasma cell"))
filter_meta <- subset(filter_meta, microsatellite_status %in% c("MSI", "MSI-H", "MSS"))

filter_meta$microsatellite_status <- gsub("MSI-H", "MSI", filter_meta$microsatellite_status)
filter_meta$atlas_cell_type_fine <- as.character(filter_meta$atlas_cell_type_fine)

prop <- prop.table(table(filter_meta$iso_celltype, filter_meta$microsatellite_status), margin = 2) * 100
df <- as.data.frame.matrix(prop)
df$celltype <- rownames(df)

# wide → long
df_long <- df %>%
  pivot_longer(cols = c("MSI", "MSS"),
               names_to = "Group",
               values_to = "Proportion")

df_long$major_type <- ifelse(grepl("^Plasma", df_long$celltype), "Plasma cell", "B cell")
df_long$isotype <- str_extract(df_long$celltype, "IG[A-Z0-9]+|unknown")

iso_colors <- c(
  "#6699CC", 
  "#F2AD4E", 
  "#74B482", 
  "#ED7D7D", 
  "#A28ECC", 
  "#81C7D4", 
  "#F5D05B", 
  "#D99694", 
  "#90BE6D", 
  "#A5A5A5"  
)

p1 <- ggplot(df_long, aes(x = Group, y = Proportion, fill = celltype)) +
  geom_flow(aes(alluvium = celltype), alpha = .5, color = "black",
            curve_type = "sigmoid", width = .5) +
  geom_col(width = .6, color = "black") +
  scale_fill_manual(values = iso_colors) +
  coord_flip() + 
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() + theme_classic()

p1

p1 <- ggplot(df_long, aes(x = Group, y = Proportion, fill = celltype)) +
  geom_flow(
    aes(alluvium = celltype),
    alpha = .5,
    color = "black",
    curve_type = "sigmoid",
    width = .5
  ) +
  geom_col(
    width = .6,
    color = "black"
  ) +
  scale_fill_manual(values = iso_colors) +
  scale_y_break(
    c(50, 75),   
    scales = 0.2,
    space = 0.1
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  coord_flip() +
  
  cowplot::theme_minimal_hgrid() +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.line.x.top = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = NULL
  )
p1


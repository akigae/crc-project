library(Seurat)
library(schard)
library(ggalluvial)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

# load Xenium data
xenium_obj <- readRDS("~/projects/CRC_xenium/xenium_data/xenium_obj.rds")
b_plasma <- subset(xenium_obj, Major_cell_type == "Plasma")

### MAX CHAIN approach
ix <- which(rownames(b_plasma@assays$SCT) %in% c("IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHG4"))
x <- as.data.frame(t(as.matrix(b_plasma@assays[["SCT"]]@counts[ix,])))
colnames(x) <- rownames(b_plasma@assays[["SCT"]])[ix]
x$max.chain <- colnames(x)[max.col(x[,1:5])]
b_plasma$final_isotype_R <- x$max.chain

ISO_MAIN <- c("IGHG1","IGHG2","IGHG3","IGHG4", "IGHM")
cts <- b_plasma@assays$SCT@counts

main_present <- intersect(ISO_MAIN, rownames(cts))
mat <- t(as.matrix(cts[main_present, , drop=FALSE]))

# MAX gene/class per cell
max_chain <- colnames(mat)[max.col(mat, ties.method="first")]

ntie <- apply(mat, 1, function(v) sum(v == max(v)))
max_chain[ntie > 1] <- "unknown"
b_plasma$final_isotype_Amerge <- max_chain

new_col <- paste(b_plasma$Major_cell_type, b_plasma$final_isotype_Amerge, sep="_")
b_plasma$iso_celltype <- new_col

# Proportion bar plot----
all_meta <- b_plasma@meta.data
all_meta$iso_celltype <- as.character(all_meta$iso_celltype)
prop <- prop.table(table(all_meta$iso_celltype, all_meta$slide_id), margin = 2) * 100
df <- as.data.frame.matrix(prop)
df$celltype <- rownames(df)

# wide → long
df_long <- df %>%
  pivot_longer(cols = c("MSI", "MSS"),
               names_to = "Group",
               values_to = "Proportion")

# color palette
iso_colors <- c(
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

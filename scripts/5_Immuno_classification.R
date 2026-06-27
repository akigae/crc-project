library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyr)

merged_obj <- readRDS("~/projects/CRC_xenium/xenium_data/xenium_obj.rds")
meta_data <- merged_obj@meta.data

#### Calculate the area size for each sample ####

tissue1_cells <-read.csv("~/projects/CRC_xenium/SPATA/area/S24-21659_border_cells_stats.csv", 
                         skip = 2, header = TRUE)
tissue1_cells$Cell.ID <- paste0(tissue1_cells$Cell.ID, "_1")

# Subset the data based on the selected cell IDs
gene_exp_data_tissue1 <- merged_obj[, tissue1_cells$Cell.ID]

# Get coordinates
fov_1 <- GetTissueCoordinates(gene_exp_data_tissue1, which = "centroids")
head(fov_1)

fov_1 <- fov_1[is.finite(fov_1$x) & is.finite(fov_1$y), ]

# Compute the convex hull
h <- chull(fov_1$x, fov_1$y)
poly <- fov_1[c(h, h[1]), c("x","y")]

# Calculate area size (shoelace formula)
area_um2 <- 0.5 * abs(sum(poly$x[-1]*poly$y[-nrow(poly)] - poly$x[-nrow(poly)]*poly$y[-1]))
area_mm2 <- area_um2 / 1e6

plot(
  fov_1$x, fov_1$y,
  pch = 16, cex = 0.4,
  xlab = "x", ylab = "y",
  asp = 1,
  main = sprintf("Convex hull (area = %.2f mm^2)", area_mm2)
)

# Visualizing
lines(poly$x, poly$y, col = "red", lwd = 2)
points(poly$x, poly$y, pch = 19, cex = 0.8, col = "red")

# calculating CD8T cell number per mm2----
# S24-21659(Interface) : 7.396217
S24_21659_border_df <- table(gene_exp_data_tissue1$final_annotation_ver4)/7.396217
S24_21659_border_df <- S24_21659_border_df[c("CD8T_Eff_GNLY", "CD8T_Eff_GZMK", "CD8T_Naive", "CD8T_Trm")]

## Repeat the above steps for each sample to generate per-sample result objects ##

# Merge results
df_list <- list(
  S24_21659_border = S24_21659_border_df,
  S24_21659_tumor  = S24_21659_tumor_df,
  S24_22134_border = S24_22134_border_df,
  S24_22134_tumor  = S24_22134_tumor_df,
  S24_23828_border = S24_23828_border_df,
  S24_23828_tumor  = S24_23828_tumor_df,
  S25_04910_border = S25_04910_border_df,
  S25_04910_tumor  = S25_04910_tumor_df,
  S25_08052_border = S25_08052_border_df,
  S25_08052_tumor  = S25_08052_tumor_df,
  S25_14704_border = S25_14704_border_df,
  S25_14704_tumor  = S25_14704_tumor_df,
  S25_15675_border = S25_15675_border_df,
  S25_15675_tumor  = S25_15675_tumor_df,
  S25_16689_border = S25_16689_border_df,
  S25_16689_tumor  = S25_16689_tumor_df,
  S25_16883_border = S25_16883_border_df,
  S25_16883_tumor  = S25_16883_tumor_df,
  S25_17181_border = S25_17181_border_df,
  S25_17181_tumor  = S25_17181_tumor_df
)

result_df <- bind_rows(
  lapply(names(df_list), function(x) {
    as.data.frame(t(df_list[[x]])) %>%
      mutate(sample_id = x, .before = 1)
  })
)

#### Statistical analysis ####

result_df　<- readRDS("~/projects/CRC_xenium/xenium_data/CD8T_density.rds")

sample_total_df <- result_df %>%
  group_by(sample_id) %>%
  summarise(
    total_value = sum(Freq, na.rm = TRUE),
    .groups = "drop"
  )

sample_total_df <- sample_total_df %>%
  mutate(
    region = ifelse(grepl("border", sample_id), "border", "tumor"),
    patient_id = gsub("_(border|tumor)$", "", sample_id)
  )

msi_samples <- c(
  "S24_21659", "S24_22134", "S24_23828",
  "S25_04910", "S25_08052"
)

mss_samples <- c(
  "S25_14704", "S25_15675", "S25_16689",
  "S25_16883", "S25_17181"
)

sample_total_df <- sample_total_df %>%
  mutate(
    status = case_when(
      patient_id %in% msi_samples ~ "MSI",
      patient_id %in% mss_samples ~ "MSS",
      TRUE ~ NA_character_
    )
  )

sample_total_df

# visualizing----
border_df <- subset(sample_total_df, region == "border")
tumor_df <- subset(sample_total_df, region == "tumor")
MSI_df <- subset(sample_total_df, status == "MSI")
MSS_df <- subset(sample_total_df, status == "MSS")

MSI_wide <- MSI_df %>%
  select(patient_id, region, total_value) %>%
  pivot_wider(
    names_from = region,
    values_from = total_value
  )

MSS_wide <- MSS_df %>%
  select(patient_id, region, total_value) %>%
  pivot_wider(
    names_from = region,
    values_from = total_value
  )

border_p <- t.test(total_value ~ status, data = border_df)
tumor_p <- t.test(total_value ~ status, data = tumor_df)
MSI_p <- t.test(MSI_wide$border, MSI_wide$tumor,paired = TRUE)
MSS_p <- t.test(MSS_wide$border, MSS_wide$tumor, data = MSS_wide, paired = TRUE)

sample_total_df <- sample_total_df %>%
  mutate(
    region = factor(region, levels = c("border", "tumor"))
  )
sample_total_df

y_max <- max(sample_total_df$total_value, na.rm = TRUE) * 1.2

p <- ggplot(
  sample_total_df,
  aes(x = region, y = total_value, group = patient_id)) +
  geom_hline(
    yintercept = 200,
    linetype = "dashed",
    linewidth = 0.5,
    color = "gray40"
  ) +
  geom_line(
    color = "gray20",
    linewidth = 0.5,
    alpha = 0.7
  ) +
  geom_point(
    aes(color = status),
    size = 3,
    alpha = 0.9
  ) +
  scale_color_manual(
    values = c("MSI" = "#E54B4B", "MSS" = "#4272D7")
  ) + 
  scale_x_discrete(
    expand = expansion(mult = c(0.3, 0.3))
  ) + 
  coord_cartesian(ylim = c(0, y_max)) + 
  labs(
    x = "Region",
    y = "CD8T density (cells/mm²)",
    color = NULL,
    title = "CD8T density"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text        = element_text(color = "black", face = "bold"),
    axis.title       = element_text(color = "black", face = "bold"),
    axis.line        = element_line(linewidth = 0.8, color = "black"),
    axis.ticks       = element_line(linewidth = 0.8, color = "black"),
    legend.position  = "right",
    legend.text      = element_text(size = 12, face = "bold"),
    plot.title       = element_text(size = 16, face = "bold", color = "black", margin = margin(b = 10))
  )

p

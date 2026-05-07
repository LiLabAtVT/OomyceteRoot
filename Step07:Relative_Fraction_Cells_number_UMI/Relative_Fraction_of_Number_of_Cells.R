# ===============================================================
# Relative fraction of cells per cluster and condition
# ===============================================================

# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plyr)

# Define input and output paths
seurat_path <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/Patho_Nonpatho_integration_07072025.Rds"
out_dir <- "/projects/intro2gds/Razan_2026/PLant_Pathogen_2026/Step_70_Relative_Fraction_Cells_number_UMI/Script"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load Seurat object
seurat_obj <- readRDS(seurat_path)

# Remove duplicated metadata columns if any
dup_cols <- duplicated(colnames(seurat_obj@meta.data))
if (any(dup_cols)) {
  cat("Removing duplicated metadata columns:\n")
  print(colnames(seurat_obj@meta.data)[dup_cols])
  seurat_obj@meta.data <- seurat_obj@meta.data[, !dup_cols]
}

# Assign simplified condition labels
seurat_obj$condition <- ifelse(
  grepl("Neg", seurat_obj$orig.ident, ignore.case = TRUE),
  "Non_Infected",
  "Infected"
)

# Define cluster labels
cluster_labels <- c(
  '0'  = "0:Meristem_1",
  '1'  = "1:Initial_cells",
  '2'  = "2:Meristem_2",
  '3'  = "3:Companion_cell",
  '4'  = "4:Atrichoblast",
  '5'  = "5:Epidermis/Cortex",
  '6'  = "6:Procambium",
  '7'  = "7:Rootcap",
  '8'  = "8:Cortex",
  '9'  = "9:G2M_phase",
  '10' = "10:Endodermis",
  '11' = "11:Trichoblast"
)

# Add cluster labels to metadata
seurat_obj$cluster_label <- plyr::revalue(
  as.character(seurat_obj$seurat_clusters),
  replace = cluster_labels,
  warn_missing = FALSE
)

# Compute relative fraction of cells per cluster and condition
cell_fraction_df <- seurat_obj@meta.data %>%
  dplyr::group_by(seurat_clusters, cluster_label, condition) %>%
  dplyr::summarise(num_cells = dplyr::n(), .groups = "drop")

# Ensure all clusters and both conditions are represented
all_combos <- expand.grid(
  seurat_clusters = as.character(0:11),
  cluster_label   = unname(cluster_labels),
  condition       = c("Non_Infected", "Infected"),
  stringsAsFactors = FALSE
)

cell_fraction_df <- dplyr::full_join(
  all_combos,
  cell_fraction_df,
  by = c("seurat_clusters", "cluster_label", "condition")
) %>%
  dplyr::mutate(num_cells = tidyr::replace_na(num_cells, 0)) %>%
  dplyr::group_by(seurat_clusters, cluster_label) %>%
  dplyr::mutate(relative_fraction = num_cells / sum(num_cells)) %>%
  dplyr::ungroup()

# Order x-axis by cluster IDs
cell_fraction_df$cluster_label <- factor(
  cell_fraction_df$cluster_label,
  levels = unname(cluster_labels)
)

# Output file
pdf_file <- file.path(out_dir, "cluster_relative_fraction_stacked_by_condition_withIDs.pdf")

# Plot
pdf(pdf_file, width = 15, height = 10)

p <- ggplot(cell_fraction_df, aes(x = cluster_label, y = relative_fraction, fill = condition)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.85) +
  scale_fill_manual(values = c("Non_Infected" = "#f7f7f7", "Infected" = "#e41a1c")) +
  labs(
    x = "Cluster",
    y = "Relative fraction of number of cells",
    fill = "Condition"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),

    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 24,
      margin = margin(t = 2)
    ),
    axis.text.y = element_text(size = 24),
    axis.title.x = element_text(size = 28, margin = margin(t = 10)),
    axis.title.y = element_text(size = 26),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)
  )

print(p)
dev.off()

cat("Stacked bar plot saved to:\n", pdf_file, "\n")

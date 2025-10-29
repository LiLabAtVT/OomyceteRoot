# ===============================================================
# Relative fraction of cells per cluster and condition (with cluster IDs on X-axis)
# ===============================================================

# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

# Define input and output paths
seurat_path <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/Patho_Nonpatho_integration_07072025.Rds"
out_dir <- "/projects/songli_lab/RazanPlantPath_2025/Script_New_Cluster_ID/OutPut"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load Seurat object
seurat_obj <- readRDS(seurat_path)

# 🧹 Remove duplicated metadata columns (if any)
dup_cols <- duplicated(colnames(seurat_obj@meta.data))
if (any(dup_cols)) {
  cat("⚠️ Removing duplicated metadata columns:\n")
  print(colnames(seurat_obj@meta.data)[dup_cols])
  seurat_obj@meta.data <- seurat_obj@meta.data[, !dup_cols]
}

# ✅ Assign simplified condition labels
seurat_obj$condition <- ifelse(grepl("Neg", seurat_obj$orig.ident, ignore.case = TRUE),
                               "Non_Infected", "Infected")

# ✅ Define cluster labels (with IDs)
cluster_labels <- c(
  '0' = "0:Meristem_1",
  '1' = "1:Meristem_2", 
  '2' = "2:Meristem_3",
  '3'= "3:Pericycle", 
  '4'= "4:Atrichoblast", 
  '5' = "5:Rootcap_1", 
  '6' = "6:Xylem", 
  '7'= "7:Rootcap_2", 
  '8' = "8:Cortex", 
  '9' = "9:Meristem_4", 
  '10'= "10:Endodermis", 
  '11' = "11:Trichoblast"
)

# Add cluster labels to metadata
seurat_obj$cluster_label <- plyr::revalue(as.character(seurat_obj$seurat_clusters),
                                          replace = cluster_labels, warn_missing = FALSE)

# 📊 Compute relative fraction of cells per cluster and condition
cell_fraction_df <- seurat_obj@meta.data %>%
  group_by(seurat_clusters, cluster_label, condition) %>%
  summarise(num_cells = n(), .groups = "drop")

# 🧩 Ensure all clusters (0–11) and both conditions are represented
all_combos <- expand.grid(
  seurat_clusters = as.character(0:11),
  cluster_label   = unname(cluster_labels),
  condition       = c("Non_Infected", "Infected")
)

cell_fraction_df <- full_join(all_combos, cell_fraction_df,
                              by = c("seurat_clusters", "cluster_label", "condition")) %>%
  mutate(num_cells = replace_na(num_cells, 0)) %>%
  group_by(seurat_clusters, cluster_label) %>%
  mutate(relative_fraction = num_cells / sum(num_cells)) %>%
  ungroup()

# Order x-axis by cluster IDs
cell_fraction_df$cluster_label <- factor(cell_fraction_df$cluster_label, levels = unname(cluster_labels))

# 🎨 Plot: stacked relative fraction per cluster (x-axis shows cluster ID + label)
pdf(file.path(out_dir, "cluster_relative_fraction_stacked_by_condition_withIDs.pdf"), width = 10, height = 5)
ggplot(cell_fraction_df, aes(x = cluster_label, y = relative_fraction, fill = condition)) +
  geom_bar(stat = "identity", color = "black", size = 0.3, width = 0.85) +
  scale_fill_manual(values = c("Non_Infected" = "#f7f7f7", "Infected" = "#e41a1c")) +
  labs(x = "Cluster ID and Annotation",
       y = "Relative fraction of number of cells",
       fill = "Condition") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )
dev.off()

cat("✅ Stacked bar plot saved to:\n",
    file.path(out_dir, "cluster_relative_fraction_stacked_by_condition_withIDs.pdf"), "\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(readr)
  library(forcats)
  library(tidyr)
})

# ===============================================================
# Paths
# ===============================================================
seurat_obj <- readRDS("/projects/intro2gds/Razan_2026/PLant_Pathogen_2026/Step_50_Dot_Plot_Heat_map/dotplot_Cluster.rds")

out_dir <- "/projects/intro2gds/Razan_2026/PLant_Pathogen_2026/Step_50_Dot_Plot_Heat_map"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

markers_csv  <- file.path(out_dir, "all_markers_heatmap.csv")
top_markers_csv <- file.path(out_dir, "top_markers_heatmap.csv")
heatmap_pdf  <- file.path(out_dir, "top_marker_heatmap_04162026.pdf")

# ===============================================================
# Define cluster labels and order
# ===============================================================
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

# ------------------ Set cluster identities ------------------
if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  warning("`seurat_clusters` not in meta.data; using current Idents().")
  orig_clusters <- as.character(Idents(seurat_obj))
} else {
  Idents(seurat_obj) <- seurat_obj$seurat_clusters
  orig_clusters <- as.character(Idents(seurat_obj))
}

# Check mapping covers all clusters
missing_levels <- setdiff(unique(orig_clusters), names(cluster_ordered_labels))
if (length(missing_levels) > 0) {
  stop("These clusters have no label mapping: ", paste(missing_levels, collapse = ", "))
}

# Apply mapping
seurat_obj$renamed_clusters <- dplyr::recode(orig_clusters, !!!cluster_ordered_labels)

# Lock desired order
desired_levels <- unname(cluster_ordered_labels)
seurat_obj$renamed_clusters <- factor(seurat_obj$renamed_clusters, levels = desired_levels)

# Use renamed clusters as identities
Idents(seurat_obj) <- seurat_obj$renamed_clusters

# ------------------ Downsample per cluster ------------------
# Helps keep the heatmap visually cleaner
downsample_per_cluster <- 250
seurat_small <- subset(seurat_obj, downsample = downsample_per_cluster)
Idents(seurat_small) <- seurat_small$renamed_clusters

# ------------------ Find markers on full object ------------------
markers <- FindAllMarkers(
  object = seurat_obj,
  only.pos = TRUE,
  min.pct = 0.01,
  logfc.threshold = 0.25
)

write.csv(markers, markers_csv, row.names = FALSE)

# ------------------ Select top markers per cluster ------------------
# Fewer genes per cluster = more readable heatmap
top_n_genes <- 6

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = top_n_genes, with_ties = FALSE) %>%
  ungroup() %>%
  filter(avg_log2FC > 0.5) %>%
  mutate(
    cluster = as.character(cluster),
    cluster_name = factor(cluster, levels = names(cluster_ordered_labels), labels = desired_levels)
  )

# Fallback if filtering is too strict
if (nrow(top_markers) == 0) {
  warning("No markers passed avg_log2FC > 0.5 filter. Using top markers without extra filtering.")
  top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = top_n_genes, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      cluster = as.character(cluster),
      cluster_name = factor(cluster, levels = names(cluster_ordered_labels), labels = desired_levels)
    )
}

write.csv(top_markers, top_markers_csv, row.names = FALSE)

# ------------------ Sort genes ------------------
sorted_genes <- top_markers %>%
  arrange(cluster_name, desc(avg_log2FC)) %>%
  pull(gene) %>%
  unique()

# Keep only genes present in the object
features_to_plot <- intersect(sorted_genes, rownames(seurat_small))
if (length(features_to_plot) == 0) {
  stop("None of the selected marker genes are present in seurat_small.")
}

cat("Number of genes selected for heatmap:", length(features_to_plot), "\n")

# ------------------ Save RDS file for heatmap plotting only ------------------
heatmap_rds_file <- file.path(out_dir, "heatmap_plot_object_downsampled.rds")

saveRDS(seurat_small, heatmap_rds_file)

cat("Heatmap plotting RDS saved to:\n", heatmap_rds_file, "\n")

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
seurat_obj <- readRDS("/projects/intro2gds/Razan_2026/PLant_Pathogen_2026/Step_50_Dot_Plot_Heat_map/heatmap_trend.rds")

out_dir <- "/projects/intro2gds/Razan_2026/PLant_Pathogen_2026/Step_50_Dot_Plot_Heat_map"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

markers_csv      <- file.path(out_dir, "all_markers_heatmap_analysis.csv")
top_markers_csv  <- file.path(out_dir, "top_markers_heatmap_analysis.csv")
heatmap_pdf      <- file.path(out_dir, "Heatmap_Top_Marker_Genes_by_Cluster_04192026_FIXED.pdf")

# ===============================================================
# Cluster labels
# ===============================================================
cluster_labels <- c(
  "0"  = "0:Meristem_1",
  "1"  = "1:Initial_cells",
  "2"  = "2:Meristem_2",
  "3"  = "3:Companion_cell",
  "4"  = "4:Atrichoblast",
  "5"  = "5:Epidermis/Cortex",
  "6"  = "6:Procambium",
  "7"  = "7:Rootcap",
  "8"  = "8:Cortex",
  "9"  = "9:G2M_phase",
  "10" = "10:Endodermis",
  "11" = "11:Trichoblast"
)

desired_levels <- unname(cluster_labels)

# ===============================================================
# Set identities
# ===============================================================
if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  Idents(seurat_obj) <- seurat_obj$seurat_clusters
}

orig_clusters <- as.character(Idents(seurat_obj))

cat("Original cluster IDs found in object:\n")
print(sort(unique(orig_clusters)))
cat("\n")

missing_levels <- setdiff(unique(orig_clusters), names(cluster_labels))
if (length(missing_levels) > 0) {
  stop("Missing cluster labels for: ", paste(missing_levels, collapse = ", "))
}

seurat_obj$renamed_clusters <- unname(cluster_labels[orig_clusters])

if (any(is.na(seurat_obj$renamed_clusters))) {
  bad_clusters <- unique(orig_clusters[is.na(seurat_obj$renamed_clusters)])
  stop("Some clusters could not be mapped: ", paste(bad_clusters, collapse = ", "))
}

seurat_obj$renamed_clusters <- factor(
  seurat_obj$renamed_clusters,
  levels = desired_levels
)

Idents(seurat_obj) <- seurat_obj$renamed_clusters

cat("Renamed cluster counts:\n")
print(table(Idents(seurat_obj), useNA = "ifany"))
cat("\n")

# ===============================================================
# Find markers
# ===============================================================
markers <- FindAllMarkers(
  object = seurat_obj,
  only.pos = TRUE,
  min.pct = 0.01,
  logfc.threshold = 0.25
)

write.csv(markers, markers_csv, row.names = FALSE)

# ===============================================================
# Select top markers
# ===============================================================
top_n_genes <- 4

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = top_n_genes, with_ties = FALSE) %>%
  ungroup() %>%
  filter(avg_log2FC > 0.5) %>%
  mutate(cluster_name = factor(as.character(cluster), levels = desired_levels))

if (nrow(top_markers) == 0) {
  warning("No strong markers found — relaxing filter")
  top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = top_n_genes, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(cluster_name = factor(as.character(cluster), levels = desired_levels))
}

write.csv(top_markers, top_markers_csv, row.names = FALSE)

# ===============================================================
# Prepare genes
# ===============================================================
sorted_genes <- top_markers %>%
  arrange(cluster_name, desc(avg_log2FC)) %>%
  pull(gene) %>%
  unique()

features_to_plot <- intersect(sorted_genes, rownames(seurat_obj))

if (length(features_to_plot) == 0) {
  stop("No genes found in Seurat object.")
}

cat("Genes selected:", length(features_to_plot), "\n")

# ===============================================================
# Scale data
# ===============================================================
seurat_obj <- ScaleData(
  seurat_obj,
  features = rownames(seurat_obj),
  verbose = FALSE
)

scaled_mat <- GetAssayData(seurat_obj, layer = "scale.data")
features_to_plot <- intersect(features_to_plot, rownames(scaled_mat))

if (length(features_to_plot) == 0) {
  stop("No selected genes found in scale.data.")
}

# ===============================================================
# Downsample cells per cluster
# ===============================================================
set.seed(123)

cells_keep <- unlist(lapply(desired_levels, function(cl) {
  cl_cells <- WhichCells(seurat_obj, idents = cl)

  if (length(cl_cells) == 0) {
    return(character(0))
  }

  if (length(cl_cells) > 80) {
    sample(cl_cells, 80)
  } else {
    cl_cells
  }
}))

cells_keep <- unique(cells_keep)
seurat_plot <- subset(seurat_obj, cells = cells_keep)

# Reset levels after subsetting
seurat_plot$renamed_clusters <- factor(
  as.character(seurat_plot$renamed_clusters),
  levels = desired_levels
)

Idents(seurat_plot) <- seurat_plot$renamed_clusters
Idents(seurat_plot) <- factor(Idents(seurat_plot), levels = desired_levels)

cat("Cluster counts in plotting object:\n")
print(table(Idents(seurat_plot), useNA = "ifany"))
cat("\n")

# ===============================================================
# Keep only genes present in plotting object
# ===============================================================
features_to_plot <- intersect(features_to_plot, rownames(seurat_plot))

if (length(features_to_plot) == 0) {
  stop("No genes found in plotting object.")
}

# ===============================================================
# Heatmap
# ===============================================================
heatmap_plot <- DoHeatmap(
  object     = seurat_plot,
  features   = features_to_plot,
  raster     = TRUE,
  draw.lines = FALSE,
  size       = 5
) +
  ggtitle("Top Marker Genes by Cluster") +
  theme(
    plot.title = element_text(size = 26, face = "bold", hjust = 0.5),

    # remove gray block = hide cell barcode labels
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),

    # gene labels
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.y = element_blank(),

    # cleaner legends
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),

    panel.grid = element_blank(),
    panel.background = element_blank()
  ) +
  guides(
    fill = guide_legend(ncol = 1, override.aes = list(size = 6)),
    color = guide_colorbar(order = 2)
  )

# ===============================================================
# Save PDF
# ===============================================================
pdf(heatmap_pdf, width = 20, height = 13)
print(heatmap_plot)
dev.off()

cat("Heatmap saved to:", heatmap_pdf, "\n")
cat("Marker CSV saved:", markers_csv, "\n")
cat("Top marker CSV saved:", top_markers_csv, "\n")

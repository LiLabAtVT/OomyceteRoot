# ------------------ Libraries ------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# ------------------ Paths ------------------
seurat_path <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/Patho_Nonpatho_integration_07072025.Rds"
out_dir     <- "/projects/songli_lab/RazanPlantPath_2025/Script_New_Cluster_ID/OutPut"
markers_csv <- file.path(out_dir, "Cluster_Markers_Heatmap_07072025.csv")
heatmap_pdf <- file.path(out_dir, "Clustered_Marker_Heatmap_Trend_07072025.pdf")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------ Load Seurat object ------------------
seurat_obj <- readRDS(seurat_path)

# ------------------ Custom cluster labels in desired order ------------------
cluster_ordered_labels <- c(
  '0'  = "0:Meristem_1",
  '1'  = "1:Meristem_2",
  '2'  = "2:Meristem_3",
  '3'  = "3:Pericycle",
  '4'  = "4:Atrichoblast",
  '5'  = "5:Rootcap_1",
  '6'  = "6:Xylem",
  '7'  = "7:Rootcap_2",
  '8'  = "8:Cortex",
  '9'  = "9:Meristem_4",
  '10' = "10:Endodermis",
  '11' = "11:Trichoblast"
)

# Make sure Idents are cluster IDs
if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  warning("`seurat_clusters` not in meta.data; using current Idents().")
} else {
  Idents(seurat_obj) <- seurat_obj$seurat_clusters
}

# Check mapping covers all clusters
orig_clusters <- as.character(Idents(seurat_obj))
missing_levels <- setdiff(unique(orig_clusters), names(cluster_ordered_labels))
if (length(missing_levels) > 0) {
  stop("These clusters have no label mapping: ", paste(missing_levels, collapse = ", "))
}

# Apply mapping (no plyr, avoids namespace issues)
seurat_obj$renamed_clusters <- dplyr::recode(orig_clusters, !!!cluster_ordered_labels)

# lock order
desired_levels <- unname(cluster_ordered_labels)
seurat_obj$renamed_clusters <- factor(seurat_obj$renamed_clusters, levels = desired_levels)

# Use renamed clusters as identities
Idents(seurat_obj) <- seurat_obj$renamed_clusters

# ------------------ (Optional but recommended) Downsample per cluster ------------------
# This prevents the x-axis from filling with thousands of cell barcodes.
# Adjust 'downsample_per_cluster' to taste (e.g., 150–400).
downsample_per_cluster <- 250
seurat_small <- subset(seurat_obj, downsample = downsample_per_cluster)
Idents(seurat_small) <- seurat_small$renamed_clusters  # keep same labels/order

# ------------------ Find markers on the full object (better power) ------------------
markers <- FindAllMarkers(
  object = seurat_obj,
  only.pos = TRUE,
  min.pct = 0.01,
  logfc.threshold = 0.25
)
write.csv(markers, markers_csv, row.names = FALSE)

# ------------------ Top 10 per cluster ------------------
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    cluster_name = dplyr::recode(as.character(cluster), !!!cluster_ordered_labels),
    cluster_name = factor(cluster_name, levels = desired_levels)
  )

# Sort genes by your cluster trend then by log2FC
sorted_genes <- top_markers %>%
  arrange(cluster_name, desc(avg_log2FC)) %>%
  pull(gene) %>%
  unique()

# Keep only genes present
features_to_plot <- intersect(sorted_genes, rownames(seurat_small))
if (length(features_to_plot) == 0) stop("None of the selected marker genes are present.")

# ------------------ Heatmap ------------------
# Key changes vs the “garbled” plot:
#   - raster=TRUE: fast & crisp for many cells
#   - draw.lines=FALSE: removes white separators
#   - remove x-axis text/ticks (these are *cell* labels, not clusters)
#   - big fonts (14), clean theme
heatmap_plot <- DoHeatmap(
  seurat_small,
  features    = features_to_plot,
  group.by    = "renamed_clusters",
  raster      = TRUE,
  draw.lines  = FALSE,
  size        = 3
) +
  ggtitle("Top Marker Genes by Cluster (Trend Ordered)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
    # Hide the cell barcodes on x-axis (these cause the stair-step text)
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    # Make gene labels readable
    axis.text.y  = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    panel.grid   = element_blank() # cleaner background (no grid)
  )

# ------------------ Save PDF ------------------
pdf(heatmap_pdf, width = 12, height = 10)
print(heatmap_plot)
dev.off()

cat("✅ Heatmap saved to:\n", heatmap_pdf, "\n")
cat("✅ Full marker table saved to:\n", markers_csv, "\n")

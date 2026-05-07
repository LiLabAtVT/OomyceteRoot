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

markers_csv     <- file.path(out_dir, "all_markers_dotplot_analysis.csv")
top_markers_csv <- file.path(out_dir, "top_markers_dotplot_analysis.csv")
dotplot_pdf     <- file.path(out_dir, "DotPlot_Top_Marker_Genes_by_Cluster_04152026.pdf")

# ===============================================================
# Cluster labels
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

# ===============================================================
# Set identities
# ===============================================================
if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  Idents(seurat_obj) <- seurat_obj$seurat_clusters
}

orig_clusters <- as.character(Idents(seurat_obj))

missing_levels <- setdiff(unique(orig_clusters), names(cluster_labels))
if (length(missing_levels) > 0) {
  stop("Missing cluster labels for: ", paste(missing_levels, collapse = ", "))
}

seurat_obj$renamed_clusters <- dplyr::recode(orig_clusters, !!!cluster_labels)

desired_levels <- unname(cluster_labels)
seurat_obj$renamed_clusters <- factor(seurat_obj$renamed_clusters, levels = desired_levels)

Idents(seurat_obj) <- seurat_obj$renamed_clusters

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
    ungroup()
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

cat("✅ Genes selected:", length(features_to_plot), "\n")

# ===============================================================
# Dot plot
# ===============================================================
dot_plot <- DotPlot(
  object = seurat_obj,
  features = features_to_plot,
  group.by = "renamed_clusters"
) +
  theme_minimal(base_size = 28) +
  ggtitle("Top Marker Genes by Cluster") +

  # tighten spacing
  scale_y_discrete(expand = expansion(mult = c(0.02, 0.05))) +

  theme(
    

    # Gene names
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 25,
      face = "bold",
      margin = margin(t = -12)
    ),

    # Cluster labels
    axis.text.y = element_text(
      size = 30
      # face = "bold"
    ),

    axis.title.x = element_blank(),
    axis.title.y = element_blank(),

    legend.title = element_text(
      size = 26
      # face = "bold"
    ),

    legend.text = element_text(size = 22),

    panel.grid = element_blank()
  )

# Save PDF
pdf(dotplot_pdf, width = 40, height = 15)
print(dot_plot)
dev.off()

cat("Dot plot saved to:", dotplot_pdf, "\n")
cat("Marker CSV saved:", markers_csv, "\n")
cat("Top marker CSV saved:", top_markers_csv, "\n")

# --- Load required libraries ---
library(Seurat)
library(dplyr)
library(ggplot2)

# --- Define input and output paths ---
seurat_obj_path <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/Patho_Nonpatho_integration_07072025.Rds"
output_dir <- "/projects/songli_lab/RazanPlantPath_2025/Script_New_Cluster_ID/OutPut"

# Define output path for UMAP PDF and new Seurat object with cluster labels
output_umap_pdf <- file.path(output_dir, "UMAP_cluster_labels_ID_07072025.pdf")
output_seurat_path <- sub("\\.Rds$", "_Cluster_ID_07072025.Rds", seurat_obj_path)

# --- Ensure output directory exists ---
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# --- Load Seurat object ---
seurat_obj <- readRDS(seurat_obj_path)

# --- Define cluster labels (match cluster IDs to biological interpretation) ---
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

# --- Assign readable labels to the Seurat object ---
seurat_obj$cluster_labels <- factor(
  seurat_obj$seurat_clusters,
  levels = names(cluster_labels),
  labels = cluster_labels
)

# --- Save updated Seurat object ---
saveRDS(seurat_obj, output_seurat_path)
cat("✅ Seurat object with cluster labels saved to:\n", output_seurat_path, "\n")

# --- Create UMAP plot grouped by cluster with labels and custom legend ---
umap_plot <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by  = "seurat_clusters",
  label     = TRUE,
  repel     = TRUE
) +
  labs(title = "UMAP") +
  theme_classic(base_size = 14) +   # sets global font size & removes grids
  theme(
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_blank()
  ) +
  scale_color_manual(
    values = scales::hue_pal()(length(cluster_labels)),
    labels = cluster_labels
  )

# --- Save the UMAP plot as PDF ---
pdf(output_umap_pdf, width = 9, height = 6.5)
print(umap_plot)
dev.off()

cat("✅ UMAP plot saved to:\n", output_umap_pdf, "\n")

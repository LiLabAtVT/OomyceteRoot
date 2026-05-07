# ===============================================================
# Relative fraction of cells per cluster and condition
# Save Seurat object before generating plot
# ===============================================================

# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plyr)

# ===============================================================
# Define input and output paths
# ===============================================================
seurat_path <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/Patho_Nonpatho_integration_07072025.Rds"

out_dir <- "/projects/intro2gds/Razan_2026/PLant_Pathogen_2026/Step_50_Dot_Plot_Heat_map"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Save processed Seurat object here
seurat_save_path <- "/projects/intro2gds/Razan_2026/PLant_Pathogen_2026/Step_50_Dot_Plot_Heat_map/heatmap_trend.rds"
dir.create(dirname(seurat_save_path), showWarnings = FALSE, recursive = TRUE)

# ===============================================================
# Load Seurat object
# ===============================================================
seurat_obj <- readRDS(seurat_path)

# ===============================================================
# Remove duplicated metadata columns (if any)
# ===============================================================
dup_cols <- duplicated(colnames(seurat_obj@meta.data))
if (any(dup_cols)) {
  cat("⚠️ Removing duplicated metadata columns:\n")
  print(colnames(seurat_obj@meta.data)[dup_cols])
  seurat_obj@meta.data <- seurat_obj@meta.data[, !dup_cols]
}

# ===============================================================
# Assign simplified condition labels
# ===============================================================
seurat_obj$condition <- ifelse(
  grepl("Neg", seurat_obj$orig.ident, ignore.case = TRUE),
  "Non_Infected",
  "Infected"
)

# ===============================================================
# Define cluster labels (with IDs)
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
# Add cluster labels to metadata
# ===============================================================
seurat_obj$cluster_label <- plyr::revalue(
  as.character(seurat_obj$seurat_clusters),
  replace = cluster_labels,
  warn_missing = FALSE
)

# ===============================================================
# Save Seurat object 
# ===============================================================
saveRDS(seurat_obj, seurat_save_path)
cat("✅ Seurat object saved to:\n", seurat_save_path, "\n")

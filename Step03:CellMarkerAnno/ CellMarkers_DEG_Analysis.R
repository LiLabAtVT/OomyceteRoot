# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Define file paths
seurat_obj_path <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/Patho_Nonpatho_integration_07072025.Rds"
output_all_markers <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/09292025/10022025_GO_DEG/all_markers_07072025.csv"
output_deg_markers <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/09292025/10022025_GO_DEG/deg_markers_07072025.csv"

# Load Seurat object
seurat_obj <- readRDS(seurat_obj_path)


# Identify **all marker genes** for clusters
all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.25)

# Save **all marker genes** to CSV
write.csv(all_markers, output_all_markers, row.names = FALSE)
cat("✅ All marker genes saved to:", output_all_markers, "\n")

# Filter **DEG markers** (Significant markers with p-value < 0.05)
deg_markers <- all_markers %>% filter(p_val_adj < 0.05)

# Save **DEG markers** to CSV
write.csv(deg_markers, output_deg_markers, row.names = FALSE)
cat("✅ DEG marker genes saved to:", output_deg_markers, "\n")



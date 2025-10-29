# =================== Paths ===================
seurat_path <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/Patho_Nonpatho_integration_07072025.Rds"
out_dir     <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/09292025"

# =================== Libraries ===================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(plyr)
})

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =================== Load Seurat object ===================
obj <- readRDS(seurat_path)

# Ensure UMAP exists
if (!("umap" %in% tolower(names(obj@reductions)))) {
  message("No UMAP found; computing PCA + UMAP...")
  if (!("pca" %in% tolower(names(obj@reductions)))) {
    obj <- RunPCA(obj, verbose = FALSE)
  }
  obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
}

reductions_lower <- tolower(names(obj@reductions))
reduction_to_use <- if ("umap" %in% reductions_lower) "umap" else names(obj@reductions)[1]
message(sprintf("Using reduction: %s", reduction_to_use))

# =================== Cluster labels ===================
cluster_labels <- c(
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

# =================== Assign Idents ===================
if (!all(names(cluster_labels) %in% levels(Idents(obj)))) {
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    Idents(obj) <- factor(obj$seurat_clusters)
  }
}
obj$celltype_label <- plyr::revalue(as.character(Idents(obj)), replace = cluster_labels, warn_missing = FALSE)
obj$celltype_label <- factor(obj$celltype_label, levels = unname(cluster_labels))
Idents(obj) <- obj$celltype_label

# =================== TAIR Marker Sets ===================
marker_map <- list(
  "0: Meristematic_1" = c("AT3G20840","AT1G51190","AT4G37490"),
  "1: Meristem_2"      = c("AT3G11260","AT1G74500","AT3G18060"),
  "2: Meristem_2"      = c("AT5G03150","AT3G21670","AT1G09750"),
  "3: Pericycle"         = c("AT1G79430","AT1G22710","AT1G26880","AT3G01680"),
  "4: Atrichoblast"       = c("AT1G79840","AT5G14750"),
  "5: Root_cap_1"       = c("AT1G79580","AT1G33280","AT4G10350","AT1G51190"),
  "6: Xylem"           = c("AT5G62380","AT1G71930","AT5G17420","AT2G17530","AT4G35350"),
  "7: Root_cap_2"     = c("AT1G79580","AT1G33280","AT4G10350"),
  "8: Cortex"         = c("AT1G09750","AT1G62500","AT5G03150"),
  "9: Meristem_4"     = c("AT1G71930","AT5G17420","AT1G20850"),
  "10: Endodermis"    = c("AT2G36100","AT3G54220","AT5G57660"),
  "11: Trichoblasts"  = c("AT1G27740","AT1G66470","AT1G12560","AT1G63650")
)

# =================== Resolve Features ===================
all_feats <- rownames(obj)
resolve_present <- function(cands, features) cands[cands %in% features]

resolved_markers <- lapply(names(marker_map), function(lbl) {
  present <- resolve_present(marker_map[[lbl]], all_feats)
  tibble(cluster_label = lbl, feature = present)
}) %>% bind_rows()

# Save marker index
readr::write_csv(resolved_markers, file.path(out_dir, "FeaturePlot_marker_index.csv"))

# =================== Plot Function ===================
plot_cluster_markers <- function(object, cluster_label, features, reduction = "umap") {
  if (length(features) == 0) return(NULL)
  Idents(object) <- object$celltype_label
  p <- FeaturePlot(
    object, features = features, reduction = reduction,
    order = TRUE, pt.size = 0.5, ncol = min(4, length(features))
  ) +
    ggplot2::labs(title = paste0("Markers: ", cluster_label)) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12))
  p
}

# =================== Generate FeaturePlots ===================
message("Generating FeaturePlots (no split)...")
for (lbl in names(marker_map)) {
  feats <- resolved_markers %>% filter(cluster_label == lbl) %>% pull(feature)
  for (gene in feats) {
    p <- FeaturePlot(obj, features = gene, reduction = reduction_to_use,
                     order = TRUE, pt.size = 0.5) +
         ggplot2::labs(title = paste0(lbl, " – ", gene)) +
         ggplot2::theme(plot.title = element_text(face = "bold", size = 12))
    outfile <- file.path(out_dir,
                         paste0("FeaturePlot_", gsub("[^A-Za-z0-9_]+","_", lbl),
                                "_", gene, ".png"))
    ggsave(outfile, plot = p, width = 5, height = 5, dpi = 300)
  }
}

message("✅ FeaturePlots saved successfully.")

# =================== Paths ===================
seurat_path <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/Patho_Nonpatho_integration_07072025.Rds"
out_dir     <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/09292025/FeaturePlots_MarkergenWithlegend_V1"

# =================== Libraries ===================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.At.tair.db)
})

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =================== Load Seurat object ===================
obj <- readRDS(seurat_path)

# Ensure UMAP exists
if (!("umap" %in% tolower(names(obj@reductions)))) {
  message("No UMAP found; computing PCA + UMAP...")
  if (!("pca" %in% tolower(names(obj@reductions)))) obj <- RunPCA(obj, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
}
red_names <- names(obj@reductions)
reduction_to_use <- if ("umap" %in% tolower(red_names)) red_names[grep("umap", tolower(red_names))[1]] else red_names[1]
message(sprintf("Using reduction: %s", reduction_to_use))

# =================== TAIR Marker Sets ===================
marker_map <- list(
  "0: Meristematic_1" = c("AT3G20840","AT1G51190","AT4G37490"),
  "1: Unknown_1"      = c("AT3G11260","AT1G74500","AT3G18060"),
  "2: Unknown_2"      = c("AT5G03150","AT3G21670","AT1G09750"),
  "3: Phloem"         = c("AT1G79430","AT1G22710","AT1G26880","AT3G01680"),
  "4: Non_hair"       = c("AT1G79840","AT5G14750"),
  "5: Root_cap"       = c("AT1G79580","AT1G33280","AT4G10350","AT1G51190"),
  "6: Xylem_1"        = c("AT5G62380","AT1G71930","AT5G17420","AT2G17530","AT4G35350"),
  "7: Root_cap_2"     = c("AT1G79580","AT1G33280","AT4G10350"),
  "8: Cortex"         = c("AT1G09750","AT1G62500","AT5G03150"),
  "9: Xylem_2"        = c("AT1G71930","AT5G17420","AT1G20850"),
  "10: Endodermis"    = c("AT2G36100","AT3G54220","AT5G57660"),
  "11: Trichoblasts"  = c("AT1G27740","AT1G66470","AT1G12560","AT1G63650")
)

# =================== TAIR -> SYMBOL for nice titles ===================
all_tair <- unique(unlist(marker_map))
ann <- suppressMessages(AnnotationDbi::select(
  org.At.tair.db, keys = all_tair, keytype = "TAIR", columns = c("SYMBOL")
)) %>% dplyr::distinct(TAIR, SYMBOL)
tair2sym <- setNames(ann$SYMBOL, ann$TAIR)

pretty_title <- function(tair) {
  sym <- tair2sym[[tair]]
  if (!is.null(sym) && !is.na(sym) && nchar(sym) > 0) paste0(sym, " (", tair, ")") else tair
}

# =================== Gather features to plot (present only) ===================
all_feats <- rownames(obj)
features_requested <- unique(unlist(marker_map))
found <- intersect(features_requested, all_feats)
if (length(found) == 0) stop("❌ None of the requested marker genes are present in the object rownames.")

# ---- one-panel FeaturePlots (no split), consistent 0–3 scale ----
mincut <- 0
maxcut <- 3

# Save **separate files** per gene (PNG and PDF)
for (g in found) {
  p <- FeaturePlot(
    obj,
    features   = g,
    reduction  = reduction_to_use,
    order      = TRUE,
    min.cutoff = mincut,   # numeric caps avoid quantile/NA issues
    max.cutoff = maxcut,
    combine    = TRUE,
    keep.scale = "feature"
  ) +
    scale_color_gradient(
      limits = c(mincut, maxcut),
      breaks = c(0, 1, 2, 3),
      low = "grey90", high = "blue4",
      name = "Expression"
    ) +
    guides(color = guide_colorbar(
      title = "Expression",
      title.position = "top",
      barwidth = 0.6,
      barheight = 4
    )) +
    labs(title = pretty_title(g)) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position  = "right",
      legend.title     = element_text(size = 10),
      legend.text      = element_text(size = 9),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid       = element_blank()
    )

  # filenames
  base <- file.path(out_dir, paste0("FeaturePlot_", g, "_capped0to3_NoSplit"))
  ggsave(paste0(base, ".png"), plot = p, width = 7.5, height = 6.5, dpi = 300)
  ggsave(paste0(base, ".pdf"), plot = p, width = 7.5, height = 6.5)
  message("✅ Saved: ", paste0(base, ".png"))
}

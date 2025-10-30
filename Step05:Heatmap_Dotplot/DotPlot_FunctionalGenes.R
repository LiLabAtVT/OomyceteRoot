# ================= DotPlot + Gene Annotation Table (robust, 1:many-safe) =================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(AnnotationDbi)
  library(org.At.tair.db)
})

# -------- Paths --------
seurat_obj_path    <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/Patho_Nonpatho_integration_07072025.Rds"
deg_markers_path   <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/deg_markers_07072025.csv"
output_dotplot_pdf <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/DotPlot_Cell_Type_Markers_Cluster_ID-07072025.pdf"
output_annot_csv   <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/DotPlot_Cell_Type_Markers_Annotations_07072025.csv"

dir.create(dirname(output_dotplot_pdf), recursive = TRUE, showWarnings = FALSE)

# -------- Load Seurat object --------
obj <- readRDS(seurat_obj_path)

# If rownames have prefixes, strip them once
rownames(obj) <- gsub("^Arabidopsis_thaliana_gene:|^gene:", "", rownames(obj))

# -------- Cluster labels (x-axis order) --------
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
obj$cluster_labels <- factor(
  obj$seurat_clusters,
  levels = names(cluster_labels),
  labels = cluster_labels
)

# -------- Load DEG markers (preserve names exactly as in file) --------
deg <- readr::read_csv(deg_markers_path, show_col_types = FALSE)

# Case-insensitive column resolver
get_col <- function(df, want) {
  nl <- tolower(names(df))
  idx <- match(tolower(want), nl)
  if (is.na(idx)) stop("Missing required column: ", want)
  names(df)[idx]
}

col_cluster   <- get_col(deg, "cluster")
col_gene      <- get_col(deg, "gene")
col_avglog2fc <- get_col(deg, "avg_log2FC")  # as in your file

# Build a working frame with canonical names
deg2 <- deg %>%
  dplyr::transmute(
    cluster    = .data[[col_cluster]],
    gene       = gsub("^gene:", "", .data[[col_gene]]),
    avg_log2FC = suppressWarnings(as.numeric(.data[[col_avglog2fc]]))
  )

# -------- Top N per cluster --------
TOP_N <- 5
top_tbl <- deg2 %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = TOP_N, with_ties = FALSE) %>%
  dplyr::ungroup()

# Keep only genes present in the object
valid_genes <- intersect(unique(top_tbl$gene), rownames(obj))
missing_in_obj <- setdiff(unique(top_tbl$gene), valid_genes)
if (length(valid_genes) == 0) stop("No selected marker genes are present in the Seurat object.")
if (length(missing_in_obj)) message("⚠️ Not found in object (skipped): ", paste(missing_in_obj, collapse = ", "))

# -------- Annotation: TAIR -> SYMBOL, GENENAME (collapse 1:many) --------
anno_raw <- suppressWarnings(
  AnnotationDbi::select(
    org.At.tair.db,
    keys    = unique(valid_genes),
    keytype = "TAIR",
    columns = c("SYMBOL","GENENAME")
  )
)

# Handle possible 1:many mappings: collapse to one row per TAIR
anno_collapsed <- anno_raw %>%
  dplyr::filter(!is.na(TAIR)) %>%
  dplyr::group_by(TAIR) %>%
  dplyr::summarise(
    SYMBOL   = paste(unique(na.omit(SYMBOL)), collapse = "; "),
    GENENAME = paste(unique(na.omit(GENENAME)), collapse = "; "),
    .groups = "drop"
  )

# -------- Build annotation table (keep 'gene' until after the join) --------
annot_out <- top_tbl %>%
  dplyr::filter(gene %in% valid_genes) %>%
  dplyr::mutate(cluster_label = dplyr::recode(as.character(cluster), !!!cluster_labels)) %>%
  dplyr::left_join(anno_collapsed, by = c("gene" = "TAIR")) %>%
  dplyr::mutate(TAIR = gene) %>%                     # now safely create TAIR column
  dplyr::select(cluster, cluster_label, TAIR, SYMBOL, GENENAME, avg_log2FC) %>%
  dplyr::arrange(suppressWarnings(as.integer(as.character(cluster))), dplyr::desc(avg_log2FC))

readr::write_csv(annot_out, output_annot_csv)
message("✅ Annotation table: ", output_annot_csv)

# -------- DotPlot (readable fonts/order; no grid) --------
# Lock gene order to table order (cluster then log2FC)
ordered_genes <- unique(annot_out$TAIR)

p <- DotPlot(
  obj,
  features  = ordered_genes,
  group.by  = "cluster_labels",
  cols      = c("lightgrey", "darkred"),
  dot.scale = 6
) +
  ggplot2::scale_color_gradient(low = "blue", high = "red") +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::labs(
    title = "Cell Type-Specific Marker Genes",
    x = "Marker Genes", y = "Clusters"
  ) +
  ggplot2::theme(
    plot.title  = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid  = ggplot2::element_blank()
  )

pdf(output_dotplot_pdf, width = 12, height = 8)
print(p)
dev.off()
message("✅ DotPlot: ", output_dotplot_pdf)

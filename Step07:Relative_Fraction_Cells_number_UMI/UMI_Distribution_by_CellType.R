# ===============================================================
# Dissemination of Total UMI by Cell Type AND Condition (auto-detect condition)
# ===============================================================

# ---------- Paths ----------
seurat_path <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/Patho_Nonpatho_integration_07072025.Rds"
out_dir     <- "/projects/songli_lab/RazanPlantPath_2025/Script_New_Cluster_ID/OutPut"

# ---------- Libraries ----------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- Load ----------
obj <- readRDS(seurat_path)
message("✅ Seurat object loaded.")

# ---------- Assign cell-type labels ----------
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

if (!all(names(cluster_labels) %in% levels(Idents(obj)))) {
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    Idents(obj) <- as.factor(obj$seurat_clusters)
  }
}

id_chr <- as.character(Idents(obj))
mapped <- ifelse(id_chr %in% names(cluster_labels), cluster_labels[id_chr], id_chr)
obj$celltype_label <- factor(mapped, levels = unname(cluster_labels))
Idents(obj) <- obj$celltype_label
message("✅ Cell-type labels assigned.")

# ---------- Detect UMI column ----------
meta_cols <- colnames(obj@meta.data)
if ("nCount_SCT" %in% meta_cols) {
  umi_col <- "nCount_SCT"
} else if ("nCount_RNA" %in% meta_cols) {
  umi_col <- "nCount_RNA"
} else {
  stop("❌ Could not find total UMI column (nCount_RNA or nCount_SCT).")
}
message(sprintf("📦 Using %s for total UMI.", umi_col))

# ---------- Auto-detect / create `condition` ----------
infer_condition <- function(df) {
  # candidate columns to inspect, in priority order
  cands <- c("condition","Condition","group","Group","status","Status",
             "infection","Infection","treatment","Treatment",
             "orig.ident","orig_ident","sample","Sample","sample_id","dataset","Dataset")
  cands <- cands[cands %in% colnames(df)]
  if (length(cands) == 0) return(NULL)

  # pick the first with >1 unique value
  for (cn in cands) {
    vals <- as.character(df[[cn]])
    if (length(unique(vals)) > 1) {
      return(list(col = cn, vals = vals))
    }
  }
  return(NULL)
}

ic <- infer_condition(obj@meta.data)

if (is.null(ic)) {
  stop("❌ Could not infer a condition column. Please tell me which metadata column encodes Infected vs Non_infected.")
} else {
  src_col <- ic$col
  raw_vals <- ic$vals
  message(sprintf("🔎 Inferring condition from '%s'.", src_col))

  # map common patterns to Infected / Non_infected
  to_infected <- str_detect(raw_vals, regex("pos|infect|inf|treated|p24|pos24|poshpi|pos24hpi", ignore_case = TRUE))
  to_noninf   <- str_detect(raw_vals, regex("neg|mock|control|ctrl|non|untreated|neg24|neg24hpi|neghpi", ignore_case = TRUE))

  cond <- ifelse(to_infected, "Infected",
           ifelse(to_noninf, "Non_infected", NA))

  # if still NA, try a fallback: two most common values → map by name heuristics
  if (any(is.na(cond))) {
    uniq <- unique(raw_vals)
    # if exactly 2 groups, label them deterministically by regex on the group name
    if (length(uniq) == 2) {
      g1 <- uniq[1]; g2 <- uniq[2]
      f <- function(x) {
        if (str_detect(x, regex("pos|infect|inf|treated", ignore_case = TRUE))) "Infected" else
          if (str_detect(x, regex("neg|mock|control|ctrl|non|untreated", ignore_case = TRUE))) "Non_infected" else x
      }
      map1 <- f(g1); map2 <- f(g2)
      # if still unmapped, assign alphabetically for reproducibility
      if (map1 == g1 && map2 == g2) {
        maps <- setNames(c("Infected","Non_infected"), sort(uniq))
        cond <- unname(maps[raw_vals])
      } else {
        maps <- setNames(c(map1, map2), c(g1, g2))
        cond <- unname(maps[raw_vals])
      }
    }
  }

  # final check
  if (any(is.na(cond))) {
    msg <- paste0(
      "❌ Could not confidently map all values in '", src_col, "' to conditions.\n",
      "Unique values detected: ", paste(unique(raw_vals), collapse = ", "), "\n",
      "Please tell me which column is the condition (and value mapping), or rename/create obj$condition yourself."
    )
    stop(msg)
  }

  obj$condition <- factor(cond, levels = c("Non_infected","Infected"))
  message("✅ condition column created (levels: Non_infected, Infected).")
  # Optional: write a mapping table for transparency
  map_tbl <- tibble(!!src_col := raw_vals, condition = obj$condition)
  write_csv(distinct(map_tbl), file.path(out_dir, "Condition_Inference_Mapping.csv"))
  message("💾 Saved Condition_Inference_Mapping.csv")
}

# ---------- Summaries by cell type AND condition ----------
umi_summary <- obj@meta.data %>%
  dplyr::group_by(celltype_label, condition) %>%
  dplyr::summarise(
    median_UMI = stats::median(.data[[umi_col]]),
    mean_UMI   = base::mean(.data[[umi_col]]),
    sd_UMI     = stats::sd(.data[[umi_col]]),
    min_UMI    = base::min(.data[[umi_col]]),
    max_UMI    = base::max(.data[[umi_col]]),
    n_cells    = dplyr::n(),
    .groups    = "drop"
  ) %>%
  dplyr::arrange(celltype_label, condition)

write_csv(umi_summary, file.path(out_dir, "UMI_Distribution_by_CellType_and_Condition.csv"))
message("💾 Saved UMI_Distribution_by_CellType_and_Condition.csv")

# ---------- Violin (split by condition) ----------
vln <- VlnPlot(
  obj,
  features = umi_col,
  group.by = "celltype_label",
  split.by = "condition",
  pt.size = 0
) +
  ggplot2::labs(
    title = "Total UMI per Cell Type (Infected vs Non_infected)",
    x = "Cell Type",
    y = paste0("Total UMI (", umi_col, ")")
  ) +
  ggplot2::theme_bw(base_size = 14) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
  )

ggsave(file.path(out_dir, "VlnPlot_UMI_by_CellType_and_Condition.png"),
       plot = vln, width = 12, height = 6, dpi = 300)
message("🎨 Saved VlnPlot_UMI_by_CellType_and_Condition.png")

# ---------- Boxplot (log10) ----------
meta_df <- obj@meta.data
p_box <- ggplot2::ggplot(meta_df, ggplot2::aes(
  x = celltype_label,
  y = .data[[umi_col]],
  fill = condition
)) +
  ggplot2::geom_boxplot(outlier.size = 0.4, alpha = 0.8) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    title = "Dissemination of Total UMI by Cell Type and Condition",
    x = "Cell Type",
    y = paste0("Total UMI (log10 ", umi_col, ")")
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
  )

ggsave(file.path(out_dir, "BoxPlot_UMI_by_CellType_and_Condition.png"),
       plot = p_box, width = 12, height = 6, dpi = 300)
message("🎨 Saved BoxPlot_UMI_by_CellType_and_Condition.png")

message("✅ Done. Outputs written to: ", out_dir)


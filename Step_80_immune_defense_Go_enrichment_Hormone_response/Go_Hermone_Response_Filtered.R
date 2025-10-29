# ===============================================================
# GO Enrichment: Hormone Response (Infected vs Non-Infected)
# OLD layout + remove GO terms where -log10(p.adjust) < 2 (p.adjust > 0.01)
# ===============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(clusterProfiler)
  library(org.At.tair.db)
  library(ggplot2)
  library(stringr)
  library(readr)
  library(forcats)
})

# ---------- Paths ----------
seurat_path <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/Patho_Nonpatho_integration_07072025.Rds"
output_dir  <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/09292025"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

out_csv <- file.path(output_dir, "GO_Enrichment_Hormone_Response_filtered.csv")
out_pdf <- file.path(output_dir, "GO_Enrichment_Hormone_Response_filtered.pdf")

# ---------- Load Seurat object ----------
obj <- readRDS(seurat_path)

if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
  stop("`seurat_clusters` not found in meta.data.")
}
Idents(obj) <- obj$seurat_clusters

# Define infection condition
obj$condition <- ifelse(grepl("Neg", obj$orig.ident, ignore.case = TRUE),
                        "Non_Infected", "Infected")

# ---------- Cluster ID → label map ----------
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
cluster_x_levels <- unname(cluster_labels)

clusters <- as.character(sort(unique(obj$seurat_clusters)))
all_go <- list()

# ===============================================================
# Loop through clusters — Compare Infected vs Non-Infected
# ===============================================================
for (cl in clusters) {
  cat("\n🚀 Running GO enrichment for cluster", cl, "\n")

  cl_cells <- WhichCells(obj, idents = cl)
  if (!length(cl_cells)) next

  sub <- subset(obj, cells = cl_cells)
  sub$condition <- ifelse(grepl("Neg", sub$orig.ident, ignore.case = TRUE),
                          "Non_Infected", "Infected")
  Idents(sub) <- sub$condition

  deg <- FindMarkers(
    sub,
    ident.1 = "Non_Infected",
    ident.2 = "Infected",
    only.pos = FALSE,
    min.pct = 0.01,
    logfc.threshold = 0.25
  )
  deg$gene <- rownames(deg)

  up_infected    <- deg %>% filter(p_val_adj < 0.05 & avg_log2FC > 0)
  up_noninfected <- deg %>% filter(p_val_adj < 0.05 & avg_log2FC < 0)

  # --- Infected ---
  if (nrow(up_infected) >= 10) {
    conv_inf <- suppressWarnings(bitr(up_infected$gene, fromType="TAIR", toType="ENTREZID", OrgDb=org.At.tair.db))
    if (!is.null(conv_inf) && nrow(conv_inf)) {
      ego_inf <- enrichGO(conv_inf$ENTREZID, OrgDb=org.At.tair.db, keyType="ENTREZID",
                          ont="BP", pAdjustMethod="BH", qvalueCutoff=0.3, readable=TRUE)
      df_inf <- as.data.frame(ego_inf)
      if (nrow(df_inf)) {
        df_inf$group <- "Infected"
        df_inf$cluster <- cl
        all_go[[paste0("Infected_", cl)]] <- df_inf
      }
    }
  }

  # --- Non-Infected ---
  if (nrow(up_noninfected) >= 10) {
    conv_non <- suppressWarnings(bitr(up_noninfected$gene, fromType="TAIR", toType="ENTREZID", OrgDb=org.At.tair.db))
    if (!is.null(conv_non) && nrow(conv_non)) {
      ego_non <- enrichGO(conv_non$ENTREZID, OrgDb=org.At.tair.db, keyType="ENTREZID",
                          ont="BP", pAdjustMethod="BH", qvalueCutoff=0.3, readable=TRUE)
      df_non <- as.data.frame(ego_non)
      if (nrow(df_non)) {
        df_non$group <- "Non_Infected"
        df_non$cluster <- cl
        all_go[[paste0("NonInfected_", cl)]] <- df_non
      }
    }
  }
}

if (!length(all_go)) stop("❌ No enriched GO terms found across clusters.")

go_all <- bind_rows(all_go)
readr::write_csv(go_all, out_csv)
cat("✅ GO results saved to:", out_csv, "\n")

# ===============================================================
# Filter for hormone-related terms & remove p.adjust > 0.01
# ===============================================================
hormone_pattern <- paste0(
  "(?i)",
  "(",
    "(response|signaling|signalling).*\\b(hormone|auxin|abscisic|ABA|jasmonic|JA|salicylic|SA|ethylene|ET|gibberellin|GA|cytokinin|CK|brassinosteroid|BR|strigolactone)\\b",
    "|",
    "\\b(hormone|auxin|abscisic|ABA|jasmonic|JA|salicylic|SA|ethylene|ET|gibberellin|GA|cytokinin|CK|brassinosteroid|BR|strigolactone)\\b.*(response|signaling|signalling)",
  ")"
)

hormone_terms <- go_all %>%
  filter(str_detect(Description, regex(hormone_pattern))) %>%
  filter(!is.na(p.adjust) & p.adjust <= 0.01) %>%  # remove -log10(p)<2
  group_by(cluster, group) %>%
  slice_min(p.adjust, n = 10, with_ties = FALSE) %>%
  ungroup()

if (!nrow(hormone_terms)) stop("❌ No hormone-response GO terms found (after p ≤ 0.01 filter).")

# ---------- Keep OLD order and cluster labeling ----------
desc_levels <- unique(hormone_terms$Description)
hormone_terms$Description <- factor(hormone_terms$Description, levels = rev(desc_levels))

hormone_terms$cluster <- as.character(hormone_terms$cluster)
hormone_terms$cluster_label <- cluster_labels[hormone_terms$cluster]
hormone_terms$cluster_label[is.na(hormone_terms$cluster_label)] <- hormone_terms$cluster
hormone_terms$cluster_label <- factor(hormone_terms$cluster_label, levels = cluster_x_levels)

# ===============================================================
# Plot: no triangles; just circles; remove all points with p>0.01
# ===============================================================
p <- ggplot(hormone_terms, aes(
  x = cluster_label,
  y = Description,
  color = group,
  size = -log10(p.adjust)
)) +
  geom_point(alpha = 0.85, shape = 16) +  # only circles
  scale_color_manual(values = c("Infected" = "firebrick3", "Non_Infected" = "steelblue3")) +
  theme_minimal(base_size = 12) +
  labs(
    title = "GO Enrichment: Hormone Response (Infected vs Non–Infected, p ≤ 0.01)",
    x = "Cluster (ID + Label)",
    y = "GO Biological Process",
    size = expression(-log[10]*"(adj. p-value)"),
    color = "Condition"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5, face = "bold")
  )

ggsave(out_pdf, plot = p, width = 15, height = 11)
cat("✅ Hormone-response dotplot (filtered, p ≤ 0.01, circles only) saved to:", out_pdf, "\n")

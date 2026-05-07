# ===============================================================
# GO Enrichment: Hormone Response (Infected vs Non-Infected)
# Save filtered results as RDS for faster plotting later
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
  library(tidyr)
})

# ---------- Paths ----------
seurat_path <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/Patho_Nonpatho_integration_07072025.Rds"
output_dir  <- "/projects/intro2gds/Razan_2026/PLant_Pathogen_2026/Step_80_immune_defense_Go_enrichment_Hormone_response"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

out_csv_full     <- file.path(output_dir, "GO_Enrichment_Hormone_Response_07072025.csv")
out_csv_filtered <- file.path(output_dir, "GO_Enrichment_Hormone_Response_Filtered_GT3clusters_07072025.csv")
out_rds_filtered <- file.path(output_dir, "GO_Enrichment_Hormone_Response_Filtered_GT3clusters_07072025.rds")
out_pdf          <- file.path(output_dir, "GO_Enrichment_Hormone_Response_Filtered_GT3clusters_07072025.pdf")

# ---------- Load Seurat object ----------
obj <- readRDS(seurat_path)
obj$condition <- ifelse(grepl("Neg", obj$orig.ident, ignore.case = TRUE),
                        "Non_Infected", "Infected")
Idents(obj) <- obj$seurat_clusters

# ---------- Define cluster ID → label map ----------
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

all_cluster_ids <- names(cluster_labels)
all_cluster_levels <- unname(cluster_labels)

clusters <- levels(factor(as.character(obj$seurat_clusters), levels = all_cluster_ids))
all_go <- list()

# ===============================================================
# Loop through clusters — Compare Infected vs Non-Infected
# ===============================================================
for (cl in clusters) {
  cat("\nRunning GO enrichment for cluster", cl, "\n")

  sub <- subset(obj, idents = cl)
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

  if (nrow(up_infected) >= 10) {
    conv_inf <- bitr(
      up_infected$gene,
      fromType = "TAIR",
      toType   = "ENTREZID",
      OrgDb    = org.At.tair.db
    )

    if (!is.null(conv_inf) && nrow(conv_inf) > 0) {
      ego_inf <- enrichGO(
        gene          = conv_inf$ENTREZID,
        OrgDb         = org.At.tair.db,
        ont           = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff  = 0.3,
        readable      = TRUE
      )

      df_inf <- as.data.frame(ego_inf)
      if (nrow(df_inf) > 0) {
        df_inf$group <- "Infected"
        df_inf$cluster <- cl
        all_go[[paste0("Infected_", cl)]] <- df_inf
      }
    }
  }

  if (nrow(up_noninfected) >= 10) {
    conv_non <- bitr(
      up_noninfected$gene,
      fromType = "TAIR",
      toType   = "ENTREZID",
      OrgDb    = org.At.tair.db
    )

    if (!is.null(conv_non) && nrow(conv_non) > 0) {
      ego_non <- enrichGO(
        gene          = conv_non$ENTREZID,
        OrgDb         = org.At.tair.db,
        ont           = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff  = 0.3,
        readable      = TRUE
      )

      df_non <- as.data.frame(ego_non)
      if (nrow(df_non) > 0) {
        df_non$group <- "Non_Infected"
        df_non$cluster <- cl
        all_go[[paste0("NonInfected_", cl)]] <- df_non
      }
    }
  }
}

if (length(all_go) == 0) stop("No enriched GO terms found across clusters.")

go_all <- bind_rows(all_go)
write_csv(go_all, out_csv_full)
cat("Full GO results saved to:", out_csv_full, "\n")

# ---------- Focus on hormone-response/signaling terms ----------
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
  group_by(cluster, group) %>%
  slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>%
  ungroup()

if (nrow(hormone_terms) == 0) stop("No hormone-response GO terms found.")

# ---------- Filter GO terms found in > 3 clusters ----------
term_cluster_freq <- hormone_terms %>%
  distinct(Description, cluster) %>%
  count(Description, name = "n_clusters")

hormone_terms_filtered <- hormone_terms %>%
  left_join(term_cluster_freq, by = "Description") %>%
  filter(n_clusters > 3)

if (nrow(hormone_terms_filtered) == 0) {
  stop("After filtering GO terms present in <= 3 clusters, no hormone terms remain.")
}

# ---------- Add labels before saving RDS ----------
hormone_terms_filtered$cluster <- as.character(hormone_terms_filtered$cluster)
hormone_terms_filtered$cluster_label <- cluster_labels[hormone_terms_filtered$cluster]
hormone_terms_filtered$cluster_label[is.na(hormone_terms_filtered$cluster_label)] <- hormone_terms_filtered$cluster
hormone_terms_filtered$cluster_label <- factor(
  hormone_terms_filtered$cluster_label,
  levels = all_cluster_levels
)

write_csv(hormone_terms_filtered, out_csv_filtered)
saveRDS(hormone_terms_filtered, out_rds_filtered)

cat("Filtered hormone CSV saved to:", out_csv_filtered, "\n")
cat("Filtered hormone RDS saved to:", out_rds_filtered, "\n")

# ---------- Check remaining terms per cluster ----------
cluster_summary <- hormone_terms_filtered %>%
  distinct(cluster, Description) %>%
  count(cluster, name = "n_terms_remaining")

cat("\nHormone GO terms remaining per cluster after filtering:\n")
print(cluster_summary)

missing_clusters <- setdiff(all_cluster_ids, unique(as.character(hormone_terms_filtered$cluster)))
if (length(missing_clusters) > 0) {
  cat("\nClusters with zero remaining hormone GO terms after filtering:\n")
  print(missing_clusters)
}

# ---------- Plot ----------
plot_df <- hormone_terms_filtered
plot_df$Description <- str_wrap(plot_df$Description, width = 60)

p <- ggplot(
  plot_df,
  aes(
    x = cluster_label,
    y = fct_reorder(Description, -log10(p.adjust)),
    color = group,
    shape = group,
    size = -log10(p.adjust)
  )
) +
  geom_point(alpha = 0.85) +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(values = c("Infected" = "firebrick3", "Non_Infected" = "steelblue3")) +
  theme_minimal(base_size = 12) +
  labs(
    title = "GO Enrichment: Hormone Response (terms in >3 clusters)",
    x = "Cluster (ID + Label)",
    y = "GO Biological Process",
    size = "-log10(adj. p-value)",
    color = "Condition",
    shape = "Condition"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 14)
  )

ggsave(out_pdf, plot = p, width = 15, height = 11)
cat("Filtered hormone-response dotplot saved to:", out_pdf, "\n")

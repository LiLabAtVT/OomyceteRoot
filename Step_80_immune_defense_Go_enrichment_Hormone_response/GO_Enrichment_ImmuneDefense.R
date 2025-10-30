# ===============================================================
# GO Enrichment: Defense / Oomycete Response (Infected vs Non-Infected)
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

out_csv <- file.path(output_dir, "GO_Enrichment_Defense_Oomycete_07072025_p01.csv")
out_pdf <- file.path(output_dir, "GO_Enrichment_Defense_Oomycete_DotPlot_07072025_p01.pdf")

# ---------- Load Seurat object ----------
obj <- readRDS(seurat_path)
Idents(obj) <- obj$seurat_clusters

# Add infection status labels
obj$condition <- ifelse(grepl("Neg", obj$orig.ident, ignore.case = TRUE),
                        "Non_Infected", "Infected")

# ---------- Define cluster ID → label map ----------
cluster_labels <- c(
  '0' = "0:Meristem_1",
  '1' = "1:Meristem_2", 
  '2' = "2:Meristem_3",
  '3' = "3:Pericycle", 
  '4' = "4:Atrichoblast", 
  '5' = "5:Rootcap_1", 
  '6' = "6:Xylem", 
  '7' = "7:Rootcap_2", 
  '8' = "8:Cortex", 
  '9' = "9:Meristem_4", 
  '10' = "10:Endodermis", 
  '11' = "11:Trichoblast"
)
clusters <- sort(unique(as.character(obj$seurat_clusters)))
all_go <- list()

# ===============================================================
# Loop through clusters — Compare Infected vs Non-Infected
# ===============================================================
for (cl in clusters) {
  cat("\n🚀 Running GO enrichment for cluster", cl, "\n")

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

  # --- Infected upregulated ---
  if (nrow(up_infected) >= 10) {
    conv_inf <- suppressWarnings(bitr(up_infected$gene, fromType = "TAIR", toType = "ENTREZID", OrgDb = org.At.tair.db))
    conv_inf <- conv_inf %>% filter(!is.na(ENTREZID))
    if (nrow(conv_inf) > 0) {
      ego_inf <- enrichGO(
        gene          = unique(conv_inf$ENTREZID),
        OrgDb         = org.At.tair.db,
        keyType       = "ENTREZID",
        ont           = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff  = 0.3,
        readable      = TRUE
      )
      df_inf <- as.data.frame(ego_inf)
      if (nrow(df_inf) > 0) {
        df_inf$group   <- "Infected"
        df_inf$cluster <- cl
        all_go[[paste0("Infected_", cl)]] <- df_inf
      }
    }
  }

  # --- Non-infected upregulated ---
  if (nrow(up_noninfected) >= 10) {
    conv_non <- suppressWarnings(bitr(up_noninfected$gene, fromType = "TAIR", toType = "ENTREZID", OrgDb = org.At.tair.db))
    conv_non <- conv_non %>% filter(!is.na(ENTREZID))
    if (nrow(conv_non) > 0) {
      ego_non <- enrichGO(
        gene          = unique(conv_non$ENTREZID),
        OrgDb         = org.At.tair.db,
        keyType       = "ENTREZID",
        ont           = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff  = 0.3,
        readable      = TRUE
      )
      df_non <- as.data.frame(ego_non)
      if (nrow(df_non) > 0) {
        df_non$group   <- "Non_Infected"
        df_non$cluster <- cl
        all_go[[paste0("NonInfected_", cl)]] <- df_non
      }
    }
  }
}

if (length(all_go) == 0) stop("❌ No enriched GO terms found across clusters.")

go_all <- bind_rows(all_go)
write_csv(go_all, out_csv)
cat("✅ GO results saved to:", out_csv, "\n")

# ===============================================================
# Filter: defense/immune/oomycete terms and keep only p.adjust ≤ 0.01
# ===============================================================
defense_terms <- go_all %>%
  filter(str_detect(
    Description,
    regex("defense|immune|oomycete|oxidative|pathogen|jasmonic|salicylic|fungus|hypersensitive", ignore_case = TRUE)
  )) %>%
  filter(!is.na(p.adjust) & p.adjust <= 0.01) %>%     # remove dots where -log10(p) < 2
  group_by(cluster, group) %>%
  slice_min(p.adjust, n = 10, with_ties = FALSE) %>%
  ungroup()

if (nrow(defense_terms) == 0) stop("❌ No defense-related GO terms remain after p.adjust ≤ 0.01 filter.")

# ---------- Add cluster labels ----------
defense_terms$cluster <- as.character(defense_terms$cluster)
defense_terms$cluster_label <- cluster_labels[defense_terms$cluster]
defense_terms$cluster_label[is.na(defense_terms$cluster_label)] <- defense_terms$cluster
defense_terms$cluster_label <- factor(defense_terms$cluster_label, levels = cluster_labels)

# ---------- Plot ----------
defense_terms$Description <- str_wrap(defense_terms$Description, width = 60)

p <- ggplot(defense_terms, aes(
  x = cluster_label,
  y = forcats::fct_reorder(Description, -log10(p.adjust)),
  color = group,
  size  = -log10(p.adjust)
)) +
  geom_point(alpha = 0.9) +
  scale_color_manual(values = c("Infected" = "firebrick3", "Non_Infected" = "steelblue3")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "GO Enrichment: Defense & Oomycete Response (Infected vs Non-Infected)",
    x = "Cluster (ID + Label)",
    y = "GO Biological Process",
    size  = expression(-log[10]*"(adj. p-value)"),
    color = "Condition"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 16)
  )

ggsave(out_pdf, plot = p, width = 15, height = 11)
cat("✅ GO Defense/Oomycete dotplot (p.adjust ≤ 0.01) saved to:", out_pdf, "\n")

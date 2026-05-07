# ===============================================================
# GO Enrichment per Cluster
# ===============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(clusterProfiler)
  library(org.At.tair.db)
  library(stringr)
})

# ---------- Paths ----------
seurat_path <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/Patho_Nonpatho_integration_07072025.Rds"
output_dir  <- "/projects/intro2gds/Razan_2026/PLant_Pathogen_2026/Step_40_GOTerm/Script"
rds_out     <- file.path(output_dir, "top_terms_combined_filtered.rds")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------- Set minimum number of clusters ----------
min_clusters <- 3

# ---------- Load Seurat object ----------
seurat_obj <- readRDS(seurat_path)
if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  Idents(seurat_obj) <- seurat_obj$seurat_clusters
}
cluster_ids <- sort(unique(Idents(seurat_obj)))

# ---------- Cluster ID → Label map ----------
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

# ---------- Initialize ----------
all_cluster_go <- list()

# ===============================================================
# Loop over clusters: DEGs + GO enrichment
# ===============================================================
for (clust in cluster_ids) {
  cat("\n🚀 Cluster", clust, "\n")

  cluster_cells <- WhichCells(seurat_obj, idents = clust)
  cluster_obj   <- subset(seurat_obj, cells = cluster_cells)

  if (!"orig.ident" %in% colnames(cluster_obj@meta.data)) {
    cat("❌ 'orig.ident' missing, skipping cluster", clust, "\n")
    next
  }

  cluster_obj$infection_status <- ifelse(
    grepl("Neg", cluster_obj$orig.ident, ignore.case = TRUE),
    "Non_Infected",
    "Infected"
  )
  Idents(cluster_obj) <- cluster_obj$infection_status

  de_markers <- FindMarkers(
    cluster_obj,
    ident.1 = "Non_Infected",
    ident.2 = "Infected",
    only.pos = FALSE,
    min.pct = 0.01,
    logfc.threshold = 0.25
  )

  de_markers$gene <- rownames(de_markers)

  signif_markers <- de_markers %>%
    filter(p_val_adj < 0.1) %>%
    mutate(
      cluster = as.character(clust),
      direction = ifelse(avg_log2FC > 0, "Up_Infected", "Up_NonInfected")
    )

  if (nrow(signif_markers) < 10) {
    cat("⚠️ Few DEGs, skipping GO for cluster", clust, "\n")
    next
  }

  gene_conversion <- suppressWarnings(
    bitr(
      signif_markers$gene,
      fromType = "TAIR",
      toType   = "ENTREZID",
      OrgDb    = org.At.tair.db
    )
  )

  merged_degs <- merge(signif_markers, gene_conversion, by.x = "gene", by.y = "TAIR")

  infected_up <- merged_degs %>%
    filter(avg_log2FC > 0) %>%
    pull(ENTREZID) %>%
    unique()

  non_infected_up <- merged_degs %>%
    filter(avg_log2FC < 0) %>%
    pull(ENTREZID) %>%
    unique()

  go_inf <- go_noninf <- data.frame()

  if (length(infected_up) > 0) {
    go_infected <- enrichGO(
      gene          = infected_up,
      OrgDb         = org.At.tair.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.3,
      readable      = TRUE
    )

    go_inf <- as.data.frame(go_infected)

    if (nrow(go_inf) > 0) {
      go_inf$group   <- "Infected"
      go_inf$cluster <- as.character(clust)
    }
  }

  if (length(non_infected_up) > 0) {
    go_non_infected <- enrichGO(
      gene          = non_infected_up,
      OrgDb         = org.At.tair.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.3,
      readable      = TRUE
    )

    go_noninf <- as.data.frame(go_non_infected)

    if (nrow(go_noninf) > 0) {
      go_noninf$group   <- "Non_Infected"
      go_noninf$cluster <- as.character(clust)
    }
  }

  cluster_go_df <- bind_rows(go_inf, go_noninf)

  if (nrow(cluster_go_df) > 0) {
    all_cluster_go[[as.character(clust)]] <- cluster_go_df
    cat("✅ GO terms added for cluster", clust, "\n")
  } else {
    cat("❌ No enriched terms for cluster", clust, "\n")
  }
}

if (!length(all_cluster_go)) {
  stop("❌ No GO enrichment results found.")
}

go_all_clusters <- bind_rows(all_cluster_go)

# ===============================================================
# Helper functions
# ===============================================================
safe_map_cluster_labels <- function(clusters_chr, cluster_labels) {
  labs <- cluster_labels[as.character(clusters_chr)]
  labs[is.na(labs)] <- paste0(clusters_chr[is.na(labs)], ": ", clusters_chr[is.na(labs)])
  labs
}

filter_terms_by_cluster_frequency <- function(df, min_clusters = 3) {
  df %>%
    group_by(Description) %>%
    mutate(n_clusters_with_term = n_distinct(cluster)) %>%
    ungroup() %>%
    filter(n_clusters_with_term >= min_clusters)
}

# ===============================================================
# Build combined filtered top_terms object only
# ===============================================================
top_terms <- go_all_clusters %>%
  group_by(cluster, group) %>%
  slice_min(p.adjust, n = 10, with_ties = FALSE) %>%
  ungroup()

top_terms <- filter_terms_by_cluster_frequency(top_terms, min_clusters = min_clusters)

if (nrow(top_terms) == 0) {
  stop("⚠️ No GO terms remain after filtering in combined object.")
}

top_terms$cluster <- as.character(top_terms$cluster)
top_terms$cluster_label <- safe_map_cluster_labels(top_terms$cluster, cluster_labels)
top_terms$cluster_label <- factor(top_terms$cluster_label, levels = unname(cluster_labels))
top_terms$Description <- str_wrap(top_terms$Description, width = 60)

term_order <- top_terms %>%
  group_by(Description) %>%
  summarise(
    n_clusters_with_term = max(n_clusters_with_term),
    best_p = min(p.adjust),
    .groups = "drop"
  ) %>%
  arrange(n_clusters_with_term, -log10(best_p)) %>%
  pull(Description)

top_terms$Description <- factor(top_terms$Description, levels = term_order)

# ---------- Save ONLY RDS ----------
saveRDS(top_terms, rds_out)

cat("✅ RDS saved to:\n", rds_out, "\n")

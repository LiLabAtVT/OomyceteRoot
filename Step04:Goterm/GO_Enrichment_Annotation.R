# ===============================================================
# GO Enrichment per Cluster with Cluster IDs & Labels on X-axis
# ===============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(clusterProfiler)
  library(org.At.tair.db)
  library(ggplot2)
  library(forcats)
  library(stringr)
  library(readr)
})

# ---------- Paths ----------
seurat_path    <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/Razan2025/OutPut_2025/07072025/Patho_Nonpatho_integration_07072025.Rds"
output_dir     <- "/projects/songli_lab/RazanPlantPath_2025/Script_New_Cluster_ID/OutPut"
go_table_path  <- file.path(output_dir, "GO_Enrichment_All_Clusters_07072025.csv")
pdf_infected   <- file.path(output_dir, "GO_Enrichment_Infected_07072025_ClusterID.pdf")
pdf_noninfected<- file.path(output_dir, "GO_Enrichment_NonInfected_07072025_ClusterID.pdf")
pdf_combined   <- file.path(output_dir, "GO_Enrichment_Combined_07072025_ClusterID.pdf")

markers_all_csv <- file.path(output_dir, "Signif_Markers_PerCluster_07072025_ClusterID.csv")
markers_dir     <- file.path(output_dir, "Markers_PerCluster")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(markers_dir)) dir.create(markers_dir, recursive = TRUE)

# ---------- Load Seurat object ----------
seurat_obj <- readRDS(seurat_path)
if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  Idents(seurat_obj) <- seurat_obj$seurat_clusters
}
cluster_ids <- sort(unique(Idents(seurat_obj)))

# ---------- Cluster ID → Label map ----------
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

# ---------- Initialize lists ----------
all_cluster_go <- list()
all_signif_markers <- list()

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

  # Infection label
  cluster_obj$infection_status <- ifelse(grepl("Neg", cluster_obj$orig.ident, ignore.case = TRUE),
                                         "Non_Infected", "Infected")
  Idents(cluster_obj) <- cluster_obj$infection_status

  # --- Differential expression ---
  de_markers <- FindMarkers(cluster_obj,
                            ident.1 = "Non_Infected",
                            ident.2 = "Infected",
                            only.pos = FALSE,
                            min.pct = 0.01,
                            logfc.threshold = 0.25)
  de_markers$gene <- rownames(de_markers)

  signif_markers <- de_markers %>%
    filter(p_val_adj < 0.1) %>%
    mutate(cluster = as.character(clust),
           direction = ifelse(avg_log2FC > 0, "Up_Infected", "Up_NonInfected"))

  all_signif_markers[[as.character(clust)]] <- signif_markers
  readr::write_csv(signif_markers, file.path(markers_dir,
                   paste0("Signif_Markers_cluster_", clust, "_07072025.csv")))

  if (nrow(signif_markers) < 10) {
    cat("⚠️ Few DEGs, skipping GO for cluster", clust, "\n")
    next
  }

  # --- TAIR → ENTREZ mapping ---
  gene_conversion <- suppressWarnings(
    bitr(signif_markers$gene, fromType = "TAIR", toType = "ENTREZID", OrgDb = org.At.tair.db)
  )
  merged_degs <- merge(signif_markers, gene_conversion, by.x = "gene", by.y = "TAIR")

  infected_up     <- merged_degs %>% filter(avg_log2FC > 0) %>% pull(ENTREZID) %>% unique()
  non_infected_up <- merged_degs %>% filter(avg_log2FC < 0) %>% pull(ENTREZID) %>% unique()

  go_inf <- go_noninf <- data.frame()

  if (length(infected_up) > 0) {
    go_infected <- enrichGO(gene = infected_up, OrgDb = org.At.tair.db,
                            keyType = "ENTREZID", ont = "BP",
                            pAdjustMethod = "BH", qvalueCutoff = 0.3, readable = TRUE)
    go_inf <- as.data.frame(go_infected)
    if (nrow(go_inf) > 0) {
      go_inf$group   <- "Infected"
      go_inf$cluster <- as.character(clust)
    }
  }

  if (length(non_infected_up) > 0) {
    go_non_infected <- enrichGO(gene = non_infected_up, OrgDb = org.At.tair.db,
                                keyType = "ENTREZID", ont = "BP",
                                pAdjustMethod = "BH", qvalueCutoff = 0.3, readable = TRUE)
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

# ---------- Save markers & GO ----------
if (length(all_signif_markers)) {
  all_markers <- bind_rows(all_signif_markers)
  write_csv(all_markers, markers_all_csv)
  cat("✅ Combined significant markers saved.\n")
}

if (length(all_cluster_go)) {
  go_all_clusters <- bind_rows(all_cluster_go)
  write_csv(go_all_clusters, go_table_path)
  cat("✅ GO enrichment table saved.\n")
} else {
  stop("❌ No GO enrichment results found.")
}

# ===============================================================
# ---- Helper: Safe map cluster IDs → Labels
# ===============================================================
safe_map_cluster_labels <- function(clusters_chr, cluster_labels) {
  labs <- cluster_labels[as.character(clusters_chr)]
  labs[is.na(labs)] <- paste0(clusters_chr[is.na(labs)], ": ", clusters_chr[is.na(labs)])
  labs
}

# ===============================================================
# ---- PLOT FUNCTIONS (Cluster ID + Label on X-axis) ----
# ===============================================================
plot_go_dotplot_single <- function(data, group_name, output_file, cluster_labels) {
  top_terms <- data %>%
    filter(group == group_name) %>%
    group_by(cluster, group) %>%
    slice_min(p.adjust, n = 10, with_ties = FALSE) %>%
    ungroup()

  top_terms$cluster <- as.character(top_terms$cluster)
  top_terms$cluster_label <- safe_map_cluster_labels(top_terms$cluster, cluster_labels)
  top_terms$cluster_label <- factor(top_terms$cluster_label, levels = unname(cluster_labels))
  top_terms$Description <- str_wrap(top_terms$Description, width = 60)

  p <- ggplot(top_terms, aes(
    x = cluster_label,
    y = fct_reorder(Description, -log10(p.adjust)),
    color = cluster_label,
    size = -log10(p.adjust)
  )) +
    geom_point() +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("GO Enrichment:", group_name),
      x = "Cluster (ID + Label)",
      y = "GO Biological Process",
      size = "-log10(adj. p-value)",
      color = "Cluster"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9)
    )

  ggsave(output_file, plot = p, width = 14, height = 12)
  cat("✅ GO dotplot saved to:", output_file, "\n")
}

plot_go_dotplot_combined <- function(data, output_file, cluster_labels) {
  top_terms <- data %>%
    group_by(cluster, group) %>%
    slice_min(p.adjust, n = 10, with_ties = FALSE) %>%
    ungroup()

  top_terms$cluster <- as.character(top_terms$cluster)
  top_terms$cluster_label <- safe_map_cluster_labels(top_terms$cluster, cluster_labels)
  top_terms$cluster_label <- factor(top_terms$cluster_label, levels = unname(cluster_labels))
  top_terms$Description <- str_wrap(top_terms$Description, width = 60)

  p <- ggplot(top_terms, aes(
    x = cluster_label,
    y = fct_reorder(Description, -log10(p.adjust)),
    color = group,
    shape = group,
    size = -log10(p.adjust)
  )) +
    geom_point(alpha = 0.85) +
    theme_minimal(base_size = 12) +
    labs(
      title = "GO Enrichment: Infected vs Non-Infected (Combined)",
      x = "Cluster (ID + Label)",
      y = "GO Biological Process",
      size = "-log10(adj. p-value)",
      color = "Condition",
      shape = "Condition"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9)
    )

  ggsave(output_file, plot = p, width = 16, height = 12)
  cat("✅ Combined GO dotplot saved to:", output_file, "\n")
}

# ===============================================================
# ---- Generate the three GO dotplots ----
# ===============================================================
go_all_clusters <- read_csv(go_table_path, show_col_types = FALSE)

plot_go_dotplot_single(go_all_clusters, "Infected",     pdf_infected,   cluster_labels)
plot_go_dotplot_single(go_all_clusters, "Non_Infected", pdf_noninfected, cluster_labels)
plot_go_dotplot_combined(go_all_clusters,               pdf_combined,   cluster_labels)

cat("🎉 All GO dotplots generated successfully.\n")

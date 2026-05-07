suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(readr)
})

# ===============================================================
# Paths
# ===============================================================
top_terms_path <- "/projects/intro2gds/Razan_2026/PLant_Pathogen_2026/Step_40_GOTerm/Script/top_terms_combined_filtered.rds"
out_dir        <- "/projects/intro2gds/Razan_2026/PLant_Pathogen_2026/Step_40_GOTerm/Script"

pdf_infected    <- file.path(out_dir, "GO_dotplot_Infected.pdf")
pdf_noninfected <- file.path(out_dir, "GO_dotplot_Non_Infected.pdf")
pdf_combined    <- file.path(out_dir, "GO_dotplot_Combined.pdf")
csv_output      <- file.path(out_dir, "GO_terms_filtered_output.csv")

# minimum number of clusters a GO term must appear in
min_clusters <- 3

# ===============================================================
# Read data
# ===============================================================
top_terms <- readRDS(top_terms_path)

# Make sure required columns exist
required_cols <- c("cluster_label", "Description", "group", "p.adjust")
missing_cols <- setdiff(required_cols, colnames(top_terms))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ===============================================================
# Filter GO terms:
# keep terms present in >= min_clusters unique clusters
# ===============================================================
top_terms_filtered <- top_terms %>%
  group_by(Description) %>%
  mutate(n_clusters = n_distinct(cluster_label)) %>%
  ungroup() %>%
  filter(n_clusters >= min_clusters)

# Save filtered CSV
write_csv(top_terms_filtered, csv_output)
cat("✅ Filtered CSV saved to:", csv_output, "\n")

# ===============================================================
# Plot function for single condition
# ===============================================================
plot_go_dotplot_single <- function(df, condition_name, output_file) {

  df_sub <- df %>%
    filter(group == condition_name)

  if (nrow(df_sub) == 0) {
    cat("⚠ No rows found for condition:", condition_name, "\n")
    return(NULL)
  }

  point_color <- if (condition_name == "Infected") "firebrick3" else "steelblue3"

  p <- ggplot(
    df_sub,
    aes(
      x = cluster_label,
      y = Description,
      size = -log10(p.adjust)
    )
  ) +
    geom_point(color = point_color, alpha = 0.85) +
    scale_size(range = c(4, 12)) +
    theme_minimal() +
    labs(
      title = paste0("GO Enrichment: ", condition_name),
      x = "Cluster",
      y = "GO Biological Process",
      size = expression(-log[10] * "(adj. p-value)")
    ) +
    guides(
      size = guide_legend(
        order = 1,
        override.aes = list(alpha = 1)
      )
    ) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 28
      ),
      axis.text.y = element_text(size = 28),
      axis.title.x = element_text(
        size = 32,
        face = "bold"
      ),
      axis.title.y = element_text(
        size = 32,
        face = "bold"
      ),
      plot.title = element_text(
        size = 34,
        face = "bold",
        hjust = 0.5
      ),
      legend.title = element_text(
        size = 28,
        face = "bold"
      ),
      legend.text = element_text(size = 26),
      panel.grid.minor = element_blank()
    )

  ggsave(output_file, plot = p, width = 24, height = 18)
  cat("✅ Single-condition GO dotplot saved to:", output_file, "\n")
}

# ===============================================================
# Plot function for combined conditions
# ===============================================================
plot_go_dotplot_combined <- function(df, output_file) {

  p <- ggplot(
    df,
    aes(
      x = cluster_label,
      y = Description,
      color = group,
      shape = group,
      size  = -log10(p.adjust)
    )
  ) +
    geom_point(alpha = 0.85) +
    scale_size(range = c(4, 12)) +
    scale_color_manual(
      values = c(
        "Infected" = "firebrick3",
        "Non_Infected" = "steelblue3"
      )
    ) +
    guides(
      color = guide_legend(
        order = 1,
        override.aes = list(size = 8, alpha = 1)
      ),
      shape = guide_legend(
        order = 1,
        override.aes = list(size = 8, alpha = 1)
      ),
      size = guide_legend(
        order = 2,
        override.aes = list(alpha = 1)
      )
    ) +
    theme_minimal(base_size = 12) +
    labs(
      title = "GO Enrichment: Infected vs Non-Infected",
      x = "Cluster",
      y = "GO Biological Process",
      size = expression(-log[10] * "(adj. p-value)"),
      color = "Condition",
      shape = "Condition"
    ) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 28
      ),
      axis.text.y = element_text(size = 28),
      axis.title.x = element_text(
        size = 32,
        face = "bold"
      ),
      axis.title.y = element_text(
        size = 32,
        face = "bold"
      ),
      plot.title = element_text(
        size = 34,
        face = "bold",
        hjust = 0.5
      ),
      legend.title = element_text(
        size = 28,
        face = "bold"
      ),
      legend.text = element_text(size = 26),
      panel.grid.minor = element_blank()
    )

  ggsave(output_file, plot = p, width = 26, height = 20)
  cat("✅ Combined GO dotplot saved to:", output_file, "\n")
}

# ===============================================================
# Generate outputs
# ===============================================================
plot_go_dotplot_single(
  top_terms_filtered,
  "Infected",
  pdf_infected
)

plot_go_dotplot_single(
  top_terms_filtered,
  "Non_Infected",
  pdf_noninfected
)

plot_go_dotplot_combined(
  top_terms_filtered,
  pdf_combined
)

cat("🎉 All GO dotplots and CSV generated successfully.\n")


# 🧬 **Single-nucleus transcriptome analysis of *Arabidopsis thaliana* roots infected with *Phytophthora capsici*** #

## **Introduction** ##

Understanding how plant roots coordinate immune responses at the cellular level is fundamental for dissecting host–pathogen interactions. In this study, we employed single-nucleus RNA sequencing (snRNA-seq) to profile the *Arabidopsis thaliana* root transcriptome following infection with the oomycete pathogen *Phytophthora capsici*.

By capturing transcriptional changes at 24 hours post-infection (24 hpi), we generated a high-resolution cellular atlas that reveals the spatial and functional dynamics of plant immune responses. Using Seurat-based integration, we analyzed over 26,000 nuclei and identified 13 transcriptionally distinct cell clusters representing major root tissues and cell types.

Comparative analyses between infected and mock samples highlighted the activation of key defense and hormone signaling pathways, including salicylic acid (SA) and ethylene (ET) responses, alongside oxidative stress mechanisms. Conversely, mock-treated roots were enriched for metabolic and translational processes.

Gene Ontology (GO) enrichment and marker gene expression profiling further revealed that immune responses to *P. capsici* are cell-type specific and involve distinct sets of defense-related genes. Together, this study provides a comprehensive reference for host transcriptional reprogramming during early stages of oomycete infection and establishes a framework for exploring plant–pathogen interactions at single-cell resolution.

This repository contains all scripts and workflows used for the analysis of single-nucleus RNA-seq (snRNA-seq) data from *Arabidopsis thaliana* roots infected with *Phytophthora capsici* (24 hours post infection).
The pipeline integrates host–pathogen mapping, reference construction, data integration, clustering, differential expression, and Gene Ontology (GO) enrichment.

## **🧩 Requirements** ##

To reproduce this workflow, install the following tools and R packages with the specified versions (or newer):

Tool / Package	Version	Description

<img width="814" height="446" alt="Screenshot 2025-10-30 at 10 38 45 AM" src="https://github.com/user-attachments/assets/2f1ede4e-6123-4ab1-9c8d-bb10390f1555" />


# **🧰 Installation** #
```bash 
# CRAN packages
install.packages(c("ggplot2", "dplyr", "cowplot"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Seurat", "clusterProfiler", "org.At.tair.db", "enrichplot"))
```

**📂 Repository Structure**
```bash
├── 1_ReadMapping/
│   ├── Cellranger_Filtred_mkgtf.sh
│   ├── Cellranger_mkref_Multispecies.sh
│   ├── Cellranger_mkref_ATH.sh
│   ├── CellRanger_Pos24hpi_1_scRNAseq_Multi_Genome.sh
│   └── README.md
│
├── 2_Integration_UMAP/
│   ├── Integration_Seurat_SCT.R
│   ├── UMAP_ClusterID.R
│   └── README.md
│
├── 3_CellMarkers_DEG/
│   ├── CellMarkers_DEG_Analysis.R
│   └── README.md
│
├── 4_GO_Annotation/
│   ├── GO_Enrichment_Annotation.R
│   ├── GO_Enrichment_ImmuneDefense.R
│   ├── GO_Enrichment_HormoneResponse.R
│   └── README.md
│
├── 5_Visualization/
│   ├── Heatmap_Markers.R
│   ├── DotPlot_FunctionalGenes.R
│   ├── FeaturePlot_CellMarkers.R
│   ├── FeaturePlot_SelectedMarkers.R
│   └── README.md
│
├── 6_QC_Summary/
│   ├── RelativeFraction_CellNumber.R
│   ├── UMI_Distribution.R
│   └── README.md
│
└── main_README.md  ← (this file)
``` 
## **🧬 Pipeline Overview** ##

### **Step 1 - Read Mapping & Reference Preparation** ###

Scripts:
```bash
Cellranger_Filtred_mkgtf.sh, Cellranger_mkref_Multispecies.sh, Cellranger_mkref_ATH.sh, CellRanger_Pos24hpi_1_scRNAseq_Multi_Genome.sh
```
Goal: build host-only and multi-species references and align reads with Cell Ranger.
Outputs: multi-species reference (refdata-Arabidopsis_Pcapsici_MultiRef/),
Arabidopsis-only reference (refdata-Arabidopsis_thaliana.TAIR10.59/), BAM files, filtered matrices.

### **Step 2 - Data Integration & UMAP Clustering** ###

Scripts:
```bash
Integration_Seurat_SCT.R, UMAP_ClusterID.R
```
Integrates infected (Pos) and non-infected (Neg) samples using Seurat v5 SCTransform workflow, performs PCA, clustering, and UMAP visualization.
Outputs: Integrated_Seurat_Object.rds, UMAP_Cluster_Labelled.png, Cluster_Metadata.csv

### **Step 3 — Marker Gene & Differential Expression** ###

Script:
```bash
CellMarkers_DEG_Analysis.R
```
Identifies cluster-specific marker genes and DEGs between infected and control samples.
Outputs:
Cluster_Markers_All.csv, Infected_vs_NonInfected_DEG.csv, Volcano_DEG_Pos_vs_Neg.png, Top10_Markers_Heatmap.png

### **Step 4 — GO Enrichment and Functional Annotation** ###

Scripts:
```bash
GO_Enrichment_Annotation.R, GO_Enrichment_ImmuneDefense.R, GO_Enrichment_HormoneResponse.R
```
Performs:

🧠 global GO enrichment (BP / MF / CC)

🦠 immune & defense pathway analysis

🌿 hormone response enrichment (SA / JA / ET)

Outputs:
GO_Enrichment_Results.csv, ImmuneDefense_HormoneResponse_GO.csv, GO_Enrichment_Dotplot.png

### **Step 5 — Visualization** ###

Scripts:
```bash
Heatmap_Markers.R, DotPlot_FunctionalGenes.R, FeaturePlot_CellMarkers.R, FeaturePlot_SelectedMarkers.R
```
Generates visual summaries for marker expression and functional genes.
Outputs:

Script	Visualization
Heatmap_Markers.R	Expression of top markers per cluster
DotPlot_FunctionalGenes.R	Defense/hormone gene expression across clusters
FeaturePlot_CellMarkers.R	UMAP of cluster-defining genes
FeaturePlot_SelectedMarkers.R	UMAPs of selected functional marker genes

 ### **Step 6 — Quality Control and Summary Metrics** ###

Scripts:
```bash
RelativeFraction_CellNumber.R, UMI_Distribution.R
```
Quantifies cluster abundance and evaluates sequencing depth (UMIs).
Outputs:
Cluster_RelativeFraction_Barplot.png, UMI_Distribution_Violin.png, Cluster_CellNumber_RelativeFraction.csv

 ### **Step 7 — Focused Enrichment (Immune & Hormone Pathways)** ###

Scripts:
```bash
GO_Enrichment_ImmuneDefense.R, GO_Enrichment_HormoneResponse.R
```
Highlights GO terms such as:

Defense response to oomycete (GO:0002239)

Response to salicylic acid (GO:0009751)

Response to jasmonic acid (GO:0009753)

Ethylene-activated signaling pathway (GO:0009873)

Outputs:
GO_ImmuneDefense_Dotplot.png, GO_HormoneResponse_Dotplot.png

## 🧠 Interpretation Summary ##

| **Step**                 | **Biological Insight**                                                     |
| ------------------------ | -------------------------------------------------------------------------- |
| Integration & Clustering | Defines transcriptionally distinct root cell populations.                  |
| Marker Analysis          | Identifies genes defining clusters and infection-induced changes.          |
| GO Enrichment            | Highlights biological processes for defense, hormone, and stress pathways. |
| Visualization            | Displays expression landscapes of tissue-specific and immune genes.        |
| QC Metrics               | Ensures balanced representation and sequencing consistency.                |



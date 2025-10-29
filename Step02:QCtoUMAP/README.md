**🧬 Integration and UMAP Visualization Scripts**

This folder contains the R scripts used for integrating single-nucleus RNA-seq datasets and visualizing the integrated data with UMAP embeddings and cluster identities.

⚙️ 1️⃣ Integration Script 

**Purpose**

This script performs data integration across multiple samples (infected vs. non-infected) using the Seurat v5 SCTransform workflow.
It aligns single-nucleus expression profiles to remove batch effects and create a unified object for downstream clustering and visualization.

🧭 2️⃣ UMAP and Cluster ID Script 

**Purpose**

This script performs dimensional reduction and clustering visualization of the integrated Seurat object generated in Step 1.
It identifies transcriptionally distinct clusters and annotates them based on known Arabidopsis cell-type markers.

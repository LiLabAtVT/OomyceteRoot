# Integration of Pathogen and Non-pathogen scRNA-seq data

library(Seurat)
library(sctransform)

# Load raw 10X count matrices
Neg1 <- Read10X(data.dir = "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Cell_Ranger_Analysis_Multiple_Genome_10X_10_22_24/10X_Combined_Pathogenome_Mapping_Neg24hpi_1_scRNAseq_10_22_24/outs/filtered_feature_bc_matrix")
Neg2 <- Read10X(data.dir = "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Cell_Ranger_Analysis_Multiple_Genome_10X_10_22_24/10x_Combined_Pathogene_Mapping_Neg24hpi_2_scRNAseq_10_22_24/outs/filtered_feature_bc_matrix")
Pos1 <- Read10X(data.dir = "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/10X_Pos24hpi_Bam_file_Barcode_File/Pos24hpi_1/11_1_24_10X/plant_fastq_output_20241101_170110/10x_Combined_Pathogene_Mapping_Pos24hpi_1_scRNAseq_10_28_24_0_1_22GH3VLT4/combined_Files_Fastq/Mapping_Pos24hpi_1_scRNAseq_11_2_24_ATH_Read_Only/outs/filtered_feature_bc_matrix")
Pos2 <- Read10X(data.dir = "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/10X_Pos24hpi_Bam_file_Barcode_File/Pos24hpi_2/11_1_24_10X/plant_fastq_output_20241101_215749/10X_Combined_Pathogene_Mapping_Pos24hpi_2_scRNAseq_10_28_24_0_1_22GH3VLT4/Combined_Files_Fastq/Mapping_Pos24hpi_2_scRNAseq_11_2_24_ATH_Read_Only/outs/filtered_feature_bc_matrix")

rownames(Neg1) <- gsub("Arabidopsis_thaliana_gene:", "", rownames(Neg1))
rownames(Neg2) <- gsub("Arabidopsis_thaliana_gene:", "", rownames(Neg2))
rownames(Pos1) <- gsub("gene:", "", rownames(Pos1))
rownames(Pos2) <- gsub("gene:", "", rownames(Pos2))

process_seurat_object <- function(counts_data, project_name,
                                  feature_upper, nfeatures,
                                  pca_dims, clustering_resolution) {
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts_data, project = project_name, 
                                   min.cells = 3, min.features = 200)
  
  # Compute mitochondrial percentage for Arabidopsis
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^ATMG")

seurat_obj[["percent.mt"]] <- 0

  # Filter based on quality control metrics
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 100 & 
                         nFeature_RNA < feature_upper & 
                         percent.mt < 5)
                         
    
  # Normalize data
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  # Identify variable features
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = nfeatures)
  # Scale and run PCA
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  # Cluster and visualize
  seurat_obj <- FindNeighbors(seurat_obj, dims = pca_dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = clustering_resolution)
  seurat_obj <- RunUMAP(seurat_obj, dims = pca_dims)
  
  # Save object
  saveRDS(seurat_obj, paste0("/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Script/OMG_Nina_Analysis/OutPut/", project_name, ".rds"))
  
  return(seurat_obj)
}


Neg24hpi_1 <- process_seurat_object(Neg1, project_name = "Neg_hpi1", 
                                    feature_upper = 6000, nfeatures = 2000,
                                    pca_dims = 1:20, clustering_resolution = 0.22)
Neg24hpi_2 <- process_seurat_object(Neg2, project_name = "Neg_hpi2", 
                                    feature_upper = 5000, nfeatures = 2000,
                                    pca_dims = 1:20, clustering_resolution = 0.22)
Pos24hpi_1 <- process_seurat_object(Pos1, project_name = "Pos_hpi1", 
                                    feature_upper = 3500, nfeatures = 2000,
                                    pca_dims = 1:20, clustering_resolution = 0.22)
Pos24hpi_2 <- process_seurat_object(Pos2, project_name = "Pos_hpi2", 
                                    feature_upper = 2000, nfeatures = 2000,
                                    pca_dims = 1:20, clustering_resolution = 0.22)


crossSamples <- SplitObject(merge(Neg24hpi_1, y= c(Neg24hpi_2, Pos24hpi_1, Pos24hpi_2)) , split.by = "orig.ident")
for (i in names(crossSamples)) {
  crossSamples[[i]] <- SCTransform(crossSamples[[i]], verbose = FALSE)
}
features <- SelectIntegrationFeatures(object.list = crossSamples, nfeatures = 10000)
print(length(features))
crossSamples <- PrepSCTIntegration(object.list = crossSamples, anchor.features = features)
crossSamples.anchors <- FindIntegrationAnchors(object.list = crossSamples, normalization.method = "SCT", anchor.features = features)
crossSamples.combine <- IntegrateData(anchorset = crossSamples.anchors, normalization.method = "SCT")
crossSamples.combine <- RunPCA(crossSamples.combine, verbose = FALSE)
crossSamples.combine <- RunUMAP(crossSamples.combine, reduction = "pca", dims = 1:30)
crossSamples.combine <- FindNeighbors(crossSamples.combine, reduction = "pca", dims = 1:30)
crossSamples.combine <- FindClusters(crossSamples.combine, resolution = 0.3) #org =0.3

saveRDS(crossSamples.combine, "./Patho_Nonpatho_integration_07072025.Rds")



# Save metadata in current directory
write.csv(crossSamples.combine@meta.data, 
          file = "./Meta_Patho_Nonpatho_Integrated_07072025.csv", 
          row.names = TRUE)

# Find and save marker genes in current directory
markers_combined <- FindAllMarkers(crossSamples.combine, 
                                   only.pos = TRUE, 
                                   min.pct = 0.1, 
                                   logfc.threshold = 0.25)

write.csv(markers_combined, 
          file = "./Markers_Patho_Nonpatho_Integrated_07072025.csv", 
          row.names = FALSE)

pdf("./Patho_Nonpatho_07072025.pdf", width = 8, height = 5)
DimPlot(crossSamples.combine, reduction = "umap", group.by = "orig.ident", cols = c("lightblue","blue", "#FFC073","#FF8C00"))
dev.off()

pdf("./Patho_Nonpatho_sepColor_07072025.pdf", width = 8, height = 5)
DimPlot(crossSamples.combine, reduction = "umap", split.by = "orig.ident", group.by = "orig.ident", cols = c("lightblue","blue", "#FFC073","#FF8C00"))
dev.off()

pdf("./Patho_Nonpatho_lable_07072025.pdf", width = 8, height = 5)
DimPlot(crossSamples.combine, reduction = "umap", label = TRUE)
dev.off()


pdf("./Patho_Nonpatho_sep_07072025.pdf", width = 8, height = 5)
DimPlot(crossSamples.combine, reduction = "umap", split.by = "orig.ident")
dev.off()

#!/bin/bash
#SBATCH --job-name=Cellranger_Analysis_Pos24hpi_1  # Job name
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1                      # Number of tasks (MPI processes)
#SBATCH --cpus-per-task=64               # Number of CPU cores per task
#SBATCH --mem=128G                       # Memory per node (in GB)
#SBATCH --time=24:00:00                 # Time limit (hh:mm:ss)
#SBATCH --mail-type=ALL                      # Send email at beginning and end of job
#SBATCH --mail-user=arazan@vt.edu       # Email address for job notifications
#SBATCH --output=Cellranger_Pos24hpi_1_ATH_Read_Only.out      # Standard output log
#SBATCH --account=introtogds            # Replace with your valid account

echo "Starting"

# Run Cell Ranger count command
/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Refdata_Arab/cellranger-8.0.1/cellranger count \
    --id=Mapping_Pos24hpi_1_scRNAseq_11_2_24_ATH_Read_Only \
    --transcriptome=/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Refdata_Arab/ATH_genome_10_3_24 \
    --fastqs=/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/10X_Pos24hpi_Bam_file_Barcode_File/Pos24hpi_1/11_1_24_10X/plant_fastq_output_20241101_170110/10x_Combined_Pathogene_Mapping_Pos24hpi_1_scRNAseq_10_28_24_0_1_22GH3VLT4/combined_Files_Fastq \
    --sample=Pos24hpi1 \
    --force-cells=8000 \
    --create-bam=true

echo "Finished"
date

exit;


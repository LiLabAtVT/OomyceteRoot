#!/bin/bash
#SBATCH --job-name=Cellranger_Analysis_Neg24hpi_2  # Job name
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1                      # Number of tasks (MPI processes)
#SBATCH --cpus-per-task=64               # Number of CPU cores per task
#SBATCH --mem=128G                       # Memory per node (in GB)
#SBATCH --time=24:00:00                 # Time limit (hh:mm:ss)
#SBATCH --mail-type=ALL                      # Send email at beginning and end of job
#SBATCH --mail-user=arazan@vt.edu       # Email address for job notifications
#SBATCH --output=Cellranger_Neg24hpi_2_10_28_24.out      # Standard output log
#SBATCH --account=introtogds            # Replace with your valid account

echo "Starting"

# Run Cell Ranger count command
/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Refdata_Arab/cellranger-8.0.1/cellranger count \
    --id=10x_Combined_Pathogene_Mapping_Neg24hpi_2_scRNAseq_10_28_24 \
    --transcriptome=/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Refdata_Arab/ATH_genome_10_3_24 \ 
    --fastqs=/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/scRNA_Seq_Arab_Data/usftp21.novogene.com/01.RawData/Neg24hpi_2 \
    --sample=Neg24hpi_2_CKDL240032443-1A_22GH3VLT4 \
    --force-cells=8000 \
    --create-bam=true


echo "Finished"
date

exit;


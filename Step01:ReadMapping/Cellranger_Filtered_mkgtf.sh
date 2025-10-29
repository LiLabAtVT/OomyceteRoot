#!/bin/bash
#SBATCH --job-name=cellranger_filtered_mkgtf           # Job name
#SBATCH --nodes=1                             # Number of nodes
#SBATCH --ntasks=1                            # Number of tasks (MPI processes)
#SBATCH --cpus-per-task=4                     # Number of CPU cores per task
#SBATCH --mem=16G                             # Memory per node (in GB)
#SBATCH --time=01:00:00                       # Time limit (hh:mm:ss)
#SBATCH --mail-user=arazan@vt.edu             # Email address for job notifications
#SBATCH --mail-type=ALL                       # Send email at beginning and end of job
#SBATCH --output=cellranger_filtered_mkgtf_%j.out      # Standard output log
#SBATCH --account=introtogds                  # Replace with your valid account

echo "Starting cellranger mkgtf"
date

# Full path to the cellranger executable
CELLRANGER_PATH="/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Refdata_Arab/cellranger-8.0.1"  # Replace with the actual path

# Set paths for GTF files
GTF_FILE="Arabidopsis_thaliana.TAIR10.59.gtf"
FILTERED_GTF_FILE="Arabidopsis_thaliana.TAIR10.59.gtf.filtered.gtf"

# Run cellranger mkgtf
$CELLRANGER_PATH/cellranger mkgtf ${GTF_FILE} ${FILTERED_GTF_FILE} --attribute=gene_biotype:protein_coding

if [ $? -ne 0 ]; then
    echo "Error: cellranger mkgtf command failed."
    exit 1
fi

echo "cellranger mkgtf completed successfully"
date

exit 0

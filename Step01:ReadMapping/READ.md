


🌍 **Build Multi-Species Reference Genome**

To analyze host–pathogen single-nucleus RNA-seq data, a combined reference genome must be created so that reads can be mapped to both *Arabidopsis thaliana* (host) and *Phytophthora capsici* (pathogen).

**Prepare the genome and annotation files**
🧩 **Step 1: Download Pathogen Genome**

Obtain *the Phytophthora capsici* genome from NCBI Datasets.

Download both of the following files:

FASTA file — contains the genomic sequence.

GTF file — contains gene annotations.

Example genome used in this project:

GCA_016618375.1 – Phytophthora capsici genome

Save both files in the designated reference directory before building the combined or species-specific reference.

🌿 **Step 2: Download Host Genome**

Download the host plant genome from Ensembl Plants.

Obtain both of the following files:

FASTA file — genomic DNA sequence.

GTF file — gene annotation file.

Example genome used in this project:

Arabidopsis thaliana – Ensembl Plants Release 59

Store both files in your reference directory to prepare for single- or multi-species reference building.


**2.Filter the GTF file**

Before building the reference, filter the GTF annotation file to include only protein-coding genes, following 10x Genomics recommendations:
```bash
Cellranger_Filtred_mkgtf.sh
```

This produces a streamlined annotation file that improves mapping accuracy and reduces reference size.

**3. Build the multi-species reference**

Combine both genomes and their filtered annotation files into one reference using cellranger mkref:

```bash
/Cellranger_mkref_Multispecies.sh
```
***4. Output***

A new directory will be created:

refdata-Arabidopsis_Pcapsici_MultiRef/
This folder contains the combined indices and annotations required for read alignment and UMI counting.

**🧫 Step 5: Run Cell Ranger Count for Infected Samples using Multi-Species Reference Genome**

This step performs read alignment and UMI quantification for each single-nucleus RNA-seq sample using the multi-species reference genome created in Step 4.
The same reference allows simultaneous mapping to both Arabidopsis thaliana (host) and Phytophthora capsici (pathogen).

**1. Sample naming convention**

Neg → Non-infected (control) samples

Pos → Infected samples (24 h post infection)

Example samples in this project:

Neg24hpi_1, Neg24hpi_2 → non-infected roots

Pos24hpi_1, Pos24hpi_2 → infected roots

**2. Run Cell Ranger Count**

Use the following SLURM script for the infected sample sample:
CellRanger_Pos24hpi_1_scRNAseq_Multi_Genome.sh
CellRanger_Pos24hpi_2_scRNAseq_Multi_Genome.sh

This script executes the cellranger count command, specifying the appropriate sample FASTQ path and the multi-reference genome generated in Step 4.



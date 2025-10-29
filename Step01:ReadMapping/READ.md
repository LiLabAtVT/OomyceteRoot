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

**Output**

A new directory will be created:

refdata-Arabidopsis_Pcapsici_MultiRef/

This step creates a Cell Ranger–compatible Host-Pathogen reference directory, which will be used in downstream quantification (Step 5).

**4. Build the *Arabidopsis* Reference (ATH)**

Use the *Arabidopsis thaliana* reference genome files — the FASTA (sequence) and GTF (annotation) — to build a single-species reference genome.
This reference will be used later for:

The filtered FASTQ files generated from the infected samples, and

The original FASTQ files from the non-infected (Neg) samples.

Run the following script to build the reference:

Cellranger_mkref_ATH.sh

This step creates a Cell Ranger–compatible Arabidopsis-only reference directory, which will be used in downstream quantification (Step 8).

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

**🧬 Step 6: Separate Host vs Pathogen Barcodes (Infected samples: Pos24hpi_1, Pos24hpi_2)**

After cellranger count (Step 5), each infected sample directory contains:

outs/filtered_feature_bc_matrix/ — matrix for downstream analysis,

outs/*.bam — alignments (used later),

outs/gem_classification_<SAMPLE>.csv — per-barcode species call (Arabidopsis thaliana, Pcap, or Multiplet). 
10x Genomics
+2
10x Genomics
+2

Do this for both infected samples to create organism-specific barcode lists.

Commands (run from the folder containing the gem_classification_*.csv files)

```bash
# list the infected (Pos) samples you want to process
SAMPLES=("Pos24hpi_1" "Pos24hpi_2")

for S in "${SAMPLES[@]}"; do
  CSV="gem_classification_${S}.csv"

  # Host: Arabidopsis thaliana
  awk -F',' 'NR > 1 && $4 == "Arabidopsis_thaliana" { print $1 }' \
    "${CSV}" > "${S}_arabidopsis_barcodes.txt"

  # Pathogen: Phytophthora capsici
  awk -F',' 'NR > 1 && $4 == "Pcap" { print $1 }' \
    "${CSV}" > "${S}_pcap_barcodes.txt"

  # Multiplets: mixed droplets
  awk -F',' 'NR > 1 && $4 == "Multiplet" { print $1 }' \
    "${CSV}" > "${S}_multiplet_barcodes.txt"

  echo "Wrote ${S}_arabidopsis_barcodes.txt, ${S}_pcap_barcodes.txt, ${S}_multiplet_barcodes.txt"
done
```

What this does

-F',' treats the CSV as comma-separated,

NR > 1 skips the header,

$4 selects the call column to filter Arabidopsis_thaliana, Pcap, or Multiplet,

$1 prints the barcode (GEM) ID.

Outputs created per sample
```bash
Pos24hpi_1_arabidopsis_barcodes.txt
Pos24hpi_1_pcap_barcodes.txt
Pos24hpi_1_multiplet_barcodes.txt
Pos24hpi_2_arabidopsis_barcodes.txt
Pos24hpi_2_pcap_barcodes.txt
Pos24hpi_2_multiplet_barcodes.txt
```

**Notes**
GEM classification and the gem_classification.csv file are standard outputs for mixed-species (“barnyard”) analyses in Cell Ranger. 

**🧪 Step 7: Generate Organism-Specific BAM and FASTQ Files**

Goal
Use the barcodes extracted in Step 6 to filter the original Cell Ranger BAM files and create separate BAM / FASTQ files for host (Arabidopsis thaliana) and pathogen (Phytophthora capsici) reads.
This enables downstream analysis of each organism individually.

**1. Install Required Tools**

Download the official 10x Genomics utilities (Linux binaries):

bamtofastq – https://github.com/10XGenomics/bamtofastq/releases

subset-bam – https://github.com/10XGenomics/subset-bam/releases

Add both executables to your $PATH.

**2. Filter BAM Files by Barcode**

Run the following script for each infected sample (Pos24hpi_1, Pos24hpi_2):

```bash
/Filter_and_convert_10X.sh
```

💡 subset-bam retains only the reads whose cell barcodes match the ones listed in each file produced in Step 6.

💡 Each conversion produces a set of *_R1.fastq.gz, *_R2.fastq.gz, and *_I1.fastq.gz files for downstream processing.

**3. Outputs**

Per infected sample you will obtain:

Pos24hpi_1_Arabidopsis_filtered.bam
Pos24hpi_1_Pcap_filtered.bam
Pos24hpi_2_Arabidopsis_filtered.bam
Pos24hpi_2_Pcap_filtered.bam
Pos24hpi_1_Arabidopsis_fastq/
Pos24hpi_1_Pcap_fastq/
Pos24hpi_2_Arabidopsis_fastq/
Pos24hpi_2_Pcap_fastq/

**🌱 Step 8: Quantify Using Arabidopsis-Only Reference**

We now quantify reads only against the Arabidopsis thaliana reference (built in Step 3).

Infected (Pos) samples: use the filtered FASTQ directories produced in Step 7 (host-only reads).

Non-infected (Neg) samples: use the original FASTQ directories from sequencing (all reads are host).

```bash 
/CellRanger_Neg24hi_1_scRNAseq_ATH_Genome.sh
/CellRanger_Neg24hi_2_scRNAseq_ATH_Genome.sh
/CellRanger_Pos24hi_1_scRNAseq_ATH_Genome.sh
/CellRanger_Pos24hi_2_scRNAseq_ATH_Genome.sh
```
**📁 Step 9: Prepare Filtered Matrices for Downstream Analysis**

After running Cell Ranger count (Step 8), each sample output directory contains a folder:

outs/filtered_feature_bc_matrix/


This folder includes the three core files required for downstream analysis in Seurat or other single-cell analysis packages:

File	Description
barcodes.tsv.gz	      List of cell barcodes corresponding to high-quality nuclei.
features.tsv.gz	      Feature (gene) annotation table linking gene IDs and names.
matrix.mtx.gz	        Sparse expression matrix storing counts for each gene × barcode pair.

**Usage**
These three files are used as input for loading count matrices into R


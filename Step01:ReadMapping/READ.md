**🌍 ReadMapping**

To analyze host–pathogen single-nucleus RNA-seq data, two reference genomes are required:

A multi-species reference for simultaneous mapping of host (Arabidopsis thaliana) and pathogen (Phytophthora capsici) reads.

An Arabidopsis-only reference (ATH) for host-specific quantification and downstream analysis.

**🧩 Step 1: Download Pathogen Genome**

Obtain the Phytophthora capsici genome from NCBI Datasets
.

Download both of the following files:

FASTA file — contains the genomic sequence.

GTF file — contains gene annotations.

Example genome used in this project:

GCA_016618375.1 – Phytophthora capsici genome


Save both files in the designated reference directory before building the combined or species-specific reference.

**🌿 Step 2: Download Host Genome**

Download the Arabidopsis thaliana genome from Ensembl Plants
.

Obtain both of the following files:

FASTA file — genomic DNA sequence.

GTF file — gene annotation file.

Example genome used in this project:

Arabidopsis thaliana – Ensembl Plants Release 59


Store both files in your reference directory to prepare for single- or multi-species reference building.

**🧬 Step 3: Filter the GTF File**

Before building the reference, filter the Arabidopsis GTF annotation file to include only protein-coding genes, following 10x Genomics recommendations.

Run the script:
```bash
Cellranger_Filtred_mkgtf.sh
```

This step produces a simplified annotation file that improves mapping accuracy and reduces reference size.

**🧫 Step 4: Build the Multi-Species Reference Genome**

Combine both genomes and their filtered annotation files into one reference using Cell Ranger’s mkref function.

Run the script:
```bash
Cellranger_mkref_Multispecies.sh
```

Output:

refdata-Arabidopsis_Pcapsici_MultiRef/


This step creates a Cell Ranger–compatible Host–Pathogen reference directory, which will be used for infected samples during initial quantification (Step 5).

**🌱 Step 5: Build the Arabidopsis-Only Reference (ATH)**

Use the Arabidopsis thaliana reference genome files — the FASTA (sequence) and GTF (annotation) — to build a single-species reference genome.

This reference will be used later for:

The filtered FASTQ files generated from infected samples (host-only reads).

The original FASTQ files from non-infected (Neg) samples.

Run the script:
```bash
Cellranger_mkref_ATH.sh
```

This creates a Cell Ranger–compatible Arabidopsis-only reference directory, used for downstream quantification (Step 8).

**🧪 Step 6: Run Cell Ranger Count (Infected Samples, Multi-Species Reference)**

Perform read alignment and UMI quantification using the multi-species reference built in Step 4.

Sample Naming Convention

**Neg** → Non-infected (control) samples

**Pos** → Infected samples (24 h post infection)

Example samples:

Neg24hpi_1, Neg24hpi_2 → non-infected roots  
Pos24hpi_1, Pos24hpi_2 → infected roots

Run Scripts
```bash
CellRanger_Pos24hpi_1_scRNAseq_Multi_Genome.sh
CellRanger_Pos24hpi_2_scRNAseq_Multi_Genome.sh
```

Each script executes cellranger count, specifying the sample FASTQ path and the multi-species reference.

**🧬 Step 7: Separate Host vs Pathogen Barcodes (Infected Samples)**

After running Cell Ranger (Step 6), each infected sample directory contains:

outs/filtered_feature_bc_matrix/ — filtered matrix

outs/*.bam — aligned reads

outs/gem_classification_<SAMPLE>.csv — per-barcode species calls (Arabidopsis thaliana, Pcap, or Multiplet)

Run the following commands to extract barcodes by species:

```bash
SAMPLES=("Pos24hpi_1" "Pos24hpi_2")

for S in "${SAMPLES[@]}"; do
  CSV="gem_classification_${S}.csv"

  # Host: Arabidopsis thaliana
  awk -F',' 'NR > 1 && $4 == "Arabidopsis_thaliana" { print $1 }' \
    "${CSV}" > "${S}_arabidopsis_barcodes.txt"

  # Pathogen: Phytophthora capsici
  awk -F',' 'NR > 1 && $4 == "Pcap" { print $1 }' \
    "${CSV}" > "${S}_pcap_barcodes.txt"

  # Multiplets
  awk -F',' 'NR > 1 && $4 == "Multiplet" { print $1 }' \
    "${CSV}" > "${S}_multiplet_barcodes.txt"
done
```

Output files per sample:

Pos24hpi_1_arabidopsis_barcodes.txt
Pos24hpi_1_pcap_barcodes.txt
Pos24hpi_1_multiplet_barcodes.txt
Pos24hpi_2_arabidopsis_barcodes.txt
Pos24hpi_2_pcap_barcodes.txt
Pos24hpi_2_multiplet_barcodes.txt

**🧪 Step 8: Generate Organism-Specific BAM and FASTQ Files**

Use the barcode lists (Step 7) to filter the original BAMs and generate organism-specific FASTQ files.

Run:
```bash
Filter_and_convert_10X.sh
```

Tools required:

subset-bam

bamtofastq

Output directories:

Pos24hpi_1_Arabidopsis_fastq/
Pos24hpi_1_Pcap_fastq/
Pos24hpi_2_Arabidopsis_fastq/
Pos24hpi_2_Pcap_fastq/

**🌾 Step 9: Quantify Using Arabidopsis-Only Reference**

Quantify host reads only against the Arabidopsis reference created in Step 5.

Infected samples (Pos): use filtered FASTQs (*_Arabidopsis_fastq/)

Non-infected samples (Neg): use original FASTQs from sequencing

Run the scripts:
```bash
CellRanger_Neg24hpi_1_scRNAseq_ATH_Genome.sh
CellRanger_Neg24hpi_2_scRNAseq_ATH_Genome.sh
CellRanger_Pos24hpi_1_scRNAseq_ATH_Genome.sh
CellRanger_Pos24hpi_2_scRNAseq_ATH_Genome.sh
```

Each script generates filtered_feature_bc_matrix/ for downstream Seurat integration.

**📁 Step 10: Prepare Filtered Matrices for Downstream Analysis**

Each outs/filtered_feature_bc_matrix/ directory includes the three files required for downstream analysis:

File	                  Description
barcodes.tsv.gz	        List of high-quality cell barcodes.
features.tsv.gz	        Gene annotation table linking IDs and gene names.
matrix.mtx.gz	          Sparse count matrix of gene × barcode expression values.

**Usage** These three files are used as input for loading count matrices into R

# Comparison of LNCaP Gene Expression Data from PMID: 36167836 Using RNA-Seq Tools 

# RNA-seq Analysis Workflow

This project provides a reproducible RNA-seq analysis workflow using **Nextflow** and **R**. It covers all steps from reference genome preparation to differential expression analysis and functional enrichment visualization.

The workflow was developed for an exercise project. The aim is to compare various treatments of LNCaP cells and to perform functional pathway analysis for the clustered genes.


## Table of Contents

- [Background](#background)
- [Overview](#overview)
- [Experimental Design](#experimental-design)
- [Data and Reference](#data-and-reference)
- [Project Structure](#project-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Workflow](#workflow)
- [Design Choices](#design-choices)
- [Usage](#usage)
- [Outputs](#outputs)
- [Resource Management](#resource-management)

## Background

- Phytogenics as a plant-based feed additive in animal nutrition to support animal growth and health
- Alternative to antibiotic growth promoters due to antibiotic resistance concern
- The paper tackles RNA-Sequencing of liver tissue and ileum to investigate the transcriptome of piglets fed an essential oil-based phytogenic for 28 days

## Overview

The workflow performs the following tasks:

1. **Reference Genome Setup**: Downloads and organizes genome FASTA and GTF annotation files.
2. **Read Alignment & Quantification**: Maps RNA-seq reads to the reference genome and counts reads per gene.
3. **Differential Expression Analysis**: Identifies genes with significant expression changes across conditions.
4. **Functional Enrichment**: Generates barplots for pathways or GO terms enriched in DE genes.
5. **Visualization**: Produces publication-ready figures for downstream analysis.

The original paper investigated the effect of oxidized oil and phytogenic feed additives in pigs. In this project, the focus is on the RNA-seq part of the study and on checking how far the published analysis can be reproduced with a Nextflow workflow.

## Experimental Design

The dataset contains paired-end RNA-seq samples from pig liver (**Sus scrofa**). The samples are divided into four treatment groups:

| Group | Description | Number of samples |
| :--- | :--- | ---: |
| `CTR` | Control diet | 6 |
| `PFA` | Phytogenic feed additive in standard diet | 6 |
| `OO` | Standard diet with oxidized oil | 6 |
| `PFA_OO` | Phytogenic feed additive and oxidized oil | 6 |

In total, 24 samples are used. The sample information is stored in `metadata/SraRunTable.csv`. The treatment information is parsed from the SRA metadata and used as condition information for the downstream analysis.



delete after README file finished

## Data and Reference

The sequencing data comes from the NCBI SRA dataset used in the pig liver transcriptomics study. The FASTQ files are not included in this repository because they are too large. They should be placed in the `data/` folder or linked there as symbolic links.

Expected input file pattern:

```text
data/*_{1,2}.fastq.gz
```

The workflow uses the **Sus scrofa Sscrofa11.1** reference from Ensembl release 111:

- Genome FASTA: `Sus_scrofa.Sscrofa11.1.dna.toplevel.fa`
- Annotation GTF: `Sus_scrofa.Sscrofa11.1.111.gtf`

Both reference files are downloaded automatically by the workflow and are then used to build the STAR index.


## Project Structure

```text
workflow_design_group4/
├── main.nf               # Nextflow pipeline orchestration
├── nextflow.config       # Pipeline parameters, experimental design, and containers
├── modules/              # DSL2 process modules
├── metadata/             # Primary Project Logic & Provenance
│   ├── SraRunTable.csv   # The 24 selected NCBI samples
│   └── README.md         # Internal documentation of metadata export (see below)
├── data/                 # Symbolic links to .fastq.gz files (Ignored by Git)
├── results/              # Output directory (ignored by Git)
├── work/                 # Nextflow work directory (ignored by Git)
└── .gitignore            # Git exclusion rules
```

## Requirements

- [Nextflow](https://www.nextflow.io/) (≥21.04)
- [R](https://www.r-project.org/) (≥4.0) with required packages (`DESeq2`, `clusterProfiler`, etc.)
- Docker, Singularity, or Apptainer for containerized processes
- Sufficient storage for raw reads and reference genomes


## Aim of Analysis

## Questions from Slides
Where is the metadata?: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159583 \
Snakemake pipeline (Where is it ?): The pipeline was not provided. \
What’s the issue in term of reproducibility?: As the pipeline was not available, recreating the entire analysis with the knowledge of the parameters (e.g. STAR alignment) will not be possible entirely. \

## Design Choices

### No Separate Trimming Step

The workflow does not include a separate read trimming step. The raw reads are checked with FastQC first, and the analysis keeps the reads close to the original downloaded SRA data. This also avoids adding another preprocessing step that could make the reproduction of the published analysis harder to compare.

For this dataset, the main limiting factor was not adapter trimming, but handling the large paired-end RNA-seq files and keeping enough disk space during alignment. If a FastQC report showed strong adapter contamination or very low-quality read ends, a trimming module could be added before STAR alignment.

### Counting with STAR GeneCounts

The first version of the workflow used STAR for alignment and HTSeq for counting. In this version, the separate HTSeq step was removed and STAR is run with `--quantMode GeneCounts`.

This was done for two reasons. First, STAR can create gene-level count files directly from the same alignment step, using the same GTF annotation. Second, it avoids writing and storing large sorted BAM files for every sample. This is important on the server because the RNA-seq files are large and the available home directory quota is limited.

The strandedness is controlled in `nextflow.config`:

```groovy
counts = [
    strandedness : "reverse"
]
```

### Enrichment Gene ID Mapping

The functional enrichment step uses `clusterProfiler` and `org.Ss.eg.db`, which require Entrez gene IDs. The clustered genes from the DESeq2/DEGreport analysis are Ensembl gene IDs, so a mapping step is needed before running GO enrichment.

BioMart is a common way to convert Ensembl IDs to Entrez IDs and was used in the first enrichment approach. During testing, the BioMart connection and cache handling were not always stable on the server. Therefore, the current workflow downloads the NCBI `gene2ensembl` table and filters it for pig (`taxid 9823`). This gives a simple Ensembl-to-Entrez mapping that is easier to reproduce in the workflow run.

The enrichment itself is still performed with `clusterProfiler` and `org.Ss.eg.db`. The NCBI table is only used for the ID conversion step before the GO enrichment.

## Usage

Run the workflow with:

```bash
nextflow run main.nf
```

For continuing an interrupted run, use:

```bash
nextflow run main.nf -resume
```

The most important parameters can be changed in `nextflow.config`. For example, the enrichment cluster is set here:

## Outputs

```text
results/
├── reference/          # Genomic FASTA, GTF, and STAR Index
├── fastqc/             # HTML & Zip reports for raw read quality
├── counts/             # Gene-level counts (raw_counts_matrix.txt)
├── stats/              # Unified alignment percentages CSV
├── plots_test/         # Recreated Figure 2A (Clustering boxplots)
└── enrichment/         # Recreated Figure 2B (Functional enrichment barplots)
```

Important output files:

- `results/counts/raw_counts_matrix.txt`
- `results/stats/alignment_stats.csv`
- `results/plots_test/Figure_2A_Recreated.pdf`
- `results/enrichment/Figure_2B_Barplot.pdf`

The downstream part of the analysis produces the recreated clustering plot for Figure 2A and a GO enrichment barplot for Figure 2B. In the tested run, the enrichment results include terms related to lipid oxidation, fatty acid oxidation, organic acid metabolism, and mitochondrial structures.

## Resource Management

To handle large RNA-seq datasets (60M+ reads), this pipeline is optimized to minimize disk usage in the user's home directory and bypass strict storage quotas (e.g., 100GB limits).

### Scratch Space (SSD)
By default, the heavy I/O step (**STAR_ALIGN**) is configured to use "scratch" space. STAR writes compact gene-count tables directly, avoiding large sorted BAM files and reducing pressure on strict storage quotas.

* **Pros:** Much lower disk usage because the workflow no longer keeps per-sample alignment BAMs.
* **Cons:** Requires enough temporary space for STAR's own alignment files while each sample is running.

### How to Toggle
You can control this behavior using the `--use_scratch` parameter without modifying the code:

| Command | Description |
| :--- | :--- |
| `nextflow run main.nf` | **Default:** Uses scratch space for STAR. |
| `nextflow run main.nf --use_scratch false` | Writes all temporary files to the `work/` directory (HDD). |

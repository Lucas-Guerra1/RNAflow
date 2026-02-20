# RNAflow — Automated RNA-seq Analysis Pipeline

A step-by-step Bash pipeline for RNA-seq data analysis on Linux, integrating widely used bioinformatics tools. Users provide FASTQ files and a reference genome; the pipeline handles quality control, genome indexing, read alignment, and BAM indexing.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Directory Structure](#directory-structure)
- [Usage](#usage)
- [Output Files](#output-files)
- [Visualizing Results in IGV](#visualizing-results-in-igv)
- [Limitations](#limitations)
- [License](#license)

---

## Overview

RNAflow automates the core steps of RNA-seq analysis:

```
FASTQ files → QC (FastQC) → Aggregated QC (MultiQC) → Genome Index (STAR) → Alignment (STAR) → Mapping QC (MultiQC) → BAM Indexing (Samtools) → IGV
```

---

## Features

- Quality control of raw reads with FastQC
- Aggregated QC reports with MultiQC
- Genome index generation with STAR
- Paired-end read alignment with STAR
- BAM file indexing with Samtools for downstream analysis and IGV visualization
- Clear directory structure to organize all inputs and outputs

---

## Requirements

- Linux-based operating system (Ubuntu, CentOS, or equivalent)
- [Conda](https://docs.conda.io/en/latest/miniconda.html) (Miniconda or Anaconda)
- The following tools (see [Installation](#installation)):
  - STAR ≥ 2.7
  - FastQC
  - MultiQC
  - Samtools

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/your-username/RNAflow.git
cd RNAflow
```

### 2. Activate Conda

```bash
eval "$(~/anaconda3/bin/conda shell.bash hook)"
```

> **Note:** Adjust the path above to match your Conda installation (e.g., `~/miniconda3/bin/conda`).

### 3. Install dependencies via Conda

```bash
conda install -c bioconda star fastqc multiqc samtools
```

#### Alternatively, install STAR manually

```bash
wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.10a.tar.gz
tar -xzf 2.7.10a.tar.gz
cd STAR-2.7.10a/bin/Linux_x86_64
export PATH=$PATH:$(pwd)
source ~/.bashrc
```

---

## Directory Structure

Create the following structure before running the pipeline:

```
Mapeamento/
├── Dados/
│   ├── ArquivosParaAnalise/     # Place your FASTQ files here (.fastq, .fastq.gz, .fastqsanger)
│   ├── SequenciadeReferencia/   # Place your reference FASTA (.fa) and annotation (.gtf) here
│   └── DadosdoStar/             # STAR genome index will be generated here
└── MapeamentosFeitos/           # Alignment outputs (BAM files) will be saved here
```

Run the commands below to create this structure:

```bash
mkdir -p Mapeamento/Dados/ArquivosParaAnalise
mkdir -p Mapeamento/Dados/SequenciadeReferencia
mkdir -p Mapeamento/Dados/DadosdoStar
mkdir -p Mapeamento/MapeamentosFeitos
cd Mapeamento
```

---

## Usage

### Step 1 — Quality control with FastQC

```bash
fastqc Dados/ArquivosParaAnalise/*
```

Expected output: one `.html` report per FASTQ file in the `ArquivosParaAnalise/` directory.

### Step 2 — Aggregate QC reports with MultiQC

```bash
multiqc Dados/ArquivosParaAnalise/*fastqc*
```

Expected output: a single `multiqc_report.html` aggregating all FastQC results.

### Step 3 — Generate STAR genome index

Create and run the indexing script:

```bash
cat > generate_DadosdoStar.sh << 'EOF'
#!/bin/bash
STAR \
  --runThreadN 4 \
  --limitGenomeGenerateRAM 2000000000 \
  --runMode genomeGenerate \
  --genomeDir Dados/DadosdoStar \
  --genomeFastaFiles Dados/SequenciadeReferencia/SequenciadeReferencia.fa \
  --sjdbGTFfile Dados/SequenciadeReferencia/SequenciadeReferencia.gtf \
  --genomeSAindexNbases 12
EOF
chmod +x generate_DadosdoStar.sh
./generate_DadosdoStar.sh
```

> **Note:** `--genomeSAindexNbases 12` is suitable for small genomes. For large genomes (e.g., human), use the default value (14) or omit this parameter.

### Step 4 — Align reads with STAR

Create and run the alignment script:

```bash
cat > generate_MapearRNAStar.sh << 'EOF'
#!/bin/bash
STAR \
  --runThreadN 4 \
  --readFilesIn Dados/ArquivosParaAnalise/Sample_R1.fastq Dados/ArquivosParaAnalise/Sample_R2.fastq \
  --genomeDir Dados/DadosdoStar \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix MapeamentosFeitos/Sample
EOF
chmod +x generate_MapearRNAStar.sh
./generate_MapearRNAStar.sh
```

> **Note:** Replace `Sample_R1.fastq` and `Sample_R2.fastq` with your actual file names. For gzip-compressed files, add `--readFilesCommand zcat`.

### Step 5 — Aggregate mapping statistics with MultiQC

```bash
multiqc MapeamentosFeitos/*Log.final.out
```

### Step 6 — Index BAM files with Samtools

```bash
samtools index MapeamentosFeitos/SampleAligned.sortedByCoord.out.bam
```

> Replace the filename with the actual BAM file produced by STAR.

---

## Output Files

| File | Description |
|---|---|
| `*_fastqc.html` | Per-sample FastQC quality report |
| `multiqc_report.html` | Aggregated QC report |
| `*Aligned.sortedByCoord.out.bam` | Coordinate-sorted BAM file |
| `*Aligned.sortedByCoord.out.bam.bai` | BAM index for IGV |
| `*Log.final.out` | STAR mapping statistics summary |
| `*SJ.out.tab` | Detected splice junctions |

---

## Visualizing Results in IGV

1. Install IGV: `sudo apt install igv` or download from [igv.org](https://igv.org)
2. Open IGV: `igv`
3. Load the reference genome: **Genomes → Load Genome from File** → select your `.fa` file
4. Load BAM files: **File → Load from File** → select `.bam` files from `MapeamentosFeitos/`
5. Optionally load the annotation: **File → Load from File** → select your `.gtf` file
6. Navigate to a gene of interest using the search bar

---

## Limitations

- Designed for Linux only; not tested on macOS or Windows (WSL)
- Processes one sample at a time; for batch processing, wrap Step 4 in a loop
- `--limitGenomeGenerateRAM 2000000000` limits RAM to ~2 GB for genome indexing; increase this value for large genomes
- Only standard paired-end Illumina reads are covered; stranded library options must be set manually

---

## License

MIT License — feel free to use, modify, and distribute with attribution.

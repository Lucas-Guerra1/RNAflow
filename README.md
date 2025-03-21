# RNAflow
Automated RNA-seq Analysis Pipeline
Overview

This project aims to create an automated pipeline for RNA-seq data analysis, integrating widely used bioinformatics tools available for Linux (e.g., FastQC, Cutadapt, MultiQC, STAR, etc.). The goal is to provide a seamless workflow where users only need to supply their FASTQ files, and the pipeline will handle the entire analysis process, saving results in the project directory.
Features

    Automated Workflow: From raw FASTQ files to processed results with minimal user input.

    Integrated Tools:

        FastQC: Quality control for raw sequencing data.

        Cutadapt: Adapter trimming and quality filtering.

        STAR: Alignment of reads to a reference genome.

        MultiQC: Aggregate and visualize quality control reports.

    Customizable: Easily modify parameters and add additional tools as needed.

    User-Friendly: Designed for both beginners and experienced bioinformaticians.

Requirements

    Linux-based operating system (Ubuntu, CentOS, etc.).

    Python 3.x.

    Required bioinformatics tools installed (see Installation for details).

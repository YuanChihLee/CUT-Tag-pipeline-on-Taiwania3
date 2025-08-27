# CUT&Tag Data Analysis Pipeline on Taiwania3

This repository contains an automated bioinformatics pipeline for processing CUT&Tag sequencing data. It takes raw sequencing files (FASTQ) and performs all essential steps, including quality control, read alignment, peak calling, and read quantification.

This pipeline is designed to be run on a high-performance computing (HPC) cluster, specifically **Taiwania3**, using the SLURM workload manager.

## 🚀 Overview

The pipeline automates the entire CUT&Tag analysis workflow, from start to finish. It is broken down into two main stages:

* **Stage 1: Pre-processing & Peak Calling** - Cleans up raw sequencing data, aligns reads to both the mouse (primary) and yeast (spike-in) genomes, and calls individual peaks. It also generates three types of BigWig files for data visualization.
* **Stage 2: Peak Merging & Quantification** - Identifies reproducible peaks across biological replicates and quantifies the number of reads within these consensus peaks for downstream differential analysis.

---

## 🔧 Prerequisites

Before you can run the pipeline, you need to ensure the following software is installed and available in your environment on the Taiwania3 cluster. This pipeline assumes you have access to the `module` and `conda` systems on the cluster.

* **Core Tools:**
    * [Trimmomatic](https://www.usadellab.org/cms/?page=trimmomatic)
    * [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
    * [SAMtools](http://www.htslib.org/)
    * [BEDTools](https://bedtools.readthedocs.io/en/latest/)
    * [Picard](https://broadinstitute.github.io/picard/)
    * [UCSC Utilities](https://genome.ucsc.edu/utilities.html) (`bedGraphToBigWig`)
    * [MACS3](https://github.com/macs3-project/MACS)
    * [SEACR](https://seacr.readthedocs.io/en/stable/)
    * [Subread (`featureCounts`)](http://subread.sourceforge.net/)
* **Parallelization:**
    * [GNU Parallel](https://www.gnu.org/software/parallel/)
* **Anaconda/Miniconda:**
    * You must create a Conda environment named **`featurecounts_env`** and install `featureCounts` into it.

---

## ⚙️ Setup

Follow these steps to configure the pipeline for your data.

### 1. File Structure

Organize your project directory with the following structure. The pipeline will create all the sub-directories automatically.
proj_root/
├── 00_raw_data/
│   └── *.fastq.gz   # Your raw sequencing files
├── CUTandTag_pkg/
│   ├── index/       # Genome indices (mouse/yeast)
│   ├── mm39.chrom.sizes
│   ├── mm39.excluderanges_2.bed
│   ├── gencode.vM36.annotation.gtf
│   ├── SEACR_1.3.sh
│   ├── default_adaptor.fa
│   └── ...
└── filelist.txt     # A list of your samples (see below)

### 2. `filelist.txt`

This file is crucial for the pipeline to identify and process your samples. It should be a tab-separated text file with the following columns:

`SampleID <tab> i5_Index <tab> i7_Index`

Example:
WT_H3K27me3_Rep1    i5_1    i7_1
WT_H3K27me3_Rep2    i5_2    i7_2
KO_H3K4me3_Rep1     i5_3    i7_3

### 3. `config.sh`

Create a file named **`config.sh`** in your main project directory. This file stores all customizable paths and settings, making the main script reusable. Copy the content below and update the paths to match your own file locations on the Taiwania3 cluster.

```bash
#!/bin/bash

# Main directory for your project
PROJ_ROOT="/work/j120885731/CUTandTAG/TS250710005"
# Directory containing genome indices, blacklist files, and other resources
PKG_DIR="/work/j120885731/CUTandTAG/CUTandTag_pkg"
# The file listing your samples
FILELIST="$PROJ_ROOT/filelist.txt"

# Path to MACS3 executable (change this if it's not in your home directory)
export MACS3_PATH="$HOME/.local/bin/macs3"

# Number of parallel jobs for peak counting. This should be a small number
# to prevent overloading the cluster (e.g., 5).
NUM_PARALLEL_JOBS=5

# Export variables to make them available to the main script
export PROJ_ROOT PKG_DIR FILELIST NUM_PARALLEL_JOBS
```

💻 Usage
Once you have set up your file structure and edited config.sh, simply submit the main pipeline script to the SLURM scheduler:
```bash
sbatch run_cutntag_pipeline.sh
```

📁 Output Files
All output files are organized into a clean directory structure within your proj_root folder:
```
01_trimmed_fastq/ : Cleaned FASTQ files.

02_alignment_bam/: Final, clean BAM files ready for analysis.

03_browser_tracks/:

raw_bw/

cpm_normalized_bw/

spike-in_normalized_bw/

04_peaks/:

SEACR_individual/

MACS3_individual/

MACS3_replicate_merged/

MACS3_mark_consensus/

05_quantification_counts/: Final *.txt count tables.

logs/: All log files and QC summaries.

QC_summary.txt: A summary file of key quality control metrics for each sample.
```

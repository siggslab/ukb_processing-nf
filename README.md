# UK Biobank Variant Processing, Filtering, and Data Analysis Workflow

## Overview

This repository provides the full workflow and scripts for the study **“Population Prevalence and Presentation of Rare Germline and Somatic Variants in Genes Associated with Dominant Inborn Errors of Immunity”**.

The workflow is implemented primarily in **Nextflow** and supports high-throughput processing, variant annotation, and downstream filtering and analysis of UK Biobank (UKB) whole-exome/whole-genome variant call files (VCFs).

The repository includes:
- Shell scripts for raw data extraction.
- Python scripts for sample sheet creation, post-annotation filtering, and downstream statistical analysis.
- Nextflow pipelines for variant annotation.
- Example submission scripts for running on HPC clusters (Gadi).

The main workflow consists of two stages:
- **Variant Processing & Annotation** (green in the diagram below).
- **Post-Annotation Filtering & Analysis** (blue in the diagram below).

![Workflow Diagram](https://github.com/user-attachments/assets/23070333-2977-4f25-9a84-e4ea8e08d8d7)

---

## Contents

- [`bin/`](./bin) — Shell scripts for data extraction and running Nextflow.
- [`scripts/`](./scripts) — Python scripts for sample sheet generation, post-annotation filtering, and data analysis.
- [`main.nf`](./main.nf) — Nextflow pipeline for VCF processing and annotation.
- [`results/`](./results) — Output directory for annotated VCFs and final result files.
- [`samplesheets/`](./samplesheets) — Generated sample sheets for batch processing.

---

## Prerequisites

Before running this workflow, ensure you have:
- Access to **UK Biobank raw VCFs** stored on your HPC filesystem.
- Access to a high-performance computing (HPC) cluster (this workflow was developed and tested on **Gadi**).
- Installed:
  - `Nextflow` (recommended version ≥ 21.10.0)
  - `Python 3.8+`
  - `bcftools`, `tar`, and other basic Linux utilities.
  - All required Python packages (see `requirements.txt` if available).

---

## User Guide

### 1. Extract UKB VCFs from Tar Archives

Raw UKB variant data is assumed to be delivered as tar archives distributed across multiple directories on **Gadi**.  
Use the extraction script to decompress and organize VCFs into batches.

**Usage:**
```bash
bash bin/run_extract_vcf.sh <source_directory> <start_index> <end_index>
```

- **source_directory:** Path containing the tar archives.
- **start_index**, **end_index:** Define the batch range to extract (e.g., 1–20).

Each batch will have:
- Its own subdirectory under `ukb_vcfs/`.
- Symlinks for easier management.

**Note:** Due to storage and compute constraints, extraction must be done in manageable batch sizes.

---

### 2. Generate Sample Sheets for Nextflow

After extraction, use the Python script to create Nextflow-compatible sample sheets for each batch.

**Example:**
```bash
for i in {1..20}; do
  python3 scripts/create_samplesheet.py ukb_vcfs/$i samplesheets/${i}_samplesheet.csv
done
```

- Each sample sheet includes up to 9,999 VCFs.
- If a batch contains more than 5,000 VCFs, the sample sheet is automatically split into parts (`part1`, `part2`).

---

### 3. Run Nextflow Annotation Pipeline

Use the provided HPC submission script (`bin/run_nf_gadi.sh`) to launch the Nextflow pipeline on each batch sample sheet.

**Example:**
```bash
for samplesheet in samplesheets/11_samplesheet_part*.csv; do
  ID=$(basename "$samplesheet" | sed -E 's/([0-9]+)_samplesheet_(part[0-9]+)\.csv/\1_\2/')
  qsub -v ID="$ID",samplesheet="$(realpath "$samplesheet")" bin/run_nf_gadi.sh
done
```

- The pipeline will produce processed and annotated VCFs stored in:
  ```
  results/<batchID>_part1/
  results/<batchID>_part2/
  ```

---

### 4. Post-Annotation Filtering

After processing all batches, combine the resulting TSV outputs from all parts into a single dataset.

Then run:
```bash
python3 scripts/post_annotation_filtering.py <combined_tsv> <output_directory>
```

This script:
- Filters annotated variants.
- Identifies participants with likely pathogenic/pathogenic (LP/P) variants.
- Produces final stratified output files for downstream analysis.

---

### 5. Downstream Analysis and Plot Generation

To reproduce the analyses and figures in the paper, use:
```bash
python3 scripts/data_analysis.py <filtered_data> <output_directory>
```

This script:
- Generates summary statistics and plots.
- Reproduces tables and figures reported in the manuscript.
- Outputs results in a format ready for supplementary material.

---

## Tools and Dependencies

The workflow uses:
- **Nextflow** for scalable, reproducible workflows.
- **Python 3** for sample sheet creation, filtering, and analysis.
- **bcftools** for VCF processing.
- **qsub** (PBS job submission) for HPC execution.
- Standard GNU/Linux utilities.

---

## Additional Notes

- **HPC Compatibility:** Scripts assume usage on the Gadi HPC cluster but can be adapted for other clusters supporting PBS or SLURM.
- **Storage Requirements:** Extracting and processing all UKB VCFs requires substantial storage (>50 TB recommended).
- **Compute Requirements:** It is strongly recommended to parallelize batches to leverage HPC resources efficiently.

---

## Acknowledgements

This workflow was developed as part of the study *“Population Prevalence and Presentation of Rare Germline and Somatic Variants in Genes Associated with Dominant Inborn Errors of Immunity”*. Please cite the corresponding publication when using this code.

**Contact:**  
For questions or contributions, please raise an issue or contact the corresponding author.

---

## Citation

If you use this repository or any of its components, please cite:  
**[Insert your citation here once published]**


#!/bin/bash
#PBS -q normal
#PBS -P np30
#PBS -l ncpus=6
#PBS -l mem=124G
#PBS -l walltime=15:00:00
#PBS -l storage=gdata/np30+scratch/np30+gdata/tn36+scratch/tn36
#PBS -l wd

set -euo pipefail

module load nextflow
module load singularity

if [ -d "$ID" ]; then
    echo "Directory $ID already exists, skipping creation and file copying."
else
    mkdir "$ID"
    cp -r main.nf nextflow.config config modules trace iuis.bed "$ID/"
fi

cd $ID

# Run the Nextflow pipeline with the specified samplesheet
nextflow run main.nf -profile gadi --vcf_csv $samplesheet && rm -r work


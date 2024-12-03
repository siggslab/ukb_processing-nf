#!/bin/bash
#PBS -q normal
#PBS -P tn36
#PBS -l ncpus=4
#PBS -l mem=16G
#PBS -l walltime=1:00:00
#PBS -l storage=gdata/tn36+scratch/tn36
#PBS -l wd

set -euo pipefail

module load nextflow
module load singularity

nextflow run main.nf -profile gadi
#!/bin/bash
#PBS -q normal
#PBS -P np30
#PBS -l ncpus=4
#PBS -l mem=64G
#PBS -l walltime=5:00:00
#PBS -l storage=gdata/np30+scratch/np30
#PBS -l wd

set -euo pipefail

module load python3

python3 scripts/extract_all_vcfs.py
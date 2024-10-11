#!/bin/sh
#SBATCH --time=2:00:00
#SBATCH --mem=12G
#SBATCH --partition=caslake
#SBATCH --account=pi-lbarreiro

cd /project/lbarreiro/USERS/daniel/chrom_pred/chrom_pred/01_prepare_inputs
module load R
Rscript 01b_format_dosages.R

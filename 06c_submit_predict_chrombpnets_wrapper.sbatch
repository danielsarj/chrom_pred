#!/bin/sh
#SBATCH --time=36:00:00
#SBATCH --mem=8G
#SBATCH --partition=caslake
#SBATCH --account=pi-lbarreiro

cd /project/lbarreiro/USERS/daniel/deeplearn_pred/github
bash 06b_submit_predict_chrombpnets.sh

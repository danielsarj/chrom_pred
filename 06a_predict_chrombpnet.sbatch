#!/bin/sh
#SBATCH --time=36:00:00
#SBATCH --mem=64G
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --account=pi-lbarreiro

module load python
module load cudnn/8.3.1
conda init
conda activate chrombpnet

Conditions=("NI" "Flu")

ID=$1
for COND in ${Conditions[@]}; do
        for i in {0..4}; do

         chrombpnet pred_bw -cmb /project/lbarreiro/USERS/daniel/deeplearn_pred/DIY/models/${COND}_fold${i}/models/chrombpnet_nobias.h5 -r /project/lbarreiro/USERS/daniel/deeplearn_pred/DIY/ATACseq_peaks.bed -g /project/lbarreiro/USERS/daniel/katieData/WGS/custom_refs/${ID}_sequence.fa -c /project/lbarreiro/USERS/daniel/katieData/WGS/chr_length.txt -op /project/lbarreiro/USERS/daniel/deeplearn_pred/DIY/prediction/chrombpnet_nobias_model/${COND}_fold${i}_${ID}

        chrombpnet pred_bw -cm /project/lbarreiro/USERS/daniel/deeplearn_pred/DIY/models/${COND}_fold${i}/models/chrombpnet.h5 -r /project/lbarreiro/USERS/daniel/deeplearn_pred/DIY/ATACseq_peaks.bed -g /project/lbarreiro/USERS/daniel/katieData/WGS/custom_refs/${ID}_sequence.fa -c /project/lbarreiro/USERS/daniel/katieData/WGS/chr_length.txt -op /project/lbarreiro/USERS/daniel/deeplearn_pred/DIY/prediction/chrombpnetbias_model/${COND}_fold${i}_${ID}

        done
done

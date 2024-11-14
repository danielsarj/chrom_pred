#!/bin/bash

Conditions=("NI")
for COND in ${Conditions[@]}; do
        for i in {0..4}; do

        SCRIPT_FILE=${COND}"_fold"${i}"_chrombpnetmodel.sbatch"
        cat > "$SCRIPT_FILE" << EOF
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
cd /project/lbarreiro/USERS/daniel/deeplearn_pred/DIY

chrombpnet pipeline -g hg19_genome.fa -c chr_length.txt -ibam AF04_${COND}.sorted.dup.auto.bam -o models/${COND}_fold${i}/ -d "ATAC" -p ATACseq_peaks.bed -n nonpeaks_fold${i}_negatives.bed -fl splits/fold_${i}.json -b bias/${COND}_fold${i}/models/${COND}_fold${i}_bias.h5

EOF
        sbatch $SCRIPT_FILE
        done
done


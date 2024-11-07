cd /project/lbarreiro/USERS/daniel/deeplearn_pred/DIY

for i in {0..4}
do
	echo $i
	chrombpnet prep nonpeaks -g hg19_genome.fa -o nonpeaks_fold$i -p ATACseq_peaks.bed -c chr_length.txt -fl splits/fold_$i.json
done


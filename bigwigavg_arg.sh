#!/bin/bash
cd /project/lbarreiro/USERS/daniel/deeplearn_pred/DIY/prediction

IDs=("AF04" "AF06" "AF10" "AF12" "AF14" "AF16" "AF18" "AF20" "AF22" "AF24" "AF26" "AF28" "AF30" "AF34" "AF36" "AF38" "EU03" "EU05" "EU07" "EU09" "EU13" "EU15" "EU17" "EU19" "EU21" "EU25" "EU27" "EU29" "EU33" "EU37" "EU39" "EU41" "EU43" "EU47")
Conditions=("NI" "Flu")

for i in ${IDs[@]}; do
for c in ${Conditions[@]}; do
for f in {0..4}; do

echo ${c}"_"${f}"_"${i}
awk '{print $1, $2, $3, $1":"$2":"$3}' ${c}_fold${f}_${i}_chrombpnet_nobias_preds.bed > ${c}_fold${f}_${i}_chrombpnet_nobias_preds.short.bed
../../../SOFTWARE/bigWigAverageOverBed ${c}_fold${f}_${i}_chrombpnet_nobias.bw ${c}_fold${f}_${i}_chrombpnet_nobias_preds.short.bed ${c}_fold${f}_${i}_bwavg.txt

done
done
done

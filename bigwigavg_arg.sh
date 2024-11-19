#!/bin/bash
cd /project/lbarreiro/USERS/daniel/deeplearn_pred/DIY/prediction

IDs=("AF04" "AF06" "AF08")
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

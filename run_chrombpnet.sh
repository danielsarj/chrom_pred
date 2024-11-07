#!/bin/bash
cd /project/lbarreiro/USERS/daniel/deeplearn_pred/chrombpnet/test_prediction/outputs

InputIDs=("AF04" "AF06" "AF08" "AF10" "AF12" "AF14" "AF16" "AF18" "AF20" "AF22" "AF24" "AF26" "AF28" "AF30" "AF34" "AF36" "AF38" "EU03" "EU05" "EU07" "EU09" "EU13" "EU15" "EU17" "EU19" "EU21" "EU25" "EU27" "EU29" "EU33" "EU37" "EU39" "EU41" "EU43" "EU47")

for ID in ${InputIDs[@]}; do

        echo ${ID}
        awk '{print $1, $2, $3, $1"_"$2"_"$3}' ${ID}_bias_preds.bed > ${ID}_bias_preds.short.bed
        ../../../../SOFTWARE/bigWigAverageOverBed ${ID}_bias.bw ${ID}_bias_preds.short.bed ${ID}_bwavg.txt

done


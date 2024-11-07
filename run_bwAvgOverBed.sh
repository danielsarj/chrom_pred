#!/bin/bash
cd /project/lbarreiro/USERS/daniel/deeplearn_pred/chrombpnet/test_prediction/outputs

InputIDs=("Epi_AF04_NI_WGS_1" "Epi_AF06_NI_WGS_1" "Epi_AF08_NI_WGS_1" "Epi_AF10_NI_WGS_1" "Epi_AF12_Flu_WGS_1" "Epi_AF14_Mock_WGS_1" "Epi_AF16_NI_WGS_1" "Epi_AF18_NI_WGS_1" "Epi_AF20_NI_WGS_1" "Epi_AF22_Flu_WGS_1" "Epi_AF24_Flu_WGS_1" "Epi_AF26_NI_WGS_1" "Epi_AF28_NI_WGS_1" "Epi_AF30_NI_WGS_1" "Epi_AF34_NI_WGS_1" "Epi_AF36_Flu_WGS_1" "Epi_AF38_NI_WGS_1" "Epi_EU03_NI_WGS_1" "Epi_EU05_NI_WGS_1" "Epi_EU07_NI_WGS_1" "Epi_EU09_Mock_WGS_1" "Epi_EU13_NI_WGS_1" "Epi_EU15_NI_WGS_1" "Epi_EU17_NI_WGS_1" "Epi_EU19_NI_WGS_1" "Epi_EU21_Flu_WGS_1" "Epi_EU25_NI_WGS_1" "Epi_EU27_NI_WGS_1" "Epi_EU29_NI_WGS_1" "Epi_EU33_Flu_WGS_1" "Epi_EU37_Flu_WGS_1" "Epi_EU39_NI_WGS_1" "Epi_EU41_NI_WGS_1" "Epi_EU43_NI_WGS_1" "Epi_EU47_NI_WGS_1")
for ID in ${InputIDs[@]}; do
	
	echo ${ID}
	awk '{print $1, $2, $3, $1":"$2":"$3}' ${ID}_bias_preds.bed > ${ID}_bias_preds.short.bed
	../../../../SOFTWARE/bigWigAverageOverBed ${ID}_bias.bw ${ID}_bias_preds.short.bed ${ID}_bwavg.txt

done

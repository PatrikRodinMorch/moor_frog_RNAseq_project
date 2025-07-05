#!/bin/bash


tissue="liver
muscle"

pop="across
P1
P2
P3
P4
P5"

for m in $tissue;do
	for p in $pop;do

		/home/prm/Skrivbord/Final_RNAseq_small_scale_results/leafcutter/leafcutter/leafviz/prepare_results.R ../${m}/${p}_intron_clustering_perind_numers.counts.FINAL.gz ./${m}/${p}/leafcutter_ds_cluster_significance.txt ./${m}/${p}/leafcutter_ds_effect_sizes.txt temporaria_annotation_code -m ../${m}/${p}_${m}_metadata_final.txt -o ${m}/${p}/${p}.RData
	done

done




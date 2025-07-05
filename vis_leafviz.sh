#!/bin/bash


tissue="liver
muscle"

pop="P1
P2
P3
P4
P5"

for m in $tissue;do
	for p in $pop;do

		../leafcutter/leafviz/run_leafviz.R /home/prm/Skrivbord/Final_RNAseq_small_scale_results/leafcutter/results/${m}/${p}/${p}.RData
	done

done




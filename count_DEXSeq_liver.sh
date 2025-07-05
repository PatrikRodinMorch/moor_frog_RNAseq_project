#!/bin/bash -l


module load bioinfo-tools
module load R/4.1.1
module load R_packages/4.1.1

python /sw/apps/R_packages/4.1.1/rackham/DEXSeq/python_scripts/dexseq_count.py -p yes -r name -s reverse -f bam GCF_905171775.1_aRanTem1.1_genomic.gff /crex/proj/snic2019-30-54/reference_guided/liver/$1 /crex/proj/snic2019-30-54/reference_guided/liver/htseq/dexseq/liver/$1.txt

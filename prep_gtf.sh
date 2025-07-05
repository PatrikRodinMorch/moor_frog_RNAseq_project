#!/bin/bash -l


module load bioinfo-tools
module load R/4.1.1
module load R_packages/4.1.1
module load htseq

python /sw/apps/R_packages/4.1.1/rackham/DEXSeq/python_scripts/dexseq_prepare_annotation.py /crex/proj/snic2019-30-54/reference_guided/GCF_905171775.1_aRanTem1.1_genomic.gtf GCF_905171775.1_aRanTem1.1_genomic.gff

#!/bin/bash -l


module load bioinfo-tools
module load samtools/1.5

/sw/apps/bioinfo/GATK/4.1.1.0/rackham/gatk --java-options "-Xmx100G" SplitNCigarReads \
-R /crex/proj/naiss2023-23-494/patrik/reference_guided/GCF_905171775.1_aRanTem1.1_genomic.fna \
--create-output-bam-index false \
-I $1.duplMarked_RG.bam \
-O $1.duplMarked_RG_splitreads.bam

